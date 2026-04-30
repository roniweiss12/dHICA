import random
import tensorflow as tf
from tqdm import tqdm
import numpy as np
import os
import json
import functools

from model_dHICA import dHICA
import time
from operator import itemgetter
import subprocess
from optparse import Option, OptionParser
#os.environ["CUDA_VISIBLE_DEVICES"] = '2'
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"

def get_metadata(data_path):
    # Keys:
    # num_targets, train_seqs, valid_seqs, test_seqs, seq_length,
    # pool_width, crop_bp, target_length
    path = os.path.join(data_path, 'statistics.json')
    with tf.io.gfile.GFile(path, 'r') as f:
        return json.load(f)

def tfrecord_files(data_path, subset):
    # Sort the values by int(*).
    return sorted(tf.io.gfile.glob(os.path.join(
        data_path, 'tfrecords', f'{subset}-*.tfr'
    )), key=lambda x: int(x.split('-')[-1].split('.')[0]))


def deserialize(serialized_example, metadata):
    feature_map = {
        'sequence': tf.io.FixedLenFeature([], tf.string),
        'atac-seq': tf.io.FixedLenFeature([], tf.string),
        'start-end': tf.io.FixedLenFeature([], tf.string)
    }
    example = tf.io.parse_example(serialized_example, feature_map)

    sequence = tf.io.decode_raw(example['sequence'], tf.bool)
    sequence = tf.reshape(sequence, (metadata['seq_length'], 4))
    sequence = tf.cast(sequence, tf.float32)

    atac = tf.io.decode_raw(example['atac-seq'], tf.float16)
    atac = tf.reshape(atac, (metadata['atac_length'], metadata['num_atacseq']))
    atac = tf.cast(atac, tf.float32)

    start_end = tf.io.decode_raw(example['start-end'], tf.int32)
    start_end = tf.reshape(start_end, (1, 3))
    start_end = tf.cast(start_end, tf.int32)

    return {
        'sequence': sequence,
        'atac-seq': atac,
        'start-end': start_end
    }

def get_dataset(data_path, subset, num_threads=8):
    metadata = get_metadata(data_path)
    dataset = tf.data.TFRecordDataset(tfrecord_files(data_path, subset),
                                      compression_type='ZLIB',
                                      num_parallel_reads=num_threads)
    # print(tfrecord_files(organism, subset))
    dataset = dataset.map(functools.partial(deserialize, metadata=metadata),
                          num_parallel_calls=num_threads)
    return dataset


def write_bedGraph(results, out_path, chr_length_human, pre_out, is_chr22=False):
    from operator import itemgetter
    import numpy as np
    import os

    # Open all files once
    file_handles = {}
    buffers = {}

    for histone in histone_list:
        path = os.path.join(out_path, f"{pre_out}-{histone}.bedgraph")
        file_handles[histone] = open(path, 'w')
        buffers[histone] = []

    BUFFER_LIMIT = 100000  # tune if needed

    def flush(histone):
        if buffers[histone]:
            file_handles[histone].writelines(buffers[histone])
            buffers[histone] = []

    for chr in results:
        chr_result = sorted(results[chr], key=itemgetter('start'))

        if not chr_result:
            continue

        for j, histone in enumerate(histone_list):

            # Initial padding
            if is_chr22:
                if chr == 'chr22' and chr_result[0]['start'] > 0:
                    buffers[histone].append(f"{chr}\t0\t{chr_result[0]['start']}\t0\n")
            else:
                if chr_result[0]['start'] > 0:
                    buffers[histone].append(f"{chr}\t0\t{chr_result[0]['start']}\t0\n")

            last_end = 0

            for item in chr_result:
                start_base = item['start']
                end_base = item['end']
                preds = item['predicted'][:, j]

                starts = start_base + np.arange(896) * 128
                ends = starts + 128

                if is_chr22:
                    # simpler path
                    for s, e, v in zip(starts, ends, preds):
                        buffers[histone].append(f"{chr}\t{s}\t{e}\t{v}\n")
                else:
                    if start_base >= last_end:
                        for s, e, v in zip(starts, ends, preds):
                            buffers[histone].append(f"{chr}\t{s}\t{e}\t{v}\n")
                    else:
                        gap_h = last_end - start_base
                        h_start = gap_h // 128

                        # overlap correction line
                        buffers[histone].append(
                            f"{chr}\t{last_end}\t{start_base + 128 * (h_start+1)}\t{preds[h_start]}\n"
                        )

                        for i in range(h_start + 1, 896):
                            s = starts[i]
                            e = ends[i]
                            v = preds[i]
                            buffers[histone].append(f"{chr}\t{s}\t{e}\t{v}\n")

                last_end = end_base

                # flush periodically
                if len(buffers[histone]) >= BUFFER_LIMIT:
                    flush(histone)

            # Final padding
            try:
                if chr_result[-1]['end'] < chr_length_human[chr]:
                    buffers[histone].append(
                        f"{chr}\t{chr_result[-1]['end']}\t{chr_length_human[chr]}\t0\n"
                    )
            except:
                print('an error', chr)

            flush(histone)

    # Close all files
    for histone in histone_list:
        flush(histone)
        file_handles[histone].close()


def evaluate_model(model, dataset, chr_length_human, ID_to_chr_dict, max_steps=None):

    def predict(batch):
        return model(batch['sequence'], batch['atac-seq'], is_training=False)['human']

    @tf.function
    def distributed_predict_step(dist_inputs):
        predicted = mirrored_strategy.run(predict, args=(dist_inputs,))
        return predicted
    # len_dataset = 0

    results = {}
    for chr in chr_length_human:
        results[chr] = []

    for i, batch in tqdm(enumerate(dataset)):
        if max_steps is not None and i > max_steps:
            break
        predicted = distributed_predict_step(batch)
        predicted = predicted.numpy()
        start_end = batch['start-end'].numpy()

        for i in range(len(predicted[0])):
            result = {}
            # len_dataset += 1
            try:
                chr_id = start_end[i][0][0]
            except:
                continue

            chr = ID_to_chr_dict[str(chr_id)]
            result['start'] = start_end[i][0][1] + 40960

            result['end'] = start_end[i][0][2] - 40960
            result['predicted'] = predicted[i]
            results[chr].append(result)

    return results

mirrored_strategy = tf.distribute.MirroredStrategy()

def main():
    usage = 'usage: %prog [options] <fasta_file> <targets_file>'
    parser = OptionParser(usage)
    parser.add_option('-d', dest='dataset',
      default='data_set',
      help='Dataset path [Default: %default]')
    parser.add_option('--model_type', dest='modelType',
      default='R', type='str',
      help='model you choosed R:roseq-only, D:DNA-only, RD:roseq+DNA [Default: %default]')
    parser.add_option('--idx', dest='index',
      default='idx_set',
      help='index file of bedGrapgToBigwig [Default: %default]')
    parser.add_option('--idx_chr', dest='indexChr',
      default='idx_set',
      help='index file of chromosome [Default: %default]')
    parser.add_option('-o', dest='outpath',
      default='outpath',
      help='Output path [Default: %default]')
    parser.add_option('-m', dest='model_path',
      default='model_path',
      help='Model Path(contain R, D, RD models) [Default: %default]')
    parser.add_option('--op', dest='output_pre',
      default='out',
      help='Prefix for output file [Default: %default]')

    (options, args) = parser.parse_args()

    # chr dict
    chr_length_file = options.indexChr
    chr_length_human = make_length_dict(chr_length_file)
    ID_to_chr_dict = make_chr_id(chr_length_file)

    # make output dir
    out_path = options.outpath
    if not(os.path.isdir(out_path)):
        os.mkdir(out_path)

    # set batchsize
    batch_size_replica = 16
    global_batch_size = batch_size_replica * mirrored_strategy.num_replicas_in_sync

    # choose model
    with mirrored_strategy.scope():
        model = dHICA(channels=768,
                             num_heads=8,
                             num_transformer_layers=11,
                             pooling_type='max')
        weight_path = os.path.join(options.model_path, 'model.ckpt')
        model.load_weights(weight_path)

    # load dataset
    dataset = (get_dataset(options.dataset, 'train', num_threads=8).batch(global_batch_size, drop_remainder=False).prefetch(tf.data.AUTOTUNE))
    distributed_dataset = mirrored_strategy.experimental_distribute_dataset(dataset)

    # predict
    results = evaluate_model(model,
                            distributed_dataset,
                            chr_length_human=chr_length_human,
                            ID_to_chr_dict=ID_to_chr_dict,
                            max_steps=100000)

    pre_out = options.output_pre + '-' + options.modelType

    # write results to file.bedgraph
    #print(results)
    write_bedGraph(results, out_path, chr_length_human, pre_out, is_chr22=False)
    
    # .bedgraph to .bw
    for histone in histone_list:
        bedgraph_path = os.path.join(out_path, pre_out + '-' + histone + '.bedgraph')

        bedgraph_path_sorted = os.path.join(out_path, pre_out + '-' + histone + '_sorted.bedgraph')
        bw_path = os.path.join(out_path, pre_out + '-' + histone + '.bw')
        hg19_idx = options.index
        cmd_bedSort = 'sort-bed ' + bedgraph_path + ' > ' + bedgraph_path_sorted
        p = subprocess.Popen(cmd_bedSort, shell=True)
        p.wait()

        cmd = ['bedGraphToBigWig', bedgraph_path_sorted, hg19_idx, bw_path]
        subprocess.call(cmd)

        cmd_rm = ['rm', '-f', bedgraph_path]
        subprocess.call(cmd_rm)

        cmd_rm = ['rm', '-f', bedgraph_path_sorted]
        subprocess.call(cmd_rm)


def make_length_dict(path):
  length_dict = {}
  for line in open(path):
    a = line.split()
    length_dict[a[0]] = int(a[2])
  return length_dict


def make_chr_id(path):
  id_dict = {}
  for line in open(path):
    a = line.split()
    id_dict[a[4]] = a[0]
  return id_dict


# sort-bed H3k4me1.bedgraph > H3k4me1_sorted.bedgraph
# bedGraphToBigWig H3k79me2.bedgraph ../../hg19/hg19.fa.fai H3k79me2.bw
histone_list = ['H3K122ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K9ac', 'H3K9me3', 'H4K20me1']




if __name__ == '__main__':
    main()
