dHICA
===============

a deep transformer-based model enables accurate histone imputation from chromatin accessibility

Cloud Computing Service:
-------------------------
We provide a computational gateway to run dHICA on GPU server, the users don't need to install any software, only upload the bigWig files and wait for the results, it is simple and easy. Please click the link to try this site:

https://dreg.dnasequence.org/


![Hi](https://github.com/Danko-Lab/dREG/raw/master/dreg-gateway.png?v=4&s=200 "dREG gateway")


Abstract:
--------
Histone modifications (HMs) are pivotal in various biological processes, including transcription, replication, and DNA repair, significantly impacting chromatin structure. These modifications underpin the molecular mechanisms of cell-type-specific gene expression and complex diseases. However, annotating HMs across different cell types solely using experimental approaches is impractical due to cost and time constraints. Herein, we present dHICA (deep histone imputation using chromatin accessibility), a novel deep learning framework that integrates DNA sequences and chromatin accessibility data to predict multiple HM tracks. Employing the transformer architecture alongside dilated convolutions, dHICA boasts an extensive receptive field and captures more cell-type-specific information. dHICA outperforms state-of-the-art baselines and achieves superior performance in cell-type-specific loci and gene elements, aligning with biological expectations. Furthermore, dHICA’s imputations hold significant potential for downstream applications, including chromatin state segmentation and elucidating the functional implications of SNPs (Single Nucleotide Polymorphisms). In conclusion, dHICA serves as a valuable tool for advancing the understanding of chromatin dynamics, offering enhanced predictive capabilities and interpretability.

Data preparation: 
==========================

dHICA takes bigWig files with single strands（ATAC-seq or DNase-seq） as the input. Performing data normalization beforehand is **NOT** recommended.


Download instruction: 
==========================
Install dHICA
------------
dHICA's source code is availiable on GitHub (https://github.com/wzhy2000/dHICA).  

Required software
-----------------
* Python 3.9
* Tensorflow 2.13.0 (https://www.tensorflow.org/)
* pyBigwig (https://github.com/deeptools/pyBigWig)


Get the dHICA models
-------------------
Pre-trained model that can be used to predict dREG scores across the genome is availiable here( https://dreg.dnasequence.org/themes/dreg/assets/file/dHICA_model.zip ).
If you are failed to download the model files, please contact us.

Usage instruction:
===================

## 1 Predict HMs with Pre-trained Models
First you need to make sure you have download the our Pre-trained Models. The type of ATAC-seq should be fold over change, and that of DNase-seq should be read-depth normalized signal.

    python predict.py -m model_path  -o output_path --atac bw_path --ref ref_genome.fa

    model_path      -- The path to the pre-trained model(DNase or ATAC).
    output_path     -- The path to save the output files.
    ref_genome.fa   -- Reference genome.(hg19 or mm10)
    bw_path         -- Read counts (not normalized) formatted as a bigWig or bw file.

After the command execution, you will obtain 10 bigWig (bw) files for different histone modifications. You can visualize and inspect these files using tools such as IGV (Integrative Genomics Viewer) or Genome Browser. Alternatively, you can utilize the Python package pyBigWig to read the signal data.

## 2 Training Your Own Model
### 1） Generate Dataset
The type of ATAC-seq should be fold over change, and that of DNase-seq should be read-depth normalized signal.

    python dHICA_data.py -l seq_length --local -o ouput_path ref_genome.fa seq.bw target_HM.txt
    
    seq_length      -- [default=131072] Sequence length.
    ouput_path      -- The path to save the output files.
    ref_genome.fa   -- Reference genome (hg19).
    seq.bw          -- The bigWig file of ATAC-seq or DNase-seq.
    target_HM.txt   -- The HMs target file.

Here is an example of target_HM.txt, please replace the 'file' with your own path for HM:
| index | identifier | file                      | clip | sum_stat | description      |
|-------|------------|---------------------------|------|----------|------------------|
| 0     | H3K122ac   | /local/w/H3K122ac.bw      | 384  | mean     | MEL-GATA-1-ER    |
| 1     | H3K4me1    | /local/w/H3K4me1.bw       | 384  | mean     | MEL-GATA-1-ER    |
| 2     | H3K4me2    | /local/w/H3K4me2.bw       | 384  | mean     | MEL-GATA-1-ER    |
| 3     | H3K4me3    | /local/w/H3K4me3.bw       | 384  | mean     | MEL-GATA-1-ER    |
| 4     | H3K27ac    | /local/w/H3K27ac.bw       | 384  | mean     | MEL-GATA-1-ER    |
| 5     | H3K27me3   | /local/w/H3K27me3.bw      | 384  | mean     | MEL-GATA-1-ER    |
| 6     | H3K36me3   | /local/w/H3K36me3.bw      | 384  | mean     | MEL-GATA-1-ER    |
| 7     | H3K9ac     | /local/w/H3K9ac.bw        | 384  | mean     | MEL-GATA-1-ER    |
| 8     | H3K9me3    | /local/w/H3K9me3.bw       | 384  | mean     | MEL-GATA-1-ER    |
| 9     | H4K20me1   | /local/w/H4K20me1.bw      | 384  | mean     | MEL-GATA-1-ER    |

### 2） Train Model 
    python train.py --data_path train_data

    train_data   -- The output path of 1) Generate Dataset

**Notice:** 
That command takes more than 8 hours (depends on size of the train datasets) to execute on NVIDA 3090 GPU. Due to very long computational time, we don't suggest to run on CPU nodes.

### 3） Calculate Performance
    python correlation.py -a dHICA.bw -b ref.bw -p resolution
    
    dHICA.bw    -- dHICA's output.
    ref.bw      -- the HMs target.
    resolution  -- a figure(128, 1000 or 10000) or a bed file.

How to cite
===================
Wen Wen, Jiaxin Zhong, Zhaoxi Zhang, Lijuan Jia, Tinyi Chu, Nating Wang, Charles G Danko, Zhong Wang, dHICA: a deep transformer-based model enables accurate histone imputation from chromatin accessibility, Briefings in Bioinformatics, Volume 25, Issue 6, November 2024, bbae459, https://doi.org/10.1093/bib/bbae459
