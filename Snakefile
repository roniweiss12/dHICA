from snakemake.io import glob_wildcards, directory

# Load the config.yaml file correctly using the `configfile` directive
configfile: "config.yaml"

# Access the configuration
SAMPLES, = glob_wildcards("input_bam/{sample}.mLb.clN.sorted.bam")

BAMDIR = "input_bam"
BWDIR = "input_bw"
PRED_OUTDIR = "output"

rule all:
    input:
        expand(f"{PRED_OUTDIR}" + "/{sample}", sample=SAMPLES)

rule index_bam:
    input:
        bam=f"{BAMDIR}" + "/{sample}.mLb.clN.sorted.bam"
    output:
        bai=f"{BAMDIR}" + "/{sample}.mLb.clN.sorted.bam.bai"
    threads: 4
    conda:
        "envs/deept.yaml"
    shell:
        r"""
        samtools index -@ {threads} {input.bam}
        """

rule bam_to_bw:
    input:
        bam=f"{BAMDIR}" + "/{sample}.mLb.clN.sorted.bam",
        bai=f"{BAMDIR}" + "/{sample}.mLb.clN.sorted.bam.bai"
    output:
        bw=f"{BWDIR}" + "/{sample}.bw"
    threads: 24
    conda:
        "envs/deept.yaml"
    shell:
        r"""
        mkdir -p {BWDIR}
        bamCoverage \
            -b {input.bam} \
            -o {output.bw} \
            --binSize 10 \
            --normalizeUsing None \
            -p {threads}
        """
rule predict:
    input:
        bw=f"{BWDIR}" + "/{sample}.bw"
    output:
        directory(f"{PRED_OUTDIR}" + "/{sample}")
    threads: 24
    conda:
        "envs/dhica.yaml"
    shell:
        r"""
        mkdir -p logs

        CONDA_ENV=$(conda info --envs | grep dhica | awk '{{print $NF}}')
        CONDA_PYTHON=$CONDA_ENV/bin/python
        SCRIPT=logs/{wildcards.sample}_predict.sh

        cat > $SCRIPT << 'HEREDOC'
#!/bin/bash
export PATH=__CONDA_ENV__/bin:$PATH
mkdir -p __OUTPUT__
__CONDA_PYTHON__ predict_code/predict.py \
    -m dna_atac/ \
    -o __OUTPUT__/ \
    --atac __INPUT_BW__ \
    --ref GRCh38_chr.fa \
    -p __THREADS__
HEREDOC

        sed -i "s|__CONDA_ENV__|$CONDA_ENV|g" $SCRIPT
        sed -i "s|__CONDA_PYTHON__|$CONDA_PYTHON|g" $SCRIPT
        sed -i "s|__OUTPUT__|{output}|g" $SCRIPT
        sed -i "s|__INPUT_BW__|{input.bw}|g" $SCRIPT
        sed -i "s|__THREADS__|{threads}|g" $SCRIPT

        chmod +x $SCRIPT

        sbatch --job-name={config[slurm][job_name]} \
               --partition={config[slurm][partition]} \
               --cpus-per-task={threads} \
               --mem={config[slurm][memory]} \
               --time={config[slurm][time]} \
               --output=logs/{wildcards.sample}.out \
               --error=logs/{wildcards.sample}.err \
               --wait \
               $SCRIPT
        """