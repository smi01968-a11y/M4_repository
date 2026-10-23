#!/bin/bash

#SBATCH --job-name=fastqc_job

#SBATCH --time=6:00:00

#SBATCH --mem=8G

#SBATCH --cpus-per-task=4

#SBATCH --output=fastqc_%j.log

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu


#before running, make sample list file

#each line contains the sample name, read group, and the lane number

#how to run

#sbatch --array=1-196%20 --export=LIST=samples.txt 2fastqc.slurm


#runs all jobs in parallel

#STEP1 = checking the quality of raw reads with fastqc

#fastQC module .12.1 load
module load fastqc/0.12.1

SAMPLE="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST} | cut -f 1)"

#location of the paired end reads that have been trimmed already
FASTQ_DIR="/scratch.global/espin164/Evan_MEE/M4"

#run fastqc on raw fastq
OUT_DIR="/scratch.global/espin164/Evan_MEE/M4/newest_fastqc_results"

#make output folder since it doesn't exist yet:
mkdir -p "$OUT_DIR"

#runs fastqc on the TRIMMED fastq data, only the paied reads

fastqc \
"$FASTQ_DIR/${SAMPLE}_R1_paired.fastq.gz" \
"$FASTQ_DIR/${SAMPLE}_R2_paired.fastq.gz" \
--outdir "$OUT_DIR" \
--threads "$SLURM_CPUS_PER_TASK"

#what an example sample looks like
#MEE00074-D-1-1-MEE00074-D-1-1-CGGACCCGGACC_R1_paired.fastq.gz

#fastqc -f fastq trimmed/${SAMPLE}_R1_paired.fastq -o newest_fastqc_results/

#fastqc -f fastq trimmed/${SAMPLE}_R2_paired.fastq -o newest_fastqc_results/

