#!/bin/bash

#SBATCH --job-name=fastqc_job

#SBATCH --time=1:00:00

#SBATCH --mem=8G

#SBATCH --cpus-per-task=4

#SBATCH --output=fastqc_%A_%a.log

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu


#Create sample list before running the fastqc analysis, named samples.txt

#Each line contains the sample name, read group, and lane number ---  tab separated

#to run:

#sbatch --array=1-196%20 --export=LIST=samples.txt 2fastqc.slurm



#fastQC module

module load fastqc/0.11.9


# location of raw FASTQ files

FASTQ_DIR="/projects/standard/stuparr/data_delivery/umgc/2026-q3/20260722_AV241910_Aviti3_Run_417_07_22_26/Stupar_Project_076-merged"



# Location where you want FastQC results

OUT_DIR="/scratch.global/espin164/Evan_MEE/M4/Fastqc_results"



# Make output directory

mkdir -p "$OUT_DIR"



# Get sample name from the appropriate line of the sample list

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | cut -f1)"



echo "Processing sample: $SAMPLE"

echo "R1: $FASTQ_DIR/${SAMPLE}_R1.fastq.gz"

echo "R2: $FASTQ_DIR/${SAMPLE}_R2.fastq.gz"



# Run FastQC

fastqc \
"$FASTQ_DIR/${SAMPLE}_R1.fastq.gz" \
"$FASTQ_DIR/${SAMPLE}_R2.fastq.gz" \
--outdir "$OUT_DIR" \
--threads "$SLURM_CPUS_PER_TASK"
