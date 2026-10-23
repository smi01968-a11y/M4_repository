#!/bin/bash
#SBATCH --job-name=coverage_array
#SBATCH --output=coverage_%A_%a.out
#SBATCH --error=coverage_%A_%a.err
#SBATCH --array=1-196
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=smi01968@umn.edu



module load samtools/1.21

# Get list of MQ20.bam files
#assumes you are already in the directory with the MQ20.bam files, and no other bam files
BAM_LIST=(*.bam)

# Select BAM for this task
BAM=${BAM_LIST[$SLURM_ARRAY_TASK_ID-1]}

# Output file per sample
OUT=${BAM%.bam}_coverage2.txt

#average depth across the bases with coverage
#mean covered depth
#minimum base quality of 20
samtools depth -Q 20 -G 1024 "$BAM" | \
awk '{sum += $3; n++} END {print file, sum/n}' \
> "$OUT"

#combine results after this with cat *_coverage2.txt > coverage_summary2.txt
