#!/bin/bash

#SBATCH --job-name=filtermq

#SBATCH --output=filtermq_%A_%a.out

#SBATCH --error=filtermq_%A_%a.err

#SBATCH --time=8:00:00

#SBATCH --mem=32G

#SBATCH --cpus-per-task=2

#SBATCH --array=1-196     #adjust number of samples

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu

#to run this script:
#sbatch filtermq.slurm

module load samtools/1.21

which samtools



#directories

DEDUP_DIR="/scratch.global/espin164/Evan_MEE/M4/deduped"

FILTER_DIR="/scratch.global/espin164/Evan_MEE/M4/filtered"



#Create output directory because it doesn’t exist
mkdir -p "$FILTER_DIR"



#get sample info from MEE_test_list.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

SAMPLE=$(echo "$LINE" | cut -f1)


#Input/output files

INPUT_BAM="${DEDUP_DIR}/${SAMPLE}.dedup.bam"

FILTERED_BAM="${FILTER_DIR}/${SAMPLE}.MQ20.bam"

FLAGSTAT_OUT="${FILTER_DIR}/${SAMPLE}.flagstat.txt"



echo "Filtering sample: $SAMPLE"

echo "Input: $INPUT_BAM"



#1: Filter BAM by mapping quality ≥ 20
#In samtools/1.21, -q is not base quality. It is specified as:  -q, --min-MQ INT           ...have mapping quality >= INT
#-q is the mapping quality threshold

samtools view -@ ${SLURM_CPUS_PER_TASK} -b -q 20 "$INPUT_BAM" -o "$FILTERED_BAM"



#2: Index filtered BAM
samtools index "$FILTERED_BAM"



#3: Generate alignment summary
samtools flagstat "$FILTERED_BAM" > "$FLAGSTAT_OUT"



echo "Finished sample: $SAMPLE"


