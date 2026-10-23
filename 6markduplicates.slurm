#!/bin/bash

#SBATCH --job-name=postrgs

#SBATCH --output=postrgs_%A_%a.out

#SBATCH --error=postrgs_%A_%a.err

#SBATCH --time=8:00:00

#SBATCH --mem=32G

#SBATCH --cpus-per-task=4

#SBATCH --array=1-196     #adjust number of samples

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu



#to run this code:

#sbatch postrgs.slurm

module load samtools/1.21

module load gatk/4.4.0

which samtools
which gatk

#directory of bam files with read groups added

RG_ADDED="/scratch.global/espin164/Evan_MEE/M4/afteralignments"

#output directory to be created
DEDUP_DIR="/scratch.global/espin164/Evan_MEE/M4/deduped"


#mkdir for dedup because it doesn't exist yet
mkdir -p "$DEDUP_DIR"


#Get line of MEE_test_list.txt that matches this array task

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

SAMPLE=$(echo "$LINE" | cut -f1)



#Define input and output bams
INPUT_BAM="${RG_ADDED}/${SAMPLE}.RG.bam"

SORTED_BAM="${DEDUP_DIR}/${SAMPLE}.sorted.bam"

DEDUP_BAM="${DEDUP_DIR}/${SAMPLE}.dedup.bam"

METRICS="${DEDUP_DIR}/${SAMPLE}.dedup.metrics.txt"



echo "Processing sample: $SAMPLE"

echo "Input BAM: $INPUT_BAM"


#the bam files with read groups added are in the rgadded directory. Now we need to mark duplicates, sort, index


#sort bam by coordinate
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "$SORTED_BAM" "$INPUT_BAM"

#then index
samtools index "$SORTED_BAM"

#mark duplicates with gatk MarkDuplicates
gatk MarkDuplicates \
-I "$SORTED_BAM" \
-O "$DEDUP_BAM" \
-M "$METRICS"

#index deduplicated bam

samtools index "$DEDUP_BAM"

echo "finished sample $SAMPLE"
