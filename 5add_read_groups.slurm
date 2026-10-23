#!/bin/bash

#SBATCH --job-name=post_alignment

#SBATCH --output=post_align_%A_%a.out

#SBATCH --error=post_align_%A_%a.err

#SBATCH --time=8:00:00

#SBATCH --mem=32G

#SBATCH --cpus-per-task=4

#SBATCH --array=1-196     #adjust number of samples

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu


#to run this code:

#sbatch postalignment.slurm

#PART1: adding read groups to the bam files using samplelist.txt

module load samtools/1.21

module load gatk/4.4.0

which samtools
which gatk


#directory of the bam files
BAM_DIR="/scratch.global/espin164/Evan_MEE/M4/alignments"
#directory to create for outputs:
RG_ADDED="/scratch.global/espin164/Evan_MEE/M4/afteralignments"

mkdir -p $RG_ADDED

#the RGsamplelist.txt looks like this:
#sample name	flowcellID	flowcellID.lane

#Get line of MEE_test_list.txt that matches this array task

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

SAMPLE=$(echo "$LINE" | cut -f1)

FCID=$(echo "$LINE" | cut -f2)

FC_LANE=$(echo "$LINE" | cut -f3)



INPUT_BAM="${BAM_DIR}/${SAMPLE}.sorted.bam"

OUTPUT_BAM="${RG_ADDED}/${SAMPLE}.RG.bam"



gatk AddOrReplaceReadGroups \
-I "$INPUT_BAM" \
-O "$OUTPUT_BAM" \
-RGID "$FCID" \
-RGLB lib1 \
-RGPL ILLUMINA \
-RGPU "$FC_LANE" \
-RGSM "$SAMPLE"


samtools index "$OUTPUT_BAM"
#indexing the new bam files


