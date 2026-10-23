#!/bin/bash

#SBATCH --job-name=bwa_align

#SBATCH --output=bwa_align_%A_%a.out

#SBATCH --error=bwa_align_%A_%a.err

#SBATCH --time=48:00:00

#SBATCH --mem=32G

#SBATCH --tmp=32G

#SBATCH --cpus-per-task=4

#SBATCH --array=1-196%20     #adjust number of samples

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu



#run 20 samples at a time to use lesscomp resources at a time
#to run this code:
#sbatch bwamem2.slurm
#48 G and 24 hours time for the bsa samples, 8 cpus per task

module load bwa/0.7.17

module load samtools/1.21

#verify correct versions are used
which bwa
which samtools



LIST=./samples.txt

# Get sample name from list

SAMPLE="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST} | cut -f1)"

echo "Task ${SLURM_ARRAY_TASK_ID} running on sample: $SAMPLE"


# Paths
#direction to reference genome .fa file (Wm82)
REF=/projects/standard/stuparr/shared/References/Gmax.Wm82.a6/assembly/Gmax_880_v6.0.fa

#directory of the paired reads post trimming
READS_DIR=/scratch.global/espin164/Evan_MEE/M4/

#output directory
OUT_DIR=/scratch.global/espin164/Evan_MEE/M4/alignments



mkdir -p $OUT_DIR

echo "Task ID: $SLURM_ARRAY_TASK_ID"

echo "Sample: $SAMPLE"

echo "Using reference: $REF"

echo "Reads R1: ${READS_DIR}/${SAMPLE}_R1_paired.fastq.gz"

echo "Reads R2: ${READS_DIR}/${SAMPLE}_R2_paired.fastq.gz"

echo "Output directory: $OUT_DIR"




echo "Task ${SLURM_ARRAY_TASK_ID} running on sample: $SAMPLE"


#align with BWA and output sorted BAM

bwa mem -t ${SLURM_CPUS_PER_TASK} \
-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
$REF \
${READS_DIR}/${SAMPLE}_R1_paired.fastq.gz \
${READS_DIR}/${SAMPLE}_R2_paired.fastq.gz \
| samtools view -@ ${SLURM_CPUS_PER_TASK} -bS - \
| samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${OUT_DIR}/${SAMPLE}.sorted.bam



#index the BAM

samtools index ${OUT_DIR}/${SAMPLE}.sorted.bam

echo "Finished sample $SAMPLE"

