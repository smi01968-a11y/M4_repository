#!/bin/bash

#SBATCH --job-name=trimmomatic

#SBATCH --time=24:00:00

#SBATCH --mem=32G

#SBATCH --cpus-per-task=4

#SBATCH --output=trim_%A_%a.out

#SBATCH --error=trim_%A_%a.err

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu


#after running qc and analyzing results with multiqc
#use this script to trim the adapter sequences, using the alladaptors.fa file
#for common adapters
#output will be paired.fastq and unpaired.fastq files that end up in the current directory

#to run this script

#ie specify number of samples in subset, along with the file containing sample info from before
#sbatch --array=1-196 --export=LIST=samples.txt trim.slurm


# load trimmomatic (version 0.39)
module load java/openjdk-8_202
module load trimmomatic/0.39
#defining the loaction of trimmomatic
TRIMJAR=/common/software/install/migrated/trimmomatic/0.39/trimmomatic-0.39.jar

FASTQDIR=/projects/standard/stuparr/data_delivery/umgc/2026-q3/20260722_AV241910_Aviti3_Run_417_07_22_26/Stupar_Project_076-merged



#make array variable(again)
SAMPLE="$(/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${LIST} | cut -f1)"

#use trimmomatic on paired end reads
#PE=pair end reads (instead of single end (SE)

java -jar $TRIMJAR PE -threads 4 \
${FASTQDIR}/${SAMPLE}_R1.fastq.gz ${FASTQDIR}/${SAMPLE}_R2.fastq.gz \
${SAMPLE}_R1_paired.fastq.gz ${SAMPLE}_R1_unpaired.fastq.gz \
${SAMPLE}_R2_paired.fastq.gz ${SAMPLE}_R2_unpaired.fastq.gz \
ILLUMINACLIP:all_adaptors.fa:2:30:10 \
LEADING:10 TRAILING:10 \
MINLEN:50

#leading trailingremoves very poor quality bases at the beginning and end of reads
#illumina clip to remove actual adapter sequence
#minlen:50 = removes any reads less than 50 bases (expecting roughly 150 bp reads)
