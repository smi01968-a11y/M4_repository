#!/bin/bash

#SBATCH --job-name=concat_vcf

#SBATCH --output=concat_vcf_%j.out

#SBATCH --error=concat_vcf_%j.err

#SBATCH --time=8:00:00

#SBATCH --mem=16G

#SBATCH --cpus-per-task=4

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu


#before running, create a file named chunks.list with the paths to each of the 20 chromosome-level
#vcfs produced by gather_chromosome vcfs step.


module load bcftools/1.21


bcftools concat \
-f chunks.list \
-Oz \
-o M4_joint.mapq20.vcf.gz

#M4_joint.mapq20.vcf.gz is the final combined vcf containing all the chromosomes (and samples)

bcftools index -t M4_joint.mapq20.vcf.gz
