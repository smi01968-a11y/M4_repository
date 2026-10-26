#!/bin/bash -l

#SBATCH --time=1:00:00

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem=32g

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu

#SBATCH --job-name=regions_check

#SBATCH --output=regions_%j.out

#SBATCH --error=regions_%j.err

#check regions
module load gatk/4.4.0

REF="/projects/standard/stuparr/shared/References/Gmax.Wm82.a6/assembly/Gmax_880_v6.0.fa"

BAM_LIST="/scratch.global/espin164/Evan_MEE/M4/filtered/bamlist.txt"

INTERVALS="/scratch.global/espin164/Evan_MEE/M4/filtered/soybean_250kb_intervals.list"

OUTDIR="/scratch.global/espin164/Evan_MEE/M4/filtered/gatk_joint_output_mapq20_parallel"

mkdir -p "$OUTDIR"

# Get the interval corresponding to this SLURM task

INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INTERVALS")

# Convert Gm01:1-2000000 into Gm01_1_2000000

# Convert Gm01:1-2000000 into Gm01_1_2000000

#REGION=$(echo "$INTERVAL" | tr ':-' '_')


while read -r INTERVAL; do
REGION=$(echo "$INTERVAL" | tr ':-' '_')
echo "${REGION}.mapq20.vcf.gz"
done < "$INTERVALS" | sort > expected_vcfs.txt

ls "$OUTDIR"/*.vcf.gz | sed 's#.*/##' | sort > actual_vcfs.txt
