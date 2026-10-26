#!/bin/bash -l

#SBATCH --time=12:00:00

#SBATCH --ntasks=1

#SBATCH --mem=16g

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu

#SBATCH --job-name=GATHER_CHR

#SBATCH --output=GATHER_CHR_%A_%a.out

#SBATCH --error=GATHER_CHR_%A_%a.err

#SBATCH --array=1-20

#takes each vcf genomic chunk and compiles them by chromosome

module load gatk/4.4.0

#specify output directory (already should exist)
OUTDIR="/scratch.global/espin164/Evan_MEE/M4/filtered/gatk_joint_output_mapq20_parallel"

#specify directory for gathering chromosome vcfs
GATHERDIR="${OUTDIR}/chromosome_vcfs"

mkdir -p "$GATHERDIR"


CHR=$(printf "Gm%02d" "$SLURM_ARRAY_TASK_ID")

#output vcf name per chromosome
OUTVCF="${GATHERDIR}/${CHR}.mapq20.vcf.gz"


echo "Chromosome: $CHR"

echo "Output: $OUTVCF"


#get all chunks for this chromosome and sort by genomic start position

find "$OUTDIR" -maxdepth 1 \
-name "${CHR}_*.mapq20.vcf.gz" |
sort -t_ -k2,2n \

> "${GATHERDIR}/${CHR}_vcfs.list"


#count chunks

N=$(wc -l < "${GATHERDIR}/${CHR}_vcfs.list")

echo "Found $N chunks for $CHR"


if [ "$N" -eq 0 ]; then
echo "ERROR: No VCF chunks found"
exit 1
fi



#build GatherVcfs command
INPUT_ARGS=""


while read -r VCF; do
INPUT_ARGS+=" -I $VCF"
done < "${GATHERDIR}/${CHR}_vcfs.list"

#gather the VCFs
gatk GatherVcfs \
$INPUT_ARGS \
-O "$OUTVCF"


#index

gatk IndexFeatureFile \
-I "$OUTVCF"


echo "Finished $CHR"
