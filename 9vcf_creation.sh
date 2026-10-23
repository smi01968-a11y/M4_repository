#!/bin/bash -l

#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=smi01968@umn.edu
#SBATCH --job-name=VCF_JOINT
#SBATCH --output=VCFJOINT_%j.out
#SBATCH --error=VCFJOINT_%j.err
#SBATCH --array=1-4055

#haplotype vcf mode calling. All samples need to go in as input simultaneously.
#all samples called together in a single run.

#operates on intervals of 250kb from the genome at a time.
#generate chunks from the reference genome index file
#awk '{
#    chrom=$1
#    chrlen=$2
#
#    for (start=1; start<=chrlen; start+=250000) {
#        end=start+249999
#        if (end>chrlen) end=chrlen
#
#        print chrom ":" start "-" end
#    }
#}' ${REF}.fai > soybean_250kb_intervals.list


module load gatk/4.4.0

REF="/projects/standard/stuparr/shared/References/Gmax.Wm82.a6/assembly/Gmax_880_v6.0.fa"

#list of the names of the filtered bam files
BAM_LIST="/scratch.global/espin164/Evan_MEE/M4/filtered/bamlist.txt"

#list of the 250kb intervals previously created, from reference genome
INTERVALS="/scratch.global/espin164/Evan_MEE/M4/filtered/soybean_250kb_intervals.list"

#output directory name
OUTDIR="/scratch.global/espin164/Evan_MEE/M4/filtered/gatk_joint_output_mapq20_parallel"

mkdir -p "$OUTDIR"

#get the interval corresponding to this SLURM task
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INTERVALS")

#convert Gm01:1-2000000 into Gm01_1_2000000
REGION=$(echo "$INTERVAL" | tr ':-' '_')

OUTVCF="${OUTDIR}/${REGION}.mapq20.vcf.gz"

#build the -I arguments for all 196 BAMs
INPUT_ARGS=""

while read -r BAM; do
INPUT_ARGS+=" -I $BAM"
done < "$BAM_LIST"

echo "Processing interval: $INTERVAL"
echo "Output: $OUTVCF"

#gatk haplotype caller vcf mode
gatk --java-options "-Xmx24g -Djava.io.tmpdir=/scratch.global/$SLURM_JOB_ID" \
HaplotypeCaller \
-R "$REF" \
$INPUT_ARGS \
-L "$INTERVAL" \
--native-pair-hmm-threads 8 \
-O "$OUTVCF"

gatk IndexFeatureFile -I "$OUTVCF"

