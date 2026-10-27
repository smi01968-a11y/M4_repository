#!/bin/bash
#SBATCH --job-name=vep
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=smi01968@umn.edu

#using Variant Effect Predictor to generate annotations for indel and snp variants

#export BCFTOOLS_PLUGINS=/projects/standard/morrellp/public/Software/bcftools-1.21/plugins/
module load htslib/1.16-gcc-8.2.0-x236tnl

WORKING_DIR="/users/1/smi01968/M4/"

mkdir -p ${WORKING_DIR}/{vep,tmp,cache}

################################################################################
# Get VeP singularity image
module load apptainer

export SINGULARITY_TMPDIR="${WORKING_DIR}/tmp"
export SINGULARITY_CACHEDIR="${WORKING_DIR}/cache"

if [ ! -f ${WORKING_DIR}/ensembl-vep-r_111.0.sif ]; then
echo "Building VeP singularity image..."
singularity build -F \
${WORKING_DIR}/ensembl-vep-r_111.0.sif \
docker://ensemblorg/ensembl-vep:release_111.0
fi

######################
# Add path to .gff file
GFF="/projects/standard/stuparr/shared/References/Gmax.Wm82.a6/annotation/Gmax_880_Wm82.a6.v1.gene_exons.gff3.gz"
GFF_SORTED="${WORKING_DIR}/Gmax_Wm82.a6.gene_exons.sorted.gff3"
GFF_SORTED_ZIPPED="${WORKING_DIR}/Gmax_Wm82.a6.gene_exons.sorted.gff3.gz"
(zgrep "^#" ${GFF}; zgrep -v "^#" ${GFF} | sort -k1,1 -k4,4n) > ${GFF_SORTED}
bgzip -c ${GFF_SORTED} > ${GFF_SORTED_ZIPPED}
tabix -f -p gff -C ${GFF_SORTED_ZIPPED}

# Add path to reference fasta
FASTA="/projects/standard/stuparr/shared/References/Gmax.Wm82.a6/assembly/Gmax_880_v6.0.fa"

FASTA_DIR="/projects/standard/stuparr/shared/References/Gmax.Wm82.a6/assembly"
GFF_DIR="/projects/standard/stuparr/shared/References/Gmax.Wm82.a6/annotation"

SPECIES="glycine_max"

################################################################################

# Variant sets should be either 'all', 'deletions', 'insertions', 'snps', 'common', or 'rare'
#specify snps for SNP VEP, and all for indels
#or use all for VEP of all variants
VARIANT_SET="all"

# Add path to .vcf file
#VCF="${WORKING_DIR}/temp/singletons/exclusive/vcf/WTP0579-40-SNPs.recode.vcf"
#bgzip -c ${VCF} > ${VCF}.gz
#VCF_ZIP="${VCF}.gz"

for i in /users/1/smi01968/M4/Filtered_VCF/M4_family_singletonINDELS.vcf.gz ; do
fname=$(basename ${i} .vcf.gz)
echo $fname
singularity run \
--bind ${WORKING_DIR}/:${WORKING_DIR}/ \
--bind ${FASTA_DIR}:${FASTA_DIR} \
--bind ${GFF_DIR}:${GFF_DIR} \
${WORKING_DIR}/ensembl-vep-r_111.0.sif \
vep \
-i "${i}" \
--gff "${GFF_SORTED_ZIPPED}" \
--fasta "${FASTA}" \
--total_length \
--check_svs \
--verbose \
--format vcf \
--force \
--fork 20 \
--everything \
--species glycine_max \
--warning_file "${WORKING_DIR}/vep/${fname}-VeP-warning.txt" \
-o "${WORKING_DIR}/vep/${fname}-VeP.txt"
done
