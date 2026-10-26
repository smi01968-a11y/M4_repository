#!/bin/bash -l

#SBATCH --time=2:00:00

#SBATCH --ntasks=24

#SBATCH --tmp=100g

#SBATCH --mem=100g

#SBATCH --mail-type=ALL

#SBATCH --mail-user=smi01968@umn.edu

#SBATCH --job-name=VCF_hardfiltering

#SBATCH --output=Part5_%j.out

#SBATCH --error=Part5_%j.err



#----------------------------------

### GATK - Hard Filtering for snps and InDels

#-----------------------------------


#Loading required modules and file paths

module load python3

module load picard/2.25.6

module load gatk/4.4.0



#path for the Wm82.a6 reference genome

REF="$MSIPROJECT/shared/References/Gmax.Wm82.a6/assembly/Gmax_880_v6.0.fa"

#path to where outputs will end up:
Filtered_VCF="/scratch.global/espin164/Evan_MEE/M4/filtered/gatk_joint_output_mapq20_parallel/Filtered_VCF"

mkdir -p $Filtered_VCF


##########################################################

#Subset SNP and INDEL variants

###########################################################################

#the vcf.gz has got to be in the same directory in which this script is being run
#select SNPs
gatk SelectVariants \
-R $REF \
-V M4_joint.mapq20.vcf.gz \
--select-type-to-include SNP \
-O M4_cohort_SNP.vcf.gz


echo "SNP variants selected"



#Select Indels

gatk SelectVariants \
-R $REF \
-V M4_joint.mapq20.vcf.gz \
--select-type-to-include INDEL \
-O M4_cohort_INDEL.vcf.gz



echo "Indel variants selected"


######################################################################################

#Hard filtering SNPs

##################################################################

#the decision to maintain the defaults was made by checking the VCF features in R

gatk VariantFiltration \
-V M4_cohort_SNP.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL40" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O $Filtered_VCF/M4_filtered_SNP.vcf.gz



if [ $? -ne 0 ]; then
echo >&2 "ERROR: SNP filtering failed"
exit 1
fi
#checks to make sure variant filtration works - if it didn't then ends script

echo "DONE: SNP filtering"


#####################################################################

#Select passed SNP variants

#################################################################################

gatk SelectVariants \
-R $REF \
-V $Filtered_VCF/M4_filtered_SNP.vcf.gz \
--exclude-filtered \
-O $Filtered_VCF/M4_filtered_pass_SNP.vcf.gz

if [ $? -ne 0 ]; then
echo >&2 "ERROR: Selecting SNPs failed"
exit 1
fi

echo "DONE: Selecting SNPs"


########################################################

#Hard filtering INDELs

############################################################


#default parameters implemented other than QUAL and MQ, which are 40 rather than 30
#mechanism for production of indels by EMS or ENU is not well understood, therefore stricter parameters may be justified here
gatk VariantFiltration \
-V M4_cohort_INDEL.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 40.0" --filter-name "QUAL40" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O $Filtered_VCF/M4_filtered_INDEL.vcf.gz


if [ $? -ne 0 ]; then
echo >&2 "ERROR: INDEL filtering failed"
exit 1
fi

echo "DONE: INDEL filtering"



####################################

#Select passed INDEL variants

#########################################


gatk SelectVariants \
-R $REF \
-V $Filtered_VCF/M4_filtered_INDEL.vcf.gz \
--exclude-filtered \
-O $Filtered_VCF/M4_filtered_pass_INDEL.vcf.gz


if [ $? -ne 0 ]; then
echo >&2 "ERROR: Selecting INDELs failed"
exit 1
fi

echo "DONE: Selecting INDELs"



#######################################################################

#Merge filtered SNPs and INDELs
#may be useful

###############################################################

cd $Filtered_VCF/

java -jar $PICARD_DIR/picard.jar MergeVcfs \
I=M4_filtered_pass_SNP.vcf.gz \
I=M4_filtered_pass_INDEL.vcf.gz \
O=M4_filtered_pass_all_variants.vcf.gz


if [ $? -ne 0 ]; then
echo >&2 "ERROR: Merging variants failed"
exit 1
fi

echo "DONE: Merging variants"
