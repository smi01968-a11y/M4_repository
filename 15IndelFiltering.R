#Script for filtering INDELS

#before running this script, normalize the vcf and restrict to biallelic indels

# Reference genome
# REF="/projects/standard/stuparr/shared/References/Gmax.Wm82.a6/assembly/Gmax_880_v6.0.fa"
# 
# INPUT="M4_filtered_INDEL.vcf.gz"
# NORMALIZED="M4_indels_normalized.vcf.gz"
# OUTPUT="M4_indels_biallelic.vcf.gz"
# 
# # 1. Normalize indels against the reference genome
# bcftools norm \
# -f "$REF" \
# -Oz \
# -o "$NORMALIZED" \
# "$INPUT"
# 
# bcftools index -t "$NORMALIZED"
# 
# # 2. Retain biallelic variants
# bcftools view \
# -m 2 \
# -M 2 \
# -Oz \
# -o "$OUTPUT" \
# "$NORMALIZED"
# 
# bcftools index -t "$OUTPUT"


#already hard filtered, with standard parameters other than
#setting QUAL<40 and MQ<40
#also, alreaady normalized/left aligned the indels relative to reference genome
#and removed multiallelic Indels


setwd("~/M4/Filtered_VCF")

#packages
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vcfR)
#vcfR version 1.16.0
sessionInfo()


OG<- read.vcfR("M4_indels_biallelic.vcf.gz")
#M4_indels_biallelic.vcf.gz


#inspect with 
# head(OG@fix)
# head(OG@gt)



#double check that all multiallelic InDel sites were actually removed
multiallelic_count<-grepl(",", OG@fix[, "ALT"])
#search for commas in the alt field, indicating more than one ALT allele
sum(multiallelic_count)
#if there are no multiallelic sites, then output sum should be zero
#sum=0 


#inspect
colnames(OG@gt)

all_samples<- colnames(OG@gt)[-1]
#drop FORMAT
length(all_samples)
#196
print(all_samples)




#subsetting to just the non-outcrossed 171 samples

#create lists of sample names to be kept together (171 + 4 parents)
#suspected outcrossed samples removed
m4_samples<- c(
  "MEE00001-C-1-1-MEE00001-C-1-1-GCCATAGCCATA",
  "MEE00001-D-1-1-MEE00001-D-1-1-GCCATAGCCATA",
  "MEE00002-C-1-1-MEE00002-C-1-1-GCCATAGCCATA",
  "MEE00002-D-1-1-MEE00002-D-1-1-GCCATAGCCATA",
  "MEE00003-C-1-1-MEE00003-C-1-1-GCCATAGCCATA",
  "MEE00003-E-1-1-MEE00003-E-1-1-GCCATAGCCATA",
  "MEE00004-C-1-1-MEE00004-C-1-1-GCCATAGCCATA",
  "MEE00004-D-1-1-MEE00004-D-1-1-GCCATAGCCATA",
  "MEE00005-E-1-1-MEE00005-E-1-1-TCTGGTTCTGGT",
  "MEE00005-F-1-1-MEE00005-F-1-1-TCTGGTTCTGGT",
  "MEE00006-C-1-1-MEE00006-C-1-1-TCTGGTTCTGGT",
  "MEE00006-D-1-1-MEE00006-D-1-1-TCTGGTTCTGGT",
  "MEE00006-E-1-1-MEE00006-E-1-1-TCTGGTTCTGGT",
  "MEE00006-F-1-1-MEE00006-F-1-1-TGGCTTTGGCTT",
  "MEE00007-C-1-1-MEE00007-C-1-1-TGGCTTTGGCTT",
  "MEE00007-D-1-1-MEE00007-D-1-1-TGGCTTTGGCTT",
  "MEE00007-E-1-1-MEE00007-E-1-1-TGGCTTTGGCTT",
  "MEE00010-E-1-1-MEE00010-E-1-1-TGGCTTTGGCTT",
  "MEE00011-C-1-1-MEE00011-C-1-1-TGGCTTTGGCTT",
  "MEE00012-A-1-1-MEE00012-A-1-1-TGGCTTTGGCTT",
  "MEE00014-A-1-1-MEE00014-A-1-1-AACACTAACACT",
  "MEE00015-A-1-1-MEE00015-A-1-1-AACACTAACACT",
  "MEE00016-A-1-1-MEE00016-A-1-1-AACACTAACACT",
  "MEE00018-A-1-1-MEE00018-A-1-1-AACACTAACACT",
  "MEE00018-B-1-1-MEE00018-B-1-1-AACACTAACACT",
  "MEE00018-C-1-1-MEE00018-C-1-1-AACACTAACACT",
  "MEE00018-D-1-1-MEE00018-D-1-1-AACACTAACACT",
  "MEE00018-E-1-1-MEE00018-E-1-1-AACACTAACACT",
  "MEE00019-A-1-1-MEE00019-A-1-1-TCAGGATCAGGA",
  "MEE00021-B-1-1-MEE00021-B-1-1-TCAGGATCAGGA",
  "MEE00022-C-1-1-MEE00022-C-1-1-TCAGGATCAGGA",
  "MEE00023-A-1-1-MEE00023-A-1-1-TCAGGATCAGGA",
  "MEE00023-C-1-1-MEE00023-C-1-1-TCAGGATCAGGA",
  "MEE00023-D-1-1-MEE00023-D-1-1-TCAGGATCAGGA",
  "MEE00024-B-1-1-MEE00024-B-1-1-TCAGGATCAGGA",
  "MEE00024-C-1-1-MEE00024-C-1-1-TCAGGATCAGGA",
  "MEE00024-D-1-1-MEE00024-D-1-1-TGGTCCTGGTCC",
  "MEE00024-E-1-1-MEE00024-E-1-1-TGGTCCTGGTCC",
  "MEE00025-B-1-1-MEE00025-B-1-1-TGGTCCTGGTCC",
  "MEE00025-C-1-1-MEE00025-C-1-1-TGGTCCTGGTCC",
  "MEE00025-E-1-1-MEE00025-E-1-1-TGGTCCTGGTCC",
  "MEE00026-B-1-1-MEE00026-B-1-1-TGGTCCTGGTCC",
  "MEE00028-A-1-1-MEE00028-A-1-1-TGGTCCTGGTCC",
  "MEE00029-A-1-1-MEE00029-A-1-1-CAACTGCAACTG",
  "MEE00029-B-1-1-MEE00029-B-1-1-CAACTGCAACTG",
  "MEE00030-B-1-1-MEE00030-B-1-1-CAACTGCAACTG",
  "MEE00030-C-1-1-MEE00030-C-1-1-CAACTGCAACTG",
  "MEE00030-E-1-1-MEE00030-E-1-1-CAACTGCAACTG",
  "MEE00032-A-1-1-MEE00032-A-1-1-CAACTGCAACTG",
  "MEE00032-B-1-1-MEE00032-B-1-1-CAACTGCAACTG",
  "MEE00032-D-1-1-MEE00032-D-1-1-CGGACCCGGACC",
  "MEE00033-A-1-1-MEE00033-A-1-1-CGGACCCGGACC",
  "MEE00033-E-1-1-MEE00033-E-1-1-CGGACCCGGACC",
  "MEE00034-A-1-1-MEE00034-A-1-1-CGGACCCGGACC",
  "MEE00034-B-1-1-MEE00034-B-1-1-CGGACCCGGACC",
  "MEE00034-C-1-1-MEE00034-C-1-1-CGGACCCGGACC",
  "MEE00034-D-1-1-MEE00034-D-1-1-CGGACCCGGACC",
  "MEE00035-B-1-1-MEE00035-B-1-1-CGGACCCGGACC",
  "MEE00037-D-1-1-MEE00037-D-1-1-ATCGAGATCGAG",
  "MEE00038-C-1-1-MEE00038-C-1-1-ATCGAGATCGAG",
  "MEE00038-D-1-1-MEE00038-D-1-1-ATCGAGATCGAG",
  "MEE00039-A-1-1-MEE00039-A-1-1-ATCGAGATCGAG",
  "MEE00039-B-1-1-MEE00039-B-1-1-ATCGAGATCGAG",
  "MEE00039-C-1-1-MEE00039-C-1-1-ATCGAGATCGAG",
  "MEE00039-D-1-1-MEE00039-D-1-1-ATCGAGATCGAG",
  "MEE00039-E-1-1-MEE00039-E-1-1-ATGGTGATGGTG",
  "MEE00040-A-1-1-MEE00040-A-1-1-ATGGTGATGGTG",
  "MEE00040-C-1-1-MEE00040-C-1-1-ATGGTGATGGTG",
  "MEE00041-A-1-1-MEE00041-A-1-1-ATGGTGATGGTG",
  "MEE00041-C-1-1-MEE00041-C-1-1-ATGGTGATGGTG",
  "MEE00042-B-1-1-MEE00042-B-1-1-ATGGTGATGGTG",
  "MEE00042-E-1-1-MEE00042-E-1-1-ATGGTGATGGTG",
  "MEE00043-A-1-1-MEE00043-A-1-1-CTTAAGCTTAAG",
  "MEE00043-B-1-1-MEE00043-B-1-1-CTTAAGCTTAAG",
  "MEE00043-C-1-1-MEE00043-C-1-1-CTTAAGCTTAAG",
  "MEE00044-A-1-1-MEE00044-A-1-1-CTTAAGCTTAAG",
  "MEE00044-C-1-1-MEE00044-C-1-1-CTTAAGCTTAAG",
  "MEE00044-D-1-1-MEE00044-D-1-1-CTTAAGCTTAAG",
  "MEE00044-E-1-1-MEE00044-E-1-1-CTTAAGCTTAAG",
  "MEE00046-A-1-1-MEE00046-A-1-1-GAGTGCGAGTGC",
  "MEE00048-A-1-1-MEE00048-A-1-1-GAGTGCGAGTGC",
  "MEE00048-C-1-1-MEE00048-C-1-1-GAGTGCGAGTGC",
  "MEE00048-D-1-1-MEE00048-D-1-1-GAGTGCGAGTGC",
  "MEE00049-B-1-1-MEE00049-B-1-1-GAGTGCGAGTGC",
  "MEE00049-C-1-1-MEE00049-C-1-1-GAGTGCGAGTGC",
  "MEE00050-B-1-1-MEE00050-B-1-1-GAGTGCGAGTGC",
  "MEE00052-A-1-1-MEE00052-A-1-1-GCCATAGCCATA",
  "MEE00052-D-1-1-MEE00052-D-1-1-GCCATAGCCATA",
  "MEE00053-A-1-1-MEE00053-A-1-1-GCCATAGCCATA",
  "MEE00053-D-1-1-MEE00053-D-1-1-GCCATAGCCATA",
  "MEE00054-B-1-1-MEE00054-B-1-1-GCCATAGCCATA",
  "MEE00054-C-1-1-MEE00054-C-1-1-TCTGGTTCTGGT",
  "MEE00054-D-1-1-MEE00054-D-1-1-TCTGGTTCTGGT",
  "MEE00054-E-1-1-MEE00054-E-1-1-TCTGGTTCTGGT",
  "MEE00055-A-1-1-MEE00055-A-1-1-TCTGGTTCTGGT",
  "MEE00055-B-1-1-MEE00055-B-1-1-TCTGGTTCTGGT",
  "MEE00055-C-1-1-MEE00055-C-1-1-TCTGGTTCTGGT",
  "MEE00056-A-1-1-MEE00056-A-1-1-TCTGGTTCTGGT",
  "MEE00056-B-1-1-MEE00056-B-1-1-TCTGGTTCTGGT",
  "MEE00056-C-1-1-MEE00056-C-1-1-TGGCTTTGGCTT",
  "MEE00056-D-1-1-MEE00056-D-1-1-TGGCTTTGGCTT",
  "MEE00057-A-1-1-MEE00057-A-1-1-TGGCTTTGGCTT",
  "MEE00057-B-1-1-MEE00057-B-1-1-TGGCTTTGGCTT",
  "MEE00057-C-1-1-MEE00057-C-1-1-TGGCTTTGGCTT",
  "MEE00057-D-1-1-MEE00057-D-1-1-TGGCTTTGGCTT",
  "MEE00057-E-1-1-MEE00057-E-1-1-TGGCTTTGGCTT",
  "MEE00058-C-1-1-MEE00058-C-1-1-AACACTAACACT",
  "MEE00059-A-1-1-MEE00059-A-1-1-AACACTAACACT",
  "MEE00059-B-1-1-MEE00059-B-1-1-AACACTAACACT",
  "MEE00059-C-1-1-MEE00059-C-1-1-AACACTAACACT",
  "MEE00059-D-1-1-MEE00059-D-1-1-AACACTAACACT",
  "MEE00061-A-1-1-MEE00061-A-1-1-AACACTAACACT",
  "MEE00061-B-1-1-MEE00061-B-1-1-TCAGGATCAGGA",
  "MEE00061-C-1-1-MEE00061-C-1-1-TCAGGATCAGGA",
  "MEE00062-A-1-1-MEE00062-A-1-1-TCAGGATCAGGA",
  "MEE00062-B-1-1-MEE00062-B-1-1-TCAGGATCAGGA",
  "MEE00063-A-1-1-MEE00063-A-1-1-TCAGGATCAGGA",
  "MEE00063-B-1-1-MEE00063-B-1-1-TCAGGATCAGGA",
  "MEE00063-C-1-1-MEE00063-C-1-1-TCAGGATCAGGA",
  "MEE00063-D-1-1-MEE00063-D-1-1-TCAGGATCAGGA",
  "MEE00064-B-1-1-MEE00064-B-1-1-TGGTCCTGGTCC",
  "MEE00064-C-1-1-MEE00064-C-1-1-TGGTCCTGGTCC",
  "MEE00065-A-1-1-MEE00065-A-1-1-TGGTCCTGGTCC",
  "MEE00065-B-1-1-MEE00065-B-1-1-TGGTCCTGGTCC",
  "MEE00066-B-1-1-MEE00066-B-1-1-TGGTCCTGGTCC",
  "MEE00067-A-1-1-MEE00067-A-1-1-TGGTCCTGGTCC",
  "MEE00067-B-1-1-MEE00067-B-1-1-TGGTCCTGGTCC",
  "MEE00067-C-1-1-MEE00067-C-1-1-TGGTCCTGGTCC",
  "MEE00067-D-1-1-MEE00067-D-1-1-CAACTGCAACTG",
  "MEE00067-E-1-1-MEE00067-E-1-1-CAACTGCAACTG",
  "MEE00069-B-1-1-MEE00069-B-1-1-CAACTGCAACTG",
  "MEE00069-C-1-1-MEE00069-C-1-1-CAACTGCAACTG",
  "MEE00069-D-1-1-MEE00069-D-1-1-CAACTGCAACTG",
  "MEE00070-A-1-1-MEE00070-A-1-1-CAACTGCAACTG",
  "MEE00070-C-1-1-MEE00070-C-1-1-CAACTGCAACTG",
  "MEE00070-D-1-1-MEE00070-D-1-1-CAACTGCAACTG",
  "MEE00071-A-1-1-MEE00071-A-1-1-CGGACCCGGACC",
  "MEE00071-C-1-1-MEE00071-C-1-1-CGGACCCGGACC",
  "MEE00072-B-1-1-MEE00072-B-1-1-CGGACCCGGACC",
  "MEE00073-B-1-1-MEE00073-B-1-1-CGGACCCGGACC",
  "MEE00073-D-1-1-MEE00073-D-1-1-CGGACCCGGACC",
  "MEE00074-D-1-1-MEE00074-D-1-1-CGGACCCGGACC",
  "MEE00075-A-1-1-MEE00075-A-1-1-CGGACCCGGACC",
  "MEE00075-B-1-1-MEE00075-B-1-1-ATCGAGATCGAG",
  "MEE00075-C-1-1-MEE00075-C-1-1-ATCGAGATCGAG",
  "MEE00076-A-1-1-MEE00076-A-1-1-ATCGAGATCGAG",
  "MEE00076-B-1-1-MEE00076-B-1-1-ATCGAGATCGAG",
  "MEE00076-C-1-1-MEE00076-C-1-1-ATCGAGATCGAG",
  "MEE00076-D-1-1-MEE00076-D-1-1-ATCGAGATCGAG",
  "MEE00077-C-1-1-MEE00077-C-1-1-ATCGAGATCGAG",
  "MEE00078-C-1-1-MEE00078-C-1-1-ATGGTGATGGTG",
  "MEE00078-D-1-1-MEE00078-D-1-1-ATGGTGATGGTG",
  "MEE00078-E-1-1-MEE00078-E-1-1-ATGGTGATGGTG",
  "MEE00079-B-1-1-MEE00079-B-1-1-ATGGTGATGGTG",
  "MEE00079-C-1-1-MEE00079-C-1-1-ATGGTGATGGTG",
  "MEE00079-D-1-1-MEE00079-D-1-1-ATGGTGATGGTG",
  "MEE00079-E-1-1-MEE00079-E-1-1-ATGGTGATGGTG",
  "MEE00080-A-1-1-MEE00080-A-1-1-CTTAAGCTTAAG",
  "MEE00080-B-1-1-MEE00080-B-1-1-CTTAAGCTTAAG",
  "MEE00080-C-1-1-MEE00080-C-1-1-CTTAAGCTTAAG",
  "MEE00082-B-1-1-MEE00082-B-1-1-CTTAAGCTTAAG",
  "MEE00087-A-1-1-MEE00087-A-1-1-CTTAAGCTTAAG",
  "MEE00087-B-1-1-MEE00087-B-1-1-CTTAAGCTTAAG",
  "MEE00087-C-1-1-MEE00087-C-1-1-GAGTGCGAGTGC",
  "MEE00088-A-1-1-MEE00088-A-1-1-GAGTGCGAGTGC",
  "MEE00088-B-1-1-MEE00088-B-1-1-GAGTGCGAGTGC",
  "MEE00088-D-1-1-MEE00088-D-1-1-GAGTGCGAGTGC",
  "MEE00088-E-1-1-MEE00088-E-1-1-GAGTGCGAGTGC",
  "MEE00089-A-1-1-MEE00089-A-1-1-GAGTGCGAGTGC",
  "MEE00089-B-1-1-MEE00089-B-1-1-GAGTGCGAGTGC",
  "MEE00454-D-44-2-MEE00454-D-44-2-GCCATAGCCATA",
  "MN1806CN-1_founder-prep04_Rep1-MN1806CN-1_founder-prep04_Rep1-CTTAAGCTTAAG",
  "MN1806CN-1_founder-prep04_Rep2-MN1806CN-1_founder-prep04_Rep2-TGGCTTTGGCTT",
  "MN1806CN-1_founder-prep11_Rep1-MN1806CN-1_founder-prep11_Rep1-GAGTGCGAGTGC",
  "MN1806CN-1_founder-prep11_Rep2-MN1806CN-1_founder-prep11_Rep2-ATGGTGATGGTG"
)

#including "MEE00454-D-44-2-MEE00454-D-44-2-GCCATAGCCATA" since it is one of the mutants


new454BSA <-c("MEE00454-D-44-2-MEE00454-D-44-2-GCCATAGCCATA",
              "MEE00454-D-WT_BSA_Rep1-MEE00454-D-WT_BSA_Rep1-TCTGGTTCTGGT",
              "MEE00454-D-WT_BSA_Rep2-MEE00454-D-WT_BSA_Rep2-TGGCTTTGGCTT",
              "MEE00454-D-t_BSA_Rep1-MEE00454-D-t_BSA_Rep1-TGGCTTTGGCTT",
              "MEE00454-D-t_BSA_Rep2-MEE00454-D-t_BSA_Rep2-AACACTAACACT",
              "MN1806CN-1_founder-prep04_Rep1-MN1806CN-1_founder-prep04_Rep1-CTTAAGCTTAAG",
              "MN1806CN-1_founder-prep04_Rep2-MN1806CN-1_founder-prep04_Rep2-TGGCTTTGGCTT",
              "MN1806CN-1_founder-prep11_Rep1-MN1806CN-1_founder-prep11_Rep1-GAGTGCGAGTGC",
              "MN1806CN-1_founder-prep11_Rep2-MN1806CN-1_founder-prep11_Rep2-ATGGTGATGGTG")

parentals <- c("MN1806CN-1_founder-prep04_Rep1-MN1806CN-1_founder-prep04_Rep1-CTTAAGCTTAAG",
               "MN1806CN-1_founder-prep04_Rep2-MN1806CN-1_founder-prep04_Rep2-TGGCTTTGGCTT",
               "MN1806CN-1_founder-prep11_Rep1-MN1806CN-1_founder-prep11_Rep1-GAGTGCGAGTGC",
               "MN1806CN-1_founder-prep11_Rep2-MN1806CN-1_founder-prep11_Rep2-ATGGTGATGGTG")

length(m4_samples)
#should return 175


###################################
#################################
################################
#filtering for:
#depth
#quality

m4_samples[!m4_samples %in% colnames(OG@gt)]
#verify the sample names for m4 are indeed the same as in jointvcf

#same depth and GQ that SNPs were filtered on
subset_and_classify <- function(vcf, samples, min_dp = 10, max_dp = 100, min_gq = 20) {
  
  sub <- vcf
  
  #keep FORMAT and specified samples
  sub@gt <- sub@gt[, c("FORMAT", samples)]
  
  # Extract AD and GQ
  ad <- extract.gt(sub, element = "AD")
  gq <- extract.gt(sub, element = "GQ")
  gq <- apply(gq, 2, as.numeric)
  
  #extract REF and ALT read counts
  parse_ref <- function(x) as.numeric(sapply(strsplit(x, ","), `[`, 1))
  parse_alt <- function(x) as.numeric(sapply(strsplit(x, ","), `[`, 2))
  
  ref_mat <- apply(ad, 2, parse_ref)
  alt_mat <- apply(ad, 2, parse_alt)
  
  #calculate depth from AD, should be the same as depth
  dp_ad <- ref_mat + alt_mat
  
  #mask genotypes that fail DP/GQ
  geno_fail <- (dp_ad < min_dp | dp_ad > max_dp) | (gq < min_gq)
  
  dp_ad[geno_fail]   <- NA
  gq[geno_fail] <- NA
  ref_mat[geno_fail] <- NA
  alt_mat[geno_fail] <- NA
  
  #calculate allele balance (VAF)
  ab_mat <- alt_mat / dp_ad
  
  #classify genotypes based on AB
  ab_genotype <- function(ab) {
    geno <- rep(NA_character_, length(ab))
    
    geno[ab < 0.4] <- "REF"
    geno[ab >= 0.4 & ab < 0.9] <- "HET"
    geno[ab >= 0.9] <- "HOM_ALT"
    #since min depth is 10, mutant MUST have at least 4 ALT reads for
    #a SNP (het or hom alt) to be called.
    #for a hom alt, must have at least 9 ALT reads to pass
    geno
  }
  
  #Apply genotype classification
  geno_calls <- apply(ab_mat, 2, ab_genotype)
  
  #identify samples carrying ALT
  has_alt <- geno_calls == "HET" | geno_calls == "HOM_ALT"
  has_alt[is.na(geno_calls)] <- FALSE
  
  #keep site if at least one sample has a passing ALT call
  variant_in_subset <- rowSums(has_alt) >= 1
  
  keep <- variant_in_subset
  
  message(sprintf(
    "Subset (%d samples): %d / %d sites have at least one passing ALT call",
    length(samples),
    sum(keep),
    length(keep)
  ))
  
  list(
    vcf = sub[keep, ],
    geno_calls = geno_calls[keep, , drop = FALSE],
    ab_mat = ab_mat[keep, , drop = FALSE],
    ref_mat = ref_mat[keep, , drop = FALSE],
    alt_mat = alt_mat[keep, , drop = FALSE],
    dp_mat  = dp_ad[keep, , drop = FALSE],
    keep_idx = which(keep)
  )
}

#run the function for each specified subsets
m4_result   <- subset_and_classify(OG, m4_samples)

saveRDS(m4_result, "checkpoint_m4_INDELS.rds")

bsa1_result <- subset_and_classify(OG, new454BSA)
##############################
###############################
###################################





###########################################################
#NEXT STEP:Remove parental genotype SNP sites
#####################################################


#specify parental genotype sample names again
parentals <- c("MN1806CN-1_founder-prep04_Rep1-MN1806CN-1_founder-prep04_Rep1-CTTAAGCTTAAG",
               "MN1806CN-1_founder-prep04_Rep2-MN1806CN-1_founder-prep04_Rep2-TGGCTTTGGCTT",
               "MN1806CN-1_founder-prep11_Rep1-MN1806CN-1_founder-prep11_Rep1-GAGTGCGAGTGC",
               "MN1806CN-1_founder-prep11_Rep2-MN1806CN-1_founder-prep11_Rep2-ATGGTGATGGTG")

#provide subset vcf object and parentals list
remove_parental_alt <- function(result, parentals) {
  gc <- result$geno_calls
  ab <- result$ab_mat
  
  p <- gc[, parentals, drop = FALSE]
  
  #any parental site with a het or hom_alt call is classified as an INDEL
  any_INDEL <- apply(p, 1, function(x) {
    any(x == "HET" | x == "HOM_ALT", na.rm = TRUE)
  })
  #keep sites where there is NO INDEL call for any parental sample
  keep <- !any_INDEL
  
  message(sprintf("Sites with no parental INDEL: %d / %d (removed: %d)",
                  sum(keep), length(keep), sum(!keep)))
  
  list(
    #drop = FALSE ensures proper matrix is returned regardless of if there are few samples
    vcf        = result$vcf[keep, ],
    geno_calls = gc[keep, , drop = FALSE],
    ab_mat     = ab[keep, , drop = FALSE],
    ref_mat    = result$ref_mat[keep, , drop = FALSE],
    alt_mat    = result$alt_mat[keep, , drop = FALSE],
    dp_mat     = result$dp_mat[keep, , drop = FALSE]
  )
}


#prior to parental INDEL removal:
parental_gc_before <- m4_result$geno_calls[, parentals, drop = FALSE]

parental_alt_before <- apply(parental_gc_before, 1, function(x) {
  any(x == "HET" | x == "HOM_ALT", na.rm = TRUE)
})

cat("\nParental INDEL sites before removal:",
    sum(parental_alt_before), "\n")


#remove parental INDELS
m4_noparental<- remove_parental_alt(m4_result, parentals)


#CHECK THAT PARENTAL SNPs WERE REMOVED

parental_gc_after <- m4_noparental$geno_calls[, parentals, drop = FALSE]

parental_alt_after <- apply(parental_gc_after, 1, function(x) {
  any(x == "HET" | x == "HOM_ALT", na.rm = TRUE)
})

cat("\n========================================\n")
cat("PARENTAL SNP REMOVAL CHECK\n")
cat("========================================\n")
cat("Parental SNP sites remaining:",
    sum(parental_alt_after), "\n")

if (sum(parental_alt_after) == 0) {
  cat("CHECK PASSED: No parental HET/HOM_ALT sites remain.\n")
} else {
  cat("CHECK FAILED: Parental SNP sites remain\n")
}

cat("========================================\n\n")

stopifnot(sum(parental_alt_after) == 0)

bsa1_noparental <- remove_parental_alt(bsa1_result, parentals)


#get back the vcf with parentals removed and also the genotype calls
#should be row-matched to each other, but if not these check:
stopifnot(nrow(m4_noparental$vcf@gt)== nrow(m4_noparental$geno_calls))
stopifnot(nrow(bsa1_noparental$vcf@gt)== nrow(bsa1_noparental$geno_calls))


#checkpoint
saveRDS(m4_noparental, "m4_noparentalINDELS.rds")
saveRDS(bsa1_noparental, "bsa1_noparentalINDELS.rds")



####################################
#####################################
########################################
#calling singleton INDELS with sibling awareness

#from this point on we are only working with the non-bsa data
gc<- m4_noparental$geno_calls
ab <- m4_noparental$ab_mat



non_parental_samples <- setdiff(colnames(gc), parentals)
#setdiff removes the parental samples from the gc all samples
gc_np <- gc[, non_parental_samples]
#keep only the non parental samples at all sites for geno call and allele balance
ab_np <- ab[, non_parental_samples]



#subset new matrices to nonparental samples
ref_np <- m4_noparental$ref_mat[, non_parental_samples]
alt_np <- m4_noparental$alt_mat[, non_parental_samples]
dp_np  <- m4_noparental$dp_mat[, non_parental_samples]

stopifnot(identical(colnames(ref_np), colnames(gc_np)),
          identical(colnames(alt_np), colnames(gc_np)),
          identical(colnames(dp_np), colnames(gc_np)))




#classify + assign family IDs 
#Family ID = the MEE##### number ONLY
#code only keeps the leading digits after MEE

get_family_id <- function(sample_names) {
  sub("^(MEE[0-9]+).*", "\\1", sample_names)
}
#gets the sample name: grabs MEE and all following DIGITS, and removes everything
#after the digits (only keeps first grouping)
#so no letters/following -1 is included, just MEE#####
#Use this to identify samples of the same family

#create two-column df, one with non-parental sample names, one with family names
sample_info<- data.frame(sample=non_parental_samples, stringsAsFactors = FALSE)
#builds a data frame. Each row=each non-parental sample. 
sample_info$family<- get_family_id(sample_info$sample)
#has two columns, one with full sample name, and one with the family

fam_check<- sample_info %>% count(family, sort = TRUE) %>% filter(n > 1)
#checks which families have more than one member.
#groups by the family column, sorts from most to least, identifies families with >1
message("Families with >1 member:")
print(fam_check)
#identifies families with more than one sample, any time n>1 prints it
#example:
# family n
# MEE00053 2
# MEE00112 2

#M4 dataset shouldn't have more than 5 siblings per family
max_expected_siblings <- 5

#warning if it identifies a family with more than 5 siblings in it
#also a warning if less than 2 since fam_check should only have families with siblings
if (any(fam_check$n > max_expected_siblings) || any(fam_check$n < 2)) {
  warning("Unexpected family grouping detected")}



#define the mutant samples and families
mutant_samples<-sample_info$sample
#full list of non parental samples

mutant_families<-unique(sample_info$family)
#list of unique families

#confident (HET/HOM_ALT) calls.
confident_alt <- gc_np == "HET" | gc_np == "HOM_ALT"
#produces TRUE when either HET or HOM_ALT at each site for each sample
confident_alt[is.na(confident_alt)] <- FALSE
#NAs become FALSE


family_map <- setNames(sample_info$family, sample_info$sample)
#creates vector where the names are sample names, while the values are family IDs
#ex: family_map["MEE00053A_S12"] gives back "MEE00053"


#loop gives back a matrix with one row per site, one column per family
#and TRUE/FALSE for whether or not the family contains a confident ALT call at site
family_has_alt <- sapply(mutant_families, function(fam) {
  cols <- mutant_samples[family_map[mutant_samples] == fam]
  #figures out which samples belong to each current family being processed
  #pulls out the sample name(s) belonging to just that family
  if (length(cols) == 1) {
    #If there's only one sample in the family, then just grabs that sample's column
    confident_alt[, cols]
  } else {
    apply(confident_alt[, cols, drop = FALSE], 1, any)
    #if there's more than one sample, then subsets to the columns for those siblings
    #For each SNP site: checks whether there are any TRUE's
    #if one of them is a TRUE, then the whole family counts as an ALT carrier
    #drop=FALSE prevents conversion into a plain vector 
  }
})


colnames(family_has_alt) <- mutant_families
#makes sure names of columns are the mutant families

#counts how many families carry the ALT
family_alt_count <- rowSums(family_has_alt)

#family-level singleton definition-only one family can have an ALT call
family_singleton_sites <- family_alt_count == 1

# CHECKPOINT
saveRDS(
  list(
    family_has_alt = family_has_alt,
    family_alt_count = family_alt_count,
    family_singleton_sites = family_singleton_sites
  ),
  "checkpoint_singletonINDELS.rds"
)



message(sprintf("Family-level singleton sites (M4): %d", sum(family_singleton_sites)))


#identify which carrier family carries each singleton
#carrier family is a vector with one family ID per singleton site
carrier_family <- apply(family_has_alt[family_singleton_sites, , drop = FALSE], 1, function(x) {
  colnames(family_has_alt)[which(x)]
})
#subsets to the singleton rows, then for each row, which(x) finds the column index
#where the TRUE is 
#colnames(family_has_alt)[...] translates that index back into the actual
#family ID string
#Identifies which family carries each singleton



#extract the singlton sites from the vcf
m4_family_singletons_vcf <- m4_noparental$vcf[which(family_singleton_sites), ]
#subsets vcfR into just the rows/sites that are family-level singletons

#Check family-level singleton definition
print(table(family_alt_count))

stopifnot(all(family_alt_count[family_singleton_sites] == 1))

#Confirm no parental ALT classifications remain
stopifnot(sum(parental_alt_after) == 0)

cat("Family-level singleton sites:",
    sum(family_singleton_sites), "\n")

write.vcf(m4_family_singletons_vcf,   "M4_family_singletonINDELS.vcf.gz")




#next step: Return counts for indels 
###########################################
geno_calls_final <- gc_np[family_singleton_sites, , drop = FALSE]
ab_final <- ab_np[family_singleton_sites, , drop = FALSE]
ref_final <- ref_np[family_singleton_sites, , drop = FALSE]
alt_final <- alt_np[family_singleton_sites, , drop = FALSE]
dp_final  <- dp_np[family_singleton_sites, , drop = FALSE]

#site info, row-matched to above matrices
fix_s <- getFIX(m4_family_singletons_vcf)
stopifnot(nrow(fix_s) == nrow(ab_final))

site_info <- data.frame(
  CHROM = fix_s[, "CHROM"],
  POS   = as.integer(fix_s[, "POS"]),
  REF   = fix_s[, "REF"],
  ALT   = fix_s[, "ALT"],
  stringsAsFactors = FALSE
)

n_site <- nrow(ab_final)
samp_rep <- rep(colnames(ab_final), each = n_site)

ab_final_long <- data.frame(
  site_info[rep(seq_len(n_site), times = ncol(ab_final)), ],
  Sample    = samp_rep,
  Family    = unname(family_map[samp_rep]),
  AB        = as.vector(ab_final),
  Genotype  = as.vector(geno_calls_final),
  REF_count = as.vector(ref_final),
  ALT_count = as.vector(alt_final),
  DP        = as.vector(dp_final),
  stringsAsFactors = FALSE
)
rownames(ab_final_long) <- NULL

ab_final_long <- ab_final_long[
  !is.na(ab_final_long$AB) &
    ab_final_long$Genotype %in% c("HET", "HOM_ALT"), ]

#checks
stopifnot(all(abs(ab_final_long$AB - ab_final_long$ALT_count / ab_final_long$DP) < 1e-8))
stopifnot(all(ab_final_long$DP == ab_final_long$REF_count + ab_final_long$ALT_count))

#retain the AB information if desired
saveRDS(ab_final_long, "M4_final_singleton_AB_valuesINDELS.rds")





#####################################################
###################
###############

all_samples_final<-colnames(geno_calls_final)
#get all sample names in the vcf

sample_info <- data.frame(sample = all_samples_final, stringsAsFactors = FALSE)
sample_info$family <- get_family_id(sample_info$sample)
#reconstruct the family IDs

family_map<-setNames(sample_info$family, sample_info$sample)
#creates the sample->family lookup

has_alt_final <- geno_calls_final == "HET" | geno_calls_final == "HOM_ALT"
#The ALT genotypes are designated as HET or HOM_ALT

has_alt_final[is.na(geno_calls_final)] <- FALSE
#missing/na calls count as FALSE

het_counts<- colSums(geno_calls_final == "HET", na.rm = TRUE)
#count the number of hets for each sample (column)
hom_alt_counts<-colSums(geno_calls_final == "HOM_ALT", na.rm = TRUE)
#count the number of Homalts per sample (column)
total_counts<- het_counts + hom_alt_counts
#total of alt calls

#summary table
per_sample<- data.frame(
  sample= colnames(geno_calls_final),
  HET = het_counts,
  HOM_ALT = hom_alt_counts,
  TOTAL = total_counts
) %>%
  mutate(PCT_HET = round(HET / TOTAL * 100, 1)) %>%
  arrange(desc(TOTAL))

#checks in case sums are off
#also checks for duplicated sample names
stopifnot(all(per_sample$TOTAL == per_sample$HET + per_sample$HOM_ALT))
#sum of per sample hets plus hom_alts should equal total
stopifnot(!any(duplicated(colnames(geno_calls_final))))

print(per_sample)



message(sprintf("Total HET INDELS : %d", sum(het_counts)))

message(sprintf("Total HOM_ALT INDELS : %d", sum(hom_alt_counts)))

message(sprintf("Total singleton INDELS : %d", sum(total_counts)))

message(sprintf("Average INDELS per sample : %.1f", mean(per_sample$TOTAL)))

message(sprintf("Overall %% HET INDELS : %.1f%%", sum(het_counts) / sum(total_counts) * 100))






###############################################################
#Sibling sharing summary for families with multiple siblings
#####################################################
#figures out how the siblings share their family-singletons
#converts singltons into a table for SNP spectra analysis

summarize_family_sharing <- function(fam_id, family_map, has_alt_final) {
  #takes family info plus the TRUE/FALSE matrix telling us if each sample has an ALT
  
  members <- names(family_map)[family_map == fam_id & names(family_map) %in% colnames(has_alt_final)]
  #finds the members of that family present in the vcf
  
  if (length(members) < 2) {
    message(sprintf("Family %s has <2 members in this VCF - skipping sharing summary", fam_id))
    return(invisible(NULL))
  }
  #skips the family if there aren't two members or more
  
  fam_sites <- rowSums(has_alt_final[, members, drop = FALSE]) > 0
  #finds sites where at least one sibling has an ALT
  
  sib_alt_count <- rowSums(has_alt_final[fam_sites, members, drop = FALSE])
  #counts how many siblings carry an ALT
  
  message(sprintf("\n--- Family %s (%d members: %s) ---",
                  fam_id, length(members), paste(members, collapse = ", ")))
  print(table(sib_alt_count))
  invisible(list(fam_sites = fam_sites, sib_alt_count = sib_alt_count, members = members))
}

#loop through all families with multiple siblings and give output in named list
#use lapply to loop
#fam_check$family is the list of family IDs with multiple siblings
family_sharing_results <- lapply(fam_check$family, function(fam_id) {
  summarize_family_sharing(fam_id, family_map, has_alt_final)})
#produces list with one entry per family with the family sites, the sibling_alt_count, and members information
names(family_sharing_results) <- fam_check$family
#labels each entry in the list with its family ID
#can use family_sharing_results[["MEE00053"]] to look up specific families


#combine it into one big table
#creates df with one row per singleton site in the family
#displays which family it belongs to, total number of siblings in that family, 
#and how many siblings carry the alt allele at the site
sharing_summary <- do.call(rbind, lapply(names(family_sharing_results), function(fam_id) {
  res <- family_sharing_results[[fam_id]]
  if (is.null(res)) return(NULL)
  data.frame(
    family = fam_id,
    n_siblings = length(res$members),
    sib_alt_count = res$sib_alt_count
  )
}))


sharing_summary %>%
  count(family, n_siblings, sib_alt_count) %>%
  #groups by family size and number of siblings sharing the SNP
  mutate(category = case_when(
    sib_alt_count == 1 ~ "private",
    #1 sample is a private SNP
    sib_alt_count == n_siblings ~ "shared_all",
    #all samples in family have the SNP
    TRUE ~ "shared_partial"
    #designate as "shared partial" when more than 1 but less than n share SNP
  ))


###############################################################
#############################################

#build spectra DF for indels

#get the REF and the ALT so you know what became what (can classify the INDELS)
ref_allele <- getFIX(m4_family_singletons_vcf)[, "REF"]
alt_allele <- getFIX(m4_family_singletons_vcf)[, "ALT"]

#calculate allele lengths and indel sizes
indel_spectra <- data.frame(
  REF = ref_allele,
  ALT = alt_allele,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    REF_length = nchar(REF),
    ALT_length = nchar(ALT),
    indel_size = abs(ALT_length - REF_length),
    indel_type = dplyr::case_when(
      ALT_length > REF_length ~ "Insertion",
      ALT_length < REF_length ~ "Deletion",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(indel_size > 0, indel_type != "Other")

#plot indel size distrib
library(ggplot2)

SISD<- ggplot(indel_spectra, aes(x = indel_size)) +
  geom_bar() +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    title = "Singleton Indel Size Distribution",
    x = "Indel size (bp)",
    y = "Number of indels"
  ) +
  theme_bw()

#optional: save plot as svg
library(svglite)
ggsave("SISD_M4_plot.svg", plot = SISD, width = 8, height = 5, units = "in")


#compare insertions vs deletiosn
#double bar plot for counts of insertions and deletions at each length
IvsD<- ggplot(indel_spectra, aes(x = indel_size, fill = indel_type)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Singleton Insertion and Deletion Size Distribution",
    x = "Indel size (bp)",
    y = "Number of indels",
    fill = "Indel type"
  ) +
  theme_bw()
#optional: save plot as svg
ggsave("IvsD_M4_plot.svg", plot = IvsD, width = 8, height = 5, units = "in")


#summarize size distribution for the indels
indel_spectra %>%
  dplyr::count(indel_size, indel_type, name = "count") %>%
  dplyr::arrange(indel_size, indel_type)

