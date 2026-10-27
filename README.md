# M4_bioinformatics
This is the entire set of scripts from initial resequencing data to filtered SNP and InDel variant effect predictions for each mutant line.
Each script is numbered from 1-18 in order of usage.

1-3: initial preprocessing, including fastqc analysis and trimming of adaptors

4-8: alignment of reads to Williams82 va6 soybean reference genome, and subsequent removal of duplicates and multimapped reads

9-12: Running GATK haplotype caller on all samples, conducted in 250kb genomic chunks. VCFs were gathered into their respective chromosome vcfs, then concatenated to produce one cohort VCF file

13: Hard filtering variants with GATK, following standard methods for SNPs and standard methods for INDELS other than increased QUAL and MQ thresholds

14-15: Final filtering of SNPs and Indels based on MQ and depth, and removal of control genotype natural variants and variants shared between unrelated lines

16: Variant Effect Predictor (VEP) used to generate inital variant annotations

17-18: Generation of csv files containing HIGH and MODERATE impact SNP and INDEL variants
