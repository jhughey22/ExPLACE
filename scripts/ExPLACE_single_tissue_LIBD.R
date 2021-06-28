#!/usr/bin/env Rscript
#Jordan Hughey
#Liu Group

#This Rscript runs the ExPLACE main function that comes from the './scripts/gtex_v7_nested_cv_elnet_weighted_hybrid.R' script

#Takes args from wrapper scripts in repo
argv <- commandArgs(trailingOnly = TRUE)
#source('./scripts/gtex_v7_nested_cv_elnet_weighted_hybrid_LIBD.R')
source('./scripts/gtex_v7_nested_cv_elnet_weighted_hybrid_LIBD_v4.R')

chrom <- argv[1]
prefix <- argv[2]
out_dir <- argv[3]
cis_window <- as.numeric(argv[4])
alpha <- as.numeric(argv[5])
n_folds <- as.numeric(argv[6])
n_train_test_folds <- as.numeric(argv[7])
expression_RDS <- argv[8]
geno_file <- argv[9]
gene_annot_RDS <- argv[10]
snp_annot_RDS <- argv[11]
covariates_file <- argv[12]
gene_snps_RDS <- argv[13]

#Based off presence of covariates file we call function differently
if (covariates_file == 'NA') {
  main(snp_annot_RDS, gene_annot_RDS, geno_file, expression_RDS, covariates_file = NA,
       chrom, prefix, n_folds, n_train_test_folds, cis_window, alpha, out_dir, gene_snps_RDS)
} else {
  main(snp_annot_RDS, gene_annot_RDS, geno_file, expression_RDS, covariates_file, 
       chrom, prefix, n_folds, n_train_test_folds, cis_window, alpha, out_dir, gene_snps_RDS)
}



