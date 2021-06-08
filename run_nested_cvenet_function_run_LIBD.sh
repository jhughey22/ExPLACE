#!/usr/bin/env bash
#Jordan Hughey
#Liu Group

#This script is the wrapper to run ExPLACE model training in a single chromosome
#All the options are hard coded below and can be changed for data of interest

#Tissue the gene expression model training is performed
tiss=lcl_test   #abbreviation for longer tissue name
tissue=lymphoblastoid

chrom=22
echo "running chrom${chrom}"
prefix=${tiss}_pred_nested_cv
out_dir=`pwd`  #outputs stored here
cis_window=1e6	#Window for picking snps in proximity of gene
alpha=0.5	#Elastic net mixing parameter
n_folds=10	#Inner loop cross validation
n_train_test_folds=5	#Outer loop folds
	

expression_RDS=./example_input/test_exp_matrix_chr22.RDS
geno_file=./example_input/test_genotype_matrix_chr22.txt
gene_annot_RDS=./example_input/gencode.v26lift37.annotation.gene.only.parsed.RDS
snp_annot_RDS=./example_input/test_snp_annot.chr22.RDS
#Covariates file; if there isn't one keep as NA
covariates_file=NA

#Call to the Rscript which calls the ExPLACE function
Rscript ./scripts/run_nested_cvenet_function_LIBD.R ${chrom} ${prefix} ${out_dir} ${cis_window} ${alpha} ${n_folds} ${n_train_test_folds} \
        ${expression_RDS} ${geno_file} ${gene_annot_RDS} ${snp_annot_RDS} ${covariates_file} 




