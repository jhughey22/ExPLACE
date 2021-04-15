#!/usr/bin/env bash
#Jordan Hughey
#Liu Group

tiss=lcl_test
tissue=lymphoblastoid

chrom=22
echo "running chrom${chrom}"
prefix=${tiss}_ExPLACE_nested_cv
out_dir=`pwd`
cis_window=1e6
alpha=0.5
n_folds=10
n_train_test_folds=5
	

expression_RDS=./example_input/test_exp_matrix_chr22.RDS
geno_file=./example_input/test_genotype_matrix_chr22.txt
gene_annot_RDS=./example_input/gencode.v26lift37.annotation.gene.only.parsed.RDS
snp_annot_RDS=./example_input/test_snp_annot.chr22.RDS
covariates_file=NA
gene_snps_RDS=./example_input/snp-gene_pairs/GM12878_newABC_sig_snplist_1mb.RDS

Rscript ./scripts/ExPLACE_single_tissue.R ${chrom} ${prefix} ${out_dir} ${cis_window} ${alpha} ${n_folds} ${n_train_test_folds} \
        ${expression_RDS} ${geno_file} ${gene_annot_RDS} ${snp_annot_RDS} ${covariates_file} ${gene_snps_RDS}




