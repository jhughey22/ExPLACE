#!/usr/bin/env bash
#Jordan Hughey
#Liu Group

#This script is the wrapper to run ExPLACE model training
#This wrapper submits a job to the PBS scheduler and does so in parallel for each chromsome

#All the options are hard coded below and can be changed for data of interest

#Tissue the gene expression model training is performed
tiss=lcl	#abbreviation for longer tissue name
tissue=lymphoblastoid

#loop through chromosomes; can change this to `seq 1 22` X
for i in `seq 21 22`
do
	chrom=${i}
	echo "running chrom${chrom}"
	prefix=${tiss}_ExPLACE_nested_cv
	out_dir=`pwd`		#outputs stored here
	cis_window=1e6		#Window for picking snps in proximity of gene
	alpha=0.5		#Elastic net mixing parameter
	n_folds=10		#Inner loop cross validation
	n_train_test_folds=5	#Outer loop folds
	

	expression_RDS=./example_input/test_exp_matrix_comb.RDS
	geno_file=./example_input/test_genotype_matrix_chr${chrom}.txt
	gene_annot_RDS=./example_input/gencode.v26lift37.annotation.gene.only.parsed.RDS
	snp_annot_RDS=./example_input/test_snp_annot.chr${chrom}.RDS
	#Covariates file; if there isn't one keep as NA
	covariates_file=NA
	#gene-snps pairs file; precomputed for 14 tissues in repo
	gene_snps_RDS=./example_input/snp-gene_pairs/GM12878_newABC_sig_snplist_1mb.RDS	        

	#qsub execution using all options above
	qsub -d ${out_dir} -N ${tiss}_ExP_chr${chrom} -v chrom="${chrom}",prefix="${prefix}",out_dir="${out_dir}",cis_window="${cis_window}",alpha="${alpha}",n_folds="${n_folds}",n_train_test_folds="${n_train_test_folds}",expression_RDS="${expression_RDS}",geno_file="${geno_file}",gene_annot_RDS="${gene_annot_RDS}",snp_annot_RDS="${snp_annot_RDS}",covariates_file="${covariates_file}",gene_snps_RDS="${gene_snps_RDS}" ./scripts/ExPLACE_single_tissue.pbs 

done

