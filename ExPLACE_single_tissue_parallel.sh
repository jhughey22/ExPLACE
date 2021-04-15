#!/usr/bin/env bash
#Jordan Hughey
#Liu Group

tiss=lcl
tissue=lymphoblastoid

for i in `seq 21 22`
do
	chrom=${i}
	echo "running chrom${chrom}"
	prefix=${tiss}_ExPLACE_nested_cv
	out_dir=`pwd`
	cis_window=1e6
	alpha=0.5
	n_folds=10
	n_train_test_folds=5
	

	expression_RDS=./example_input/test_exp_matrix_comb.RDS
	geno_file=./example_input/test_genotype_matrix_chr${chrom}.txt
	gene_annot_RDS=./example_input/gencode.v26lift37.annotation.gene.only.parsed.RDS
	snp_annot_RDS=./example_input/test_snp_annot.chr${chrom}.RDS
	covariates_file=NA
	gene_snps_RDS=./example_input/snp-gene_pairs/GM12878_newABC_sig_snplist_1mb.RDS	        

	qsub -d ${out_dir} -N ${tiss}_ExP_chr${chrom} -v chrom="${chrom}",prefix="${prefix}",out_dir="${out_dir}",cis_window="${cis_window}",alpha="${alpha}",n_folds="${n_folds}",n_train_test_folds="${n_train_test_folds}",expression_RDS="${expression_RDS}",geno_file="${geno_file}",gene_annot_RDS="${gene_annot_RDS}",snp_annot_RDS="${snp_annot_RDS}",covariates_file="${covariates_file}",gene_snps_RDS="${gene_snps_RDS}" ./scripts/ExPLACE_single_tissue.pbs 

done

