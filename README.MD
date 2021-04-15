<!-- ABOUT THE PROJECT -->
## ExPLACE

ExPLACE (**Ex**pression **P**rediction **L**everaging **A**ctivity and **C**ontact **E**vidence) is an integrative method that uses chromatin conformation data and enhancer activity marks to choose the genetic variants used for predicting gene expression.  Ultimately, this tool creates biologically informed gene expression prediction models for transcriptome-wide association studies. 

<p align="center">
   <img src="https://github.com/jhughey22/ExPLACE/raw/main/image/single_tissue_ExPLACE.png" width="734" height="300.7">
</p>

<!-- GETTING STARTED -->
## Prerequisites

To run ExPLACE, you need [R](https://www.r-project.org/) with several packages including

[dplyr](https://cran.r-project.org/web/packages/dplyr/index.html), [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html), and [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)

<!-- Preparing Input -->
## ExPLACE Input data

ExPLACE will train models for gene expression prediction.  There are six main inputs that you can specifiy by using a text editor to change inputs in the `ExPLACE_single_tissue_run.sh` wrapper shell script.  Some preprocessing will be needed prior to execution to format the following inputs:

1) expression_RDS - Expression matrix by tissue as a R data object.  Rows are the sample IDs and columns are the ensemble gene ids. Should look as follows:

|        | ENSG001 | ENSG002 | ENSG003 | 
| :----: | ------- | ------- | ------- | 
|  ID1   | 0.0134  | 10.047  | 0.667   | 
|  ID2   | 0.0269  | 9.8065  | 0.577   |

2) geno_file - genotype matrix by chromosome as a tab delimited file.  Rows are the genetic variant ID (ie. chr_pos_refall_effall_build) and columns are the sample IDs.  Should look as follows:

| ID                     | ID1 | ID2 | ID3 |
| :--------------------: | --- | --- | --- |
|  22_16051249_T_C_b37   | 0   | 1   | 2   |
|  22_16051453_A_C_b37   | 2   | 1   | 0   |

* note "ID" must be colnames for first field

3) gene_annot_RDS - gene annotation file. Supplied in the `example_input` folder as `gencode.v26lift37.annotation.gene.only.parsed.RDS`

4) snp_annot_RDS - snp annotation file as R data object.  Rows are the genetic variant ID and columns the following:

|                        | chr | pos      | varID               | refAllele | effectAllele | rsid |
| :--------------------: | --- | -------- | ------------------  | --------- | ------------ | ---- |
|  22_16051249_T_C_b37   | 22  | 16051249 | 22_16051249_T_C_b37 | T         | C            | rs1  |
|  22_16051453_A_C_b37   | 22  | 16051453 | 22_16051453_A_C_b37 | A         | C            | rs2  | 

* note rsid can just be `NA` since we focus on varID

5) covariates_file - covariate data for adjusting expression values.  We use PEER Factors (as described by GTEx), Sex, Top 3 genetic PCs, and Sequencing platform.  Looks like:

| ID                | ID1   | ID2   | ID3   |
| :---------------: | ----- | ----- | ----- |
|  PC1              | 0.112 | 0.233 | -0.03 |
|  PC2              | -0.23 | -0.11 | -0.08 |
|  PC3              | 0.032 | -0.08 | 0.002 |
|  PEER_Factor1     | 0.053 | -0.72 | 0.552 |
|  PEER_Factor2     | -0.13 | 0.323 | -0.04 |
|  Sex              | 1     | 2     | 2     |
|  Platform         | 2     | 1     | 1     |

* if there is no covariates file added denote input parameter in the `ExPLACE_single_tissue_run.sh` wrapper script as `NA`

6) gene_snps_RDS - This is a R data object that consists of a list where each item in the list is a gene and its contents are all the snps that pair with this gene based off activity and contact information.  We provide snp-gene RDS objects for 13 tissues in the `example_input/snp-gene_pairs/` directory.  We will be creating another repository documenting our analysis for selecting snp-gene pairs using activity and contact evidence soon.


<!-- USAGE EXAMPLES -->
## Setup and Usage Example

1) Clone this repository.
```bash
$ git clone https://github.com/jhughey22/ExPLACE
```

2) Enter ExPLACE folder.
```bash
$ cd ExPLACE
```

3) Run the ExPLACE wrapper shell script on example data
```bash
bash ExPLACE_single_tissue_run.sh 
```
This should take ~4mins on a 3GHZ computer

The output with model performance statistics is stored in the `output/summary` folder as `lcl_ExPLACE_nested_cv_chr22_model_summaries.txt`.  The first 7 fields will look like:

| gene_id            | gene_name | gene_type      | alpha | n_snps_in_window | n_snps_in_model | lambda_min_mse     |
| :----------------: | --------- | -------------- | ----- | ---------------- | --------------- | ------------------ |
|  ENSG00000070371   | CLTCL1    | protein_coding | 0.5   | 5197             | 35              | 0.0108200838909036 | 
|  ENSG00000070413   | DGCR2     | protein_coding | 0.5   | 36               | 1               | 0.365140539710332  | 

Full specification of the output fields are documented [here](https://github.com/jhughey22/ExPLACE/tree/main/example_output/output_fields.txt)

We also output the computed weights for snps in each gene model in the `output/weights` folder and also covariances in the `output/covariances` folder

Furthermore, precomputed outputs for the example data are in the `example_output` directories

#### Running in parallel using PBS scheduler

We highly recommend running ExPLACE on a HPC preferably on a job scheduling software.  We provide a script to run model training in parallel across chromosomes with the `ExPLACE_single_tissue_parallel.sh` shell script.  This script is wrapper that submits a PBS job for each chromosome.  It calls the `/scripts/ExPLACE_single_tissue.pbs` PBS script that has the specific PBS job allocation parameters.  If your cluster utilizes PBS you can run the following on the provided example data

```bash
bash ExPLACE_single_tissue_parallel.sh
```

The input parameters are the same as the `ExPLACE_single_tissue_run.sh` shell script.  The only difference is we are running the jobs using the PBS scheduler and in parallel for each chromosome

<!-- FUTURE IMPLEMENTATION -->
## Future Directions

Coming Soon!:

Documentation for using epigonomic and 3D genomic data for snp-gene pair matching

Implementation of Ens-ExPLACE (ensemble approach that adapts ExPLACE to cross-tissue imputation)

<!-- CONTACT -->
## Contact

Jordan Hughey - jordanhugheypsu@gmail.com

<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [Dajiang J. Liu](https://dajiangliu.blog/)
* Chachrit (Poom) Khunsriraksakul - [@ChachritK](https://twitter.com/ChachritK) - ckhunsr2@gmail.com

Part of the code is modified from MetaXcan https://github.com/hakyimlab/MetaXcan. We thank the authors for sharing the code.



