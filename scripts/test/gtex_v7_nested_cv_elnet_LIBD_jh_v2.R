#!/usr/bin/env Rscript
#Jordan Hughey
#Liu Group

suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))
"%&%" <- function(a,b) paste(a,b, sep = "")

#set functions

get_gene_type <- function(gene_annot, gene) {
  filter(gene_annot, gene_id == gene)$gene_type
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window)) & (pos <= (coords[2] + cis_window)))
  if (nrow(snp_info) == 0)
    return(NA)
  gt_df <- as.data.frame(gt_df)
  cis_gt <- gt_df %>% select(one_of(intersect(snp_info$varID, colnames(gt_df))))
  column_labels <- colnames(cis_gt)
  row_labels <- rownames(cis_gt)
  # Convert cis_gt to a matrix for glmnet
  cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt)) # R is such a bad language.
  colnames(cis_gt) <- column_labels
  rownames(cis_gt) <- row_labels
  cis_gt
}

get_covariates <- function(covariate_file_name, samples) {
  cov_df <- read.table(covariate_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  #cov_df_cols <- gsub("GTEX.", "GTEX-", colnames(cov_df))
  #colnames(cov_df) <- cov_df_cols
  cov_df <- cov_df[,samples] %>% t() %>% as.data.frame()
  cov_df
}

generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
  
  #This is lines to do if you want to keep fold breakdown same over methods
  #train_test_fold_ids <- sample(fold_ids[1:n_samples])
  #train_test_fold_ids_out <- './Whole_Blood_train_test_fold_ids.RDS'
  #saveRDS(train_test_fold_ids, train_test_fold_ids_out)
  
  #may need to use this to keep fold ids constant between methods
  
  #also do for cv fold ids without nested cv for fitting all data
  #cv_fold_ids <- generate_fold_ids(n_samples, n_folds)
  #cv_fold_ids_out <- './Whole_Blood_cv_fold_ids.RDS'
  #saveRDS(cv_fold_ids, cv_fold_ids_out)
}

adjust_for_covariates <- function(expression_vec, cov_df) {
  # combined_df <- cbind(expression_vec, cov_df)
  # expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  # expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  # expr_resid
  #first adjust covariates just for LIBD set study == 1
  cov_df_libd <- subset(cov_df, cov_df$Study == 1)
  cov_df_libd_final <- cov_df_libd[,'Smoking.status.n', drop=FALSE]
  expression_vec_libd_names <- intersect(names(expression_vec), rownames(cov_df_libd))
  expression_vec_libd <- expression_vec[expression_vec_libd_names]
  
  combined_df_libd <- cbind(expression_vec_libd, cov_df_libd_final)
  expr_resid_libd <- summary(lm(expression_vec_libd ~ ., data=combined_df_libd))$residuals
  expr_resid_libd <- scale(expr_resid_libd, center = TRUE, scale = TRUE)
  
  #now add the above adjusted expression vec from LIBD back to GTEx samples and do full adjustment
  expression_vec_gtex_names <- setdiff(names(expression_vec), rownames(cov_df_libd))
  expression_vec_gtex <- expression_vec[expression_vec_gtex_names]
  expression_vec_gtex_df <- as.data.frame(expression_vec_gtex)
  colnames(expression_vec_gtex_df) <- 'V1'
  expression_vec_combined_df <- rbind(expression_vec_gtex_df, expr_resid_libd) 
  
  #now cbind df but match by name
  combined_df_full <- merge(cov_df, expression_vec_combined_df, by = 0, sort = F)
  col_order <- c(1, ncol(combined_df_full), 2:(ncol(combined_df_full)-1))
  combined_df_final <- combined_df_full[, col_order]
  rownames(combined_df_final) <- combined_df_final$Row.names
  combined_df_final <- subset(combined_df_final, select = -c(Row.names, Study, Smoking.status.n))
  colnames(combined_df_final)[[1]] <- 'expression_vec'
  
  #combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df_final))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

calc_R2 <- function(y, y_pred) {
  tss <- sum(y**2)
  rss <- sum((y - y_pred)**2)
  1 - rss/tss
}

calc_corr <- function(y, y_pred) {
  sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))
}

nested_cv_elastic_net_perf <- function(x, y, n_samples, n_train_test_folds, n_k_folds, alpha) {
  # Gets performance estimates for k-fold cross-validated elastic-net models.
  # Splits data into n_train_test_folds disjoint folds, roughly equal in size,
  # and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is
  # cross validated. Then get performance measures for how the model predicts on the hold-out
  # fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis
  # is there is no correlation between prediction and observed.
  #
  # The mean and standard deviation of R^2 over all folds is then reported, and the p-values
  # are combined using Fisher's method.
  R2_folds <- rep(0, n_train_test_folds)
  corr_folds <- rep(0, n_train_test_folds)
  zscore_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  # Outer-loop split into training and test set.
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  #train_test_fold_ids <- readRDS(train_test_fold_ids_RDS)
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    x_train <- x[train_idxs, ]
    y_train <- y[train_idxs]
    x_test <- x[test_idxs, ]
    y_test <- y[test_idxs]
    # Inner-loop - split up training set for cross-validation to choose lambda.
    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
    y_pred <- tryCatch({
      # Fit model with training data.
      fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids)
      # Predict test data using model that had minimal mean-squared error in cross validation.
      predict(fit, x_test, s = 'lambda.min')},
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) rep(mean(y_train), length(y_test)))
    R2_folds[test_fold] <- calc_R2(y_test, y_pred)
    # Get p-value for correlation test between predicted y and actual y.
    # If there was no model, y_pred will have var=0, so cor.test will yield NA.
    # In that case, give a random number from uniform distribution, which is what would
    # usually happen under the null.
    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
  }
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  rho_avg <- mean(corr_folds)
  rho_se <- sd(corr_folds)
  rho_avg_squared <- rho_avg**2
  # Stouffer's method for combining z scores.
  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  list(R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}

do_covariance <- function(gene_id, cis_gt, rsids, varIDs) {
  model_gt <- cis_gt[,varIDs, drop=FALSE]
  #colnames(model_gt) <- rsids
  colnames(model_gt) <- varIDs
  geno_cov <- cov(model_gt)
  geno_cov[lower.tri(geno_cov)] <- NA
  #cov_df <- melt(geno_cov, varnames = c("rsid1", "rsid2"), na.rm = TRUE) %>%
  cov_df <- melt(geno_cov, varnames = c("varID1", "varID2"), na.rm = TRUE) %>%
    mutate(gene=gene_id) %>%
    #select(GENE=gene, RSID1=rsid1, RSID2=rsid2, VALUE=value) %>%
    select(GENE=gene, VARID1=varID1, VARID2=varID2, VALUE=value) %>%
    #arrange(GENE, RSID1, RSID2)
    arrange(GENE, VARID1, VARID2)
  cov_df
}

###################################################################
#Main function where all gets pulled together for execution
###################################################################

main <- function(snp_annot_RDS, gene_annot_RDS, geno_file, expression_RDS, 
                 covariates_file=NA, chrom, prefix, n_folds, n_train_test_folds, 
                 cis_window, alpha, out_dir, seed=NA) {
  
  #read in all files
  #exp file first
  expr_df <- readRDS(expression_RDS)
  expr_df <- as.matrix(expr_df)
  class(expr_df) <- 'numeric'
  
  #now genotype file
  gt_df <- read.table(geno_file, header = TRUE, row.names = 'ID', stringsAsFactors = FALSE)
  # Transpose genotype for glmnet
  gt_df <- t(gt_df)
  
  #gene annot file
  gene_annot <- readRDS(gene_annot_RDS)
  gene_annot <- subset(gene_annot, gene_annot$chr == chrom)
  gene_types=c('protein_coding', 'pseudogene', 'lincRNA')
  gene_annot <- subset(gene_annot, gene_annot$gene_type %in% gene_types)
  
  #read snp_annot file
  snp_annot <- readRDS(snp_annot_RDS)
  
  # Subset expression data to only include genes with gene_info
  expr_df <- expr_df[, intersect(colnames(expr_df), rownames(gene_annot))]
  
  #get sample numbers
  samples <- rownames(expr_df)
  n_samples <- length(samples)
  genes <- colnames(expr_df)
  n_genes <- length(genes)

  #get covariates file
  #covariates_df <- get_covariates(covariates_file, samples)
  #get covariates file
  if(!is.na(covariates_file)) {
    covariates_df <- get_covariates(covariates_file, samples)
  }
  
  # Set seed----
  seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
  set.seed(seed)
  
  # Prepare output data----
  model_summary_file <- out_dir %&% '/output/summary/' %&% prefix %&% '_chr' %&% chrom %&% '_model_summaries.txt'
  model_summary_cols <- c('gene_id', 'gene_name', 'gene_type', 'alpha', 'n_snps_in_window', 'n_snps_in_model', 'lambda_min_mse',
                          'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2',
                          'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',
                          'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')
  write(model_summary_cols, file = model_summary_file, ncol = 24, sep = '\t')
  
  weights_file <- out_dir %&% '/output/weights/' %&% prefix %&% '_chr' %&% chrom %&% '_weights.txt'
  weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')
  write(weights_col, file = weights_file, ncol = 6, sep = '\t')
  
  tiss_chr_summ_f <- out_dir %&% '/output/summary/' %&% prefix %&% '_chr' %&% chrom %&% '_tiss_chr_summary.txt'
  tiss_chr_summ_col <- c('n_samples', 'chrom', 'cv_seed', 'n_genes')
  tiss_chr_summ <- data.frame(n_samples, chrom, seed, n_genes)
  colnames(tiss_chr_summ) <- tiss_chr_summ_col
  write.table(tiss_chr_summ, file = tiss_chr_summ_f, quote = FALSE, row.names = FALSE, sep = '\t')
  
  covariance_file <- out_dir %&% '/output/covariances/' %&% prefix %&% '_chr' %&% chrom %&% '_covariances.txt'
  covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
  write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')
  
  # Attempt to build model for each gene----
  #i <- 2
  cat("Processing Gene \n")
  for (i in 1:n_genes) {
    cat(i, "/", n_genes, "\n")
    gene <- genes[i]
    gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene]
    gene_type <- get_gene_type(gene_annot, gene)
    coords <- get_gene_coords(gene_annot, gene)
    cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)
    if (all(is.na(cis_gt))) {
      # No snps within window for gene.
      model_summary <- c(gene, gene_name, gene_type, alpha, 0, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      write(model_summary, file = model_summary_file, append = TRUE, ncol = 24, sep = '\t')
      next
    }
    model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    if (ncol(cis_gt) >= 2) {
      expression_vec <- expr_df[,i]
      
      #if covariates file exists, adjust expression      
      if(!is.na(covariates_file)) {
        adj_expression <- adjust_for_covariates(expression_vec, covariates_df)
      } else {
        adj_expression <- expression_vec
      }
      #adj_expression <- adjust_for_covariates(expression_vec, covariates_df)
      
      #if (null_testing)
      #  adj_expression <- sample(adj_expression)
      #adj_expression <- sample(adj_expression)
      #adj_expression <- expression_vec
      perf_measures <- nested_cv_elastic_net_perf(cis_gt, adj_expression, n_samples, n_train_test_folds, n_folds, alpha)
      R2_avg <- perf_measures$R2_avg
      R2_sd <- perf_measures$R2_sd
      pval_est <- perf_measures$pval_est
      rho_avg <- perf_measures$rho_avg
      rho_se <- perf_measures$rho_se
      rho_zscore <- perf_measures$rho_zscore
      rho_avg_squared <- perf_measures$rho_avg_squared
      zscore_pval <- perf_measures$zscore_pval
      # Fit on all data
      cv_fold_ids <- generate_fold_ids(length(adj_expression), n_folds)
      #cv_fold_ids <- readRDS(all_cv_fold_ids_RDS)
      fit <- tryCatch(cv.glmnet(cis_gt, adj_expression, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE),
                      error = function(cond) {message('tryCatch - Warning'); message(geterrmessage()); list()})
      if (length(fit) > 0) {
        cv_R2_folds <- rep(0, n_folds)
        cv_corr_folds <- rep(0, n_folds)
        cv_zscore_folds <- rep(0, n_folds)
        cv_pval_folds <- rep(0, n_folds)
        best_lam_ind <- which.min(fit$cvm)
        for (j in 1:n_folds) {
          fold_idxs <- which(cv_fold_ids == j)
          adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
          cv_R2_folds[j] <- calc_R2(adj_expression[fold_idxs], adj_expr_fold_pred)
          cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, adj_expression[fold_idxs]), 0)
          cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(adj_expression[fold_idxs]) - 3) # Fisher transformation
          cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, adj_expression[fold_idxs])$p.value, runif(1))
        }
        cv_R2_avg <- mean(cv_R2_folds)
        cv_R2_sd <- sd(cv_R2_folds)
        adj_expr_pred <- predict(fit, as.matrix(cis_gt), s = 'lambda.min')
        training_R2 <- calc_R2(adj_expression, adj_expr_pred)
        
        cv_rho_avg <- mean(cv_corr_folds)
        cv_rho_se <- sd(cv_corr_folds)
        cv_rho_avg_squared <- cv_rho_avg**2
        # Stouffer's method for combining z scores.
        cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_folds)
        cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
        cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_folds, lower.tail = F)
        
        if (fit$nzero[best_lam_ind] > 0) {
          
          weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
          weighted_snps <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
          weighted_snps_info <- snp_annot %>% filter(varID %in% weighted_snps) %>% select(rsid, varID, refAllele, effectAllele)
          weighted_snps_info$gene <- gene
          weighted_snps_info <- weighted_snps_info %>%
            merge(data.frame(weights = weights, varID=weighted_snps), by = 'varID') %>%
            select(gene, rsid, varID, refAllele, effectAllele, weights)
          write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
          covariance_df <- do_covariance(gene, cis_gt, weighted_snps_info$rsid, weighted_snps_info$varID)
          write.table(covariance_df, file = covariance_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
          model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), fit$nzero[best_lam_ind], fit$lambda[best_lam_ind], R2_avg, R2_sd, cv_R2_avg, cv_R2_sd, training_R2, pval_est,
                             rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval, cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
        } else {
          model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, fit$lambda[best_lam_ind], R2_avg, R2_sd,
                             cv_R2_avg, cv_R2_sd, training_R2, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                             cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
        }
      } else {
        model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, R2_avg, R2_sd, NA, NA, NA, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                           NA, NA, NA, NA, NA, NA)
      }
    }
    write(model_summary, file = model_summary_file, append = TRUE, ncol = 24, sep = '\t')
  }
}

# #Set variables to run main function
# expression_RDS <- '../../data/intermediate/expression_phenotypes/Whole_Blood.v7.normalized_expression_tr_euro_pca.RDS'
# geno_file <- '../../data/intermediate/genotypes/GTEx_V7_DGN_merge_genotype_final_allchrs_WB_ids_pca.matrix.chr22.txt'
# gene_annot_RDS <- '../../data/intermediate/annotations/gene_annotation/gencode.v26lift37.annotation.gene.only.parsed.RDS'
# snp_annot_RDS <- '../../data/intermediate/annotations/snp_annotation/GTEx_V7_snp_annot_final_allchrs.chr22.RDS'
# train_test_fold_ids_RDS <- '../../data/intermediate/annotations/snp_annotation/Whole_Blood_train_test_fold_ids.RDS'
# all_cv_fold_ids_RDS <- '../../data/intermediate/annotations/snp_annotation/Whole_Blood_cv_fold_ids.RDS'
# 
# n_folds<-10
# n_train_test_folds<-5
# alpha<-0.5
# #null_testing<-FALSE
# cis_window<-1e6
# chrom<-22
# prefix <- 'GTEx_wb_test' %&% '_nested_cv'
# seed <- 1993
# out_dir <- '/storage/home/jmh791/group/default/private/jordan/ABC_TWAS/Enet_joblogs/wb_1mb'
# 
# 
# # main <- function(snp_annot_RDS, gene_annot_RDS, geno_file, expression_RDS, 
# #                  chrom, prefix, n_folds, n_train_test_folds,
# #                  seed=NA, cis_window, alpha, train_test_fold_ids_RDS, all_cv_fold_ids_RDS)
# 
# main(snp_annot_RDS, gene_annot_RDS, geno_file, expression_RDS, chrom, prefix, n_folds, n_train_test_folds, 
#      seed, cis_window, alpha, train_test_fold_ids_RDS, all_cv_fold_ids_RDS, out_dir)





