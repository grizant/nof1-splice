## Validate TCGA all patients through aggregation of target pathway information
## AG Schissler
## Created 22 Feb 2019
## Last modified 

##############################################################################
#### Source validation functions

source("~/Dropbox/Splice-n-of-1-pathways/Code/target_functions.R")
source("~/Dropbox/Splice-n-of-1-pathways/Code/surv_functions.R")

set.seed(44)

##load and aggregate manually
all_data <- NULL

## 1. BLCA
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## remove filtered pathways pathways
min_genes <- 15
max_genes <- 500
scores_list <- lapply(scores_list, FUN = function(tmp_data){
    tmp_data <- tmp_data[tmp_data$num_genes_annot >= min_genes & tmp_data$num_genes_annot <= max_genes,]
})
## explore gene set bias
## tmp_data <- scores_list[[sample(length(scores_list), 1)]]
## tmp_data <- tmp_data[!is.na(tmp_data$pathway_score),]
## p0 <- ggplot(tmp_data, aes(x = num_genes_annot, y = pathway_score, color = diff_splice_call))
## p0 + geom_point()
## diff_size <- tmp_data$num_genes_measured[tmp_data$diff_splice_call == 1]
## rest_size <- tmp_data$num_genes_measured[tmp_data$diff_splice_call == 0]
## summary(diff_size)
## summary(rest_size)
(cancer_pathways <- get_cancer_pathways(scores_list))
(cancer_ids <- names(cancer_pathways))
(target <- names(cancer_pathways)[grep("Bladder", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## null pathway assessment (remove filtered pathways)
effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)
convert_to_list <- function(x) {split(x, rep(1:ncol(x), each = nrow(x)))}
## X = convert_to_list(effect_mat)[[1]]
## Y = convert_to_list(fdr_mat)[[1]]
fdr_threshold <- 0.2
call_mat <- fil_index <- mapply(function(X, Y){
    as.numeric(X > 1 & Y < fdr_threshold) ## Odds ratio > 1 and fdr > 20%
}, X = convert_to_list(effect_mat), Y = convert_to_list(fdr_mat), SIMPLIFY = T)
## null target pathway
## null 9 pathways
reps = 2000
rand_mat <- replicate(n = reps, apply(call_mat, 1, function(tmp_pat){
    unname(sample(tmp_pat,1))
}))
## summary
rand_capture <- apply(rand_mat, 2, function(tmp_col) {sum(tmp_col)/length(tmp_col)})
rand_data <- round(quantile(rand_capture, c(0.5,0.9))*100, 2)
if (rand_data[1] == 0) {rand_data[1] <- 1/(reps + 1)}
## at least one of nine? model as binominal even though not independent
(rand_cancer <- round(100*(1 - dbinom(x = 0, size = length(cancer_pathways), prob = rand_data[1]/100)), 2))
## THE ABOVE CALCULATION IS NOT REPORTED!
## summarize target capture observations
tmp_target <- summarize_target_data(target_data)
## find empirical p-value for target pathway
tmp_emp1 <- sum(tmp_target$target_capture_rate/100 <= rand_capture)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
## find empirical p-value for cancer pathways
rand9_capture <- replicate(n = reps, expr =  {
    tmp_cols <- sample(1:ncol(rand_mat), 9)
    sum(apply(rand_mat[,tmp_cols], 1, function(tmp_pat){
        any(tmp_pat == 1)
    }))/nrow(rand_mat)
})
tmp_emp9 <- sum(tmp_target$cancer_capture_rate/100 <= rand9_capture)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "BLCA"
(all_data <- rbind(all_data, to_return))

################################################
## 1b. BLCA using reviewer's suggested technique
## shuffle patient-specific values
## count how many null pathways are greater than the target_capture
tmp_count1 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    rand_call1 <- rand_mat[,sample(ncol(rand_mat), 1)]
    rand_capture1 <- sum(rand_call1)/length(rand_call1)
    as.numeric(tmp_target$target_capture_rate/100 <= rand_capture1)
})
tmp_emp1 <- sum(tmp_count1)/reps
tmp_emp1
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
## count how often the null capture9 rate is greater than the obsevered cancer_capture
tmp_rand9 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    ## select 9 pathways at random
    rand_call9 <- rand_mat[,sample(ncol(rand_mat), 9)]
    ## how many patients found at least 1
    rand_capture9 <- sum(apply(rand_call9, 1, function(tmp_pat) {sum(tmp_pat) > 0 }))/nrow(rand_call9)
    as.numeric(tmp_target$cancer_capture_rate/100 <= rand_capture9)
})
tmp_emp9 <- sum(tmp_rand9)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
tmp_emp9
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "BLCA"
(all_data <- rbind(all_data, to_return))

################################################
## 2. THCA
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_THCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## remove filtered pathways pathways
min_genes <- 15
max_genes <- 500
scores_list <- lapply(scores_list, FUN = function(tmp_data){
    tmp_data <- tmp_data[tmp_data$num_genes_annot >= min_genes & tmp_data$num_genes_annot <= max_genes,]
})
(target <- names(cancer_pathways)[grep("Thyroid", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## null pathway assessment
effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)
convert_to_list <- function(x) {split(x, rep(1:ncol(x), each = nrow(x)))}
## X = convert_to_list(effect_mat)[[1]]
## Y = convert_to_list(fdr_mat)[[1]]
fdr_threshold <- 0.2
call_mat <- fil_index <- mapply(function(X, Y){
    as.numeric(X > 1 & Y < fdr_threshold) ## Odds ratio > 1 and fdr > 20%
}, X = convert_to_list(effect_mat), Y = convert_to_list(fdr_mat), SIMPLIFY = T)
## null target pathway
## null 9 pathways
rand_mat <- replicate(n = reps, apply(call_mat, 1, function(tmp_pat){
    unname(sample(tmp_pat,1))
}))
## summary
rand_capture <- apply(rand_mat, 2, function(tmp_col) {sum(tmp_col)/length(tmp_col)})
rand_data <- round(quantile(rand_capture, c(0.5,0.9))*100, 2)
if (rand_data[1] == 0) {rand_data[1] <- 1/(reps + 1)}
## at least one of nine? model as binominal even though not independent
(rand_cancer <- round(100*(1 - dbinom(x = 0, size = length(cancer_pathways), prob = rand_data[1]/100)), 2))
## summarize target capture observations
tmp_target <- summarize_target_data(target_data)
## find empirical p-value for target pathway
tmp_emp1 <- sum(tmp_target$target_capture_rate/100 <= rand_capture)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
## find empirical p-value for cancer pathways
rand9_capture <- replicate(n = reps, expr =  {
    tmp_cols <- sample(1:ncol(rand_mat), 9)
    sum(apply(rand_mat[,tmp_cols], 1, function(tmp_pat){
        any(tmp_pat == 1)
    }))/nrow(rand_mat)
})
tmp_emp9 <- sum(tmp_target$cancer_capture_rate/100 <= rand9_capture)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "THCA"
(all_data <- rbind(all_data, to_return))

################################################
## 2b. THCA using reviewer's suggested technique
## shuffle patient-specific values
## count how many null pathways are greater than the target_capture
tmp_count1 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    rand_call1 <- rand_mat[,sample(ncol(rand_mat), 1)]
    rand_capture1 <- sum(rand_call1)/length(rand_call1)
    as.numeric(tmp_target$target_capture_rate/100 <= rand_capture1)
})
tmp_emp1 <- sum(tmp_count1)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
tmp_emp1
## count how often the null capture9 rate is greater than the obsevered cancer_capture
tmp_rand9 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    ## select 9 pathways at random
    rand_call9 <- rand_mat[,sample(ncol(rand_mat), 9)]
    ## how many patients found at least 1
    rand_capture9 <- sum(apply(rand_call9, 1, function(tmp_pat) {sum(tmp_pat) > 0 }))/nrow(rand_call9)
    as.numeric(tmp_target$cancer_capture_rate/100 <= rand_capture9)
})
tmp_emp9 <- sum(tmp_rand9)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
tmp_emp9
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "THCA"
(all_data <- rbind(all_data, to_return))

## 3. UCEC
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_UCEC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## remove filtered pathways pathways
min_genes <- 15
max_genes <- 500
scores_list <- lapply(scores_list, FUN = function(tmp_data){
    tmp_data <- tmp_data[tmp_data$num_genes_annot >= min_genes & tmp_data$num_genes_annot <= max_genes,]
})
(target <- names(cancer_pathways)[grep("Endo", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## null pathway assessment
effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)
convert_to_list <- function(x) {split(x, rep(1:ncol(x), each = nrow(x)))}
## X = convert_to_list(effect_mat)[[1]]
## Y = convert_to_list(fdr_mat)[[1]]
fdr_threshold <- 0.2
call_mat <- fil_index <- mapply(function(X, Y){
    as.numeric(X > 1 & Y < fdr_threshold) ## Odds ratio > 1 and fdr > 20%
}, X = convert_to_list(effect_mat), Y = convert_to_list(fdr_mat), SIMPLIFY = T)
## null target pathway
## null 9 pathways
rand_mat <- replicate(n = reps, apply(call_mat, 1, function(tmp_pat){
    unname(sample(tmp_pat,1))
}))
## summary
rand_capture <- apply(rand_mat, 2, function(tmp_col) {sum(tmp_col)/length(tmp_col)})
rand_data <- round(quantile(rand_capture, c(0.5,0.9))*100, 2)
if (rand_data[1] == 0) {rand_data[1] <- 1/(reps + 1)}
## at least one of nine? model as binominal even though not independent
(rand_cancer <- round(100*(1 - dbinom(x = 0, size = length(cancer_pathways), prob = rand_data[1]/100)), 2))
## summarize target capture observations
tmp_target <- summarize_target_data(target_data)
## find empirical p-value for target pathway
tmp_emp1 <- sum(tmp_target$target_capture_rate/100 <= rand_capture)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
## find empirical p-value for cancer pathways
rand9_capture <- replicate(n = reps, expr =  {
    tmp_cols <- sample(1:ncol(rand_mat), 9)
    sum(apply(rand_mat[,tmp_cols], 1, function(tmp_pat){
        any(tmp_pat == 1)
    }))/nrow(rand_mat)
})
tmp_emp9 <- sum(tmp_target$cancer_capture_rate/100 <= rand9_capture)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "UCEC"
(all_data <- rbind(all_data, to_return))

################################################
## 3b. UCEC using reviewer's suggested technique
## shuffle patient-specific values
## count how many null pathways are greater than the target_capture
tmp_count1 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    rand_call1 <- rand_mat[,sample(ncol(rand_mat), 1)]
    rand_capture1 <- sum(rand_call1)/length(rand_call1)
    as.numeric(tmp_target$target_capture_rate/100 <= rand_capture1)
})
tmp_emp1 <- sum(tmp_count1)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
tmp_emp1
## count how often the null capture9 rate is greater than the obsevered cancer_capture
tmp_rand9 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    ## select 9 pathways at random
    rand_call9 <- rand_mat[,sample(ncol(rand_mat), 9)]
    ## how many patients found at least 1
    rand_capture9 <- sum(apply(rand_call9, 1, function(tmp_pat) {sum(tmp_pat) > 0 }))/nrow(rand_call9)
    as.numeric(tmp_target$cancer_capture_rate/100 <= rand_capture9)
})
tmp_emp9 <- sum(tmp_rand9)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
tmp_emp9
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "UCEC"
(all_data <- rbind(all_data, to_return))

## 4. PRAD
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_PRAD_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## remove filtered pathways pathways
min_genes <- 15
max_genes <- 500
scores_list <- lapply(scores_list, FUN = function(tmp_data){
    tmp_data <- tmp_data[tmp_data$num_genes_annot >= min_genes & tmp_data$num_genes_annot <= max_genes,]
})
(target <- names(cancer_pathways)[grep("Prostate", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## null pathway assessment
effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)
convert_to_list <- function(x) {split(x, rep(1:ncol(x), each = nrow(x)))}
## X = convert_to_list(effect_mat)[[1]]
## Y = convert_to_list(fdr_mat)[[1]]
fdr_threshold <- 0.2
call_mat <- fil_index <- mapply(function(X, Y){
    as.numeric(X > 1 & Y < fdr_threshold) ## Odds ratio > 1 and fdr > 20%
}, X = convert_to_list(effect_mat), Y = convert_to_list(fdr_mat), SIMPLIFY = T)
## null target pathway
## null 9 pathways
rand_mat <- replicate(n = reps, apply(call_mat, 1, function(tmp_pat){
    unname(sample(tmp_pat,1))
}))
## summary
rand_capture <- apply(rand_mat, 2, function(tmp_col) {sum(tmp_col)/length(tmp_col)})
rand_data <- round(quantile(rand_capture, c(0.5,0.9))*100, 2)
if (rand_data[1] == 0) {rand_data[1] <- 1/(reps + 1)}
## at least one of nine? model as binominal even though not independent
(rand_cancer <- round(100*(1 - dbinom(x = 0, size = length(cancer_pathways), prob = rand_data[1]/100)), 2))
## summarize target capture observations
tmp_target <- summarize_target_data(target_data)
## find empirical p-value for target pathway
tmp_emp1 <- sum(tmp_target$target_capture_rate/100 <= rand_capture)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
## find empirical p-value for cancer pathways
rand9_capture <- replicate(n = reps, expr =  {
    tmp_cols <- sample(1:ncol(rand_mat), 9)
    sum(apply(rand_mat[,tmp_cols], 1, function(tmp_pat){
        any(tmp_pat == 1)
    }))/nrow(rand_mat)
})
tmp_emp9 <- sum(tmp_target$cancer_capture_rate/100 <= rand9_capture)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "PRAD"
(all_data <- rbind(all_data, to_return))

################################################
## 4b. PRAD using reviewer's suggested technique
## shuffle patient-specific values
## count how many null pathways are greater than the target_capture
tmp_count1 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    rand_call1 <- rand_mat[,sample(ncol(rand_mat), 1)]
    rand_capture1 <- sum(rand_call1)/length(rand_call1)
    as.numeric(tmp_target$target_capture_rate/100 <= rand_capture1)
})
tmp_emp1 <- sum(tmp_count1)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
tmp_emp1
## count how often the null capture9 rate is greater than the obsevered cancer_capture
tmp_rand9 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    ## select 9 pathways at random
    rand_call9 <- rand_mat[,sample(ncol(rand_mat), 9)]
    ## how many patients found at least 1
    rand_capture9 <- sum(apply(rand_call9, 1, function(tmp_pat) {sum(tmp_pat) > 0 }))/nrow(rand_call9)
    as.numeric(tmp_target$cancer_capture_rate/100 <= rand_capture9)
})
tmp_emp9 <- sum(tmp_rand9)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
tmp_emp9
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "PRAD"
(all_data <- rbind(all_data, to_return))

## 5. LUSC
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## remove filtered pathways pathways
min_genes <- 15
max_genes <- 500
scores_list <- lapply(scores_list, FUN = function(tmp_data){
    tmp_data <- tmp_data[tmp_data$num_genes_annot >= min_genes & tmp_data$num_genes_annot <= max_genes,]
})
(target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## null pathway assessment
effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)
convert_to_list <- function(x) {split(x, rep(1:ncol(x), each = nrow(x)))}
## X = convert_to_list(effect_mat)[[1]]
## Y = convert_to_list(fdr_mat)[[1]]
fdr_threshold <- 0.2
call_mat <- fil_index <- mapply(function(X, Y){
    as.numeric(X > 1 & Y < fdr_threshold) ## Odds ratio > 1 and fdr > 20%
}, X = convert_to_list(effect_mat), Y = convert_to_list(fdr_mat), SIMPLIFY = T)
## null target pathway
## null 9 pathways
rand_mat <- replicate(n = reps, apply(call_mat, 1, function(tmp_pat){
    unname(sample(tmp_pat,1))
}))
## summary
rand_capture <- apply(rand_mat, 2, function(tmp_col) {sum(tmp_col)/length(tmp_col)})
rand_data <- round(quantile(rand_capture, c(0.5,0.9))*100, 2)
if (rand_data[1] == 0) {rand_data[1] <- 1/(reps + 1)}
## at least one of nine? model as binominal even though not independent
(rand_cancer <- round(100*(1 - dbinom(x = 0, size = length(cancer_pathways), prob = rand_data[1]/100)), 2))
## summarize target capture observations
tmp_target <- summarize_target_data(target_data)
## find empirical p-value for target pathway
tmp_emp1 <- sum(tmp_target$target_capture_rate/100 <= rand_capture)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
## find empirical p-value for cancer pathways
rand9_capture <- replicate(n = reps, expr =  {
    tmp_cols <- sample(1:ncol(rand_mat), 9)
    sum(apply(rand_mat[,tmp_cols], 1, function(tmp_pat){
        any(tmp_pat == 1)
    }))/nrow(rand_mat)
})
tmp_emp9 <- sum(tmp_target$cancer_capture_rate/100 <= rand9_capture)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "LUSC"
(all_data <- rbind(all_data, to_return))

################################################
## 5b. LUSC using reviewer's suggested technique
## shuffle patient-specific values
## count how many null pathways are greater than the target_capture
tmp_count1 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    rand_call1 <- rand_mat[,sample(ncol(rand_mat), 1)]
    rand_capture1 <- sum(rand_call1)/length(rand_call1)
    as.numeric(tmp_target$target_capture_rate/100 <= rand_capture1)
})
tmp_emp1 <- sum(tmp_count1)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
tmp_emp1
## count how often the null capture9 rate is greater than the obsevered cancer_capture
tmp_rand9 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    ## select 9 pathways at random
    rand_call9 <- rand_mat[,sample(ncol(rand_mat), 9)]
    ## how many patients found at least 1
    rand_capture9 <- sum(apply(rand_call9, 1, function(tmp_pat) {sum(tmp_pat) > 0 }))/nrow(rand_call9)
    as.numeric(tmp_target$cancer_capture_rate/100 <= rand_capture9)
})
tmp_emp9 <- sum(tmp_rand9)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
tmp_emp9
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "LUSC"
(all_data <- rbind(all_data, to_return))

## ## 5.2 LUSC EEv2
## load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_EEv2_KEGG_29july2017.RData")
## (target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
## target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## (to_return <- summarize_target_data(target_data))
## rownames(to_return) <- "LUSCv2"
## (all_data <- rbind(all_data, to_return))
## 

## 6. LUAD
load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUAD_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## remove filtered pathways pathways
min_genes <- 15
max_genes <- 500
scores_list <- lapply(scores_list, FUN = function(tmp_data){
    tmp_data <- tmp_data[tmp_data$num_genes_annot >= min_genes & tmp_data$num_genes_annot <= max_genes,]
})
(target <- names(cancer_pathways)[grep("Non-small cell lung", cancer_pathways)])
target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## null pathway assessment
effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)
convert_to_list <- function(x) {split(x, rep(1:ncol(x), each = nrow(x)))}
## X = convert_to_list(effect_mat)[[1]]
## Y = convert_to_list(fdr_mat)[[1]]
fdr_threshold <- 0.2
call_mat <- fil_index <- mapply(function(X, Y){
    as.numeric(X > 1 & Y < fdr_threshold) ## Odds ratio > 1 and fdr > 20%
}, X = convert_to_list(effect_mat), Y = convert_to_list(fdr_mat), SIMPLIFY = T)
## null target pathway
## null 9 pathways
rand_mat <- replicate(n = reps, apply(call_mat, 1, function(tmp_pat){
    unname(sample(tmp_pat,1))
}))
## summary
rand_capture <- apply(rand_mat, 2, function(tmp_col) {sum(tmp_col)/length(tmp_col)})
rand_data <- round(quantile(rand_capture, c(0.5,0.9))*100, 2)
if (rand_data[1] == 0) {rand_data[1] <- 1/(reps + 1)}
## at least one of nine? model as binominal even though not independent
(rand_cancer <- round(100*(1 - dbinom(x = 0, size = length(cancer_pathways), prob = rand_data[1]/100)), 2))
## summarize target capture observations
tmp_target <- summarize_target_data(target_data)
## find empirical p-value for target pathway
tmp_emp1 <- sum(tmp_target$target_capture_rate/100 <= rand_capture)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
## find empirical p-value for cancer pathways
rand9_capture <- replicate(n = reps, expr =  {
    tmp_cols <- sample(1:ncol(rand_mat), 9)
    sum(apply(rand_mat[,tmp_cols], 1, function(tmp_pat){
        any(tmp_pat == 1)
    }))/nrow(rand_mat)
})
tmp_emp9 <- sum(tmp_target$cancer_capture_rate/100 <= rand9_capture)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "LUAD"
(all_data <- rbind(all_data, to_return))

################################################
## 6b. LUAD using reviewer's suggested technique
## shuffle patient-specific values
## count how many null pathways are greater than the target_capture
tmp_count1 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    rand_call1 <- rand_mat[,sample(ncol(rand_mat), 1)]
    rand_capture1 <- sum(rand_call1)/length(rand_call1)
    as.numeric(tmp_target$target_capture_rate/100 <= rand_capture1)
})
tmp_emp1 <- sum(tmp_count1)/reps
if (tmp_emp1 == 0) {tmp_emp1 <- 1/(reps+1)}
tmp_emp1
## count how often the null capture9 rate is greater than the obsevered cancer_capture
tmp_rand9 <- replicate(n = reps, expr = {
    rand_mat <- t(apply(call_mat, MARGIN = 1, FUN = sample))
    ## select 9 pathways at random
    rand_call9 <- rand_mat[,sample(ncol(rand_mat), 9)]
    ## how many patients found at least 1
    rand_capture9 <- sum(apply(rand_call9, 1, function(tmp_pat) {sum(tmp_pat) > 0 }))/nrow(rand_call9)
    as.numeric(tmp_target$cancer_capture_rate/100 <= rand_capture9)
})
tmp_emp9 <- sum(tmp_rand9)/reps
if (tmp_emp9 == 0) {tmp_emp9 <- 1/(reps+1)}
tmp_emp9
## format and combine
(to_return <- cbind(tmp_target, rand50 = rand_data[1], rand90 = rand_data[2], emp1 = tmp_emp1, rand9 = rand_cancer, emp9 = tmp_emp9))
rownames(to_return) <- "LUAD"
(all_data <- rbind(all_data, to_return))

## Table 3
## all_data[order(all_data$target_capture_rate),]
## all_data[order(all_data$emp1),]
all_data

## ## READ (drop as this target is not a perfect match)
## load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_READ_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## (target <- names(cancer_pathways)[grep("Colorectal", cancer_pathways)])
## target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## (to_return <- summarize_target_data(target_data))
## rownames(to_return) <- "READ"
## (all_data <- rbind(all_data, to_return))
## 
## ## COAD
## load("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_COAD_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## (target <- names(cancer_pathways)[grep("Colorectal", cancer_pathways)])
## target_data <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
## (to_return <- summarize_target_data(target_data))
## rownames(to_return) <- "COAD"
## (all_data <- rbind(all_data, to_return))
## 

##### Table 1

## get clinical data

## LUSC
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## LUAD
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUAD_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## PRAD
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_PRAD_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## THCA
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_THCA_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## UCEC
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_UCEC_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

## BLCA
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_paired_clinical.RData")
nrow(clin_data)
table(clin_data$vitalstatus)

