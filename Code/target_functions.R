## R Library for the validation of Splice-N-of-1-pathways
## AG Schissler
## Created 26 Jul 2017

############################################################
## i. Code development objects (DO NOT RUN)

## load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")

############################################################
## 1. Retrieve relevant pathways

tmp_data <- scores_list[[1]]

get_cancer_pathways <- function(scores_list) {
    if (length(scores_list) > 0){
        ## 1. retrieve cancer pathways
        tmp_data <- scores_list[[1]]
        cancer_pathways <- tmp_data[grep("cancer", tmp_data$pathway_desc), ]
        cancer_ids <- rownames(cancer_pathways)
        cancer_pathways <- cancer_pathways$pathway_desc
        names(cancer_pathways) <- cancer_ids
    } else stop("0 patients pathway data entered.")
    return(cancer_pathways)
}

############################################################
## 2. Retrieve relevant pathways and mark target(s) by id

cancer_pathways <- get_cancer_pathways(scores_list)
target <- names(cancer_pathways)[grep("Bladder", cancer_pathways)]
target <- names(cancer_pathways)

mark_target_list <- function(scores_list, target) {
    if (length(scores_list) > 0){
        ## Add targets 
        to_return <-lapply(scores_list, function(tmp_data){
            tmp_data$target <- 0
            tmp_data[target, "target"] <- 1
            tmp_data
        })
    } else stop("0 patients pathway data entered.")
    return(to_return)
}

target <- names(cancer_pathways)[grep("Bladder", cancer_pathways)]
target_scores_list <- mark_target_list(scores_list, target = target)

cancer_scores_list <- mark_target_list(scores_list, target = names(cancer_pathways))

############################################################
## 3. Summarize target characteristics 

tmp_data <- scores_list[["TCGA-BT-A20R"]]
str(tmp_data)
one_sided = T
one_sided = F
fdr = 0.2
cancer_ids <- names(cancer_pathways)

get_target_info <- function(scores_list, target, cancer_ids, fdr = 0.2, one_sided = T) {
    if (length(scores_list) > 0){
        ## one patient at a time
        target_list <- lapply(scores_list, function(tmp_data){
            ## summarize all pathways
            if (one_sided) {
                ## 1. overall summary
                num_scored <- sum(!is.na(tmp_data$pathway_score))
                num_hits <- sum(tmp_data$fdr_value <= fdr & tmp_data$pathway_score > 1)
                ## 2. summarize target pathway
                target_data <- tmp_data[target, ]
                target_capture <- target_data$fdr_value <= fdr & target_data$pathway_score > 1
                target_rank <- which(rownames(tmp_data) == target)
                ## 3. summarize cancer pathways
                cancer_data <- tmp_data[cancer_ids, ]
                cancer_capture <- any(cancer_data$fdr_value <= fdr & cancer_data$pathway_score > 1)
                ## 4. return summarize as a data frame
                to_return <- data.frame(num_scored, num_hits, cancer_capture, target_capture, target_rank)
            } else {
                num_scored <- sum(!is.na(tmp_data$pathway_score))
                num_hits <- sum(tmp_data$fdr_value <= fdr)
                ## 2. summarize target pathway
                target_data <- tmp_data[target, ]
                target_capture <- target_data$fdr_value <= fdr
                ## rank now based on fdr
                target_rank <- which(rownames(tmp_data[order(tmp_data$fdr_value), ]) == target)
                ## 3. summarize cancer pathways
                cancer_data <- tmp_data[cancer_ids, ]
                cancer_capture <- any(cancer_data$fdr_value <= fdr)
                ## 4. return summarize as a data frame
                to_return <- data.frame(num_scored, num_hits, cancer_capture, target_capture, target_rank)
            }
            return(to_return)
        })
        ## aggregate
        target_data <- do.call("rbind", target_list)
    } else stop("0 patients pathway data entered.")
    return(target_data)
}

target <- names(cancer_pathways)[grep("Bladder", cancer_pathways)]
one_sided = T
fdr = 0.2
cancer_ids <- names(cancer_pathways)


## one_sided
one_sided <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.05, one_sided = T)
sum(one_sided$cancer_capture)/nrow(one_sided)
sum(one_sided$target_capture)/nrow(one_sided)

one_sided <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.1, one_sided = T)
sum(one_sided$cancer_capture)/nrow(one_sided)
sum(one_sided$target_capture)/nrow(one_sided)

one_sided <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = T)
sum(one_sided$cancer_capture)/nrow(one_sided)
sum(one_sided$target_capture)/nrow(one_sided)

one_sided <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.5, one_sided = T)
sum(one_sided$cancer_capture)/nrow(one_sided)
sum(one_sided$target_capture)/nrow(one_sided)

## two sided
two_sided <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.05, one_sided = F)
sum(two_sided$cancer_capture)/nrow(two_sided)
sum(two_sided$target_capture)/nrow(two_sided)

two_sided <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.1, one_sided = F)
sum(two_sided$cancer_capture)/nrow(two_sided)
sum(two_sided$target_capture)/nrow(two_sided)

two_sided <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.2, one_sided = F)
sum(two_sided$cancer_capture)/nrow(two_sided)
sum(two_sided$target_capture)/nrow(two_sided)

two_sided <- get_target_info(scores_list, target = target, cancer_ids = cancer_ids, fdr = 0.5, one_sided = F)
sum(two_sided$cancer_capture)/nrow(two_sided)
sum(two_sided$target_capture)/nrow(two_sided)

