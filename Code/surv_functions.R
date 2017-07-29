## R Library for the validation of Splice-N-of-1-pathways
## AG Schissler
## Created 26 Jul 2017

############################################################
## i. Code development objects (DO NOT RUN)

## load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_paired_clinical.RData")
## load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_paired_clinical.RData")
## ## load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## patients <- names(scores_list)
## clin_data <- clin_data[rownames(clin_data) %in% patients,]
## 

############################################################
## 1. Aggregate pathway scores and fdr

## type = "pathway_score"
## type = "fdr_value"
## tmp_data = scores_list[[1]]
## remove_missing = T
 
compile_scores <- function(scores_list, type = c("pathway_score", "fdr_value"), remove_missing = T) {
    if (length(scores_list) > 0) {
        tmp_index <- which(names(scores_list[[1]]) == type)
    } else stop("Scores list is empty")
    value_list <- lapply(scores_list, function(tmp_data){
        ## retrieve score
        tmp_value <- tmp_data[, tmp_index]
        ## add pathway names
        names(tmp_value) <- rownames(tmp_data)
        ## reorder to allow for aggregation
        tmp_value[order(names(tmp_value))]
    })
    ## aggregate into a matrix
    value_mat <- do.call("rbind", value_list)
    ## remove NAs if desired
    if (remove_missing) {
        missing_logic <- apply(value_mat, 2, function(tmp_pathway) {
            any(is.na(tmp_pathway))
        })
        ## table(missing_logic) 
        value_mat <- value_mat[, !missing_logic]
        print(paste("Removed", sum(missing_logic), "unscored pathways out of", length(missing_logic)))
    }
    return(value_mat)
}

## str(compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = T))

############################################################
## 1.2 Filter to only pathways that are significant

## must not pre filter NA!

filter_value_mat <- function(effect_mat, fdr_mat, type = c("pathway_score", "fdr_value"), fdr_threshold = 0.2, one_sided = T, remove_missing = T) {
    ## 1. Find pathways enriched for at least one patient
    ## tmp_val <- fdr_mat[,1]
    ## X <- effect_mat[,1]
    ## Y <- fdr_mat[,1]
    if (one_sided) {
        ## only count enriched pathways
        ## convert to lists for mapply
        convert_to_list <- function(x) {split(x, rep(1:ncol(x), each = nrow(x)))}
        fil_index <- unlist(mapply(function(X, Y){
            any(X > 1 & Y < fdr_threshold) ## Odds ratio > 1 and fdr > 20%
        }, X = convert_to_list(effect_mat), Y = convert_to_list(fdr_mat), SIMPLIFY = F))
        ## sum(fil_index)
        names(fil_index) <- colnames(effect_mat)
        ##tmp_id <- names(fil_index[fil_index][1])
        ## any(effect_mat[,tmp_id] > 1)
        ## sort(fdr_mat[,tmp_id])
        ## any(fdr_mat[,tmp_id] < fdr_threshold)
    } else {
        ## count both under and over enriched pathways
        fil_index <- apply(fdr_mat, 2, function(tmp_val){
            any(tmp_val < fdr_threshold)
        })
    }
    ## 2. return data type specified
    if (type == "pathway_score") {
        value_mat <- effect_mat[,fil_index]
        ## filter missing pathways to pass to cluster functions
        if (remove_missing) {
            missing_logic <- apply(value_mat, 2, function(tmp_pathway) {
                any(is.na(tmp_pathway))
            })
        }
        ## table(missing_logic) 
        value_mat <- value_mat[, !missing_logic]
        print(paste("Removed", sum(missing_logic), "unscored pathways out of", length(missing_logic), "pathways that met thresholds"))
        ## return fdr values
    } else {
        value_mat <- fdr_mat[,fil_index]
    }
    return(value_mat)
}

## effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
## fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)
## fdr_threshold = 0.2
## one_sided = T
## type = "pathway_score"
## type = "fdr_value"

## fil_effect_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "pathway_score")
## fil_fdr_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "fdr_value")
## str(fil_fdr_mat)

############################################################
## 2. Produce survival objects

produce_clin_surv <- function(clin_data) {
    clin_data$status <- ifelse(clin_data$vitalstatus == "Alive", 0, 1)
    ## table(clin_data$status)
    clin_data$time <- clin_data$daystolastfollowup
    ## then insert the survival time
    clin_data$time[is.na(clin_data$time)] <- clin_data$daystodeath[!is.na(clin_data$daystodeath)]
    ## create the survival object
    clin_surv <- survival::Surv(time = clin_data$time, event = clin_data$status)
    ## head(clin_data[,c(2,3,4)], 10)
    ## head(clin_surv, 10)
    return(clin_surv)
}

## produce_clin_surv(clin_data = clin_data)

############################################################
## 3. Cluster patients based on pathway values (perhaps also on one pathway later)

## type = "sc"
## type = "pam"
## num_clusters = 2

cluster_pat <- function(value_mat, type = "pam", num_clusters = 2, min_nj = 4, max_clusters = 5) {
    if (type == "sc") {
        ## Spectral clustering (SC) using the eigengap heuristic
        ## only works for multiple pathways
        if ( !is.vector(value_mat) ) {
            ## 1. Configure settings for SC for matrix
            max_centers <- ceiling(nrow(value_mat) / min_nj)
            neighbors  <- ceiling(log(nrow(value_mat)) + 1)
            if (is.null(num_clusters)) {
                ## 2. Compute eigenvalues for max number of centers to use eignegap heuritic
                sc_max <- kknn::specClust(data = value_mat, centers = max_centers, nn = neighbors, method = "random-walk")
                ## ggplot2::qplot(factor(1:length(sc_max$eigenvalue)), sc_max$eigenvalue)
                ## 3. Find differences in eigenvalues
                eigen_diff <- diff(sc_max$eigenvalue)[-1]
                diff_data <- data.frame(m = 2:(max_centers-1), eigen_diff)
                ## find the largest gap with min_nj met (or take m = 1)
                diff_data <- diff_data[order(diff_data$eigen_diff, decreasing = T), ]
                ## add terminal clustering
                diff_data <- rbind(diff_data, c(m = 1, 0))
                ## 4. pick best clustering no more than the maximum and at least 2 clusters
                best_m <- diff_data[diff_data$m <= max_clusters & diff_data$m > 1, "m"][1]
            } else best_m <- num_clusters
            ## best_m <- 3
            ## 5. cluster with best number of clusters
            sc_best <- kknn::specClust(data = value_mat, centers = best_m, nn = neighbors, method = "random-walk")
            ## add patient names
            names(sc_best$cluster) <- rownames(value_mat)
            to_return <- sc_best$cluster
        } else stop("Spectral clustering must have more than one pathway score.")
    }
    if (type == "pam") {
        if (!is.null(num_clusters)) {
            pam_clusters <- cluster::pam(x = value_mat, k = num_clusters)
            to_return <- pam_clusters$clustering
        } else stop("kmeans-type clustering requires a number of clusters")
    }
    ## Return clusters for the patients as a vector
    return(to_return)
}

## min_nj <- 4
## max_clusters = 5
## value_mat = value_mat[,1]
## 
## (my_clusters <- cluster_pat(value_mat = value_mat, type = "sc"))
## 
## ## cluster on only significant pathways
## (fil_effect_clusters <- cluster_pat(value_mat = fil_effect_mat, type = "sc"))
## (fil_fdr_clusters <- cluster_pat(value_mat = fil_fdr_mat, type = "sc"))
## 
## effect_mat <- compile_scores(scores_list, type = "pathway_score", remove_missing = T)
## ## try clustering on one pathway (doesn't work)
## (one_effect_clusters <- cluster_pat(value_mat = effect_mat[,1], type = "pam", num_clusters = 2))

############################################################
## 4. Fit survival curve based on a clustering

## clusters = fil_effect_clusters
## clusters <- top_clust
## clusters = fil_fdr_clusters
## plot = F

fit_surv <- function(clusters, clin_data, plot = F) {
    ## 1. remove any patients that are not in the clinical data
    clusters <- clusters[names(clusters) %in% rownames(clin_data)]
    ## 2. produce survival object
    clin_surv <- produce_clin_surv(clin_data = clin_data)
    ## 3. add cluster to clinical and status
    clin_data$status <- ifelse(clin_data$vitalstatus == "Alive", 0, 1)
    clin_data$time <- clin_data$daystolastfollowup
    ## then insert the survival time
    clin_data$time[is.na(clin_data$time)] <- clin_data$daystodeath[!is.na(clin_data$daystodeath)]
    clin_data$cluster <- clusters
    ## 3. Optionally plot survival curves
    if (plot) {
        require(cowplot)
        sf_survfit <- survival::survfit(clin_surv ~ cluster, data = clin_data)
        sf_tidy = broom::tidy(sf_survfit)
        mx = max(sf_tidy$n.censor)
        p0 <- ggplot2::ggplot(sf_tidy, aes(time, estimate, fill = strata)) + 
            geom_line() +
            geom_point(aes(shape = as.factor(n.censor)), size = 3) + 
            scale_shape_manual(values=c(NA, 1:mx))+
            ## geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.25) + 
            xlab("days") + 
            ylab("Proportion Survival") +
            cowplot::background_grid(major = "xy", minor = "none") +
            ## theme_bw() +
            theme(legend.position = "none")
    }
    ## 4. Log-rank test
    surv_diff <- survival::survdiff(survival::Surv(time, status) ~ cluster, data = clin_data, rho = 0)
    ## 5. return plot also if desired
    if (plot) {
        to_return <- list(surv = surv_diff, plot = p0)
    } else {
        to_return <- surv_diff
    }
    return(to_return)
}

## fit_surv(clin_data = clin_data, clusters = fil_effect_clusters, plot = T)
## fit_surv(clin_data = clin_data, clusters = fil_fdr_clusters, plot = T)
## fit_surv(clin_data = clin_data, clusters = one_effect_clusters, plot = T)

############################################################
## 4.2 Extract p-value from survival object

get_surv_pvalue <- function(surv_diff) {
    num_clusters <- length(dimnames(surv_diff$n)$groups)
    df <- num_clusters - 1
    pchisq(surv_diff$chisq, df = df, lower.tail = F)
}

## (surv_diff <- fit_surv(clin_data = clin_data, clusters = fil_effect_clusters, plot = F))
## get_surv_pvalue(surv_diff)
## surv_diff <- fit_surv(clin_data = clin_data, clusters = one_effect_clusters, plot = F)

############################################################
## 5. Wrappers to fit multiple survival curves, one per pathway (may want to generalize later)

get_pathway_pvalue <- function(value_mat, clin_data, type, num_clusters, ...){
    ## cluster one pathway at a time
    pathway_clusters <- apply(value_mat, 2, cluster_pat, type = "pam", num_clusters = num_clusters)
    ## fit survival objects and store in a list
    pathway_surv <- apply(pathway_clusters, 2, fit_surv, clin_data = clin_data, plot = F)
    ## now find the p-values
    return(unlist(lapply(pathway_surv, get_surv_pvalue)))
}

############################################################
## 6. Comparison to standards

## helper function

## tmp_genes <- sample( unique(iso_data$geneSymbol),  size = num_genes )

clust_surv_pvalue <- function(tmp_genes, iso_data, clin_data, ...) {
    tmp_mat <- as.matrix(t(iso_data[iso_data$geneSymbol %in% tmp_genes, -(1:2)]))
    tmp_clust <- cluster_pat(tmp_mat)
    ## fit survival objects on observed data
    tmp_surv <- fit_surv(clusters = tmp_clust, clin_data, plot = F)
    ## tmp_surv$plot
    return(get_surv_pvalue(tmp_surv))
}

get_empirical_pvalue <- function(pathway_genes, iso_data, clin_data, reps = 2000, ...){
    ## subset to isoforms in gene set
    num_genes <- length(unique(iso_data$geneSymbol[(iso_data$geneSymbol %in% pathway_genes)]))
    ### get observed p-value
    obs_pvalue <- clust_surv_pvalue(pathway_genes, iso_data, clin_data)
    ## now sample of same size
    ## find random genes
    rand_genes_mat <- replicate(n = reps, sample( unique(iso_data$geneSymbol),  size = num_genes ))
    ## compute p-values
    ## system.time(rand_pvalue_mat <- apply(rand_genes, 2, clust_surv_pvalue, iso_data, clin_data)) ## 12 seconds
    rand_pvalue <- apply(rand_genes_mat, 2, clust_surv_pvalue, iso_data = iso_data, clin_data = clin_data) ## 12 seconds for 65 genes
    ## hist(rand_pvalue_mat)
    ## get proportion smaller
    empirical_pvalue <- sum(rand_pvalue <= obs_pvalue)/reps
    ## now find the p-values
    return(empirical_pvalue)
}


