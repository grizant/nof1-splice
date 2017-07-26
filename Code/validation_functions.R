## R Library for the validation of Splice-N-of-1-pathways
## AG Schissler
## Created 26 Jul 2017

############################################################
## i. Code development objects (DO NOT RUN)

load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_paired_clinical.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
patients <- names(scores_list)
clin_data <- clin_data[rownames(clin_data) %in% patients,]

############################################################
## 1. Aggregate pathway scores and fdr

type = "pathway_score"
type = "fdr_value"
tmp_data = scores_list[[1]]
remove_missing = T

## 1.1 compute Hellinger distance for isoforms of a gene given pair of measurements
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
## 2. Isoform-to-gene functions

## X is a genewise list
## remove_DEGs: perform a kMEn like clustering step to classify responsive transcripts

## transform isoform-level expression to gene level
transform_iso_gene <- function(X, method = c("delta", "hellinger"), remove_DEGs = F, ...) {
    
    ## restructure into a list of genes
    gene_list <- X
        
    ## remove DEGs (responsive transcripts) if desired via a clustering step
    if (remove_DEGs) {
        ## find log fold change at the gene level
        logFC_gene_vec <- unlist(lapply(gene_list, function(tmp_gene) {
            ## sum the isoform level to find the transcript-level expression
            log2(sum(tmp_gene[,2]) + 1) - log2(sum(tmp_gene[,1]) + 1)
        }))
        ## apply k-means to cluster (first step in kMEn) and remove DEG cluster
        tmp_cluster <- kmeans(x = abs(logFC_gene_vec), centers = 2, nstart = 100)
        ## figure out which cluster corresponds to non-DEGs
        clust_index <- which.min(tmp_cluster$centers)
        genes_to_keep <- names(tmp_cluster$cluster)[tmp_cluster$cluster == clust_index]
        gene_list <- gene_list[genes_to_keep]
    }
    
    ## transform via specified method
    if (method == "hellinger") {
        ## first compute the relative expression within a gene
        rel_list <- lapply(gene_list, function(tmp_gene) {
            apply(tmp_gene, 2, function(tmp_col) {
                ## check that there is no division by 0 
                if (sum(tmp_col) > 0) {
                    tmp_col/sum(tmp_col)
                } else tmp_col
            })  
        })
        ## then compute the Hellinger distance
        to_return <- unlist(lapply(rel_list, compute_hellinger))
        ## remove missing distances (due to non-expression in either sample)
        to_return <- to_return[!is.na(to_return)]
    }

    ## Yves method of maximum difference in expression
    if (method == "delta") {
        to_return <- unlist(lapply(gene_list, function(tmp_gene){
            ## find the largest difference 
            ## first calculate log fold change
            tmp_iso_fc <- log2(tmp_gene[,2] + 1) - log2(tmp_gene[,1] + 1)
            ## calculate and store the largest fold change
            tmp_iso_fc[which.max(abs(tmp_iso_fc))]
        }))
    }

    ## return the gene-level metrics
    return(to_return)
}

############################################################
## 3. Gene-to-pathway functions

## transform gene-level expression to pathway level
transform_gene_pathway <- function(gene_dist, annot_file, desc_file = NULL, method = c("EE", "avg", "fet"), genes_range = c(15, 500), ...) {
    ### 1. Structure gene set definitions
    ## read in gene set (pathway) annotation
    annot_data <- read.delim2(file = annot_file)
    ## restructure into a list of pathways
    annot_list <- split(annot_data, annot_data$path_id)
    ## store the number of genes annotated
    annot_lengths <- unlist(lapply(annot_list, function(tmp_set) {dim(tmp_set)[1]}))
    ## count genes actually measured
    measured_lengths <- unlist(lapply(annot_list, function(tmp_set) {
        tmp_genes <- as.character(tmp_set$symbol)
        sum(tmp_genes %in% names(gene_dist))
    }))

### 2 Score pathways via the desired method
    ## 2.1 Score by average
    if (method == "avg") {
        ###  2.1.1 Find the average distance for each pathway
        ## tmp_set <- annot_list[["GO:0000003"]]
        avg_dist <- unlist(lapply(annot_list, function(tmp_set) {
            tmp_genes <- as.character(tmp_set$symbol)
            measured_genes <- tmp_genes[tmp_genes %in% names(gene_dist)]
            ## impose genes range
            if (length(measured_genes) < genes_range[1] | length(measured_genes) > genes_range[2]) {
                to_return <- NA
            } else {
                to_return <- mean(gene_dist[measured_genes])
            }
            return(to_return)
        }))
        ## qplot(avg_dist)
        ## qplot(avg_iso_count, avg_dist)
        ## cor.test(avg_iso_count, avg_dist)
        ## center the data
        ## avg_dist <- avg_dist - mean(avg_dist)
        #### 2.1.2 Apply local FDR of Efron
        ## pdf(file = paste0("~/Dropbox/Splice-n-of-1-pathways/Preliminary_figures/Local_FDR_GO-BP_Hel_Avg_Iso30_expressed_pathwayfilter_TPM_patient_", tmp_pat,".pdf"))
        ## tmp_locfdr <- locfdr::locfdr(zz = avg_dist[!is.na(avg_dist)])
        ## dev.off()
        tmp_locfdr <- locfdr::locfdr(zz = avg_dist[!is.na(avg_dist)], plot = 0)
        #### 2.1.3 Find interesting pathways (more differentially splicing)
        ## first find the upper threshold
        tmp_upper <- tmp_locfdr$z.2[2]
        ## then indicate the hits the hits
        ## tmp_hits <- names(avg_dist)[avg_dist >= tmp_upper]
        tmp_call <- rep(0, length(avg_dist))
        tmp_call[avg_dist >= tmp_upper] <- 1
        ## 2.1.4 format the output
        to_return <- data.frame(pathway_score = avg_dist, direction = NA, fdr_value = NA, num_genes_annot = annot_lengths, num_genes_measured = measured_lengths, row.names = names(avg_dist), upper_fdr_threshold = tmp_upper, diff_splice_call = tmp_call, stringsAsFactors = F)
        ## update fdr_value
        to_return[which(!is.na(avg_dist)), "fdr_value"] <- tmp_locfdr$fdr
        summary(to_return$fdr_value)
        ## sort
        to_return <- to_return[order(to_return$pathway_score, decreasing = T),]
    }


### 2.2 FET after gene-wise local FDR
    if (method == "EE") {
        ##  2.2.1 Find alternatively spliced genes across whole transcriptome
        ## hist(gene_dist, breaks = 120)
        gene_clusters <- kmeans(x = gene_dist, 2)
        ## gene_clusters$centers
        ## remove dist = 1 (or near it) or 0 due to fitting difficulties.
        ## ... give 1 a positive call and 0 negative call manually later
        ##tolerance <- .Machine$double.eps^0.5
        ##dist_one <- which(abs(gene_dist - 1) < tolerance)
        ##dist_zero <- which(abs(gene_dist) < tolerance)
        ##remove_index <- c(dist_one, dist_zero)
        #### now filter
        ##transformed_dist <- gene_dist[-remove_index]
        #### transformed_dist <- ecdf(transformed_dist)(transformed_dist)
        ##transformed_dist <- qnorm(transformed_dist)
        ##summary(transformed_dist)
        #### replace Inf (corresponding to 1) with max + esp
        #### eps <- 0.001
        #### transformed_dist[is.infinite(transformed_dist)] <-
        #### max(transformed_dist[!is.infinite(transformed_dist)], na.rm = T) + eps
        #### now transform
        #### tail(sort(transformed_dist))
        ##hist(transformed_dist, 120)
        #### summary(transformed_dist)
        ##tmp_locfdr <- locfdr::locfdr(zz = transformed_dist, bre = 120, df = 7, pct = 0, nulltype = 1, plot = 4,
        ##                             mlests = c(mean(transformed_dist), sd(transformed_dist)))
        #### retrieve hits
        ##tmp_upper <- tmp_locfdr$z.2[2]
        ##alt_genes <- names(which(transformed_dist > tmp_upper))
        ##alt_genes <- c(alt_genes, names(dist_one))
        #### Apply calls to gene dist despite filtering.
        ## This may be anti-conservative/conservative depending on the ones and zeros.
        
        dist_data <- data.frame(gene_dist, call = 0)
        ## dist_data[transformed_dist > tmp_upper, "call"] <- 1
        ## determine which cluster has larger center
        alt_center <- which.max(gene_clusters$centers)
        dist_data[names(gene_clusters$cluster), "call"] <- ifelse(gene_clusters$cluster == alt_center, 1, 0)
        ## table(dist_data$call)
        ## qplot(x = gene_dist, fill = factor(call), data = dist_data, bins = 120)
        
        ## 2.2.2 enrichment
        ## tmp_set <- annot_list[[sample(1:length(annot_list), 1)]]
        odds_ratio <- unlist(lapply(annot_list, function(tmp_set) {
            tmp_genes <- as.character(tmp_set$symbol)
            measured_genes <- tmp_genes[tmp_genes %in% names(gene_dist)]
            ## genes outside pathway
            not_genes <- rownames(dist_data)[!(rownames(dist_data) %in% measured_genes)]
            ## impose genes range
            if (length(measured_genes) >= genes_range[1] & length(measured_genes) <= genes_range[2]) {
                ## calculate counts
                (x11 <- sum(dist_data[measured_genes, "call"]))
                (x21 <- sum(dist_data[measured_genes, "call"] == 0))
                (x12 <- sum(dist_data[not_genes, "call"]))
                (x22 <- sum(dist_data[not_genes, "call"] == 0))
                ## sum(x11, x21, x12, x22) == nrow(dist_data)
                ## compute FET
                ## tmp_fet <- fisher.test(x = matrix( c(x11, x21, x12, x22), nrow = 2, ncol = 2))
                ## to_return <- tmp_fet$estimate
                ## names(to_return) <- NULL
                ## compute odds ratio manually
                (to_return <-  (x11/x21)/(x12/x22))
            } else to_return <- NA
            return(to_return)
        }))
        
        ## hist(odds_ratio, breaks = 20)
        ## summary(odds_ratio)
        
        ## 2.2.3 local FDR at the pathway level, since we need to control here
        ## 1. first remove unscored pathways
        fil_odds <- odds_ratio[!is.na(odds_ratio)]
        ## 2. then log transform to handle extreme values and to normalize
        # (Chinn 2000, A simple method for converting an odds ratio to effect size for use in meta-analysis
        ## fil_odds <- log(fil_odds)/1.81 
        ## 3. filter any -Inf due to 0
        ## fil_odds <- fil_odds[!is.infinite(fil_odds)]
        ## detect outlying pathways above prior to fitting
        outliers <- boxplot.stats(fil_odds)$out
        if (length(outliers) > 0) fil_odds <- fil_odds[!(names(fil_odds) %in% names(outliers))]
        ## hist(fil_odds, breaks = 20)

        ## tmp_locfdr <- locfdr::locfdr(zz = fil_odds, bre = ceiling(length(fil_odds)/8), df = 4, pct = 0, pct0 = 1/4, nulltype = 2, plot = 1, mlests = c(mean(fil_odds), sd(fil_odds)))
        suppressWarnings(tmp_locfdr <- locfdr::locfdr(zz = fil_odds, bre = ceiling(length(fil_odds)/8), df = 4, pct = 0,
                                                      pct0 = 1/4, nulltype = 2, plot = 0, mlests = c(1, sd(fil_odds))))
        (tmp_upper <- tmp_locfdr$z.2[2])
        ## then indicate the hits the hits
        ## tmp_hits <- names(avg_dist)[avg_dist >= tmp_upper]
        tmp_call <- rep(0, length(odds_ratio))
        names(tmp_call) <- names(odds_ratio)
        if (!is.na(tmp_upper)) {
            odds_ratio[which(odds_ratio >= tmp_upper)]
            tmp_call[which(odds_ratio >= tmp_upper)] <- 1
        } else {
            if (length(outliers) > 0) {
                tmp_call[names(outliers)] <- 1
            }
        }
        ## table(tmp_call)

        ## 2.2.4 format the output
        to_return <- data.frame(pathway_score = odds_ratio, direction = NA, fdr_value = 1, num_genes_annot = annot_lengths, num_genes_measured = measured_lengths, row.names = names(odds_ratio), upper_fdr_threshold = tmp_upper, diff_splice_call = tmp_call, stringsAsFactors = F)

        ## str(tmp_locfdr)
        ## update fdr_value
        to_return[names(fil_odds), "fdr_value"] <- tmp_locfdr$fdr
        ## and give p-value of 0 for outliers
        if (length(outliers) > 0) {to_return[names(outliers), "fdr_value"] <- 0}
        ## head(sort(to_return$fdr_value), 50)
        ## sort
        to_return <- to_return[order(to_return$pathway_score, decreasing = T),]
    }
    
    
    ### add pathway descriptions
    to_return$pathway_desc <- "N/A"
    if (!is.null(desc_file)) {
        ## read in
        desc_data <- read.delim2(file = desc_file, stringsAsFactors = F)
        rownames(desc_data) <- desc_data$path_id
        ## check that scored pathways have annotation
        if (any(rownames(to_return) %in% desc_data$path_id)) {
            ## add to the output
            to_return$pathway_desc <- desc_data[rownames(to_return), "description"]
        } else {
            warning("Pathway ids in description data do not match ids in gene set definitions.")
        }
    }
    ## return the pathway-level metrics
    return(to_return)
}


## gene_dist <- example_hel
## gene_dist <- example_delta
## annot_file <- "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt"
## 
## system.time(example_avg <- transform_gene_pathway(gene_dist = example_hel, annot_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_description.txt", method = "avg")) ## 6.76 seconds
## 
## str(example_avg)
## head(example_avg)
## 
## 
## qplot(x = num_genes_measured, y = pathway_score, data = example_avg)

############################################################
## 4. Wrapper for the entire isoform -> gene -> pathway

## filter based on 

## transform gene-level expression to pathway level
transform_iso_pathway <- function(iso_data, annot_file, desc_file = NULL, pathway_method = c("EE", "avg", "fet"), gene_method = c("delta", "hellinger"), remove_DEGs = F, iso_range = c(2,30), genes_range = c(15,500), ...) {
    ## i. pre-processing
    ## restructure into a list of genes
    gene_list <- split(iso_data[,-1], iso_data[,1])

    ## isoform distribution
    ## tmp_gene <- gene_list[[1]]
##     num_iso <- unlist(lapply(gene_list, function(tmp_gene) {
##         ## retrive number of isoforms
##         dim(tmp_gene)[1]
##     }))
##     summary(num_iso)
##     sum(num_iso == 1)/length(num_iso) ## 28% only have one form
##     hist(log(num_iso))
##     
    ## impose min/max number of isoforms and at least one alternative isoform
    filter_logic <- unlist(lapply(gene_list, function(tmp_gene) {
        ## retrive number of isoforms
        tmp_num_iso <- dim(tmp_gene)[1]
        tmp_num_iso < iso_range[1] | tmp_num_iso > iso_range[2]
    }))

    genes_to_keep <- names(filter_logic)[!filter_logic]
    gene_list <- gene_list[genes_to_keep]

    ## 1. Find genewise distances
    gene_dist <- transform_iso_gene(X = gene_list, method = gene_method, ...)
 
    ## 2. Compute pathway-level metrics
    transform_gene_pathway(gene_dist = gene_dist, annot_file = annot_file,
                           desc_file = desc_file, method = pathway_method,
                           genes_range = genes_range,...)
}


## iso_data <- iso_kegg_list[["TCGA.BH.A1FM"]]
## iso_data <- iso_kegg_list[[91]]
## annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt"
## desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt"
## pathway_method <- "EE"
## gene_method <- "hellinger"
## iso_range = c(2,30)
## genes_range = c(15,500)
## 
## iso_data <- iso_kegg_list[[sample(1:length(iso_kegg_list),1)]]
## system.time(example_avg <- transform_iso_pathway(iso_data = iso_data, annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", pathway_method = "EE", gene_method = "hellinger")) ## 2 seconds
## 
## table(example_avg$diff_splice_call)
## head(example_avg, sum(example_avg$diff_splice_call) + 1)
## example_avg[is.na(example_avg$pathway_score), "diff_splice_call"]
## example_avg[!is.na(example_avg$pathway_score), "diff_splice_call"]
## qplot(x = pathway_score, data = example_avg)
## head(sort(example_avg$num_genes_measured))
## head(example_avg[order(example_avg$num_genes_measured),])
## sum(example_avg$num_genes_measured >= 15)
## example_avg[grep("cancer", x = example_avg$pathway_desc),]
## 
