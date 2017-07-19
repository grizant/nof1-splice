## R Library for the development of Splice-N-of-1-pathways
## AG Schissler
## Created 12 Apr 2017
## Last modified 13 Apr 2017

############################################################
## i. Code development objects (DO NOT RUN)

## load("~/Dropbox/Splice-n-of-1-pathways/Data/example_iso_tpm_data.RData")
## tmp_pat <- "TCGA.A7.A0CE"
## X <- rel_list[[tmp_symbol]]
## X <- example_iso_data
## tmp_gene <- gene_list[[sample(1:length(gene_list), 1)]]
## tmp_gene <- rel_list[[sample(1:length(rel_list), 1)]]
## tmp_symbol <- names(which(tmp_logic))[1]
## system.time(example_hel <- transform_iso_gene(example_iso_data, method = "hellinger"))
## system.time(example_delta <- transform_iso_gene(example_iso_data, method = "delta"))
## system.time(example_delta <- transform_iso_gene(example_iso_data, method = "delta", remove_DEGs = T))
## str(example_delta)
## summary(example_delta)
## qplot(example_delta)

############################################################
## 1. Isoform distance functions

## 1.1 compute Hellinger distance for isoforms of a gene given pair of measurements
compute_hellinger <- function(X) {
    ## if there is more than one isoform compute the hellinger distance
    if (length(dim(X)) == 2) {
        ## check that the gene is expressed in both samples
        if (colSums(X)[1] == 0 | colSums(X)[2] == 0) {
            ## only alternatively spliced genes are of interest
            ## hel_dist <- 0
            hel_dist <- NA
        } else {
            hel_dist <- sqrt(0.5*sum((sqrt(X[,1]) - sqrt(X[,2]))^2))
        }
    } else hel_dist <- NA
    ## return the value
    return(hel_dist)
}

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
        hist(gene_dist, breaks = 120)
        ## remove dist = 1 (or near it) or 0 due to fitting difficulties.
        ## ... give 1 a positive call and 0 negative call manually later
        tolerance <- .Machine$double.eps^0.5
        dist_one <- which(abs(gene_dist - 1) < tolerance)
        dist_zero <- which(abs(gene_dist) < tolerance)
        remove_index <- c(dist_one, dist_zero)
        ## now filter
        transformed_dist <- gene_dist[-remove_index]
        transformed_dist <- ecdf(transformed_dist)(transformed_dist)
        transformed_dist <- qnorm(transformed_dist)
        ## replace Inf (corresponding to 1) with max + esp
        eps <- 0.001
        transformed_dist[is.infinite(transformed_dist)] <- max(transformed_dist[!is.infinite(transformed_dist)], na.rm = T) + eps
        ## now transform
        ## tail(sort(transformed_dist))
        ## hist(transformed_dist, 120)
        ## summary(transformed_dist)
        tmp_locfdr <- locfdr::locfdr(zz = transformed_dist, bre = 120, df = 7, pct = 0, plot = 4)
        tmp_locfdr <- locfdr::locfdr(zz = gene_dist, bre = 120, df = 7, pct = 0, pct0 = 1/4, plot = 0)
        tmp_upper <- tmp_locfdr$z.2[2]
        ## all calls to gene dist (this may be anti-conservative/conservative depending on the ones and zeros)
        dist_data <- data.frame(gene_dist, call = 0)
        dist_data[gene_dist > tmp_upper, "call"] <- 1
        ## table(dist_data$call)
        
        ## 2.2.2 enrichment
        ## tmp_set <- annot_list[[6]]
        odds_ratio <- unlist(lapply(annot_list, function(tmp_set) {
            tmp_genes <- as.character(tmp_set$symbol)
            measured_genes <- tmp_genes[tmp_genes %in% names(gene_dist)]
            ## genes outside pathway
            not_genes <- rownames(dist_data)[!(rownames(dist_data) %in% measured_genes)]
            ## impose genes range
            if (length(measured_genes) >= genes_range[1] & length(measured_genes) <= genes_range[2]) {
                ## calculate counts
                x11 <- sum(dist_data[measured_genes, "call"])
                x21 <- sum(dist_data[measured_genes, "call"] == 0)
                x12 <- sum(dist_data[not_genes, "call"])
                x22 <- sum(dist_data[not_genes, "call"] == 0)
                ## sum(x11, x21, x12, x22) == nrow(dist_data)
                ## compute FET
                tmp_fet <- fisher.test(x = matrix( c(x11, x21, x12, x22), nrow = 2, ncol = 2))
                to_return <- tmp_fet$estimate
                names(to_return) <- NULL
            } else to_return <- NA
            return(to_return)
        }))
        
        ## qplot(odds_ratio)
        ## summary(odds_ratio)
        
        ## 2.2.3 local FDR at the pathway level
        fil_odds <- odds_ratio[!is.na(odds_ratio)] ## use log odds??
        tmp_locfdr <- locfdr::locfdr(zz = fil_odds, bre = ceiling(length(fil_odds)/8), df = 7, pct = 0, pct0 = 1/4, plot = 0)
        tmp_upper <- tmp_locfdr$z.2[2]
        ## then indicate the hits the hits
        ## tmp_hits <- names(avg_dist)[avg_dist >= tmp_upper]
        tmp_call <- rep(0, length(odds_ratio))
        tmp_call[odds_ratio >= tmp_upper] <- 1
        ## table(tmp_call)

        ## 2.2.4 format the output
        to_return <- data.frame(pathway_score = odds_ratio, direction = NA, fdr_value = 1, num_genes_annot = annot_lengths, num_genes_measured = measured_lengths, row.names = names(odds_ratio), upper_fdr_threshold = tmp_upper, diff_splice_call = tmp_call, stringsAsFactors = F)

        ## str(tmp_locfdr)
        ## update fdr_value
        to_return[names(fil_odds), "fdr_value"] <- tmp_locfdr$fdr
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

## pathway_method <- "avg"
## gene_method <- "hellinger"
## pathway_method <- "EE"
## an_iso_data <- iso_kegg_list[["TCGA.BH.A1FM"]]
## system.time(example_avg <- transform_iso_pathway(iso_data = an_iso_data, annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt", pathway_method = "EE", gene_method = "hellinger")) ## 25 seconds
##  
## str(example_avg)
## head(example_avg)
## example_avg[is.na(example_avg$pathway_score), "diff_splice_call"]
## example_avg[!is.na(example_avg$pathway_score), "diff_splice_call"]
## table(example_avg$diff_splice_call)
## qplot(x = pathway_score, data = example_avg)
## head(sort(example_avg$num_genes_measured))
## head(example_avg[order(example_avg$num_genes_measured),])
## sum(example_avg$num_genes_measured >= 15)

