## R Library for the development of Splice-N-of-1-pathways
## AG Schissler
## Created 12 Apr 2017
## Last modified 13 Apr 2017

############################################################
## i. Code development objects (DO NOT RUN)

## load("~/Dropbox/Splice-n-of-1-pathways/Data/example_iso_tpm_data.RData")
## X <- rel_list[[tmp_symbol]]
## X <- example_iso_data
## tmp_gene <- gene_list[[1]]
## tmp_gene <- rel_list[[1]]
## tmp_symbol <- names(which(tmp_logic))[1]
## system.time(example_hel <- transform_iso_gene(example_iso_data, method = "hellinger"))
## str(example_hel)
## summary(example_hel)
## qplot(example_hel)

############################################################
## 1. Isoform distance functions

## 1.1 compute Hellinger distance for isoforms of a gene given pair of measurements
compute_hellinger <- function(X) {
    ## if there is more than one isoform compute the hellinger distance
    if (length(dim(X)) == 2) {
        ## check that the gene is expressed in both samples
        if (colSums(X)[1] == 0 | colSums(X)[2] == 0) {
            ## only alternatively spliced genes are of interest
            hel_dist <- 0
        } else {
            hel_dist <- sqrt(0.5*sum((sqrt(X[,1]) - sqrt(X[,2]))^2))
        }
    } else hel_dist <- 0
    ## return the value
    return(hel_dist)
}

############################################################
## 2. Isoform-to-gene functions

## transform isoform-level expression to gene level
transform_iso_gene <- function(X, method = c("delta", "hellinger"), ...) {
    ## restructure into a list of genes
    gene_list <- split(X[,-1], X[,1])
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
    }

    ## Yves method of maximum difference in expression
    if (method == "delta") {

    }

    ## return the gene-level metrics
    return(to_return)
}

############################################################
## 3. Gene-to-pathway functions

## transform gene-level expression to pathway level
transform_gene_pathway <- function(gene_dist, annot_file, desc_file = NULL, method = c("kMEn", "avg", "fet"), ...) {
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
    ### 2. Score pathways via the desired method 
    if (method == "avg") {
        ###  2.1 Find the average distance for each pathway
        ## tmp_set <- annot_list[["GO:0000003"]]
        avg_dist <- unlist(lapply(annot_list, function(tmp_set) {
            tmp_genes <- as.character(tmp_set$symbol)
            measured_genes <- tmp_genes[tmp_genes %in% names(gene_dist)]
            mean(gene_dist[measured_genes])
        }))
        ## qplot(avg_dist)
        ## center the data
        ## avg_dist <- avg_dist - mean(avg_dist)
        #### 2.2 Apply local FDR of Efron
        ## pdf(file = paste0("~/Dropbox/Splice-n-of-1-pathways/Preliminary_figures/Local_FDR_GO-BP_Avg_Hellinger_distances_TPM_patient_", tmp_pat,".pdf"))
        tmp_locfdr <- locfdr::locfdr(zz = avg_dist, plot = 0)
        ## dev.off()
        #### 2.3 Find interesting pathways (more differentially splicing)
        ## first find the upper threshold
        tmp_upper <- tmp_locfdr$z.2[2]
        ## then indicate the hits the hits
        ## tmp_hits <- names(avg_dist)[avg_dist >= tmp_upper]
        tmp_call <- rep(0, length(avg_dist))
        tmp_call[avg_dist >= tmp_upper] <- 1
        ## 2.4 format the output
        to_return <- data.frame(pathway_score = avg_dist, direction = NA, fdr_value = tmp_locfdr$fdr, num_genes_annot = annot_lengths, num_genes_measured = measured_lengths, row.names = names(avg_dist), upper_fdr_threshold = tmp_upper, diff_splice_call = tmp_call, stringsAsFactors = F)
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

## transform gene-level expression to pathway level
transform_iso_pathway <- function(iso_data, annot_file, desc_file = NULL, pathway_method = c("kMEn", "avg", "fet"), gene_method = c("delta", "hellinger"), ...) {
    
    ## 1. Find genewise distances
    gene_dist <- transform_iso_gene(X = iso_data, method = gene_method, ...)

    ## 2. Compute pathway-level metrics
    transform_gene_pathway(gene_dist = gene_dist, annot_file = annot_file,
                           desc_file = desc_file, method = pathway_method, ...)
                           
}

## pathway_method <- "avg"
## gene_method <- "hellinger"
## 
## system.time(example_avg <- transform_iso_pathway(iso_data = example_iso_data, annot_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt", desc_file = "~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_description.txt", pathway_method = "avg", gene_method = "hellinger")) ## 12.8 seconds
##  
## str(example_avg)
## head(example_avg)
