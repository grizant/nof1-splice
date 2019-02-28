## AG Schissler
## Run N1PAS KEGG based simulation
## 28 Feb 2019

#################################################################################
## i. Set up workspace

setwd("~/Research/n1pas/sim_files")
print("Running n1pas simulation")

num.args <- 5
## add a clustering method parameter
## num.args <- 

print("Checking parameters")
args = commandArgs(TRUE)
if (length(args) != num.args) {
    stop("Not enough arguments! Need to specify dataset, pathway size (P), proportion of ASGs (pi), seed, outfile")
}

#######################################################
### 1. Read in user-specified parameters
#######################################################

## mimic the gpusim parameters

### test parameters
## args <- c("ucec", 50, 0.5, 44, "~/Research/n1pas/sim_results/test.txt")

## parse options
args.values <- args

## 1. all iso data
dataset <- args.values[1]
## 2. expressed genes in pathway size
p <- as.numeric(args.values[2])
## 3. pi, proportion of ASGs above background
pi <- as.numeric(args.values[3])
## 4. seed
seed <- as.numeric(args.values[4])
## 6. output file
outfile <- args.values[5]

## id <- paste(n,"-",d,"-","copula","1")
## today <- format(Sys.time(), "%Y-%m-%d.%H:%M:%S")
print("Arguments passed to simulation script")
print(args)

#######################################################
### 2. Read in input files, restructure, and set standard parameters
#######################################################

print("Reading iso data into memory")
## use rds next time
iso_file <- paste("~/Research/n1pas/data/", dataset, "_iso_kegg_data.RData", sep = "")
system.time(load(file=iso_file))
iso_kegg_data <- get(paste0(dataset, "_iso_kegg_data"))
rm(list = paste0(dataset, "_iso_kegg_data"))

print("Reading in pathway information and custom splice functions")

## source functions
source("~/Research/n1pas/splice_functions.R")

## set standard parameters
annot_file = "~/Research/n1pas/data/kegg_tb.txt"
desc_file = "~/Research/n1pas/data/kegg.description_tb.txt"
go_desc_file = F
pathway_method <- "EE"
gene_method <- "hellinger"
iso_range = c(2,30)
genes_range = c(15,500)
expr_threshold = 1
path_id_name = "path_id"
DEGs = NULL
pct0 = 1/4

###### Restructure into convenient data structures
## Retrieve patient IDs
pat_col <- grep("TCGA", x = names(iso_kegg_data))
patients_chr <- unique(substring(names(iso_kegg_data[pat_col]), 1, 12))

## create a empty list
iso_kegg_list <- vector(mode = "list", length =  length(patients_chr))
names(iso_kegg_list) <- patients_chr

## tmp_pat <- patients_chr[1]
for (tmp_pat in patients_chr) {
    ## retrieve gene symbols and the paired transcriptomes
    iso_kegg_list[[tmp_pat]] <- (data.frame(geneSymbol = iso_kegg_data$geneSymbol,
                                            iso_kegg_data[, grep(tmp_pat, names(iso_kegg_data))]))
    iso_kegg_list[[tmp_pat]][, "geneSymbol"] <- as.character(iso_kegg_data$geneSymbol)
}
rm(iso_kegg_data)

#################################################################################
## 3. N1PAS simulation function
#######################################################

## iso_data <- iso_kegg_list[[2]]
runSim_n1pas <- function(iso_data, p, pi){

    ## 1. filter by isoform number and expressed genes
    ## convert iso data into a genewise list
    gene_list <- split(iso_data[,-1], iso_data[,1]) ## 5757 genes
    filter_logic <- unlist(lapply(gene_list, function(tmp_gene) {
        ## retrive number of isoforms
        tmp_num_iso <- dim(tmp_gene)[1]
        tmp_num_iso < iso_range[1] | tmp_num_iso > iso_range[2]
    }))
    genes_to_keep <- names(filter_logic)[!filter_logic]
    gene_list <- gene_list[genes_to_keep] ## 4133 genes (this is NOT patient specific!!)
    ## remove low expressing genes
    if (expr_threshold > 0) {
        print("removing low expressing genes from both samples")
        ## tmp_gene <- gene_list[[1]]
        keep_logic <- unlist(lapply(gene_list, function(tmp_gene) {
            ## sum the isoform level to find the gene-level expression
            (colSums(tmp_gene)[1] >= expr_threshold) & (colSums(tmp_gene)[2] >= expr_threshold)
        }))
        num_genes_removed <- length(keep_logic) - sum(keep_logic)
        prop_genes_removed <- round(num_genes_removed/length(keep_logic),3)
        print(paste("Removed", num_genes_removed, "genes, proportion of total is", prop_genes_removed))
        gene_list <- gene_list[keep_logic]
    }

    ## 2. transform into gene level via specified method
    if (gene_method == "hellinger") {
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
        gene_dist <- unlist(lapply(rel_list, compute_hellinger, only_expressed = T))
        ## sum(is.na(to_return))
        ## remove missing distances (due to non-expression or low expression in either sample)
        gene_dist <- gene_dist[!is.na(gene_dist)]
    }

    ## 3. Shuffle gene labels and classify genes into ASG or not
    rand_gene_dist <- sample(gene_dist)
    names(rand_gene_dist) <- sample(x = names(rand_gene_dist))
    rand_cluster_fit <- kmeans(x = rand_gene_dist, 2)
    ## get the cluster labels
    asg_label <- which.max(rand_cluster_fit$centers)
    non_asg_label <- which.min(rand_cluster_fit$centers)
    num_not_asg <- table(rand_cluster_fit$cluster)[1]
    num_asg <- table(rand_cluster_fit$cluster)[2]
    ## report background ASG
    prop_asg <- unname((table(rand_cluster_fit$cluster) / length(gene_dist))[asg_label])
    print(paste("Proportion of total genes that are ASGs is", prop_genes_removed))
    asg_values <- unname(rand_gene_dist[rand_cluster_fit$cluster == asg_label])
    non_asg_values <- unname(rand_gene_dist[rand_cluster_fit$cluster == non_asg_label])

    ## 4. Filter pathways based on expressed genes and select a pathway matches specified size
    ## select p pathways, each pathway has m_p genes
    annot_data <- read.delim2(file = annot_file, stringsAsFactors = F)
    id_column <- which(names(annot_data) == path_id_name)
    annot_list <- split(annot_data$symbol, annot_data[,id_column])
    ## store the number of genes annotated
    annot_lengths <- unlist(lapply(annot_list, function(tmp_genes) {length(unique(tmp_genes))}))
    ## count expressed genes
    measured_lengths <- unlist(lapply(annot_list, function(tmp_genes) {
        tmp_genes <- as.character(tmp_genes)
        sum(tmp_genes %in% names(gene_dist))
    }))
    ## select pathway with specified expressed genes
    tmp_index <- which(measured_lengths == p)
    if (length(tmp_index) > 0){
        print("At least one pathway of target size.")
        path_id <- sample(names(tmp_index), 1)
    } else {
        stop("TARGET PATHWAY SIZE NOT FOUND, DEVELOP A FUZZY MATCH SYSTEM?")
    }

    ## 5. Inject signal into target pathway
    pi <- min(prop_asg + pi, 1)
    print(paste("Updated pi over background is", pi))
    tmp_genes <- annot_list[[path_id]]
    tmp_measured_genes <- names(gene_dist[names(gene_dist) %in% tmp_genes])
    tmp_num_asg <- round(pi * length(tmp_measured_genes))
    tmp_num_non_asg <- length(tmp_measured_genes) - tmp_num_asg
    ## Select values for genes, careful about pi = 0, 1
    ## doesn't matter which genes you select for first case,
    ## as long as they are from the same cluster
    if (tmp_num_asg > 0 & tmp_num_non_asg > 0) {
        asg_labels <- tmp_measured_genes[1:tmp_num_asg]
        non_asg_labels <- tmp_measured_genes[-(1:tmp_num_asg)]
    } else {
        if (tmp_num_asg == 0) {
            asg_labels <- NULL
            non_asg_labels <- tmp_measured_genes[1:tmp_num_non_asg]
        } else {
            asg_labels <- tmp_measured_genes[1:tmp_num_asg]
            non_asg_labels <- NULL
        }
    }
    ## fill in the rest of the gene labels
    ## find unused genes
    assigned_genes <- c(asg_labels, non_asg_labels)
    length(assigned_genes)
    unassigned_genes <- names(rand_gene_dist)[!(names(rand_gene_dist) %in% assigned_genes)]
    if ((length(unassigned_genes) + length(assigned_genes)) != length(rand_gene_dist)){
        stop("TOTAL OF ASSIGNED AND UNASSIGNED LABLES DO NOT MATCH")
    }
    assigned_genes <- c(asg_labels, non_asg_labels)
    to_add_asg <- length(asg_values) - length(asg_labels)
    to_add_asg_labels <- unassigned_genes[1:to_add_asg]
    names(asg_values) <- c(asg_labels, to_add_asg_labels)
    to_add_non_asg_labels <- unassigned_genes[(to_add_asg+1):length(unassigned_genes)]
    names(non_asg_values) <- c(non_asg_labels, to_add_non_asg_labels)
    labeled_rand_gene_dist <- c(asg_values, non_asg_values)

    ## 6. Now run N1PAS
    gene_clusters <- kmeans(x = labeled_rand_gene_dist, 2)
    dist_data <- data.frame(labeled_rand_gene_dist, call = 0)
    alt_center <- which.max(gene_clusters$centers)
    dist_data[names(gene_clusters$cluster), "call"] <- ifelse(gene_clusters$cluster == alt_center, 1, 0)
    odds_ratio <- lapply(annot_list, get_OR, dist_data, genes_range=genes_range)
    num_genes <- unlist(lapply(odds_ratio, function(tmp_id) {tmp_id$num_genes}))
    odds_ratio <- unlist(lapply(odds_ratio, function(tmp_id) {tmp_id$odds_ratio}))
    fil_odds <- odds_ratio[!is.na(odds_ratio)]
    (outliers <- boxplot.stats(fil_odds)$out)
    if (length(outliers) > 0) fil_odds <- fil_odds[!(names(fil_odds) %in% names(outliers))]
    suppressWarnings(tmp_locfdr <- locfdr::locfdr(zz = fil_odds, bre = ceiling(length(fil_odds)/8), df = 4, pct = 0,
                                                  pct0 = pct0, nulltype = 2, plot = 0, mlests = c(1, sd(fil_odds))))
    tmp_upper <- tmp_locfdr$z.2[2]
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

    ## 7. Format output
    scored_data <- data.frame(pathway_score = odds_ratio, direction = NA, fdr_value = 1, num_genes_annot = annot_lengths, num_genes_measured = measured_lengths, row.names = names(odds_ratio), upper_fdr_threshold = tmp_upper, diff_splice_call = tmp_call, stringsAsFactors = F)
    ## str(tmp_locfdr)
    ## update fdr_value
    scored_data[names(fil_odds), "fdr_value"] <- tmp_locfdr$fdr
    ## and give p-value of 0 for outliers
    if (length(outliers) > 0) {scored_data[names(outliers), "fdr_value"] <- 0}
    ## head(sort(to_return$fdr_value), 50)
    ## sort
    scored_data <- scored_data[order(scored_data$pathway_score, decreasing = T),]

    ## 8. Produce simulation output
    target_captured <- scored_data[path_id, "diff_splice_call"]
    ## caclutate fpr
    if (pi - prop_asg > 0) {
        ## false positive rate discounting target pathway
        false_positive_rate <- sum(scored_data[!(rownames(scored_data) %in% path_id), "diff_splice_call"]) / ( nrow(scored_data) - 1)
    } else {
        ## false positive rate counting target pathway
        false_positive_rate <- sum(scored_data[, "diff_splice_call"]) / ( nrow(scored_data) )
    }
    output_data <- data.frame(dataset = dataset, pathway_size = p, prop_asg = prop_asg, pi = pi - prop_asg,
                              target_captured = target_captured, false_positive_rate = false_positive_rate)
    return(output_data)
}

#######################################################
### 4. Run simulation
#######################################################
print("Simulating N1PAS")

iso_data <- iso_kegg_list[[2]]
set.seed(50)
runSim_n1pas(iso_data = iso_data, p = 50, pi = 0.2)

system.time(sim_data <- my_norta(n=n, r_mat=R_hat, margins=param_margins[[1]],
                                 paramMargins=param_margins[[2]]))

#######################################################
### 5. Check
#######################################################

## R_hat[1:5,1:5]
## cor(sim_data)[1:5,1:5]

#######################################################
### 6. Writing files
#######################################################
print("Save simulated data")
system.time(MASS::write.matrix(sim_data, outfile))
