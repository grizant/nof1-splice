## Effect size simulation script
## AG Schissler
## Created 25 Feb 2019

##############################################################################
#### 1. Begin with example TCGA + KEGG data set to develop simulation

## load TCGA BLCA TPM isoform data
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/blca_iso_kegg_data.RData") ## NEW ISOFORM 25 Jul 2017

## source functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/splice_functions.R")

##############################################################################
#### 2. Restructure iso data into a patient-wise list for parallel processing

## Retrieve patient IDs
pat_col <- grep("TCGA", x = names(blca_iso_kegg_data))
patients_chr <- unique(substring(names(blca_iso_kegg_data[pat_col]), 1, 12))

## create a empty list
iso_kegg_list <- vector(mode = "list", length =  length(patients_chr))
names(iso_kegg_list) <- patients_chr

## tmp_pat <- patients_chr[1]
for (tmp_pat in patients_chr) {
    ## retrieve gene symbols and the paired transcriptomes
    iso_kegg_list[[tmp_pat]] <- (data.frame(geneSymbol = blca_iso_kegg_data$geneSymbol,
                                            blca_iso_kegg_data[, grep(tmp_pat, names(blca_iso_kegg_data))]))
    iso_kegg_list[[tmp_pat]][, "geneSymbol"] <- as.character(blca_iso_kegg_data$geneSymbol)
}

##############################################################################
#### 3. Work with an example patient

set.seed(1)
iso_data <- as.data.frame(sample(iso_kegg_list, 1))
annot_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt"
desc_file = "~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt"
go_desc_file = F
pathway_method <- "EE"
gene_method <- "hellinger"
iso_range = c(2,30)
genes_range = c(15,500)
expr_threshold = 1
path_id_name = "path_id"
DEGs = NULL
pct0 = 1/4

system.time(example_ee <- transform_iso_pathway(iso_data = iso_data, annot_file = annot_file, desc_file = desc_file, pathway_method = pathway_method, gene_method = gene_method, DEGs=DEGs, genes_range=genes_range, expr_threshold=expr_threshold, path_id_name = path_id_name, go_desc_file = go_desc_file)) 

str(example_ee)
qplot(example_ee$pathway_score)

## work through splice functions

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
    print("removing low expressing genes from NORMAL sample only")
        ## tmp_gene <- gene_list[[1]]
    keep_logic <- unlist(lapply(gene_list, function(tmp_gene) {
        ## sum the isoform level to find the gene-level expression
        ## (colSums(tmp_gene)[1] >= expr_threshold) & (colSums(tmp_gene)[2] >= expr_threshold)
        (colSums(tmp_gene)[1] >= expr_threshold)
    }))
    num_genes_removed <- length(keep_logic) - sum(keep_logic)
    prop_genes_removed <- round(num_genes_removed/length(keep_logic),3)
    print(paste("Removed", num_genes_removed, "genes, proportion of total is", prop_genes_removed)) ## 619 genes removed
    gene_list <- gene_list[keep_logic]
}



## transform via specified method
if (gene_method == "hellinger") {
    ## first compute the relative expression within a gene
    ## tmp_gene <- gene_list[[1]]
    ## tmp_col <- gene_list[[1]]
    rel_list <- lapply(gene_list, function(tmp_gene) {
        apply(tmp_gene, 2, function(tmp_col) {
            ## check that there is no division by 0
            if (sum(tmp_col) > 0) {
                tmp_col/sum(tmp_col)
            } else tmp_col
        })
    })
    ## then compute the Hellinger distance
    to_return <- unlist(lapply(rel_list, compute_hellinger, only_expressed = T))
    ## sum(is.na(to_return))
    ## remove missing distances (due to non-expression or low expression in either sample)
    to_return <- to_return[!is.na(to_return)]
    gene_dist <- to_return
}

##############################
## explore Hellinger distances
qplot(gene_dist)

## proportions of 0's and 1's to see if I should put more weight there
num_dist_genes <- length(gene_dist) ## 3449 genes
(num0 <- sum(gene_dist == 0) / length(gene_dist)) ## 8.7%
(num1 <- sum(gene_dist == 1) / length(gene_dist)) ## 0.0171

## view this as a four-component mixture with two parts for AGSs and non-ASGs
cluster_fit <- kmeans(x = gene_dist, 2)
table(cluster_fit$cluster)
num_not_asg <- table(cluster_fit$cluster)[1]
num_asg <- table(cluster_fit$cluster)[2]
table(cluster_fit$cluster) / length(gene_dist)

## qplot(x = gene_dist, fill = cluster_fit$cluster)

summary(gene_dist[names(which(cluster_fit$cluster == 2))])
summary(gene_dist[names(which(cluster_fit$cluster == 1))])

qplot(gene_dist[names(which(cluster_fit$cluster == 2))])

## ## two part component for ASG
## ## induce variation for the mixing proportion
## pi <- num1
## (sim_pi <- rnorm(1, num1, sd = 0.01))
## (sim_num1 <- round(sim_pi * num_dist_genes))
## my_a <- 10
## my_b <- 0.01
## sim_asg <- sample(c(rep(1, sim_num1), rgamma(num_asg - sim_num1, shape = my_a, scale = my_b)))
## qplot(sim_asg)
## qplot(gene_dist[names(which(cluster_fit$cluster == 2))])
## sum(sim_asg == 1)
## summary(sim_asg)
## summary(gene_dist[names(which(cluster_fit$cluster == 2))])
## 
## ## now for non ASG part
## pi <- num0
## sim_pi <- rnorm(1, num0, sd = 0.01)
## (sim_num0 <- round(sim_pi * num_dist_genes))
## my_a <- 0.5
## my_b <- 2
## sim_not_asg <- sample(c(rep(0, sim_num0), rbeta(num_not_asg - sim_num0, my_a, my_b)))
## qplot(sim_not_asg)
## qplot(gene_dist[names(which(cluster_fit$cluster == 1))])
## sum(sim_not_asg == 0)
## summary(sim_not_asg)
## summary(gene_dist[names(which(cluster_fit$cluster == 1))])
## 

##############################
## 2. pathways information

annot_data <- read.delim2(file = annot_file, stringsAsFactors = F)
id_column <- which(names(annot_data) == path_id_name)
annot_list <- split(annot_data$symbol, annot_data[,id_column])
## store the number of genes annotated
annot_lengths <- unlist(lapply(annot_list, function(tmp_genes) {length(unique(tmp_genes))}))
## count genes actually measured
measured_lengths <- unlist(lapply(annot_list, function(tmp_genes) {
    tmp_genes <- as.character(tmp_genes)
    sum(tmp_genes %in% names(gene_dist))
}))


##############################
## OR just suffle labels of genes?
measured_genes <- names(gene_dist)
rand_gene_dist <- sample(gene_dist)
names(rand_gene_dist) <- sample(x = names(rand_gene_dist))

rand_cluster_fit <- kmeans(x = rand_gene_dist, 2)
table(rand_cluster_fit$cluster)
(asg_label <- which.max(rand_cluster_fit$centers))
(non_asg_label <- which.min(rand_cluster_fit$centers))
(num_not_asg <- table(rand_cluster_fit$cluster)[1])
(num_asg <- table(rand_cluster_fit$cluster)[2])
table(rand_cluster_fit$cluster) / length(gene_dist)
(prop_asg <- unname((table(rand_cluster_fit$cluster) / length(gene_dist))[asg_label]))

asg_values <- unname(rand_gene_dist[rand_cluster_fit$cluster == asg_label])
non_asg_values <- unname(rand_gene_dist[rand_cluster_fit$cluster == non_asg_label])

## select p pathways, each pathway has m_p genes
str(annot_list)
set.seed(4444)
p <- 1
min_filter <- 15
max_filter <- 500
## select pathway with enough measured genes
num_measured_genes <- unlist(lapply(annot_list, FUN = function(tmp_genes) {
    sum(names(gene_dist) %in% tmp_genes)
}))

filter_logic <- num_measured_genes >= min_filter & num_measured_genes <= max_filter
str(annot_list[filter_logic]) ## 184 pathways
target_pathways <- sample(annot_list[filter_logic], p)

## with pi * m_p genes ASGs
## proportion above background
pi <- 0.25
(pi <- prop_asg + pi)
## reassign dist values to match specified enrichment
## tmp_genes <- target_pathways[[1]]
pathway_count <- 1
asg_labels <- NULL
non_asg_labels <- NULL


for (tmp_genes in target_pathways) {
    ## get measured genes
    tmp_measured_genes <- names(gene_dist[names(gene_dist) %in% tmp_genes])
    ## table(rand_cluster_fit$cluster[names(rand_cluster_fit$cluster) %in% tmp_genes])
    ## select genes to give values from different clusters
    tmp_num_asg <- round(pi * length(tmp_measured_genes))
    tmp_num_non_asg <- length(tmp_measured_genes) - tmp_num_asg
    ## doesn't matter which gene gets which value
    ## careful here as the ASG vs NON-ASG assignment must be consistent across pathways
    if (is.null(asg_labels)) {
        ## doesn't matter which genes you select for first case
        tmp_asg_labels <- tmp_measured_genes[1:tmp_num_asg]
        tmp_non_asg_labels <- tmp_measured_genes[-(1:tmp_num_asg)]
    } else {
        ## find genes already in the label lists
        tmp_already_asg <- asg_labels[asg_labels %in% tmp_measured_genes]
        if (!is.null(tmp_already_asg)) {
            print("Gene overlap in asg list found.")
            tmp_num_asg <- tmp_num_asg - length(tmp_already_asg)
            tmp_measured_genes <- tmp_measured_genes[!(tmp_measured_genes %in% tmp_already_asg)]
        }
        tmp_already_non_asg <- non_asg_labels[non_asg_labels %in% tmp_measured_genes]
        if (!is.null(tmp_already_non_asg)) {
            print("Gene overlap in non-asg list found.")
            tmp_num_non_asg <- tmp_num_non_asg - length(tmp_already_non_asg)
            tmp_measured_genes <- tmp_measured_genes[!(tmp_measured_genes %in% tmp_already_non_asg)]
        }
        tmp_asg_labels <- tmp_measured_genes[1:tmp_num_asg]
        tmp_non_asg_labels <- tmp_measured_genes[-(1:tmp_num_asg)]
    }

    ## aggregate labels
    asg_labels <- c(asg_labels, tmp_asg_labels)
    non_asg_labels <- c(non_asg_labels, tmp_non_asg_labels)
}

length(asg_labels)
length(asg_values)
length(non_asg_labels)

## fill in the rest of the gene labels
## find unused genes
assigned_genes <- c(asg_labels, non_asg_labels)
length(assigned_genes)
unassigned_genes <- names(rand_gene_dist)[!(names(rand_gene_dist) %in% assigned_genes)]
length(unassigned_genes)
length(unassigned_genes) + length(assigned_genes)
length(rand_gene_dist)

to_add_asg <- length(asg_values) - length(asg_labels)
to_add_asg_labels <- unassigned_genes[1:to_add_asg]
names(asg_values) <- c(asg_labels, to_add_asg_labels)

to_add_non_asg_labels <- unassigned_genes[(to_add_asg+1):length(unassigned_genes)]
names(non_asg_values) <- c(non_asg_labels, to_add_non_asg_labels)

labeled_rand_gene_dist <- c(asg_values, non_asg_values)
length(labeled_rand_gene_dist)
length(gene_dist)

### now run the method from this point
gene_clusters <- kmeans(x = labeled_rand_gene_dist, 2)
dist_data <- data.frame(labeled_rand_gene_dist, call = 0)
alt_center <- which.max(gene_clusters$centers)
dist_data[names(gene_clusters$cluster), "call"] <- ifelse(gene_clusters$cluster == alt_center, 1, 0)
odds_ratio <- lapply(annot_list, get_OR, dist_data, genes_range=genes_range)
num_genes <- unlist(lapply(odds_ratio, function(tmp_id) {tmp_id$num_genes}))
odds_ratio <- unlist(lapply(odds_ratio, function(tmp_id) {tmp_id$odds_ratio}))
qplot(odds_ratio)
## detect outlying pathways above prior to fitting
fil_odds <- odds_ratio[!is.na(odds_ratio)]
(outliers <- boxplot.stats(fil_odds)$out)
if (length(outliers) > 0) fil_odds <- fil_odds[!(names(fil_odds) %in% names(outliers))]

suppressWarnings(tmp_locfdr <- locfdr::locfdr(zz = fil_odds, bre = ceiling(length(fil_odds)/8), df = 4, pct = 0, pct0 = pct0, nulltype = 2, plot = 0, mlests = c(1, sd(fil_odds))))

tmp_locfdr <- locfdr::locfdr(zz = fil_odds, bre = ceiling(length(fil_odds)/8), df = 4, pct = 0, pct0 = pct0, nulltype = 2, plot = 1, mlests = c(1, sd(fil_odds)))

(tmp_upper <- tmp_locfdr$z.2[2])
head(sort(tmp_locfdr$fdr))
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

target_names <- names(target_pathways)
to_return[target_names,]
sum(to_return$diff_splice_call)
sum(to_return$diff_splice_call) / nrow(to_return)

## number of patients

num_lusc <- 51
num_luad <- 58
num_prad <- 52
num_thca <- 59
num_ucec <- 7
num_blca <- 19

(num_total <- num_lusc + num_luad + num_prad + num_thca + num_ucec + num_blca)

#################################################################################
## transfer files to okapi server for simulation

## transfer isoform data
system("rsync -vt ~/Dropbox/Splice-n-of-1-pathways/Data/blca_iso_kegg_data.RData aschissler@okapi.math.unr.edu:/home/aschissler/Research/n1pas/data")
system("rsync -vt ~/Dropbox/Splice-n-of-1-pathways/Data/ucec_iso_kegg_data.RData aschissler@okapi.math.unr.edu:/home/aschissler/Research/n1pas/data")
system("rsync -vt ~/Dropbox/Splice-n-of-1-pathways/Data/lusc_iso_kegg_data.RData aschissler@okapi.math.unr.edu:/home/aschissler/Research/n1pas/data")
system("rsync -vt ~/Dropbox/Splice-n-of-1-pathways/Data/luad_iso_kegg_data.RData aschissler@okapi.math.unr.edu:/home/aschissler/Research/n1pas/data")
system("rsync -vt ~/Dropbox/Splice-n-of-1-pathways/Data/prad_iso_kegg_data.RData aschissler@okapi.math.unr.edu:/home/aschissler/Research/n1pas/data")
system("rsync -vt ~/Dropbox/Splice-n-of-1-pathways/Data/thca_iso_kegg_data.RData aschissler@okapi.math.unr.edu:/home/aschissler/Research/n1pas/data")

## splice functions
system("rsync -vt ~/Dropbox/Splice-n-of-1-pathways/Code/splice_functions.R aschissler@okapi.math.unr.edu:/home/aschissler/Research/n1pas")

## KEGG information
system("rsync -vt ~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg_tb.txt aschissler@okapi.math.unr.edu:/home/aschissler/Research/n1pas/data")
system("rsync -vt ~/Dropbox/Lab-Tools/GeneSets/KEGG/kegg.description_tb.txt aschissler@okapi.math.unr.edu:/home/aschissler/Research/n1pas/data")
