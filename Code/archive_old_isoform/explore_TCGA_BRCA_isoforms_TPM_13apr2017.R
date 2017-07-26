## Explore TCGA-BRCA Isoforms on TPM scale
## AG Schissler
## Created 13 Apr 2017
## Last modified 13 Apr 2017

################################
## 1. Create environment and read in counts
setwd("~/Dropbox/Splice-n-of-1-pathways/UofUtah")

## read in FPKM Isoform data
library(data.table)
norm_data <- fread(file = "Total_TPM_normal.txt", sep = "\t", header = T, stringsAsFactors = F)
tumor_data <- fread(file = "Total_TPM_tumor.txt", sep = "\t", header = T, stringsAsFactors = F)

## check normalization
## colSums(norm_data[,-(1:4)])
## colSums(tumor_data[,-(1:4)])

## combine into one data frame
iso_data <- merge(norm_data, tumor_data)

## 113 matched pairs of patients with 4 columns of loci ids

## define the Hellinger distance (bounded between 0 and sqrt(2))
compute_hellinger <- function(X) {
    ## if there is more than one isoform compute the hellinger distance
    if (length(dim(X)) == 2) {
        sqrt(0.5*sum((sqrt(X[,1]) - sqrt(X[,2]))^2))
    } else 0
}

################################
## 2. reorganize and rename using T and N
names(iso_data) <- gsub("01[A-Z]","T", names(iso_data))
names(iso_data) <- gsub("11[A-Z]","N", names(iso_data))

################################
## 3. Load gene set definitions to subset iso_data (makes computation more manageable)
go_bp <- read.delim2("~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_filtered15-500.txt")

go_genes <- unique(go_bp$symbol)
length(go_genes) ## 11530 genes

iso_go_data <- iso_data[iso_data$Gene_symbol %in% go_genes,]
nrow(iso_go_data)
nrow(iso_data)

## free some memory
rm(iso_data)

## remove the data.table class for convenience
class(iso_go_data) <- "data.frame"

################################
## 4. explore isoform distribution within genes

## retrieve the single pair isoform expression data and save for later use
tmp_pat <- "TCGA.A7.A0CE"

i <- grep(tmp_pat, names(iso_go_data))[1]
j <- grep(tmp_pat, names(iso_go_data))[2]

example_iso_data <- iso_go_data[, c(3,i,j)]
save(example_iso_data, file = "~/Dropbox/Splice-n-of-1-pathways/Data/example_iso_tpm_data.RData")

## get distribution of the number of alternative spliced proteins
system.time(gene_go_data <- split(iso_go_data, iso_go_data$Gene_symbol)) ## 25 seconds

length(gene_go_data) ## 11,436 unique genes (gained 1 gene??)

gene_splice_count <- unlist(lapply(gene_go_data, nrow))

summary(gene_splice_count) ## some very large number of isoforms
quantile(gene_splice_count, seq(0,1,0.05)) ## some very large number of isoforms
tail(sort(gene_splice_count), 100)
qplot(gene_splice_count)
qplot(log(gene_splice_count))

sum(gene_splice_count == 1) ## 1129 genes only have one isoform
sum(gene_splice_count > 1) ## 10306 genes only have one isoform

## 4.1 choose one gene with the median number of isoforms
sum(gene_splice_count == median(gene_splice_count)) ## 724 genes equaling median
set.seed(444)
tmp_gene <- sample(names(gene_splice_count)[gene_splice_count == median(gene_splice_count)], 1)


## select the first pair of transcriptomes
i <- 5
j <- 118
tmp_gene_data <- gene_go_data[[tmp_gene]]
(tmp_tumor <- as.numeric(tmp_gene_data[[i]]))
(tmp_normal <- as.numeric(tmp_gene_data[[j]]))

## relative proportion
rel_normal <- tmp_normal/sum(tmp_normal)
rel_tumor <- tmp_tumor/sum(tmp_tumor)

## clean up a bit
rel_normal[rel_normal < 1e-04] <- 0
rel_tumor[rel_tumor < 1e-04] <- 0

rel_mat <- cbind(rel_normal, rel_tumor)
plot(rel_normal, rel_tumor)

compute_hellinger(rel_mat)

## 4.2 choose on with many isoforms
sum(gene_splice_count == max(gene_splice_count)) ## 724 genes equaling median
set.seed(44)
tmp_gene <- sample(names(gene_splice_count)[gene_splice_count == max(gene_splice_count)], 1)

tmp_gene_data <- gene_go_data[[tmp_gene]]
(tmp_tumor <- as.numeric(tmp_gene_data[[i]]))
(tmp_normal <- as.numeric(tmp_gene_data[[j]]))

## relative proportion
rel_normal <- tmp_normal/sum(tmp_normal)
rel_tumor <- tmp_tumor/sum(tmp_tumor)

## clean up a bit
rel_normal[rel_normal < 1e-04] <- 0
rel_tumor[rel_tumor < 1e-04] <- 0

rel_mat <- cbind(rel_normal, rel_tumor)
plot(rel_normal, rel_tumor)

compute_hellinger(rel_mat)

################################
## 5. explore hellinger dist distributions across genes

## for one patient
tmp_pat <- "TCGA.A7.A0CE"

i <- grep(tmp_pat, names(tmp_gene_data))[1]
j <- grep(tmp_pat, names(tmp_gene_data))[2]


## tmp_gene <- gene_go_data[[1]]
tmp_pat_data <- lapply(gene_go_data, function(tmp_gene) {
    cbind(tmp_gene[[j]], tmp_gene[[i]])
})

### find the relative alternative splicing
## tmp_gene <- tmp_pat_data[[1]]
tmp_rel_splice <- lapply(tmp_pat_data, function(tmp_gene) {
    apply(tmp_gene, 2, function(tmp_col) {
        if (sum(tmp_col) > 0) {
            tmp_col/sum(tmp_col)
        } else tmp_col
    })  
})

str(tmp_rel_splice)

tmp_hel_dist <- unlist(lapply(tmp_rel_splice, compute_hellinger))
str(tmp_hel_dist)

summary(tmp_hel_dist)

pdf("~/Dropbox/Splice-n-of-1-pathways/Preliminary_figures/Example_patient_gene-wise alterative_spliced_Hellinger_distances.pdf")
qplot(tmp_hel_dist, xlab = paste0("Gene-level Hellinger distances for 11,436 genes in filtered GO-BP for ", tmp_pat))
dev.off()
## bounded by 0 and sqrt(2)

## see if sparsity (large number of isoforms) is driving the hellinger distance

qplot(gene_splice_count, tmp_hel_dist, ylab = paste0("Gene-level Hellinger distances for 11,436 genes in filtered GO-BP for ", tmp_pat))

## not really... (that's good!)

### questions:

## Should we filter to genes with > 1 isoforms or less than 30 as in Johnson and Purdom (2017)?
## Apply empirical Bayes at the gene level or pathway level?

################################
## 6. average the alt splice differences into the pathways

go_bp_list <- split(go_bp, go_bp$path_id)

## str(go_bp_list)

### compute the average Hellinger distance
## tmp_go <- go_bp_list[[2]]

tmp_go_hel <- unlist(lapply(go_bp_list, function(tmp_go) {
    tmp_genes <- as.character(tmp_go[,2])
    tmp_genes <- tmp_genes[tmp_genes %in% names(tmp_hel_dist)]
    mean(tmp_hel_dist[tmp_genes])
}))

str(tmp_go_hel)

summary(tmp_go_hel)

qplot(tmp_go_hel)
qqnorm(tmp_go_hel) + abline(0,1)

################################
## 7. apply empirical Bayes via local FDR on the average distances
require("locfdr")

head(tmp_go_hel)

pdf(file = paste0("~/Dropbox/Splice-n-of-1-pathways/Preliminary_figures/Local_FDR_GO-BP_Hellinger_distances_patient_", tmp_pat,".pdf"))
tmp_locfdr <- locfdr(zz = tmp_go_hel)
dev.off()
attributes(tmp_locfdr)
tmp_locfdr$fp0
head(tmp_locfdr$mat)
which(tmp_locfdr$fdr < 0.2)

tmp_upper <- tmp_locfdr$z.2[2]

## find interesting pathways
tmp_hits <- names(tmp_go_hel[tmp_go_hel > tmp_upper])

go_desc <- read.table("~/Dropbox/Lab-Tools/GeneSets/GO/2015/go_bp_description.txt", header = T, sep = "\t", quote = "")
row.names(go_desc) <- go_desc$path_id

x_tab <- xtable(go_desc[tmp_hits,], caption = paste0("Alternative spliced GO-BP terms via local FDR of average Hellinger distances for patient ", tmp_pat))

print.xtable(x_tab, include.rownames = F)

plot(tmp_locfdr)

################################
## 8. apply empirical Bayes via local FDR gene-wise

tmp_gene_locfdr <- locfdr(zz = tmp_hel_dist) ## TODO: optimize this fit
attributes(tmp_gene_locfdr)

tmp_gene_locfdr$fp0
tmp_gene_upper <- tmp_gene_locfdr$z.2[2]
tmp_gene_hits <- names(tmp_hel_dist[tmp_hel_dist > tmp_gene_upper])

### Fisher's exact test for over abundance of alternatively spliced genes
## tmp_go <- go_bp_list[[2]]

tmp_go_fet <- lapply(go_bp_list, function(tmp_go) {
    ## pathway genes
    in_go <- as.character(tmp_go[,2])
    ## not pathway genes
    not_go <- names(tmp_hel_dist)[!(names(tmp_hel_dist) %in% in_go)]
    ## not alt. spliced
    not_as <- names(tmp_hel_dist)[!(names(tmp_hel_dist) %in% tmp_gene_hits)]
    ## get counts for 2x2 table
    x11 <- length(intersect(tmp_gene_hits, in_go))
    x21 <- length(intersect(not_as, in_go))
    x12 <- length(intersect(tmp_gene_hits, not_go))
    x22 <- length(intersect(not_as, not_go))
    ## check grand sum
    if ( sum(x11, x21, x12, x22) == length(tmp_hel_dist) ) {
        spliced_pathways <- matrix( c( x11, x21, x12, x22),
                                   nrow = 2,
                                   dimnames = list(AS = c("AS", "NotAS"),
                                                   Pathway = c("Pathway", "NotPathway")))
        toReturn <- fisher.test(spliced_pathways, alternative = "greater")
    } else toReturn <- NULL
    ## return value
    toReturn
})

## just p-values
tmp_go_fet_pvalue <- unlist(lapply(tmp_go_fet, function(tmp_go) tmp_go$p.value))

summary(tmp_go_fet_pvalue)
hist(tmp_go_fet_pvalue)
