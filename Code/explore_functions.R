## R Library for exploration of results
## AG Schissler
## Created 23 Jul 2018

############################################################
## i. Code development objects (DO NOT RUN)

## load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
## load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BLCA_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")


############################################################
## 1. Explore impact of the expression on Hellinger distance

load("~/Dropbox/Splice-n-of-1-pathways/Data/example_iso_tpm_data.RData")
iso_range = c(2,30)
genes_range = c(15,500)
gene_method = "hellinger"
gene_list <- split(example_iso_data[,-1], iso_data[,1])
filter_logic <- unlist(lapply(gene_list, function(tmp_gene) {
## retrive number of isoforms
    tmp_num_iso <- dim(tmp_gene)[1]
    tmp_num_iso < iso_range[1] | tmp_num_iso > iso_range[2]
}))

genes_to_keep <- names(filter_logic)[!filter_logic]
gene_list <- gene_list[genes_to_keep]
gene_dist <- transform_iso_gene(X = gene_list, method = gene_method)

## find the total gene expression from gene list
gene_expr_N <- unlist(lapply(gene_list, function(X) {sum(X[,1])}))
str(gene_expr_N)
gene_expr_T <- unlist(lapply(gene_list, function(X) {sum(X[,2])}))

## filter out names not in gene dist (due to 0s)
gene_expr_N <- gene_expr_N[names(gene_expr_N)[names(gene_expr_N) %in% names(gene_dist)]]

gene_expr_T <- gene_expr_T[names(gene_expr_T)[names(gene_expr_T) %in% names(gene_dist)]]

## check that names match
all(names(gene_expr_N) == names(gene_dist))
all(names(gene_expr_T) == names(gene_dist))
all(names(gene_expr_N) == names(gene_expr_T))

## combine into one data set
dat <- data.frame(Dist=gene_dist, Normal=log2(gene_expr_N+1), Tumor=log2(gene_expr_T+1))
str(dat)

qplot(x=Normal, y=Dist, data=dat)
qplot(x=Tumor, y=Dist, data=dat)

p0 <- ggplot(data=dat, aes(x=Normal, y=Tumor, color=Dist)) + geom_point()
ggsave("Expression colored by Hellinger distance.pdf", p0)
