## LUSC subtyping via Splice-N-of-1-pathways
## AG Schissler
## Created 27 Jul 2017
## Last modified 23 Feb 2019

############################################################
## i. Load clinical data and results

load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_paired_clinical.RData")
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_LUSC_hel_EE_Iso30_expressiod_pathwayfilter_KEGG_25july2017.RData")
patients <- names(scores_list)
clin_data <- clin_data[rownames(clin_data) %in% patients,]

## load survival analysis functions
source("~/Dropbox/Splice-n-of-1-pathways/Code/surv_functions.R")

## explore clinical
table(clin_data$vitalstatus)

############################################################
## 1. Aggregate pathway scores

effect_mat <- compile_scores(scores_list = scores_list, type = "pathway_score", remove_missing = F)
fdr_mat <- compile_scores(scores_list = scores_list, type = "fdr_value", remove_missing = F)

## filter to only those pathways found dysregulated in at least one patient
fil_effect_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "pathway_score")
fil_fdr_mat <- filter_value_mat(effect_mat = effect_mat, fdr_mat = fdr_mat, type = "fdr_value")
str(fil_effect_mat) ## 100 pathways now
str(fil_fdr_mat)

############################################################
## 2. Select pathways that produce distinct survival curves

## two clusters (not sure how to justify this other than simplicity)
clust_pvalue <- get_pathway_pvalue(value_mat = fil_effect_mat, clin_data = clin_data, type = "pam", num_clusters = 2)
hist(clust_pvalue, breaks = 15)
clust_pvalue <- sort(clust_pvalue)
head(clust_pvalue)

## explore hits
threshold <- 0.05
any(clust_pvalue < threshold)
(hits <- names(clust_pvalue)[clust_pvalue < threshold])
sum(clust_pvalue < threshold)/length(clust_pvalue)

## retrive the hit pathway description 
tmp_data <- scores_list[[1]]
tmp_data[hits,]

#### try fdr adjustments
## standard adjustment
## p.adjust(p = clust_pvalue, "fdr")

## locFDR works and gets around the permutation approach!
transformed_pvalue <- qnorm(clust_pvalue)
eps <- 0.001
transformed_pvalue[is.infinite(transformed_pvalue)] <- max(transformed_pvalue[!is.infinite(transformed_pvalue)], na.rm = T) + eps

tmp_locfdr <- locfdr::locfdr(zz = transformed_pvalue, bre = ceiling(length(clust_pvalue)/8), df = 4, pct = 0, pct0 = 1/64, nulltype = 1, plot = 1)
tmp_locfdr

## 4 hits
top_hits <- names(transformed_pvalue)[transformed_pvalue < tmp_locfdr$z.2[1]]
tmp_data[top_hits,] ## top hits are all interesting including staph infection!

######################################################
## explore top hits

### hit 1
i <- 1
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T) )

## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Good", "Poor")))
p1 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio, fill = cluster))
(p1 <- p1 + geom_boxplot())

## collect clustering across top hits
all_data <- data.frame(cluster1 = factor(top_clust, labels = c("Good", "Poor")))

### hit 2
i <- 2
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T) )

## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Good", "Poor")))
p2 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio, fill = cluster))
(p2 <- p2 + geom_boxplot())

## collect clustering across top hits
all_data[rownames(top_data), "cluster2"] <- top_data$cluster
all.equal(all_data$cluster1, all_data$cluster2)

### hit 3
i <- 3
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T) )

## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Good", "Poor")))
p3 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio, fill = cluster))
(p3 <- p3 + geom_boxplot())

## collect clustering across top hits
all_data[rownames(top_data), "cluster3"] <- top_data$cluster
all.equal(all_data$cluster1, all_data$cluster3)

### hit 4
i <- 4
(top_or <- fil_effect_mat[,hits[i]])
(top_clust <- cluster_pat(top_or, type = "pam", num_clusters = 2))
table(top_clust)
(top_surv <- fit_surv(clusters = top_clust, clin_data = clin_data, plot = T) )

## cluster 1 is better
## odds ratio within top pathway by cluster
top_data <- data.frame(odds_ratio = top_or, cluster = factor(top_clust, labels = c("Good", "Poor")))
p4 <- ggplot(data = top_data, aes(x = cluster, y = odds_ratio, fill = cluster))
(p4 <- p4 + geom_boxplot())

## collect clustering across top hits
all_data[rownames(top_data), "cluster4"] <- top_data$cluster
all.equal(all_data$cluster1, all_data$cluster4)

all_data

#####
## construct the "worse" (poor) group Venn Diagram

clust_type <- "Poor"
## tmp_clust <- all_data$cluster1
## names(tmp_clust) <- rownames(all_data)
clust_list <- lapply(all_data, function(tmp_clust) {
    rownames(all_data[tmp_clust == clust_type,])
})

## all patients in each cluster
(areas <- unlist(lapply(clust_list, length)))

## all pairwise comparisons
num_clusterings <- length(clust_list)
## tmp_index <- pairs_mat[,1]
pairs_mat <- combn(1:num_clusterings, 2)
pairs <- apply(pairs_mat, 2, function(tmp_index){
    x = unname(unlist(clust_list[tmp_index[1]]))
    y = unname(unlist(clust_list[tmp_index[2]]))
    length(intersect(x, y))
})

pairs_worse <- pairs

pairs_Jaccard_worse <- apply(pairs_mat, 2, function(tmp_index){
    x = unname(unlist(clust_list[tmp_index[1]]))
    y = unname(unlist(clust_list[tmp_index[2]]))
    overlap <- length(intersect(x, y))
    union <- length(union(x, y))
    overlap / union
})

pairs_Jaccard_worse

## all 3-way comparisons
## tmp_index <- trips_mat[,1]
trips_mat <- combn(1:num_clusterings, 3)
trips <- apply(trips_mat, 2, function(tmp_index){
    x = unname(unlist(clust_list[tmp_index[1]]))
    y = unname(unlist(clust_list[tmp_index[2]]))
    z = unname(unlist(clust_list[tmp_index[3]]))
    x_and_y = intersect(x, y)
    x_and_y_and_z = intersect(x_and_y, z)
    length(x_and_y_and_z)
})

## all 4-way comparisons
## tmp_index <- quad_mat[,1]
quad_mat <- combn(1:num_clusterings, 4)
quads <- apply(quad_mat, 2, function(tmp_index){
    x = unname(unlist(clust_list[tmp_index[1]]))
    y = unname(unlist(clust_list[tmp_index[2]]))
    z = unname(unlist(clust_list[tmp_index[3]]))
    zed = unname(unlist(clust_list[tmp_index[4]]))
    x_and_y = intersect(x, y)
    x_and_y_and_z = intersect(x_and_y, z)
    x_and_y_and_z_and_zed = intersect(x_and_y_and_z, zed)
    length(x_and_y_and_z_and_zed)
})



## venn diagram
require("VennDiagram")
?draw.quad.venn
venn.plot <- draw.quad.venn(
    area1 = areas[1],
    area2 = areas[2],
    area3 = areas[3],
    area4 = areas[4],
    n12 = pairs[1],
    n13 = pairs[2],
    n14 = pairs[3],
    n23 = pairs[4],
    n24 = pairs[5],
    n34 = pairs[6],
    n123 = trips[1],
    n124 = trips[2],
    n134 = trips[3],
    n234 = trips[4],
    n1234 = quads[1],
    category = c("Bladder cancer", "Melanoma", "Staphylococcus aureus infection", "Non-small cell lung cancer"),
    fill = c("orange", "red", "green", "blue"),
    lty = "dashed",
    cex = 2,
    cat.cex = 2,
    cat.col = c("orange", "red", "green", "blue")
)

# Writing to file
tiff(filename = "Worse_clusterings_Quad_Venn_diagram.tiff", compression = "lzw");
grid.draw(venn.plot);
dev.off();

worse_list <- list(unname(areas), pairs, trips, quads)


######################################################
## construct the "better" (Good) group Venn Diagram

clust_type <- "Good"
## tmp_clust <- all_data$cluster1
## names(tmp_clust) <- rownames(all_data)
clust_list <- lapply(all_data, function(tmp_clust) {
    rownames(all_data[tmp_clust == clust_type,])
})

## all patients in each cluster
(areas <- unlist(lapply(clust_list, length)))

## all pairwise comparisons
num_clusterings <- length(clust_list)
## tmp_index <- pairs_mat[,1]
pairs_mat <- combn(1:num_clusterings, 2)
pairs <- apply(pairs_mat, 2, function(tmp_index){
    x = unname(unlist(clust_list[tmp_index[1]]))
    y = unname(unlist(clust_list[tmp_index[2]]))
    length(intersect(x, y))
})

pairs_better <- pairs
pairs / sum(areas)

pairs_Jaccard_better <- apply(pairs_mat, 2, function(tmp_index){
    x = unname(unlist(clust_list[tmp_index[1]]))
    y = unname(unlist(clust_list[tmp_index[2]]))
    overlap <- length(intersect(x, y))
    union <- length(union(x, y))
    overlap / union
})


pairs_worse /
pairs_better
pairs_Jaccard_better
pairs_Jaccard_worse
summary(pairs_Jaccard_better)
summary(pairs_Jaccard_worse)

## all 3-way comparisons
## tmp_index <- trips_mat[,1]
trips_mat <- combn(1:num_clusterings, 3)
trips <- apply(trips_mat, 2, function(tmp_index){
    x = unname(unlist(clust_list[tmp_index[1]]))
    y = unname(unlist(clust_list[tmp_index[2]]))
    z = unname(unlist(clust_list[tmp_index[3]]))
    x_and_y = intersect(x, y)
    x_and_y_and_z = intersect(x_and_y, z)
    overlap <- length(x_and_y_and_z)
    union <- length(union(x, y))
    Jaccard <- overlap / union
    list(overlap, union, Jaccard)

})

## all 4-way comparisons
## tmp_index <- quad_mat[,1]
quad_mat <- combn(1:num_clusterings, 4)
quads <- apply(quad_mat, 2, function(tmp_index){
    x = unname(unlist(clust_list[tmp_index[1]]))
    y = unname(unlist(clust_list[tmp_index[2]]))
    z = unname(unlist(clust_list[tmp_index[3]]))
    zed = unname(unlist(clust_list[tmp_index[4]]))
    x_and_y = intersect(x, y)
    x_and_y_and_z = intersect(x_and_y, z)
    x_and_y_and_z_and_zed = intersect(x_and_y_and_z, zed)
    length(x_and_y_and_z_and_zed)
})



## venn diagram
require("VennDiagram")
?draw.quad.venn
venn.plot <- draw.quad.venn(
    area1 = areas[1],
    area2 = areas[2],
    area3 = areas[3],
    area4 = areas[4],
    n12 = pairs[1],
    n13 = pairs[2],
    n14 = pairs[3],
    n23 = pairs[4],
    n24 = pairs[5],
    n34 = pairs[6],
    n123 = trips[1],
    n124 = trips[2],
    n134 = trips[3],
    n234 = trips[4],
    n1234 = quads[1],
    category = c("Bladder cancer", "Melanoma", "Staphylococcus aureus infection", "Non-small cell lung cancer"),
    fill = c("orange", "red", "green", "blue"),
    lty = "dashed",
    cex = 2,
    cat.cex = 2,
    cat.col = c("orange", "red", "green", "blue")
)

# Writing to file
tiff(filename = "Better_clusterings_Quad_Venn_diagram.tiff", compression = "lzw");
grid.draw(venn.plot);
dev.off();

better_list <- list(unname(areas), pairs, trips, quads)

## compare and contrast two lists
unlist(lapply(worse_list, sum)) / sum(worse_list[[1]])
unlist(lapply(better_list, sum)) / sum(better_list[[1]])

unlist(lapply(worse_list, sum))
unlist(lapply(better_list, sum))

