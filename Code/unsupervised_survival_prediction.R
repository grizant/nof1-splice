## Develop unsupervised survival prediction evaluation
## AG Schissler
## Created 12 Apr 2017
## Last modified 14 Apr 2017

##############################################################################
#### 1. Load clinical data
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_paired_clinical.RData")

##############################################################################
#### 2. Load pathway scores and update clinical information
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_avg_14apr2017.RData")

## rename patients to provide consistency
rownames(clin_data) <- gsub("-", ".", rownames(clin_data))

##############################################################################
#### 3. Cluster via PAM

## structure as a matrix for clustering
## tmp_data <- scores_list[[1]]

effect_list <- lapply(scores_list, function(tmp_data){
    ## retrieve effect size
    tmp_scores <- tmp_data$pathway_score
    ## add pathway names
    names(tmp_scores) <- rownames(tmp_data)
    ## reorder to allow for aggregation
    tmp_scores[order(names(tmp_scores))]
})

## aggregate into a matrix
effect_mat <- do.call("rbind", effect_list)

## cluster the patients (do this more sophicately in future study)
tmp_k <- 2
tmp_clustering <- cluster::pam(x = effect_mat, k = tmp_k)

table(tmp_clustering$clustering)

## remove any patients that are not in the clinical data
my_clusters <- tmp_clustering$clustering[names(tmp_clustering$clustering) %in% rownames(clin_data)]

## add cluster to clinical
clin_data$cluster <- 0
clin_data[names(my_clusters), "cluster"] <- my_clusters
table(clin_data$cluster)

##############################################################################
#### 4. Explore Kaplan-Meier curves

## load GGally for nice plots
library(GGally)

##### 4.1 provide censorship info in the way the survival package expects
clin_data$status <- ifelse(clin_data$vitalstatus == "Alive", 1, 2)
table(clin_data$status)

#### 4.2 provide survival times in the way the survival package expects
## first retrieve the censorship times
clin_data$time <- clin_data$daystolastfollowup
## then insert the survival time
clin_data$time[is.na(clin_data$time)] <- clin_data$daystodeath[!is.na(clin_data$daystodeath)]

## modified from ggsurv exampl
sf_fit <- survival::survfit(Surv(time, status) ~ cluster, data = clin_data)
summary(sf_fit)
ggsurv(sf_fit)
