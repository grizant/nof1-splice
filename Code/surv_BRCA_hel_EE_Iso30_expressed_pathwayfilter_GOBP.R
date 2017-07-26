## Develop unsupervised survival prediction evaluation
## AG Schissler
## Created 25 Jul 2017

##############################################################################
#### 1. Load clinical data
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_paired_clinical.RData")

## included patients via MD criteria
## clin_data <- clin_data[rownames(clin_data) %in% patients,]

## rename patients to provide consistency
## rownames(clin_data) <- gsub("-", ".", rownames(clin_data))

##############################################################################
#### 2. Load pathway scores and update clinical information
load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_hel_EE_Iso30_expressiod_pathwayfilter_GOBP_25july2017.RData")

## patients <- readLines("~/Dropbox/N1PS-MahDist-Paper/Data/excludePatients.txt")
patients <- names(scores_list)
## patients <- readLines("~/Dropbox/N1PS-MahDist-Paper/Data/excludePatients.txt")
clin_data <- clin_data[rownames(clin_data) %in% patients,]

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

#### New code after iso/pathway filtering
## understand NA distribution
missing_logic <- apply(effect_mat, 2, function(tmp_pathway) {
    any(is.na(tmp_pathway))
})
table(missing_logic) ## 1013 out of 3485 pathways to remove

fil_effect_mat <- effect_mat[, !missing_logic]

str(fil_effect_mat)

## retrieve fdr
fdr_list <- lapply(scores_list, function(tmp_data){
    ## retrieve effect size
    tmp_scores <- tmp_data$fdr_value
    ## add pathway names
    names(tmp_scores) <- rownames(tmp_data)
    ## reorder to allow for aggregation
    tmp_scores[order(names(tmp_scores))]
})

## aggregate into a matrix
fdr_mat <- do.call("rbind", fdr_list)

#### New code after iso/pathway filtering
## understand NA distribution
missing_logic <- apply(fdr_mat, 2, function(tmp_pathway) {
    any(is.na(tmp_pathway))
})
table(missing_logic) ## only 498 pathways to remove

fil_fdr_mat <- fdr_mat[, !missing_logic]

str(fil_fdr_mat)

######### 3.2 start clustering attempts
## cluster the patients (do this more sophicately in future study)
tmp_k <- 2
tmp_clustering <- cluster::pam(x = fil_effect_mat, k = tmp_k)
## tmp_clustering <- cluster::pam(x = fil_effect_mat[,sample(1:ncol(fil_effect_mat),1)], k = tmp_k)
table(tmp_clustering$clustering)

## remove any patients that are not in the clinical data
my_clusters <- tmp_clustering$clustering[names(tmp_clustering$clustering) %in% rownames(clin_data)]

## add cluster to clinical
clin_data$cluster <- 0
clin_data[names(my_clusters), "cluster"] <- my_clusters
table(clin_data$cluster)

## cluster by fdrthe patients (do this more sophicately in future study)
tmp_k <- 2
## tmp_k <- 3
tmp_clustering <- cluster::pam(x = fil_fdr_mat, k = tmp_k)
table(tmp_clustering$clustering)

## remove any patients that are not in the clinical data
my_clusters <- tmp_clustering$clustering[names(tmp_clustering$clustering) %in% rownames(clin_data)]

## add cluster to clinical
clin_data$cluster_fdr <- 0
clin_data[names(my_clusters), "cluster_fdr"] <- my_clusters
table(clin_data$cluster_fdr)

## try spectral clustering on effect size
min_nj <- 4
max_centers <- ceiling(nrow(fil_effect_mat) / min_nj)
neighbors  <- ceiling(log(nrow(fil_effect_mat)) + 1)

## Compute eigenvalues for max number of centers to use eignegap heuritic
sc_max <- kknn::specClust(data = fil_effect_mat, centers = max_centers, nn = neighbors, method = "random-walk")

qplot(factor(1:length(sc_max$eigenvalue)), sc_max$eigenvalue)
eigen_diff <- diff(sc_max$eigenvalue)[-1]
diff_data <- data.frame(m = 2:(max_centers-1), eigen_diff)
## find the largest gap with min_nj met (or take m = 1)
diff_data <- diff_data[order(diff_data$eigen_diff, decreasing = T), ]
## add terminal clustering
diff_data <- rbind(diff_data, c(m = 1, 0))
(best_m <- diff_data$m[1])
## best_m <- 3

sc_best <- kknn::specClust(data = fil_effect_mat, centers = best_m, nn = neighbors, method = "random-walk")

names(sc_best$cluster) <- rownames(fil_effect_mat)

## remove any patients that are not in the clinical data
my_clusters <- sc_best$cluster[names(sc_best$cluster) %in% rownames(clin_data)]

## add cluster to clinical
clin_data$cluster_effect_sc <- 0
clin_data[names(my_clusters), "cluster_effect_sc"] <- my_clusters
table(clin_data$cluster_effect_sc)

##########################################
## try spectral clustering on fdr size
min_nj <- 4
max_centers <- ceiling(nrow(fil_fdr_mat) / min_nj)
neighbors  <- ceiling(log(nrow(fil_fdr_mat)) + 1)

## Compute eigenvalues for max number of centers to use eignegap heuritic
sc_max <- kknn::specClust(data = fil_fdr_mat, centers = max_centers, nn = neighbors, method = "random-walk")

qplot(factor(1:length(sc_max$eigenvalue)), sc_max$eigenvalue)
eigen_diff <- diff(sc_max$eigenvalue)[-1]
diff_data <- data.frame(m = 2:(max_centers-1), eigen_diff)
## find the largest gap with min_nj met (or take m = 1)
diff_data <- diff_data[order(diff_data$eigen_diff, decreasing = T), ]
## add terminal clustering
diff_data <- rbind(diff_data, c(m = 1, 0))
(best_m <- diff_data$m[1])
## best_m <- 3

sc_best <- kknn::specClust(data = fil_fdr_mat, centers = best_m, nn = neighbors, method = "random-walk")

names(sc_best$cluster) <- rownames(fil_fdr_mat)

## remove any patients that are not in the clinical data
my_clusters <- sc_best$cluster[names(sc_best$cluster) %in% rownames(clin_data)]

## add cluster to clinical
clin_data$cluster_fdr_sc <- 0
clin_data[names(my_clusters), "cluster_fdr_sc"] <- my_clusters
table(clin_data$cluster_fdr_sc)


##############################################################################
#### 4. Explore Kaplan-Meier curves

## https://rawgit.com/steadyfish/JustAnotherDataBlog/master/Code/PlayingAroundwithSurvivalAnalysis.html

## load ggplot and broom for nice plots
library(ggplot2)
library(broom)
library(survival)

##### 4.1 provide censorship info in the way the survival package expects
clin_data$status <- ifelse(clin_data$vitalstatus == "Alive", 0, 1)
table(clin_data$status)

#### 4.2 provide survival times in the way the survival package expects
## first retrieve the censorship times
clin_data$time <- clin_data$daystolastfollowup
## then insert the survival time
clin_data$time[is.na(clin_data$time)] <- clin_data$daystodeath[!is.na(clin_data$daystodeath)]

clin_surv <- survival::Surv(time = clin_data$time, event = clin_data$status)
head(clin_data[,c(2,3,4)], 10)
head(clin_surv, 10)

## modified from ggsurv exampl
sf_survfit <- survival::survfit(clin_surv ~ cluster, data = clin_data)
## summary(sf_survfit)

## better plot through ggplot2 and tidyr
sf_tidy = broom::tidy(sf_survfit)
glance(sf_tidy)
mx = max(sf_tidy$n.censor)
ggplot(sf_tidy, aes(time, estimate, fill = strata)) + 
    geom_line() +
    geom_point(aes(shape = as.factor(n.censor)), size = 3) + 
    scale_shape_manual(values=c(NA, 1:mx))+
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.25) + 
    xlab("days") + 
    ylab("Proportion Survival")

## log-rank test
survival::survdiff(Surv(time, status) ~ cluster, data = clin_data, rho = 0)

## try fdr-based clusters
sf_survfit <- survival::survfit(Surv(time, status) ~ cluster_fdr, data = clin_data)
## summary(sf_survfit)
sf_tidy = tidy(sf_survfit)
broom::glance(sf_tidy)
mx = max(sf_tidy$n.censor)
ggplot(sf_tidy, aes(time, estimate, fill = strata)) + 
  geom_line() +
  geom_point(aes(shape = as.factor(n.censor)), size = 3) + 
  scale_shape_manual(values=c(NA, 1:mx))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.25) + 
  xlab("days") + 
  ylab("Proportion Survival")

## log rank test
survival::survdiff(Surv(time, status) ~ cluster_fdr, data = clin_data, rho = 0)

##############################
## try sc effect size clusters
sf_survfit <- survival::survfit(Surv(time, status) ~ cluster_effect_sc, data = clin_data)
## summary(sf_survfit)
sf_tidy = tidy(sf_survfit)
glance(sf_tidy)
mx = max(sf_tidy$n.censor)
ggplot(sf_tidy, aes(time, estimate, fill = strata)) + 
  geom_line() +
  geom_point(aes(shape = as.factor(n.censor)), size = 3) + 
  scale_shape_manual(values=c(NA, 1:mx))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.25) + 
  xlab("days") + 
  ylab("Proportion Survival")

## log-rank test
survival::survdiff(Surv(time, status) ~ cluster_effect_sc, data = clin_data, rho = 0)

##############################
## try sc fdr size clusters
sf_survfit <- survival::survfit(Surv(time, status) ~ cluster_fdr_sc, data = clin_data)
## summary(sf_survfit)
sf_tidy = tidy(sf_survfit)
glance(sf_tidy)
mx = max(sf_tidy$n.censor)
ggplot(sf_tidy, aes(time, estimate, fill = strata)) + 
  geom_line() +
  geom_point(aes(shape = as.factor(n.censor)), size = 3) + 
  scale_shape_manual(values=c(NA, 1:mx))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.25) + 
  xlab("days") + 
  ylab("Proportion Survival")

## log-rank test
survival::survdiff(Surv(time, status) ~ cluster_fdr_sc, data = clin_data, rho = 0)

##############################################################################
#### 4. Explore Cox proportional hazards

## https://www.r-bloggers.com/survival-analysis-2/
## https://stats.stackexchange.com/questions/23042/interpretation-and-validation-of-a-cox-proportional-hazards-regression-model-usi
## http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Cox-Regression.pdf

# COX PH
sf1_coxph = coxph(clin_surv ~ as.factor(cluster_fdr_sc), data = clin_data)
sf1_coxph_tidy = tidy(sf1_coxph)
sf1_coxph_tidy # equivalent of print()
glance(sf1_coxph) # equivalent of summary()
# Validating Cox PH Assumptions
validate_sf1_coxph = cox.zph(sf1_coxph, transform = "km")
validate_sf1_coxph
plot(validate_sf1_coxph) ## should lie on a straight line
abline(h=0)
validate_gg = data.frame(cbind(x = validate_sf1_coxph[["x"]],
                               y = validate_sf1_coxph[["y"]]))
names(validate_gg) = c("transformed_t", "scaled_schoenfeld_residual")

ggplot(data = validate_gg, aes(x = transformed_t, y = scaled_schoenfeld_residual)) +
    geom_point() +
    geom_smooth(method = "loess") +
    abline(h = 0)


## cluster by effect size
sf_survfit <- survival::survfit(clin_surv ~ cluster, data = clin_data)
sf_coxph <- survival::coxph(clin_surv ~ factor(cluster), data = clin_data)

##############################################################################
#### 5. Try PCA

## remove patient without clinical data
to_keep <- rownames(fil_effect_mat)[(rownames(fil_effect_mat) %in% rownames(clin_data))]

## 5.1 pca on effect size (hellinger distance)
effect_pca <- prcomp(fil_effect_mat[to_keep,])
effect_pca_mat <- effect_pca$x

## add clinical data

effect_pca_dat <- data.frame(cbind(effect_pca_mat, clin_data$vitalstatus))
names(effect_pca_dat) <- gsub(rev(names(effect_pca_dat))[1], "vital", names(effect_pca_dat))

effect_pca_dat <- data.frame(cbind(effect_pca_dat, clin_data$daystodeath))
names(effect_pca_dat) <- gsub(rev(names(effect_pca_dat))[1], "Survival", names(effect_pca_dat))


## visualize first two PCs
qplot(x = PC1, y = PC2, color = factor(vital), data = effect_pca_dat)
qplot(x = PC2, y = PC3, color = factor(vital), data = effect_pca_dat)
qplot(x = PC2, y = PC3, color = Survival, data = effect_pca_dat)


##### 5.2 pca on fdr value
fdr_pca <- prcomp(t(fdr_mat))
