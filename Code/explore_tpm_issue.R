## Figure out if data are tpm or normalized or what
## AG Schissler
## 27 Jul 2018

load("~/Dropbox/Splice-n-of-1-pathways/Data/blca_iso_paired.RData")
raw_data <- blca_iso
str(raw_data)
max(raw_data[,3])
raw_sum <- colSums(raw_data[,-c(1,2)])
summary(raw_sum)

threshold <- 5
## tmp_pat <- raw_data[,3]
apply(raw_data[,-c(1,2)], 2, function(tmp_pat){
    sum(tmp_pat <= threshold) / length(tmp_pat)
})
## this remove far too many genes
      
sum(dat$Normal <= log_threshold) / nrow(dat)
sum(dat$Tumor <= log_threshold) / nrow(dat)

## not log2 or tpm (this is good)


load("~/Dropbox/Splice-n-of-1-pathways/Data/blca_iso_paired_tpm.RData")
tpm_data <- blca_iso
str(tpm_data)
max(tpm_data[,3]) 
tpm_sum <- colSums(tpm_data[,-c(1,2)])
summary(tpm_sum)
## not log2 fyi
## not tpm (this is bad...)


load(file = "~/Dropbox/Splice-n-of-1-pathways/Data/blca_iso_kegg_data.RData")
kegg_data <- blca_iso_kegg_data
str(kegg_data) ## filtered
max(kegg_data[,3])
kegg_sum <- colSums(kegg_data[,-c(1,2)])
summary(kegg_sum)
## not log2 fyi
## not tpm (this is bad...)
