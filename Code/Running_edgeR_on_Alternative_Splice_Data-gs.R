#############Alternative Splicing Script (Running edgeR) #######################

##Source the functions to run the other libraries (if you need to)
 ##source("https://bioconductor.org/biocLite.R")
 ##biocLite("edgeR")

##Run the edgeR library
#### library(edgeR)

#####Define the edgeR function ##########

require(edgeR)

##Note: for common.dispersion this is your choice. Originally it was
##1e-02, but manual recommends .4 for human, .1 for genetically identical,
##and .01 for technical replicates

## testing parameters
## data.matrix <-
## common.dispersion=.1

get_ss_edgeR_preds <- function(data.matrix, common.dispersion=.1 ){

    ##Naming the conditions of the samples
    conds <- c('untreated','treated')
    ##Round the data to make whole counts
    roundSamples<-round(data.matrix)

    ##Transforming the counts into a DGEList object so that edgeR can work
    cds = DGEList(roundSamples, group= conds)

    ##Normalize the RNA composition signals to account for differing composition of
    ##RNA for different samples. It minimizes the log2FC between genes from different
    ##samples by using trimmed mean of M values between each pair of samples
    cds <- calcNormFactors( cds )

    ##Estimates common dispersion factor 
    cds <- estimateCommonDisp( cds )
    ##Note: Since there aren't replicates this returns NA
    ## this is a hack to introduce a common dispersion
    
    cds$common.dispersion=common.dispersion
    ##Perform the exact test with the inputted dispersion parameter
    de.cmn <- exactTest( cds , dispersion = 'common' ,
                        pair = c( "untreated" , "treated" ) )
    ##Here the code both grabs the logFC and PValue counts and puts them into a 
    ##matrix AND it rounds those numbers to the 5th decimal place
    de.cmn <- round(de.cmn$table,5)
    ##Adjust the p-values by in this case 'BY' method
    de.cmn$BY <- p.adjust(de.cmn$PValue, method = 'BY')
    ##Add id names to these results in the form of a column
    de.cmn$id <- row.names(de.cmn)

    ## Return the results for differentially expressed gene analysis
    de.cmn
}

##########End define the edgeR function ############

##Load the dataset use the paired dataset before the tpm
#### blca_iso<-read.table('blca_iso.txt',header=TRUE)
#### blca_iso<-read.table('blca_iso_paired.txt',header=TRUE)

## Note blca stands for 'bladder cancer' so this is bladder cancer dataset
## Use the RData in the common cloud storage for speed and reproducibility
load('~/Dropbox/Splice-n-of-1-pathways/Data/blca_iso_paired.RData')

## 1. split by gene symbol and remove iso/gene ID
## this could be done easier by tidyverse I think.
gene_list <- split(blca_iso[,-(1:2)], blca_iso$geneSymbol)

## 2. summarize by adding isoform expression and rounding to approximate counts
## gene_dat <- gene_list[[1]]
## I prefer to round as a data pre-processing step rather than the DEG analysis
## in a modular way.
gene_dat <- do.call("rbind", lapply(gene_list, function(gene_dat) {
    apply(gene_dat, 2, function(pat_dat) {
        round(sum(pat_dat),0)
    })
}))

## 3. now transpose and add patient ID to split on
gene_dat <- as.data.frame(t(gene_dat))
gene_dat$PatID <- substr( x=rownames(gene_dat), 1, 12 )

## 4. Now DEG analysis by patient
pat_list <- split(gene_dat[,-ncol(gene_dat)], gene_dat$PatID)

## format as the edgeR function expects
## pat_dat <- pat_list[[1]]
pat_list <- lapply(pat_list, function(pat_dat) {
    as.matrix(t(pat_dat))
})

## format as the edgeR function expects
## time to see if parallelization is required
## if other datasets take no longer, see ~/Dropbox/Splice-n-of-1-pathways/Code/score_TCGA_UCEC_hel_Avg_EE_Iso30_expressed_pathwayfilter_KEGG.R as an example to parallelize. 
## pat_dat <- pat_list[[1]]
system.time(deg_list <- lapply(pat_list, get_ss_edgeR_preds))
blca_deg_list <- deg_list

## 5. Store entire list for easier processing in R
save(blca_deg_list, file="~/Dropbox/Splice-n-of-1-pathways/Data/blca_deg_list_all.RData")

## check on fresh R session
load(file="~/Dropbox/Splice-n-of-1-pathways/Data/blca_deg_list_all.RData")
