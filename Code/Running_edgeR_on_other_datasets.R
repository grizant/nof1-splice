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

#############Define the function that is the DEG analysis pathway #########

DEGAnalysis<-function(DataSet )
 {
  ## 1. split by gene symbol and remove iso/gene ID
  ## this could be done easier by tidyverse I think.
  gene_list <- split(DataSet[,-(1:2)], DataSet$geneSymbol)

   #Each gene has rows = isoforms, and columns = sample

  ## 2. summarize by adding isoform expression and rounding to approximate counts
  ## gene_dat <- gene_list[[1]]
  ## I prefer to round as a data pre-processing step rather than the DEG analysis
  ## in a modular way.
  gene_dat <- do.call("rbind", lapply(gene_list, function(gene_dat) {
      apply(gene_dat, 2, function(pat_dat) {
          round(sum(pat_dat),0)
      })
      #Adds the counts for each isoform and rounds for each sample for given gene
  }))
   #Applies the addition process to all the genes

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
  deg_list <- lapply(pat_list, get_ss_edgeR_preds)
 }
###############End of defining the DEG analysis pathway ##################

#Define your home directory
 HOME_DIR = 'C:/Users/daberasturi'

#1.Load each of the datasets into R
 #Load the LUSC dataset
 load(file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/lusc_iso_paired.RData'))

 #Load the LUAD dataset
 load(file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/luad_iso_paired.RData'))

 #Load the PRAD dataset
 load(file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/prad_iso_paired.RData'))

 #Load the THCA dataset
 load(file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/thca_iso_paired.RData'))

 #Load the UCEC dataset
 load(file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/ucec_iso_paired.RData'))

#2.Combine the datasets into one large list
 #Create the list object to hold the 5 datasets
 datasets_list<-vector('list',5)

 #Enter in one by one the datasets
  #Enter the lusc dataset
  datasets_list[[1]]<-lusc_iso

  #Enter the luad dataset
  datasets_list[[2]]<-luad_iso

  #Enter the prad dataset
  datasets_list[[3]]<-prad_iso

  #Enter the thca dataset
  datasets_list[[4]]<-thca_iso

  #Enter the ucec dataset
  datasets_list[[5]]<-ucec_iso
 
#3.Apply the DEG Analysis to each of the Datasets
  deg_lists<-lapply(datasets_list,DEGAnalysis)

  ##Note: this step takes somewhere between 5 and 10 minutes to finish running

#4.Store all the lists for easier processing in R
 #4.a)Create named versions of the 5 result lists
  ##Note: Apparently save won't accept parts of lists, only full list objects
          #so, this step is how I get around that

   #Create the lucs list
   lusc_deg_list_all<-deg_lists[[1]]

   #Create the luad list
   luad_deg_list_all<-deg_lists[[2]]

   #Create the prad list
   prad_deg_list_all<-deg_lists[[3]]

   #Create the thca list
   thca_deg_list_all<-deg_lists[[4]]

   #Create the ucec list
   ucec_deg_list_all<-deg_lists[[5]]

 #4.b)Store each individual deg list in properly named file
  #Store the lusc deg list
   save(lusc_deg_list_all, file=file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/lusc_deg_list_all.RData'))

  #Store the luad deg list
   save(luad_deg_list_all, file=file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/luad_deg_list_all.RData'))

  #Store the prad deg list
   save(prad_deg_list_all, file=file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/prad_deg_list_all.RData'))

  #Store the thca deg list
   save(thca_deg_list_all, file=file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/thca_deg_list_all.RData'))

  #Store the ucec deg list
   save(ucec_deg_list_all, file=file.path(HOME_DIR,'Dropbox/Splice-n-of-1-pathways/Data/ucec_deg_list_all.RData'))
