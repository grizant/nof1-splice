#############Alternative Splicing Script (Running edgeR) #######################

#Load all necessary libraries
library('locfdr')
library(stringi)

#Source the functions to run the other libraries (if you need to)
 #source("https://bioconductor.org/biocLite.R")
 #biocLite("edgeR")

#Run the edgeR library
 library(edgeR)

#####Define the edgeR function ##########

require(edgeR)

#Note: for common.dispersion this is your choice. Originally it was
#1e-02, but manual recommends .4 for human, .1 for genetically identical,
#and .01 for technical replicates

get_ss_edgeR_preds <- function(data.matrix, common.dispersion=.1 ){
 #Naming the conditions of the samples
  conds <- c('untreated','treated')

 #Round the data to make whole counts
 roundSamples<-round(data.matrix)

 #Transforming the counts into a DGEList object so that edgeR can work
 cds = DGEList(roundSamples, group= conds)

  #Normalize the RNA composition signals to account for differing composition of
   #RNA for different samples. It minimizes the log2FC between genes from different
   #samples by using trimmed mean of M values between each pair of samples
  cds <- calcNormFactors( cds )

  #Estimates common dispersion factor 
  cds <- estimateCommonDisp( cds )

   #Note: Since there aren't replicates this returns NA
  
  ### this is a hack to introduce a common dispersion
  cds$common.dispersion=common.dispersion
  
  #Perform the exact test with the inputted dispersion parameter
  de.cmn <- exactTest( cds , dispersion = 'common' , pair = c( "untreated" , "treated" ) )

  #Here the code both grabs the logFC and PValue counts and puts them into a 
   #matrix AND it rounds those numbers to the 5th decimal place
  de.cmn <- round(de.cmn$table,5)

  #Adjust the p-values by in this case 'BY' method
  de.cmn$PValue <- p.adjust(de.cmn$PValue, method = 'BY')
  
  #Add id names to these results in the form of a column
  de.cmn$id <- row.names(de.cmn)

  #Return the results for differentially expressed gene analysis
  de.cmn
}

#####End define the edgeR function ######

#Load the dataset use the paired dataset before the tpm
#blca_iso<-read.table('blca_iso.txt',header=TRUE)
blca_iso<-read.table('blca_iso_paired.txt',header=TRUE)

  #Note blca stands for 'bladder cancer' so this is bladder cancer dataset

#Load the KEGG data
kegg <- read.delim2("kegg_tb.txt")

#Filter to include only those genes with KEGG annotation
 #Find all the unique genes corresponding to the KEGG terms
 kegg_genes <- unique(kegg$symbol)

 #Show the number of the KEGG genes
 length(kegg_genes) ## 5879 genes (It matches :) )

 #Find the intersection of the blca genes and the KEGG annotated genes
 blca_iso_kegg_data <- blca_iso[blca_iso$geneSymbol %in% kegg_genes,]

 #Show the number of genes in the intersection
 nrow(blca_iso_kegg_data) ##This step left 18823 rows

 #Find the number of unique genes in the intersection
 kegg_measured <- length(unique(blca_iso_kegg_data$geneSymbol)) ##Equals 5757

  ## free some memory
  rm(blca_iso)

## Filter to include only those genes with between 2 and 30 isoforms
  #Set the desired isoform range
  iso_range <- c(2,30)

  #Create a list of isoforms for each gene symbol
  iso_list <- split(blca_iso_kegg_data$isoform_id, factor(blca_iso_kegg_data$geneSymbol))

  #Read the first gene's isoforms into tmp_gene
  tmp_gene <- iso_list[[1]]

  #Determine which genes have correct number of isoforms
  to_keep <- unlist(lapply(iso_list, function(tmp_gene) {
      #Find the number of isoforms in the gene
      tmp_num <- length(tmp_gene)
      
      #Return whether or not the number of isoforms is within the desired range
      tmp_num >= iso_range[1] & tmp_num <= iso_range[2]
   }))
   #End of function calls for unlist, lapply, and function()

  #Find the number of genes outside the desired isoform range 
  kegg_measured - sum(!to_keep) ## filter 1624 genes

  #Find the number of genes inside of the desired isoform range
  num_kept<-sum(to_keep) ##4133 genes leftover from filtering

  #Obtain the intersection between those genes annotated in KEGG and those with correct number of isoforms
  filtered_blca_iso<-blca_iso_kegg_data[blca_iso_kegg_data$geneSymbol %in% names(which(to_keep)),] ## 17088

   #Find the number of these genes 
   nrow(filtered_blca_iso)

#Obtain the counts for the genes 
 #Find all the unique genes
 uniqGenes<-unique(filtered_blca_iso$geneSymbol) 
 
 #Create a matrix that can contain all the counts for each gene for each patient
 geneCounts<-matrix(0,nrow=length(uniqGenes),ncol=38)

  #Make the rownames of geneCounts the names of the unique genes
  rownames(geneCounts)<-unique(filtered_blca_iso$geneSymbol)

  #Make the colnames of geneCounts the names of the samples
  colnames(geneCounts)<-colnames(filtered_blca_iso)[c(-1,-2)]

 #Prep the data so that only the counts are left
  #Get rid of the list of isoforms in blca_iso
  Temp<-filtered_blca_iso[,-1]

  #Get rid of the list of geneSymbols
  Data<-Temp[,-1]

  #Free up some of the memory
  rm(Temp)

 #Loop through all the unique genes
 for(i in 1:length(uniqGenes) )
 {
  #Find the rows which correspond to the unique gene
  aUniqGene<-Data[(uniqGenes[i] == filtered_blca_iso$geneSymbol),]

  #Sum up the column to obtain the vector of counts and input it into geneCounts
  geneCounts[i,]<-colSums(aUniqGene)
 }
 #End of loop through all the unique genes 

#Perform edgeR for each comparison and save the results
 #Find the number of comparisons to be made
 numCompare<-dim(geneCounts)[2]/2

 #Create an array to hold results for each comparison
 comparisons<-vector('list',numCompare)

  #For extra info find the number of DEGs per comparison
  numDEGs<-matrix(0,nrow=1,ncol=numCompare)

  #Loop through each comparison
  for( j in 1:numCompare )
  {
   #Run edgeR on the comparison and read into a column of the matrix
   comparisons[[j]] <- get_ss_edgeR_preds(cbind(geneCounts[,(2*j-1)], geneCounts[,(2*j)]))

   #Input the number of DEGs for the comparison for alpha = .1
   numDEGs[1,j]<-sum(comparisons[[j]]$PValue < .1)
  }
  #End of loop through each comparison

#Output the results
 #Output a space to open up the file and prime the loop 
 write(' ', file = 'blca_output.txt')

 #Loop through all the comparisons
 for(i in 1:numCompare)
 {
  #Write that this is a comparison
  write('This is data from a single comparison',file='blca_output.txt',append=TRUE)
 
  #Write a space
  write(' ',file='blca_output.txt',append=TRUE)

  #Write the results of the edgeR comparison to a file
  write.table(comparisons[[i]],file='blca_output.txt',append=TRUE,sep=' ',row.names=FALSE,col.names=TRUE)

  #Write a space
  write(' ',file='blca_output.txt',append=TRUE)
 }
 #End of loop through all the comparisons

 #Output the number of DEGs for all the comparisons
  #Output the comment explaining what the next line is
  write('This is number of DEGs for each comparison ',file='blca_output.txt',append=TRUE)
 
  #Output the actual numbers
    #Write the results of the edgeR comparison to a file
    write.table(numDEGs,file='blca_output.txt',append=TRUE,sep=' ',row.names=FALSE,col.names=FALSE)



