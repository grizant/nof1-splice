## Explore if getting isoform level TCGA is easy
## AG Schissler
## Created 20 Jul 2017

## use BLCA as a test case
######### use TCGA2STAT to easily find patient IDs
library(TCGA2STAT)

blca.rnaseq2 <- getTCGA(disease="BLCA", data.type="RNASeq2")
# tumor-normal matched profiles
blca.rnaseq2.tum.norm <- TumorNormalMatch(blca.rnaseq2$dat)
str(blca.rnaseq2.tum.norm)

patients <- colnames(blca.rnaseq2.tum.norm[[1]])

######### use Broad Firehose to retrieve

## https://github.com/mariodeng/FirebrowseR
## try using the firehouse to get isoform level
devtools::install_github("mariodeng/FirebrowseR")
library("FirebrowseR")

mRNA.Exp = Samples.mRNASeq(format = "csv",
                           gene = c("PTEN", "RUNX1"),
                           tcga_participant_barcode = c("TCGA-GF-A4EO",
                                                        "TCGA-AC-A2FG")
                           )
mRNA.Exp[, c("tcga_participant_barcode", "expression_log2", "z.score")]

### doesn't work for isoform

### try manual
library(data.table)
iso_data <- fread(file = "~/Data/TCGA/BLCA/gdac.broadinstitute.org_BLCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.Level_3.2016012800.0.0/BLCA_RSEM_isoforms_normalized_modified_to_read.txt", sep = "\t", header = T, stringsAsFactors = F)

str(iso_data)
head(iso_data$isoform_id)

writeLines(iso_data$isoform_id, "~/Data/TCGA/BLCA/gdac.broadinstitute.org_BLCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.Level_3.2016012800.0.0/ucsc_id.txt", sep = '\n')



## retrieve gene
## library(biomaRt)
## 
## GENES = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")   
## 
## getBM(attributes = c('ensembl_gene_id','hgnc_symbol','ensembl_transcript_id','refseq_mrna','ucsc','chromosome_name','transcript_start','transcript_end'), mart = GENES)
## 
## my_map <- getBM(attributes = c('hgnc_symbol','ucsc'), mart = GENES)
## head(my_map[my_map$ucsc != "",])
## my_map[my_map$ucsc == "uc011lsn.1",]
## my_map[my_map$ucsc == "uc011lsn.1",]
## 
## sum(my_map$ucsc %in% iso_data$isoform_id)
## head(my_map[my_map$ucsc %in% iso_data$isoform_id,])
