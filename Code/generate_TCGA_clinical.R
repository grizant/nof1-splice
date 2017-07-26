## Use TCGA2STAT to retrieve clinical data from TCGA BRCA
## and save RData for all KEGG-related, paired data
## AG Schissler
## Created 11 Apr 2017
## Last modified 11 Apr 2017

#### 1. Load R library to retireve data (must have Internet connection)
library(TCGA2STAT)

#### 2. Manually list desired data sets
tcga_ids <- c("THCA", "COAD", "READ", "LUSC", "LUAD", "PRAD", "BLCA", "UCEC")

#### 3. Loop through and save paired clinical data nicely

## tmp_id <- tcga_ids[3]
for (tmp_id in tcga_ids) {

    tmp_data <- getTCGA(disease = tmp_id, data.type="RNASeq2", clinical = TRUE)

    ## Find patients with paired RNA-seq data
    paired_data <- TumorNormalMatch(tmp_data$dat)
    patients <- colnames(paired_data$normal)
    print(paste("Only", length(patients), tmp_id, "patients with survival data."))
    
    ## Subset data and isolate clinical data
    clin_data <- tmp_data$clinical[rownames(tmp_data$clinical) %in% patients, ]
    clin_data <- as.data.frame(clin_data, stringsAsFactors = F)

    ## 4. Reformat nicely only concerned with survival and age
    ## first explore clinical data
    ## str(clin_data)
    ## head(clin_data)
    ## names(clin_data)
    
    ## vitalstatus is a factor (check levels!)
    if ( all(clin_data[, "vitalstatus"] == 0) ) {
        (clin_data[,"vitalstatus"] <- factor(clin_data[,"vitalstatus"], labels = c("Alive")))
    } else {  
        if ( all(clin_data[, "vitalstatus"] == 1) ) {
            (clin_data[,"vitalstatus"] <- factor(clin_data[,"vitalstatus"], labels = c("Deceased")))
        } else {
            (clin_data[,"vitalstatus"] <- factor(clin_data[,"vitalstatus"], labels = c("Alive", "Deceased")))
        }
    }
    ## summary(clin_data[,2]) 

    ## daystodeath is numeric
    (clin_data[, "daystodeath"] <- as.numeric(clin_data[, "daystodeath"]))
    ## summary(clin_data[,3])
    ## clin_data[clin_data$vitalstatus == "Alive",3]
    ## clin_data[clin_data$vitalstatus == "Deceased",3]

    ## daystolastfollowup (censor date) is numeric 
    (clin_data[, "daystolastfollowup"] <- as.numeric(clin_data[, "daystolastfollowup"]))
    ## summary(clin_data[,4])

    ## save data frame as an R object
    print(paste("Saving TCGA", tmp_id, "paired_clinical.RData", sep = "_"))
    save(clin_data,
         file = paste("~/Dropbox/Splice-n-of-1-pathways/Data/TCGA", tmp_id, "paired_clinical.RData", sep = "_"))

}
