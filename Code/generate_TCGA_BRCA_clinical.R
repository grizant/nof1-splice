## Use TCGA2STAT to retrieve clinical data from TCGA BRCA
## and save RData
## AG Schissler
## Created 11 Apr 2017
## Last modified 11 Apr 2017

#### 1. Load library and retireve data (must have Internet connection)
library(TCGA2STAT)
brca_data <- getTCGA(disease="BRCA", data.type="RNASeq2", clinical = TRUE)

#### 2. Subset to patients with paired RNA-seq data and isolate clinical data
paired_data <- TumorNormalMatch(brca_data$dat)
patients <- colnames(paired_data$normal)
clin_data <- brca_data$clinical[rownames(brca_data$clinical) %in% patients, ]
clin_data <- as.data.frame(clin_data, stringsAsFactors = F)

#### 3. Reformat nicely
## first explore clinical data
str(clin_data)
head(clin_data)
names(clin_data)

## var 1 is useless
clin_data <- clin_data[, -1]

## yearstobirth (former var 2) is numeric
clin_data[,1] <- as.numeric(clin_data[,1])
summary(clin_data[,1])

## vitalstatus is a factor (check levels!)
clin_data[,2] <- factor(clin_data[,2], labels = c("Alive", "Deceased"))
summary(clin_data[,2]) ## 68 alive and 44 deceased

## daystodeath is numeric
clin_data[,3] <- as.numeric(clin_data[,3])
summary(clin_data[,3])
clin_data[clin_data$vitalstatus == "Alive",3]
clin_data[clin_data$vitalstatus == "Deceased",3]

## daystolastfollowup (censor date) is numeric 
clin_data[,4] <- as.numeric(clin_data[,4])
summary(clin_data[,4])

## tumortissuesite is a factor (all breast though)
clin_data[,5] <- factor(clin_data[,5])

## pathologicstage is a factor
clin_data[,6] <- factor(clin_data[,6])

## pathologyTstage is a factor
clin_data[,7] <- factor(clin_data[,7])

## pathologyNstage is a factor
clin_data[,8] <- factor(clin_data[,8])

## pathologyMstage is a factor
clin_data[,9] <- factor(clin_data[,9])

## gender is a factor (only 1 male)
clin_data[,10] <- factor(clin_data[,10])

## treat dateofinitialpathologicdiagnosis as a factor
clin_data[,11] <- factor(clin_data[,11])

## daystolastknownalive is empty
clin_data <- clin_data[,-12]

## radiationtherapy is a factor
clin_data[,12] <- factor(clin_data[,12])
table(clin_data[,12], useNA = "ifany")

## histologicaltype is a factor
clin_data[,13] <- factor(clin_data[,13])
table(clin_data[,13], useNA = "ifany")

## histologicaltype is a factor
clin_data[,13] <- factor(clin_data[,13])
table(clin_data[,13], useNA = "ifany")

## numberoflymphnodes is numeric
clin_data[,14] <- as.numeric(clin_data[,14])
summary(clin_data[,14])

## race is a factor
clin_data[,15] <- factor(clin_data[,15])
table(clin_data[,15], useNA = "ifany")

## ethnicity is a factor
clin_data[,16] <- factor(clin_data[,16])
table(clin_data[,16], useNA = "ifany")

## save data frame as an R object
save(clin_data, file = "~/Dropbox/Splice-n-of-1-pathways/Data/TCGA_BRCA_paired_clinical.RData")
