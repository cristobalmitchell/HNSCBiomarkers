# Author: Cristobal Mitchell
# Course: STAT 685 700 - Dr. Moumita Karmakar
# Title: TCGA-HNSC Preprocess

# Set seed to reproduce results
set.seed(1618)

#### IMPORTS ####
# install.packages("UCSCXenaTools") # https://docs.ropensci.org/UCSCXenaTools/
# install.packages("TCGA2STAT_1.2.tar", repos = NULL, type = "source") # http://www.liuzlab.org/TCGA2STAT/ 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefilter", update=FALSE, ask=FALSE) # https://bioconductor.org/packages/release/bioc/html/genefilter.html
BiocManager::install("CNTools", update=FALSE, ask=FALSE) # http://bioconductor.org/packages/release/bioc/html/CNTools.html
BiocManager::install("miRBaseConverter", update=FALSE, ask=FALSE) # http://bioconductor.org/packages/release/bioc/html/miRBaseConverter.html

library(UCSCXenaTools)
library(ggplot2)
library(TCGA2STAT)
library(miRBaseConverter)
library(missForest)
library(glmnet)
library(tidyverse) 


#### LOAD DATA ####
# Read in data from files
clinical = read.table('HNSC_clinicalMatrix.txt', sep='\t', header=TRUE)
miRNA = read.table('miRNA_HiSeq_gene.txt', sep='\t')
survival = read.table('HNSC_survival.txt', sep='\t', header=TRUE)

# Read in data from Xena API package
# data(XenaData)
# 
# XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
#   XenaFilter(filterDatasets = "clinical") %>% 
#   XenaFilter(filterDatasets = "HNSC") -> df_todo
# XenaQuery(df_todo) %>%
#   XenaDownload() -> xe_download
# clinical = XenaPrepare(xe_download)
# 
# XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
#   XenaFilter(filterDatasets = "survival") %>% 
#   XenaFilter(filterDatasets = "HNSC") -> df_todo
# XenaQuery(df_todo) %>%
#   XenaDownload() -> xe_download
# survival = XenaPrepare(xe_download)
# 
# XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
#   XenaFilter(filterDatasets = "miRNA_HiSeq_gene") %>% 
#   XenaFilter(filterDatasets = "HNSC") -> df_todo
# XenaQuery(df_todo) %>%
#   XenaDownload() -> xe_download
# miRNA = XenaPrepare(xe_download)

#### GATHER NAMING ####
# Gather miRNA names from the accession values found in the first column of miRNA dataset
miRNAAccessions = miRNA$V1
miRNANames = miRNA_AccessionToName(miRNAAccessions,targetVersion = "v21") # v22 results in some NA target names
miRNANames[1,2] = miRNANames[1,1]

# Replace the column names
miRNA$V1 = miRNANames$TargetName

# Transpose miRNA, methyl, and HiSeq datasets so each row represent a patient
miRNA = setNames(data.frame(t(miRNA[,-1])), miRNA[,1])

# Convert MIMA columns to correct datatype for miRNA dataset
miRNA[,-1] = sapply(miRNA[,-1], as.character)
miRNA[,-1] = sapply(miRNA[,-1], as.numeric)

#### FILTER DATA ####
# Remove patients from clinical dataset where lymphovascular_invasion_present is NA
clinical = clinical[clinical$lymphovascular_invasion_present!='',]
# Remove patients who have days_to_last_followup < 30
clinical = clinical[as.numeric(clinical$days_to_last_followup)>=30,]
# Remove observations not shared with miRNA data.frame
clinical = clinical[clinical$sampleID %in% miRNA$sample,]
# Remove observations not shared with clinical data.frame
miRNA = miRNA[miRNA$sample %in% clinical$sampleID,]
# Order data.frame 
miRNA = miRNA[order(miRNA$sample),]
clinical = clinical[order(clinical$sampleID),]
# Reset row names
rownames(miRNA) = rownames(clinical) = NULL
# Remove columns with > 30% missing values from miRNA data.frame
miRNA = miRNA[,colSums(is.na(miRNA))<=(dim(miRNA)[1]*.3)]


#### LOAD MATCHED DATA #### 
# Download the tumor normal matched data for the HNSC cancer
hnsc.miRNAseq = getTCGA(disease="HNSC", data.type="miRNASeq")
hnsc.miRNAseq.tum.norm = TumorNormalMatch(hnsc.miRNAseq$dat)

# Transpose the two matched pair datasets and log2 transform the values
hnsc.primary = data.frame(hnsc.miRNAseq.tum.norm[["primary.tumor"]])
hnsc.primary$column = rownames(hnsc.primary)
hnsc.primary = hnsc.primary[,order(ncol(hnsc.primary):1)]
hnsc.primary = setNames(data.frame(t(hnsc.primary[,-1])), hnsc.primary[,1])
hnsc.primary$sample = rownames(hnsc.primary)
hnsc.primary = hnsc.primary[,order(ncol(hnsc.primary):1)]
hnsc.primary$sample <- gsub('\\.', '-', hnsc.primary$sample)
hnsc.primary[, 2:dim(hnsc.primary)[2]] = log(hnsc.primary[2:dim(hnsc.primary)[2]]+1, 2) # NOTE: Adding small value of +1 to ensure valid transformation
hnsc.normal = data.frame(hnsc.miRNAseq.tum.norm[["normal"]])
hnsc.normal$column = rownames(hnsc.normal)
hnsc.normal = hnsc.normal[,order(ncol(hnsc.normal):1)]
hnsc.normal = setNames(data.frame(t(hnsc.normal[,-1])), hnsc.normal[,1])
hnsc.normal$sample = rownames(hnsc.normal)
hnsc.normal = hnsc.normal[,order(ncol(hnsc.normal):1)]
hnsc.normal$sample <- gsub('\\.', '-', hnsc.normal$sample)
hnsc.normal[, 2:dim(hnsc.normal)[2]] = log(hnsc.normal[2:dim(hnsc.normal)[2]]+1, 2)  # NOTE: Adding small value of +1 to ensure valid transformation
rownames(hnsc.primary) = rownames(hnsc.normal) = NULL


#### IMPUTE MISSING VALUES ####
# Remove sample column (categorical variable with >53 categories)
miRNA.mis = miRNA[,-1]
# Impute missing data 
miRNA.imp = missForest(miRNA.mis,variablewise=TRUE)
# Add sample column back
miRNA.imp = cbind(miRNA$sample,miRNA.imp$ximp)
# Rename column
colnames(miRNA.imp)[1] = 'sample'
# Convert column names to lowercase
colnames(miRNA.imp) = tolower(colnames(miRNA.imp))
# Reset row names
rownames(miRNA.imp) = NULL

# Remove sample column (categorical variable with >53 categories)
hnsc.normal.mis = hnsc.normal[-1]
# Impute missing data 
hnsc.normal.imp = missForest(hnsc.normal.mis,variablewise=TRUE)
# Add sample column back
hnsc.normal.imp = cbind(hnsc.normal$sample,hnsc.normal.imp$ximp)
# Rename column
colnames(hnsc.normal.imp)[1] = 'sample'
# Convert column names to lowercase
colnames(hnsc.normal.imp) = tolower(colnames(hnsc.normal.imp))
# Reset row names
rownames(hnsc.normal.imp) = NULL

# Remove sample column (categorical variable with >53 categories)
hnsc.primary.mis = hnsc.primary[-1]
# Impute missing data 
hnsc.primary.imp = missForest(hnsc.primary.mis,variablewise=TRUE)
# Add sample column back
hnsc.primary.imp = cbind(hnsc.primary$sample,hnsc.primary.imp$ximp)
# Rename column
colnames(hnsc.primary.imp)[1] = 'sample'
# Convert column names to lowercase
colnames(hnsc.primary.imp) = tolower(colnames(hnsc.primary.imp))
# Reset row names
rownames(hnsc.primary.imp) = NULL


#### IDENTIFY DYSREGULATED MIRNA ####
# Columns not shared among three datasets 
cols = intersect(colnames(miRNA.imp),intersect(colnames(hnsc.normal.imp),colnames(hnsc.primary.imp)))

# Add status column and combine normal and primary datasets
hnsc.normal.imp$status = 0
hnsc.primary.imp$status = 1
hnsc = rbind(hnsc.normal.imp[,c(cols,'status')],hnsc.primary.imp[,c(cols,'status')])

# Iterate over shared columns to identify significant columns
sig.cols = c()
coef = c()
pval = c()
regulated = c()
for(i in cols){
  model = glm(status ~ hnsc[,i], data=hnsc, family=binomial(link="logit"))
  if(summary(model)$coefficients[2,4] < 0.05){
    sig.cols = c(sig.cols, i)
    coef = c(coef,summary(model)$coefficients[2,1])
    pval= c(pval,round(summary(model)$coefficients[2,4],5))
    if(summary(model)$coefficients[2,1]>0){
      regulated = c(regulated,'up')
    } else {
      regulated = c(regulated,'down')
    }
  }
}

# Combine significant miRNA columns, coefs, and regulations
sig.miRNA = cbind(sig.cols,coef,pval,regulated)
sig.miRNA = data.frame(sig.miRNA)
# Count up/down regulated
sig.miRNA %>% count(regulated)


#### ADD TARGET VARIABLE ####
# Remove non-significant columns from miRNA.imp dataset
miRNA.imp = miRNA.imp[ , c("sample", sig.cols)]
# Add lymphovascular_invasion_present to miRNA.imp dataset
miRNA.imp$lymphovascular_invasion_present = clinical$lymphovascular_invasion_present[match(miRNA.imp$sample, clinical$sampleID)]
# Convert lymphovascular_invasion_present to 1/0 factor
miRNA.imp$lymphovascular_invasion_present = ifelse(miRNA.imp$lymphovascular_invasion_present=='YES',1L,0L)


#### TEST TRAIN SPLIT ####
# Split miRNA.imp to train and test sets
train_index = sample(1:nrow(miRNA.imp), 0.8 * nrow(miRNA.imp))
test_index = setdiff(1:nrow(miRNA.imp), train_index)
x_train = miRNA.imp[train_index, sig.cols]
y_train = miRNA.imp[train_index, "lymphovascular_invasion_present"]
x_test = miRNA.imp[test_index, sig.cols]
y_test = miRNA.imp[test_index, "lymphovascular_invasion_present"]
train = cbind(x_train,y_train)
test = cbind(x_test,y_test)


#### OUTPUT FILES ####
# # output CSV files
# write.csv(clinical, file='clinical.csv')
# write.csv(survival, file='survival.csv')
# write.csv(miRNA.imp, file='miRNAimp.csv')
# write.csv(x_train, file='x_train.csv')
# write.csv(y_train, file='y_train.csv')
# write.csv(x_test, file='x_test.csv')
# write.csv(y_test, file='y_test.csv')
# write.csv(train, file='train.csv')
# write.csv(test, file='test.csv')
# write.csv(train_index,'train_index.csv')
# write.csv(test_index,'test_index.csv')
