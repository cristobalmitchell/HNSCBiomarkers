# Author: Cristobal Mitchell
# Course: STAT 685 700 - Dr. Moumita Karmakar
# Title: TCGA HNSC Analysis

# install.packages("UCSCXenaTools") # https://docs.ropensci.org/UCSCXenaTools/
# install.packages("TCGA2STAT_1.2.tar", repos = NULL, type = "source") # http://www.liuzlab.org/TCGA2STAT/ 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefilter") # https://bioconductor.org/packages/release/bioc/html/genefilter.html
BiocManager::install("CNTools") # http://bioconductor.org/packages/release/bioc/html/CNTools.html
BiocManager::install("miRBaseConverter") # http://bioconductor.org/packages/release/bioc/html/miRBaseConverter.html

library(UCSCXenaTools)
library(ggplot2)
library(TCGA2STAT)
library(miRBaseConverter)

# Read in data from files
clinical = read.table('HNSC_clinicalMatrix', sep='\t', header=TRUE)
miRNA = read.table('miRNA_HiSeq_gene', sep='\t')


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
#   XenaFilter(filterDatasets = "miRNA_HiSeq_gene") %>% 
#   XenaFilter(filterDatasets = "HNSC") -> df_todo
# XenaQuery(df_todo) %>%
#   XenaDownload() -> xe_download
# miRNA = XenaPrepare(xe_download)

# Gather miRNA names from the accession values found in the first column of miRNA dataset
miRNAAccessions = miRNA$V1
miRNANames = miRNA_AccessionToName(miRNAAccessions,targetVersion = "v21") # v22 results in some NA target names
miRNANames[1,2] = miRNANames[1,1]

# Replace the column names
miRNA$V1 = miRNANames$TargetName

# Transpose miRNA dataset so each row represent a patient
miRNA = setNames(data.frame(t(miRNA[,-1])), miRNA[,1])

# Convert MIMA columns to correct datatype for miRNA dataset
miRNA[,-1] = sapply(miRNA[,-1], as.character)
miRNA[,-1] = sapply(miRNA[,-1], as.numeric)

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
# Reset rownames
rownames(miRNA) = NULL
rownames(clinical) = NULL
# Remove columns with > 30% missing values from miRNA data.frame
miRNA = miRNA[,colSums(is.na(miRNA))<=(dim(miRNA)[1]*.3)]

# Download the tumor normal matched data for the HNSC cancer
hnsc.miRNAseq = getTCGA(disease="HNSC", data.type="miRNASeq")
hnsc.miRNAseq.tum.norm = TumorNormalMatch(hnsc.miRNAseq$dat)

# Transpose the two datasets and log2 transform
hnsc.primary = data.frame(hnsc.miRNAseq.tum.norm[["primary.tumor"]])
hnsc.primary$column = rownames(hnsc.primary)
hnsc.primary=hnsc.primary[,order(ncol(hnsc.primary):1)]
hnsc.primary = setNames(data.frame(t(hnsc.primary[,-1])), hnsc.primary[,1])
hnsc.primary$sample = rownames(hnsc.primary)
hnsc.primary=hnsc.primary[,order(ncol(hnsc.primary):1)]
hnsc.primary$sample <- gsub('\\.', '-', hnsc.primary$sample)
hnsc.primary[, 2:dim(hnsc.primary)[2]] = log(hnsc.primary[2:dim(hnsc.primary)[2]]+1, 2) # NOTE: Adding small value of +1 to ensure valid transformation
rownames(hnsc.primary) = NULL

hnsc.normal = data.frame(hnsc.miRNAseq.tum.norm[["normal"]])
hnsc.normal$column = rownames(hnsc.normal)
hnsc.normal=hnsc.normal[,order(ncol(hnsc.normal):1)]
hnsc.normal = setNames(data.frame(t(hnsc.normal[,-1])), hnsc.normal[,1])
hnsc.normal$sample = rownames(hnsc.normal)
hnsc.normal=hnsc.normal[,order(ncol(hnsc.normal):1)]
hnsc.normal$sample <- gsub('\\.', '-', hnsc.normal$sample)
hnsc.normal[, 2:dim(hnsc.normal)[2]] = df = log(hnsc.normal[2:dim(hnsc.normal)[2]]+1, 2)  # NOTE: Adding small value of +1 to ensure valid transformation
rownames(hnsc.normal) = NULL

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
# Reset rownames
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
# Reset rownames
rownames(hnsc.normal.imp) = NULL

# Remove columns not shared among three datasets 
cols = intersect(sort(colnames(miRNA.imp)),sort(colnames(hnsc.normal.imp)))

# Add status column and combine normal and primary datasets
hnsc.normal$status = 0
hnsc.primary$status = 1
hnsc = rbind(hnsc.normal,hnsc.primary)

# Iterate over shared columns to identify significant columns
sig.cols = c()
for(i in cols){
  model = glm(status ~ hnsc[,i], data=hnsc, family = binomial(link="logit"))
  if(summary(model)$coefficients[2,4] < 0.05){
    sig.cols = c(sig.cols, i)
  }
}

# Remove non-significant columns from miRNA.imp dataset
miRNA.imp = miRNA.imp[ , c("sample", sig.cols)]
# Add lymphovascular_invasion_present to miRNA.imp dataset
miRNA.imp$lymphovascular_invasion_present = clinical$lymphovascular_invasion_present[match(miRNA.imp$sample, clinical$sampleID)]
# Convert lymphovascular_invasion_present to 1/0 factor
miRNA.imp$lymphovascular_invasion_present = ifelse(miRNA.imp$lymphovascular_invasion_present=='YES',1L,0L)

# Subgroup based on clinical_stage
stg.i = miRNA.imp[miRNA.imp$sample %in% clinical[clinical$clinical_stage == 'Stage I', 'sampleID'],]
stg.ii = miRNA.imp[miRNA.imp$sample %in% clinical[clinical$clinical_stage == 'Stage II', 'sampleID'],]
stg.iii = miRNA.imp[miRNA.imp$sample %in% clinical[clinical$clinical_stage == 'Stage III', 'sampleID'],]
stg.iva = miRNA.imp[miRNA.imp$sample %in% clinical[clinical$clinical_stage == 'Stage IVA', 'sampleID'],]
stg.ivb = miRNA.imp[miRNA.imp$sample %in% clinical[clinical$clinical_stage == 'Stage IVB', 'sampleID'],]
stg.ivc = miRNA.imp[miRNA.imp$sample %in% clinical[clinical$clinical_stage == 'Stage IVC', 'sampleID'],]


# Combine all stage iv
stg.iv = rbind(stg.iva, stg.ivb, stg.ivc)

fit = glmnet(x = as.matrix(stg.i[,sig.cols]), y = stg.i$lymphovascular_invasion_present)
plot(fit)
coef(fit,s=0.1)

fit = glmnet(x = as.matrix(stg.ii[,sig.cols]), y = stg.ii$lymphovascular_invasion_present)
plot(fit)
coef(fit,s=0.1)

fit = glmnet(x = as.matrix(stg.iii[,sig.cols]), y = stg.iii$lymphovascular_invasion_present)
plot(fit)
coef(fit,s=0.1)

fit = glmnet(x = as.matrix(stg.iv[,sig.cols]), y = stg.iv$lymphovascular_invasion_present)
plot(fit)
coef(fit,s=0.1)

# Export data.frame to file
write.csv(clinical, file='clinical.csv')
write.csv(miRNA, file='miRNA.csv')
write.csv(hnsc.normal, file='hnsc_normal.csv')
write.csv(hnsc.primary, file='hnsc_primary.csv')

write.csv(miRNA.imp, file='miRNA_imputed.csv')
write.csv(stg.i, file='stage_i_imp_sig.csv')
write.csv(stg.ii, file='stage_ii_imp_sig.csv')
write.csv(stg.iii, file='stage_iii_imp_sig.csv')
write.csv(stg.iva, file='stage_iva_imp_sig.csv')
write.csv(stg.ivb, file='stage_ivb_imp_sig.csv')
write.csv(stg.ivc, file='stage_ivc_imp_sig.csv')

