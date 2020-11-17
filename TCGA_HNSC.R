# Author: Cristobal Mitchell
# Course: STAT 685 700 - Dr. Moumita Karmakar
# Title: TCGA-HNSC Analysis

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

# Read in data from files
clinical = read.table('HNSC_clinicalMatrix', sep='\t', header=TRUE)
miRNA = read.table('miRNA_HiSeq_gene', sep='\t')
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
#   XenaFilter(filterDatasets = "miRNA_HiSeq_gene") %>% 
#   XenaFilter(filterDatasets = "HNSC") -> df_todo
# XenaQuery(df_todo) %>%
#   XenaDownload() -> xe_download
# miRNA = XenaPrepare(xe_download)
#
# XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
#   XenaFilter(filterDatasets = "survival") %>% 
#   XenaFilter(filterDatasets = "HNSC") -> df_todo
# XenaQuery(df_todo) %>%
#   XenaDownload() -> xe_download
# survival = XenaPrepare(xe_download)

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
# Reset row names
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
# Set seed to reproduce
set.seed(1618)
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

# output miRNA.imp
write.csv(miRNA.imp, file='miRNAimp.csv')

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
  model = glm(status ~ hnsc[,i], data=hnsc, family = binomial(link="logit"))
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

# Remove non-significant columns from miRNA.imp dataset
miRNA.imp = miRNA.imp[ , c("sample", sig.cols)]
# Add lymphovascular_invasion_present to miRNA.imp dataset
miRNA.imp$lymphovascular_invasion_present = clinical$lymphovascular_invasion_present[match(miRNA.imp$sample, clinical$sampleID)]
# Convert lymphovascular_invasion_present to 1/0 factor
miRNA.imp$lymphovascular_invasion_present = ifelse(miRNA.imp$lymphovascular_invasion_present=='YES',1L,0L)

# Which patients are shared between imputed dataset and tumor normal
miRNA.imp$patient = substr(miRNA.imp$sample, 1, 12)
intersect(miRNA.imp$patient,hnsc$sample)
# Drop patient column
miRNA.imp = subset(miRNA.imp,select=-c(patient))

# Split miRNA.imp to train and test sets
train_index = sample(1:nrow(miRNA.imp), 0.8 * nrow(miRNA.imp))
test_index = setdiff(1:nrow(miRNA.imp), train_index)
x_train = miRNA.imp[train_index, sig.cols]
y_train = miRNA.imp[train_index, "lymphovascular_invasion_present"]
x_test = miRNA.imp[test_index, sig.cols]
y_test = miRNA.imp[test_index, "lymphovascular_invasion_present"]
train = cbind(x_train,y_train)
test = cbind(x_test,y_test)

# Fit a penalized lasso/glmnet with respect to LVI using train set
cv.lasso = cv.glmnet(as.matrix(x_train), y_train, alpha = 1, family = "binomial")
model = glmnet(as.matrix(x_train), y_train, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.min)
coef(model)

plot(cv.lasso)
cv.lasso$lambda.min

miRNA.simple = subset(miRNA.imp, select=-c(sample))

fit = glm(lymphovascular_invasion_present~.,data=miRNA.simple, family = "binomial")
summary(fit)
coef(fit)

# Prognostic score based on coefs from lasso glmnet
miRNA.imp$prog = (0.03238346*miRNA.imp[,'hsa-mir-378c'])

# ROC plot based on the above prog score calculation
library(ROCR)
pred <- prediction(miRNA.imp$prog, miRNA.imp$lymphovascular_invasion_present)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)

library(pROC)
plot.roc(miRNA.imp$lymphovascular_invasion_present, miRNA.imp$prog, 
         main="Confidence interval of a threshold", 
         percent=TRUE,
         ci=TRUE, 
         of="thresholds", 
         thresholds="best", 
         print.thres="best")

# Group based on prog threshold from ROC
miRNA.imp$risk = ifelse(miRNA.imp$prog >= 0.1,"high","low")
# Add OS to miRNA.imp dataset
miRNA.imp$OS = survival$OS[match(miRNA.imp$sample, survival$sample)]
# Add OS.time to miRNA.imp dataset
miRNA.imp$time = survival$OS.time[match(miRNA.imp$sample, survival$sample)]

# Fit a cox regression predictive model for overall survival on the training data set
library(survival)
library(survminer)
# Use same train/test split from previously
cox_train = miRNA.imp[train_index,]
cox_test = miRNA.imp[test_index,]
surv_object <- Surv(cox_train$time,cox_train$OS)
cox.fit = coxph(Surv(time, OS) ~ risk, data=cox_train)
summary(cox.fit)
# Plot the Kaplan-Meier estimate for both high/low risk group.
km <- survfit(surv_object ~ risk, data = cox_train)
ggsurvplot(km, data = cox_train, pval = TRUE)


# Stepwise varaiable select from train set
library(MASS)
full.model <- glm(train$y_train ~ ., data=train, family = "binomial")
summary(full.model)

nothing <- glm(train$y_train ~ 1, data=train, family = "binomial")
summary(nothing)

backwards = step(full.model,trace=0)
formula(backwards)
summary(backwards)

forwards = step(nothing, scope=list(lower=formula(nothing),upper=formula(full.model)), direction="forward",trace=0)
formula(forwards)
summary(forwards)

bothways = step(nothing, list(lower=formula(nothing),upper=formula(full.model)), direction="both",trace=0)
formula(bothways)
summary(bothways)

# Accuracy of nothing
probabilities <- nothing %>% predict(x_test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
observed.classes <- y_test
mean(predicted.classes == observed.classes)

# Accuracy of full.model
probabilities <- full.model %>% predict(x_test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
observed.classes <- y_test
mean(predicted.classes == observed.classes)

# Accuracy of forwards model
probabilities <- forwards %>% predict(x_test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
observed.classes <- y_test
mean(predicted.classes == observed.classes)

# Accuracy of backwards model
probabilities <- backwards %>% predict(x_test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
observed.classes <- y_test
mean(predicted.classes == observed.classes)

# Accuracy of bothways model
probabilities <- bothways %>% predict(x_test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
observed.classes <- y_test
mean(predicted.classes == observed.classes)

# Unsupervised PCA
miRNA.pca = prcomp(miRNA.imp[,sig.cols], center=TRUE, scale=TRUE)
summary(miRNA.pca)
biplot(miRNA.pca, groups=miRNA.imp$lymphovascular_invasion_present)
screeplot(miRNA.pca, type = "l", npcs = 15, main = "Screeplot of the first 15 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
# Based on the above screenplot looking at PC1-PC9
miRNA.pca_sub = miRNA.pca$rotation[,colnames(miRNA.pca$rotation)[1:9]]
round(miRNA.pca_sub, 4)

# Conduct logistic regression using PC

## MethylIT package wont install correctly hard coding the pcaLogisticR function
# install.packages(c("ArgumentCheck", "caret", "Epi", "e1071", "minpack.lm", "nls2", "caTools", "rmarkdown",   "RCurl"),dependencies=TRUE)
# install.packages("devtools")
# BiocManager::install("S4Vectors")
# devtools::install_git("https://github.com/genomaths/MethylIT.git")
# library(MethylIT)

pcaLogisticR <- function(formula = NULL,
                         data = NULL,
                         n.pc = 1,
                         scale = FALSE,
                         center = FALSE,
                         tol = 1e-04,
                         max.pc = NULL) {
  
  Check <- ArgumentCheck::newArgCheck()
  if (!is.null(formula) && !is(formula, "formula")) {
    ans <- paste("A formula of the form groups ~ x1 + x2 + ... (see",
                 "?pcaLogisticR or ?glm).", sep = "")
    ArgumentCheck::addError(msg = ans, argcheck = Check)
  }
  if (!is.null(formula) && is(formula, "formula")) {
    vn <- try(attr(terms(formula), "term.labels"),
              silent = TRUE)
    if (inherits(vn, "try-error")) {
      vn <- try(setdiff(colnames(data), as.character(formula)[2]),
                silent = TRUE)
    }
    if (inherits(vn, "try-error"))
      stop("* Error in the formula")
    if (length(vn) < n.pc) {
      ans <- "The number of number predictor variables greater than "
      ans1 <- "the number of principal components: "
      ArgumentCheck::addError(msg = paste0(ans,
                                           ans1, n.pc), argcheck = Check)
    }
  }
  if (is.null(formula)) {
    ans <- "A formula or grouping variable must be provided."
    ArgumentCheck::addError(msg = ans, argcheck = Check)
  }
  ArgumentCheck::finishArgCheck(Check)
  
  m = nrow(data)
  if (floor(m/3) < n.pc) {
    ans <- "The number principal components: "
    ans1 <- " must be lower than the number of individuals N/3"
    warning(paste0(ans, n.pc, ans1))
  }
  pc <- prcomp(x = data[vn], retx = TRUE, center = center,
               scale. = scale, tol = tol, rank. = max.pc)
  
  cn <- colnames(pc$x)
  if (ncol(pc$x) > n.pc) {
    ind.coord <- as.data.frame(pc$x[, seq_len(n.pc)])
    colnames(ind.coord) <- cn[seq_len(n.pc)]
  } else {
    ind.coord <- as.data.frame(pc$x)
    colnames(ind.coord) <- cn[seq_len(ncol(pc$x))]
  }
  
  resp <- as.character(formula)[2]
  res <- data[, resp]
  l <- levels(res)
  res <- as.character(res)
  res[res == l[1]] <- 0
  res[res == l[2]] <- 1
  res <- as.numeric(res)
  data <- data.frame(res, ind.coord)
  colnames(data) <- c(resp, colnames(ind.coord))
  predictors <- paste(colnames(ind.coord), collapse = " + ")
  formula <- paste0(resp, " ~ ", predictors)
  model <- suppressWarnings(glm(formula = formula,
                                family = binomial(link = "logit"), data = data))
  
  model <- structure(list(logistic = model, pca = pc,
                          reference.level = l[1], positive.level = l[2]),
                     class = "pcaLogisticR")
  return(model)
}

predict.pcaLogisticR <- function(object, ...) UseMethod("predict")
predict.pcaLogisticR <- function(object, newdata, type = c("class",
                                                           "posterior", "pca.ind.coord", "all"), ...) {
  
  if (!is(object, "pcaLogisticR")) {
    stop("* Parameter 'object' must be a model from class 'pcaLogisticR'")
  }
  type <- match.arg(type)
  ## predictor names
  vn <- rownames(object$pca$rotation)
  
  if (!is.null(newdata) && inherits(newdata, c("pDMP",
                                               "InfDiv")))
    newdata <- unlist(newdata)
  if (inherits(newdata, "GRanges")) {
    if (is.element("pos", vn)) {
      position <- function(gr) {
        chrs <- split(gr, seqnames(gr))
        gr <- lapply(chrs, function(grc) {
          x <- start(grc)
          x.min <- min(x)
          x.max <- max(x)
          if (x.min == Inf)
            x.min = 0
          if (x.max == -Inf)
            x.max = 1
          delta <- max(c(x.max - x, 1))
          return((x - x.min)/(delta))
        })
        return(unlist(gr))
      }
      newdata$pos <- position(newdata)
    }
    newdata$logP <- log10(newdata$wprob + 2.2e-308)
    newdata <- mcols(newdata)
  }
  newdata <- newdata[vn]
  newdata <- as.matrix(newdata)
  
  ## Centering and scaling new individuals
  dt.scaled <- scale(newdata, center = object$pca$center,
                     scale = object$pca$scale)
  ## Coordinates of the individuals
  coord_func <- function(ind, loadings) {
    x <- loadings * ind
    return(apply(x, 2, sum))
  }
  nc <- ncol(object$pca$x)
  loadings <- object$pca$rotation[, seq_len(nc)]
  if (nc == 1)
    loadings <- as.matrix(loadings)
  
  ind.coord <- data.frame(t(apply(dt.scaled, 1, coord_func,
                                  loadings)))
  if (nc == 1) {
    ind.coord <- as.data.frame(t(ind.coord))
    colnames(ind.coord) <- "PC1"
    row.names(ind.coord) <- NULL
  }
  
  predictClass <- function(object, dt) {
    pred <- predict(object$logistic, newdata = dt,
                    type = "response")
    PredClass <- rep(object$reference.level, nrow(newdata))
    PredClass[pred > 0.5] <- object$positive.level
    return(PredClass)
  }
  
  pred <- switch(type, posterior = predict(object$logistic,
                                           newdata = ind.coord, type = "response"),
                 class = predictClass(object = object,
                                      dt = ind.coord), pca.ind.coord = ind.coord,
                 all = list(class = predictClass(object = object,
                                                 dt = ind.coord), posterior = predict(object$logistic,
                                                                                      newdata = ind.coord, type = "response"),
                            pca.ind.coord = ind.coord))
  return(pred)
}

# Logistic regression for LVI using PC
formula <- lymphovascular_invasion_present ~ .
pca.logistic <- pcaLogisticR(formula=formula, 
             data=miRNA.simple, 
             n.pc = 9,
             scale=TRUE, 
             center=TRUE
             )
summary(pca.logistic$logistic)

## Different biplot of PC1 & PC2
# install.packages('factoextra')
library(factoextra)
fviz_pca_biplot(miRNA.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             ind.fill = miRNA.imp$lymphovascular_invasion_present, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Lymphovascular invasion") +
  ggtitle("2D PCA-plot from 28 feature dataset") +
  theme(plot.title = element_text(hjust = 0.5))

## The below block of commented code is an attempt to do Supervised PCA based on Bair et al.
## however the code does not run as expected and the example provided by the author is 
## based on simulated data - http://statweb.stanford.edu/~tibs/superpc/tutorial.html

# library(superpc)
# set.seed(464)
# x<-as.matrix(subset(miRNA.imp,select=-c(lymphovascular_invasion_present,OS,time,sample)))
# y<-miRNA.imp[,'time']
# #xtest<-x
# #ytest<-2+5*v1+ .05*rnorm(100)
# censoring.status<- miRNA.imp[,'OS']
# #censoring.status.test<- sample(c(rep(1,80),rep(0,20)))
# featurenames <- colnames(x)
# data<-list(x=x, y=y, censoring.status=censoring.status, featurenames=featurenames)
# #data.test<-list(x=xtest,y=ytest, censoring.status=censoring.status.test, featurenames= featurenames)
# train.obj<- superpc.train(data, type="survival")
# cv.obj<-superpc.cv(train.obj, data)
# superpc.plotcv(cv.obj)
# lrtest.obj<-superpc.lrtest.curv(train.obj, data,data.test)
# superpc.plot.lrtest(lrtest.obj)
# fit.cts<- superpc.predict(train.obj, data, data.test, threshold=0.7, n.components=3, prediction.type="continuous")
# superpc.fit.to.outcome(train.obj, data.test, fit.cts$v.pred)
# fit.groups<- superpc.predict(train.obj, data, data.test, threshold=0.7, n.components=1, prediction.type="discrete")
# superpc.fit.to.outcome(train.obj, data.test, fit.groups$v.pred)
# plot(survfit(Surv(data.test$y,data.test$censoring.status)~fit.groups$v.pred), col=2:3, xlab="time", ylab="Prob survival")


#### Subgroup Stratification Analysis ####
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

