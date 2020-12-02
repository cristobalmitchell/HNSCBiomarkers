# Author: Cristobal Mitchell
# Course: STAT 685 700 - Dr. Moumita Karmakar
# Title: TCGA-HNSC Analysis

# Set seed to reproduce results
set.seed(1618)

#### IMPORTS ####
library(glmnet)
library(tidyverse) 

#### READ PREPROCESSED FILES ####
# Read in previous datasets for consistent results
miRNA.imp = read.csv('miRNAimp.csv', row.names = 1, header=TRUE)
clinical = read.csv('clinical.csv', header=TRUE)
survival = read.csv('survival.csv', header=TRUE)
x_train = read.csv('x_train.csv', row.names = 1, header=TRUE)
y_train = read.csv('y_train.csv', row.names = 1, header=TRUE)
x_test = read.csv('x_test.csv', row.names = 1, header=TRUE)
y_test = read.csv('y_test.csv', row.names = 1, header=TRUE)
train = read.csv('train.csv', header=TRUE)
test = read.csv('test.csv', header=TRUE)

colnames(miRNA.imp) = gsub('\\.','-',colnames(miRNA.imp))
colnames(clinical) = gsub('\\.','-',colnames(clinical))
# colnames(survival) = gsub('\\.','-',colnames(survival))
colnames(x_train) = gsub('\\.','-',colnames(x_train))
colnames(y_train) = gsub('\\.','-',colnames(y_train))
colnames(x_test) = gsub('\\.','-',colnames(x_test))
colnames(y_test) = gsub('\\.','-',colnames(y_test))
colnames(train) = gsub('\\.','-',colnames(train))
colnames(test) = gsub('\\.','-',colnames(test))

rownames(x_train) = NULL
colnames(miRNA.imp)[2:646]
# Iterate over shared columns to identify significant columns
sig.cols = c()
coef = c()
pval = c()
regulated = c()
for(i in colnames(miRNA.imp)[2:646]){
  model = glm(lymphovascular_invasion_present ~ miRNA.imp[,i], data=miRNA.imp, family=binomial(link="logit"))
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


#### MODEL PENALIZED LASSO ####
# Fit a penalized lasso/glmnet with respect to LVI using train set
cv.lasso = cv.glmnet(as.matrix(x_train), as.double(y_train), alpha = 1, family = "binomial")
model_pl = glmnet(as.matrix(x_train), y_train, alpha = 1, family = "binomial",
               lambda = cv.lasso$lambda.min)
coef(model_pl)
plot(cv.lasso)


#### MODEL SIMPLE LOGISTIC ####
# simple LR using full data set
miRNA.simple = subset(miRNA.imp, select=-c(sample))

fit = glm(lymphovascular_invasion_present~.,data=miRNA.simple, family = "binomial")
summary(fit)
coef(fit)


#### MODEL STEPWISE BOTHWAYS ####
# Stepwise varaiable select from train set
library(MASS)
full_model <- glm(train$y_train ~ ., data=train, family = "binomial")
summary(full_model)

int_only <- glm(train$y_train ~ 1, data=train, family = "binomial")
summary(int_only)

bothways = step(int_only, list(lower=formula(int_only),upper=formula(full_model)), direction="both",trace=0)
formula(bothways)
summary(bothways)

backwards = step(full_model,trace=0)
formula(backwards)
summary(backwards)

#### PROGNOSTIC SCORE ####
# # prognostic score based on penalized lasso
# miRNA.imp$prog = (0.081126474*miRNA.imp[,'hsa-mir-378c'])-(0.070099102*miRNA.imp[,'hsa-mir-1293'])+
#   (0.005710965*miRNA.imp[,'hsa-mir-133b'])
# 
# # prognostic score based on stepwise backward
# miRNA.imp$prog = (0.08908*miRNA.imp[,'hsa-mir-206'])-(0.37072*miRNA.imp[,'hsa-mir-4326'])+
#   (0.36513*miRNA.imp[,'hsa-mir-3662'])+(0.32535*miRNA.imp[,'hsa-mir-137'])-
#   (0.40514*miRNA.imp[,'hsa-mir-1293'])
# 
# # REVISITED prognostic score based on sig stepwise backward
# miRNA.imp$prog = (0.08908*miRNA.imp[,'hsa-mir-206'])-(0.37072*miRNA.imp[,'hsa-mir-4326'])-
#   (0.40514*miRNA.imp[,'hsa-mir-1293'])

# prognostic score based on stepwise bothways
miRNA.imp$prog = -(0.29541*miRNA.imp[,'hsa-mir-1293'])+(0.29090*miRNA.imp[,'hsa-mir-378c'])+
  (0.35463*miRNA.imp[,'hsa-mir-137'])+(0.06019*miRNA.imp[,'hsa-mir-206'])

# # REVISITED prognostic score based on sig stepwise bothways
# miRNA.imp$prog = -(0.29541*miRNA.imp[,'hsa-mir-1293'])+(0.35463*miRNA.imp[,'hsa-mir-137'])

# # prognostic score based on sig cols from miRNA.imp
# miRNA.imp$prog = (0.02045796*miRNA.imp[,'hsa-mir-885-5p'])-(0.02747145*miRNA.imp[,'hsa-mir-203a-3p'])+
#   (0.04620911*miRNA.imp[,'hsa-let-7d-5p'])+(0.08784594*miRNA.imp[,'hsa-mir-30a-3p'])+
#   (0.01576196*miRNA.imp[,'hsa-mir-128-3p'])+(0.22968588*miRNA.imp[,'hsa-mir-4491'])+
#   (0.11470938*miRNA.imp[,'hsa-mir-6087'])

# prognostic score based on sig cols from miRNA.imp
miRNA.imp$prog = (0.6204*miRNA.imp[,'hsa-mir-4491'])+(0.4428*miRNA.imp[,'hsa-mir-6087'])+
  (0.9956*miRNA.imp[,'hsa-mir-885-5p'])-(1.7396*miRNA.imp[,'hsa-mir-6820-3p'])-
  (1.0011*miRNA.imp[,'hsa-mir-629-5p'])-(0.4245*miRNA.imp[,'hsa-mir-1304-3p'])+
  (0.6534*miRNA.imp[,'hsa-mir-151a-3p'])+(0.6467*miRNA.imp[,'hsa-mir-25-3p'])-
  (0.3439*miRNA.imp[,'hsa-mir-101-3p'])


#### ROC BASED THRESHOLD ####
# ROC plot based on the above prog score calculation
library(ROCR)
pred <- prediction(miRNA.imp$prog, miRNA.imp$lymphovascular_invasion_present)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)

library(pROC)
plot.roc(miRNA.imp$lymphovascular_invasion_present, miRNA.imp$prog, 
         main="Confidence interval of a threshold", 
         percent=TRUE,
         legacy.axes=TRUE,
         ci=TRUE, 
         of="thresholds", 
         thresholds="best", 
         print.thres="best")
threshold = coords(roc(miRNA.imp$lymphovascular_invasion_present, miRNA.imp$prog), "best", "threshold")[[1]]
auc(roc(miRNA.imp$lymphovascular_invasion_present, miRNA.imp$prog))


#### COX REGRESSION ####
# Group based on prog threshold from ROC
miRNA.imp$risk = ifelse(miRNA.imp$prog >= threshold,"high","low")
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

#### KAPLAN MEIER ####
# Plot the Kaplan-Meier estimate for both high/low risk group
km <- survfit(Surv(time,OS) ~ risk, data = miRNA.imp)
ggsurvplot(km, data = miRNA.imp, pval = TRUE)

km <- survfit(Surv(time,OS) ~ risk, data = cox_train)
ggsurvplot(km, data = cox_train, pval = TRUE)

km <- survfit(Surv(time,OS) ~ risk, data = cox_test)
ggsurvplot(km, data = cox_test, pval = TRUE)

