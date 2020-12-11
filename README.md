# HNSC Biomarkers
[The Cancer Geonome Atlas](https://cancergenome.nih.gov) (TCGA) "a landmark cancer genomics program, molecularly characterized over 20,000 primary cancer and matched normal samples spanning 33 cancer types. This joint effort between the National Cancer Institute and the National Human Genome Research Institute began in 2006, bringing together researchers from diverse disciplines and multiple institutions.

Over the next dozen years, TCGA generated over 2.5 petabytes of genomic, epigenomic, transcriptomic, and proteomic data. The data, which has already lead to improvements in our ability to diagnose, treat, and prevent cancer, will remain publicly available for anyone in the research community to use."

This analysis, written in R, looks at [microRNA (miRNA)](https://en.wikipedia.org/wiki/MicroRNA) sequences and clinical data for Head and Neck Squamous Cell Carcoma (HNSC) from TCGA attempting to identify statistically significant biomarkers indicating lymphovascular invasion (LVI) presence to group patients in high and low risk classes and predicting survival. With [A novel microRNA signature predicts survival in liver hepatocellular carcinoma after hepatectomy](https://doi.org/10.1038/s41598-018-26374-9) as inspiration, this analysis was completed under the supervision of [Dr. Moumita Karmakar](https://stat.tamu.edu/people/#all-k) as a capstone project for STAT 685 while completing an MS in Applied Statistics at Texas A&M University in Fall 2020. 


## Materials and Methods
### TCGA dataset
A total of 2246 miRNA expression profiles in HNSC patients along with their corresponding clinical and survival data were downloaded from the GDC TCGA data portal using [UCSC's Xena Browser](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Head%20and%20Neck%20Cancer%20(HNSC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) data hub (September 2020). Patients without lymphovascular invasion present data and follow-up less than 30 days were removed. Additionally, miRNA expression profiles with more then 30% missing values were also excluded from the analysis.

### Identification of dysregulated miRNAs in HNSC
The raw counts of miRNA expression data of 43 HNSC patients with their paired normal tissues were obtained from the TCGA dataset using the [TCGA2STAT](http://www.liuzlab.org/TCGA2STAT/) library. The miRNAs with missing values were imputed using the same missForest package mentioned previously and were then normalized using log2(x+1). The miRNA expressions shared between the TCGA dataset and the match paired datasets were subsequently evaluated using univariate logistic regression and those with p-values < 0.05 were considered to be differentiaally expressed miRNA and included for subsequent analysis. 

### Identification of miRNAs with prognostic score in HNSC
To identify significant miRNAs related to LVI a stepwise backward logisitc regression model was used, and HR > 1 or HR < 1 with P < 0.05 was used as the cutoff. The patient subclasses in each group of clinical characteristics represented non-overlapping sets.

### Definition of prognostic risk model and ROC curve analysis
The TCGA data were randomly divided into a training and testing set. A  stepwise bothways logisitc regression model was used to select significant miRNAs and develop the linear miRNA signature prognostic model. The prognostic score was calculated as follows: Prognostic score = (-0.29541 × hsa-mir-1293) + (0.29090 × hsa-mir-378c) + (0.35463 × hsa-mir-137) + (0.06019 × hsa-mir-206). The prognostic scores were then computed for all 233 patients. The optimal cutoff point of the prognostic score was obtained in ROC curve analysis for predicting the 5-year survival of the training set. Kaplan-Meier and log-rank methods were used for evaluating OS.

## Results
### Identification of differentially expressed miRNAs in HNSC. 
After filtering out unqualified cases with a follow-up of less than 30 days, a total of 233 HNSC patients (170 male/63 female) were retained for LVI anal- ysis, among which was a total of 43 patients with solid tumour normal tissues. These cases were randomly divided into a training set (n = 186) and a testing set (n = 47). There was no significant difference in clinical covariates observed between the two sets. Differentially expressed analysis in 43 pairs of primary tumour and solid tumour normal tissues indicated that 28 miRNAs were differentially expressed (log FC > 1 or log FC < −1, P < 0.05 after FDR adjustment). Among these miRNAs, 21 miRNAs were up-regulated, and 7 miRNAs were down-regulated.

### Establishment of the miRNA prognostic model. 
The common miRNAs associated with LVI were characterized with the univariate Cox regression method. Four miRNAs were selected using stepwise backward feature selection for a logistic regression model using the training set. Thereafter, we developed a miRNA prognostic model for diagnosis in the HNSC cohort. According to ROC curve, the optimum cutoff point was calculated and used for classifying the patients into high-risk and low-risk groups. This analysis showed the prognostic scores and miRNA expression distribution of all 233 patients and 43 pairs of solid tumour normal tissues, which were ranked according to the prognostic scores for the four-miRNA signature. It indicated that the patients with higher prognostic scores showed a tendency towards expression of high-risk miRNAs, whereas patients with low prognostic scores showed a tendency towards non-protective miRNA expression.

<img src="/images/ROC_prog_threshold.png" />

### Validation of the four-miRNA signature in HNSC patients. 
Using the optimum cutoff value obtained from the training set, patients were assigned to high-risk and low-risk groups. The ability of the four-miRNA signature to predict 5-year survival prognostication was examined in the testing set and in all HNSC patients. The results failed to indicate that patients in the high-risk group had poor OS, nor that patients in the low-risk group had good outcomes in the training set (P > 0.05, with the log-rank test method) and in the entire HNSC cohort (P > 0.05).

<img src="/images/Kaplan_Meier_Risk.png" />

## Discussion
Reviewing the data sets used and comparing our setup with that of the article that inspired this analysis a few issues were present that likely influenced our inconclusive results. To begin the "Novel miRNA signature..." study had a normal tumor matched data set that contained 49 patients all of which were included in the TCGA-LIHC cohort data set. This was not reflected in our setup with only 9 of the 43 patients included in our normal tumor matched data set being shared with the TCGA-HNSC cohort data set. Additionally, the study used as inspiration had common miRNA expression profiles across different TMN stages indicating a generalized prognostic formula could be generated for the entire cohort. This too was not reflected in our TCGA-HNSC cohort data set with virtually no miRNA expression profiles being shared across more than 2 clinical stages. 

Taking a slightly different path and applying the univariate logistic regression to identify significant (p-value < 0.05) miRNA expression profiles directly from the TCGA-HNSC cohort data then building a stepwise bothways logisitic regression model had the best results but were still insignificant. 

## Acknowledgements
The results shown here are in whole or part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga.

## Citations
Fu, Q., Yang, F., Xiang, T. et al. A novel microRNA signature predicts survival in liver hepatocellular carcinoma after hepatectomy. Sci Rep 8, 7933 (2018). https://doi.org/10.1038/s41598-018-26374-9

Goldman, M.J., Craft, B., Hastie, M. et al. Visualizing and interpreting cancer genomics data via the Xena platform. Nat Biotechnol (2020). https://doi.org/10.1038/s41587-020-0546-8

Ying-Wooi Wan, Genevera I. Allen, Zhandong Liu, TCGA2STAT: simple TCGA data access for integrated statistical analysis in R, Bioinformatics, Volume 32, Issue 6, 15 March 2016, Pages 952–954, https://doi.org/10.1093/bioinformatics/btv677
