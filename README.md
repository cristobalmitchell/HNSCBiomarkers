# HNSC Biomarkers
[The Cancer Geonome Atlas](https://cancergenome.nih.gov) (TCGA) "a landmark cancer genomics program, molecularly characterized over 20,000 primary cancer and matched normal samples spanning 33 cancer types. This joint effort between the National Cancer Institute and the National Human Genome Research Institute began in 2006, bringing together researchers from diverse disciplines and multiple institutions.

Over the next dozen years, TCGA generated over 2.5 petabytes of genomic, epigenomic, transcriptomic, and proteomic data. The data, which has already lead to improvements in our ability to diagnose, treat, and prevent cancer, will remain publicly available for anyone in the research community to use."

This analysis, written in R, looks at microRNA (miRNA) sequences and clinical data for Head and Neck Squamous Cell Carcoma (HNSC) from TCGA attempting to identify statistically significant biomarkers indicating lymphovascular invasion presence across clinical stages and predicting survival. With [A novel microRNA signature predicts survival in liver hepatocellular carcinoma after hepatectomy](https://doi.org/10.1038/s41598-018-26374-9) as inspiration, this analysis was completed as a capstone project for STAT 685 under the supervision of [Dr. Moumita Karmakar](https://stat.tamu.edu/people/#all-k) for a MS in Applied Statistics at Texas A&M University in 2020. 

## Materials and Methods
### TCGA dataset
A total of 2246 miRNA expression profiles in HNSC patients along with their corresponding clinical and survival data were downloaded from the GDC TCGA data portal using UCSC's Xena Browser data hub (September 2020). Patients without lymphovascular invasion present data and follow-up less than 30 days were removed. Additionally, miRNA expression profiles with more then 30% missing values were also excluded from the analysis resulting in 233 patient observations with 646 miRNA expression profiles.

### Identification of dysregulated miRNAs in HNSC

### Identification of miRNAs with prognostic score in HNSC


## Citations
Fu, Q., Yang, F., Xiang, T. et al. A novel microRNA signature predicts survival in liver hepatocellular carcinoma after hepatectomy. Sci Rep 8, 7933 (2018). https://doi.org/10.1038/s41598-018-26374-9

Goldman, M.J., Craft, B., Hastie, M. et al. Visualizing and interpreting cancer genomics data via the Xena platform. Nat Biotechnol (2020). https://doi.org/10.1038/s41587-020-0546-8
