# Bulk Transcriptomics data
## The core objective of this project is to identify genes that are significantly upregulated or downregulated in patients suffering from specific heart conditions compared to non-failing controls.
### The analysis processes transcriptomics data (Counts Per Million - CPM) and clinical metadata across four distinct patient groups (Etiologies)
NF: Non-Failing (Control group)
DCM: Dilated Cardiomyopathy
HCM: Hypertrophic Cardiomyopathy
PPCM: Peripartum Cardiomyopathy
The final output delivers an annotated, filtered, and statistically robust list of genes alongside their expression metrics (FPKM) across these medical conditions.

### General outline:
Exploratory Diagnostic Plots: Before performing statistics, visualizing raw data distributions using box plots and density plots helps verify sample uniformity and detect major outliers.

Principal Component Analysis (PCA) reduces complex high-dimensional genetic data into simple visual clusters to spot batch effects or natural clinical groupings.

Covariate Adjustment (Limma): Biological data is naturally noisy. While a crude analysis reveals general patterns, this pipeline explicitly controls for Age and Gender using a linear model. This minimizes confounding variables, isolating the true genetic signals driven purely by the heart disease.

Multiple Testing Correction: Testing thousands of genes simultaneously increases the risk of false positives. The pipeline uses adjusted p-values (q-values) to ensure statistical stringency.

Expression Noise Filtering: Low-abundance genes often represent technical background noise. To solve this mathematically, the script targets Y-chromosome genes within female samples. Because females do not have a Y chromosome, any expression detected there is pure background noise. Calculating this mean "noise" threshold creates an empirical baseline to filter out unexpressed or irrelevant genes safely.


### The workflow is divided into 6 distinct execution phases:

#### 1. Data Import & Demographic Overview
Loads the core gene expression matrix (CPM) and corresponding clinical sample profiles.
Validates data alignment to guarantee that sample IDs match between both datasets perfectly.
Generates a publication-ready patient demographic table (tbl_summary) stratified by disease etiology.

#### 2. Quality Control & Diagnostic Visualizations
Distribution Check: Transforms data structure to build comparative Boxplots and Density Plots for checking distribution similarities across phenotypes.
Dimensionality Reduction: Conducts PCA on the expression profile matrix to project samples onto PC1 and PC2, sizing points by patient age and coloring by etiology.

#### 3. Differential Expression Profiling
Constructs an ExpressionSet object compatible with Bioconductor standards.
Runs a crude analysis versus a covariate-adjusted analysis utilizing limma (Linear Models for Microarrays).
Applies empirical Bayes smoothing (eBayes) to calculate moderated t-statistics for comparisons: DCM vs Control, HCM vs Control, and PPCM vs Control.

#### 4. Biomart Gene Annotation
Establishes a connection to the Ensembl database (incorporating automatic mirror-switching logic in case of primary server downtime).
Maps raw Ensembl Gene IDs (ENSG...) to human-readable gene symbols (external_gene_name).

#### 5. Expression Conversion & Noise Filtering
Converts relative expression from CPM to FPKM (Fragments Per Kilobase million) using external exon length data.
Isolates female sample Y-chromosome data to calculate the baseline global background noise.
Slices the matrix to preserve only genes showing expression values strictly above this calculated noise ceiling.

#### 6. Alignment & Downstream Export
Aligns clinical statistics with structural gene symbols and calculated average expression metrics (Global average, DCM specific average, and Non-Failing specific average).
Exports a tidy, tab-delimited final file ready for further biological pathway analyses.

