# =============================================================================#
# MSB1005_DEGs_project.R                                          																                                       		  #
# Version: 1.0   															                                #
# Date: Dec 14, 2023											                                    #
# Author: Mostafa Hassanein, Maastricht University            														                                #
# =============================================================================#
# This R script serves as portable, annotated R script
#-----------------------------------------------------------------------------#
# Assignment:
#-----------------------------------------------------------------------------#
# Installing packages:
# List of required packages
packages <- c(
  "styler", "tidyverse", "gtsummary", "ggplot2", "cowplot", 
  "BiocManager", "Biobase", "limma", "biomaRt", "dplyr"
)

# Function to install missing packages and load all
load_packages <- function(pkg_list) {
  new_packages <- pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
  
  if (length(new_packages) > 0) {
    # Special handling for Bioconductor packages
    bioc_pkgs <- c("Biobase", "limma", "biomaRt")
    
    # Install standard packages
    std_pkgs <- new_packages[!(new_packages %in% bioc_pkgs)]
    if (length(std_pkgs) > 0) install.packages(std_pkgs)
    
    # Install Bioconductor packages
    if (any(new_packages %in% bioc_pkgs)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(new_packages[new_packages %in% bioc_pkgs], update = FALSE, ask = FALSE)
    }
  }
  
  # Load all packages
  suppressPackageStartupMessages(invisible(lapply(pkg_list, library, character.only = TRUE)))
}

# Execute
load_packages(packages)

# 1. Data import:
# a: To Import all the data files:
# After unzipping the required transcriptomics (gene expression) file and
# sample information file, import both files as two separate objects.
# First set the active working directory to the folder containing the files.

setwd("D:/Maastricht University/Period_2/Data Analysis/assignment")
gxData <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt",
  as.is = T,
  row.names = 1
)
sampleInfo <- read.csv("MAGNET_SampleData_18112022.csv",
  as.is = T,
  row.names = 1
)

# b) A publication-ready table:

library(tidyverse)
library(gtsummary)

sampleInfo %>%
  select(!c(Library.Pool, minexpr)) %>%
  tbl_summary(
    by = etiology,
    missing_text = "Missing values",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    digits = all_continuous() ~ 2
  ) %>%
  add_p() %>%
  add_q() %>%
  add_n() %>%
  add_overall() %>%
  modify_spanning_header((c("stat_1", "stat_2", "stat_3", "stat_4") ~ "**Etiology**"))


# 2) Diagnostic plots:

library(ggplot2)
library(cowplot)

# a) Data distribution boxplot:

NF_index <- which(sampleInfo$etiology == "NF")
NF_Ids <- rownames(sampleInfo)[NF_index]
gxData_NF <- gxData[, NF_Ids]
plot_NF <- gather(gxData_NF, key = "SampleID", value = "CPM")
p1 <- ggplot(plot_NF, aes(x = SampleID, y = CPM)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

PPCM_index <- which(sampleInfo$etiology == "PPCM")
PPCM_Ids <- rownames(sampleInfo)[PPCM_index]
gxData_PPCM <- gxData[, PPCM_Ids]
plot_PPCM <- gather(gxData_PPCM, key = "SampleID", value = "CPM")
p2 <- ggplot(plot_PPCM, aes(x = SampleID, y = CPM)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

DCM_index <- which(sampleInfo$etiology == "DCM")
DCM_Ids <- rownames(sampleInfo)[DCM_index]
gxData_DCM <- gxData[, DCM_Ids]
plot_DCM <- gather(gxData_DCM, key = "SampleID", value = "CPM")
p3 <- ggplot(plot_DCM, aes(x = SampleID, y = CPM)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

HCM_index <- which(sampleInfo$etiology == "HCM")
HCM_Ids <- rownames(sampleInfo)[HCM_index]
gxData_HCM <- gxData[, HCM_Ids]
plot_HCM <- gather(gxData_HCM, key = "SampleID", value = "CPM")
p4 <- ggplot(plot_HCM, aes(x = SampleID, y = CPM)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_grid(p1, p2, p3, p4, labels = c("NF", "PPCM", "DCM", "HCM"), label_size = 12)
ggsave("distrib_box_plot.png")

# a) Data distribution density plot "better for visualization than box blot":

p_1 <- ggplot(plot_NF, aes(x = CPM, color = SampleID)) +
  geom_density(alpha = 0.2) + # with density only x-axis
  theme(panel.background = element_blank()) +
  theme(legend.position = "none")

p_2 <- ggplot(plot_PPCM, aes(x = CPM, color = SampleID)) +
  geom_density(alpha = 0.2) + # with density only x-axis
  theme(panel.background = element_blank()) +
  theme(legend.position = "none")

p_3 <- ggplot(plot_DCM, aes(x = CPM, color = SampleID)) +
  geom_density(alpha = 0.2) + # with density only x-axis
  theme(panel.background = element_blank()) +
  theme(legend.position = "none")

p_4 <- ggplot(plot_HCM, aes(x = CPM, color = SampleID, )) +
  geom_density(alpha = 0.2) + # with density only x-axis
  theme(panel.background = element_blank()) +
  theme(legend.position = "none")

plot_grid(p_1, p_2, p_3, p_4,
  labels = c("NF", "PPCM", "DCM", "HCM"),
  label_size = 12, vjust = c(1.5, 1.5, 0.5, 0.5), hjust = -1
)
ggsave("distrib_density_plot.png")

# b) PCA figure:

require(pcaMethods)
library(ggplot2)
trans_gxData <- t(gxData)
pcaRes <- pca(t(gxData), nPcs = 10)

pca_plot_data <- cbind(data.frame(pcaRes@scores), sampleInfo)

ggplot(pca_plot_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(size = age, col = etiology))
ggsave("pca_plot.png")

# Here I make PCA because I have many samples, so
# I have to reduce dimensions to be able to visualize our data
# I choose PC1 and PC2 because they have the most variation of the variables.

# 3) a) Statistical analysis:

library(BiocManager)
library(Biobase)
library(limma)

# I perform a differential gene expression analysis using the limma package.
# I convert the data frame to a matrix.
# Then, create a design matrix for the linear model.
# After that, creates a contrast matrix for the linear model,
# and fits the linear model to the ExpressionSet object using the design matrix design.
# Then, I apply the contrast matrix to the fitted mode.
# Finally, I apply empirical Bayes smoothing to the contrast estimates.

gxData1 <- as.matrix(gxData)
eset <- ExpressionSet(
  assayData = gxData1,
  phenoData = AnnotatedDataFrame(sampleInfo)
)

design <- model.matrix(~ 0 + etiology, data = pData(eset))

cont.matrix <- makeContrasts(
  DCMvsControl = etiologyDCM - etiologyNF,
  HCMvsControl = etiologyHCM - etiologyNF,
  PPCMvsControl = etiologyPPCM - etiologyNF,
  levels = design
)
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrasts = cont.matrix)
fit2 <- eBayes(fit2, trend = TRUE)
results <- decideTests(fit2)
summary(results)
dgeRes <- topTable(fit2, coef = c("DCMvsControl", "HCMvsControl", "PPCMvsControl"), number = nrow(gxData))
summary(dgeRes)
dgeRes

# b) Correct for co-variates:

design3 <- model.matrix(~ 0 + etiology + gender + age, data = sampleInfo)
fit <- lmFit(gxData, design3)
cont.matrix <- makeContrasts(
  DCMvsControl = etiologyDCM - etiologyNF,
  HCMvsControl = etiologyHCM - etiologyNF,
  PPCMvsControl = etiologyPPCM - etiologyNF,
  levels = design3
)
fit2 <- contrasts.fit(fit, cont.matrix)
ebFit <- eBayes(fit2, trend = TRUE)
results <- decideTests(ebFit)
summary(results)
dgeRes2 <- topTable(ebFit, coef = c("DCMvsControl", "HCMvsControl", "PPCMvsControl"), number = nrow(gxData))
summary(dgeRes2)
stat_data <- dgeRes2

# Here I correct for the main variables that have the highest correlation
# in this data. These variables are age and gender.
# Many of the other variables have incomplete data or have low correlation.
# Moreover, I have the adjusted p.value that is adjusted for multiple testing

# 4) a) Gene annotation:

# Using biomart package to choose the database and dataset. Connecting the database and the dataset.
# Setting the attributes and fliters. Finally, build the query.


library(biomaRt)
library(dplyr)

gxData_IDs <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt")
gxData_IDs <- gxData_IDs[, 1]

ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensemble_connection <- useMart("ensembl",
  dataset = "hsapiens_gene_ensembl"
)
attrib <- listAttributes(ensemble_connection)
filters <- listFilters(ensemble_connection)

gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = gxData_IDs,
  mart = ensemble_connection
)

# 4) b) Merge dataframes:

updated_gene_symbols <- gene_symbols %>%
  remove_rownames() %>%
  column_to_rownames(var = "ensembl_gene_id")
combined_gxData <- updated_gene_symbols %>%
  rownames_to_column() %>%
  left_join(gxData %>% rownames_to_column())

# 5) Relative expression levels:

# a) Transform the data to FPKM value:

geneTotExonLengths <- read.delim("MAGNET_exonLengths.txt",
  as.is = T,
  row.names = 1
)
all(rownames(geneTotExonLengths) == rownames(gxData))
cpm2fpkm <- function(x) {
  geneTotExonLengths_kb <- geneTotExonLengths[, 1] / 1E3
  .t <- 2^(x) / geneTotExonLengths_kb
}
gxData_fpkm <- cpm2fpkm(gxData)

# b) Background (noise) level:

# Selecting the expression data of female subjects only, and build the query of them.
# Then filter these genes for the genes that only on Y chromosome.

gxData_fpkm_female <- subset(gxData_fpkm[, sampleInfo$gender == "Female"])
gxData_fpkm_female <- cbind(rownames(gxData_fpkm_female), gxData_fpkm_female)
gxData_fpkm_female_IDs <- gxData_fpkm_female[, 1]

ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensemble_connection <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attrib <- listAttributes(ensemble_connection)
filters <- listFilters(ensemble_connection)

# Y chromosome genes in female subjects:

gene_positions_female <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = gxData_fpkm_female_IDs,
  mart = ensemble_connection
)

filter_criteria <- (gene_positions_female$chromosome_name == "Y")
female_y_genes <- gene_positions_female[filter_criteria, ]
gxData_fpkm_female_y <- gxData_fpkm_female[female_y_genes$ensembl_gene_id, ]

# Average gene expression of Y chromosome genes in female subjects:

gxData_fpkm_female_y_updated <- subset(gxData_fpkm_female_y, select = -1)
mean_df <- colMeans(gxData_fpkm_female_y_updated)
mean_noise <- mean(mean_df)
gxData_fpkm_unnois <- gxData_fpkm - mean_noise
gxData_fpkm_unnois_avg <- gxData_fpkm_unnois %>%
  mutate(average = rowMeans(.))
gxData_fpkm_abov_noise <- gxData_fpkm_unnois_avg %>%
  select(average) %>%
  filter(average > 0)

# Therefore, I filtered for genes in the dataset that is
# expressed above background (noise) level.


# 6) Export the result:

library(tibble)

# Setting gene_id as the row names through converting the dataframe into a tibble and back again to a dataframe
# Using the last column which contain the average values
# Matching function to filter the original data to that only have the genes' names (20009 out of 20781)
# Repeating with different datasets (fpkm_avg, DCM, NF, and stat_data)
# Writing a txt file containing external gene names, average fpkm, DCM fpkm, NF fpkm, and statistics comparing
# the different groups including also the adjusted P-value.

tb <- as_tibble(gene_symbols)
tb <- column_to_rownames(tb, var = "ensembl_gene_id")
gene_symbols_rows <- as.data.frame(tb)
gxData_fpkm_unnois_avg_only <- gxData_fpkm_unnois_avg[, ncol(gxData_fpkm_unnois_avg), drop = FALSE]
gxData_fpkm_unnois_avg_only <- data.frame(gxData_fpkm_unnois_avg_only[match(
  row.names(gene_symbols_rows),
  row.names(gxData_fpkm_unnois_avg_only)
), ])
colnames(gxData_fpkm_unnois_avg_only)[1] <- "fpkm_average"
gxData_fpkm_annotated <- gene_symbols_rows %>% mutate(gxData_fpkm_unnois_avg_only = gxData_fpkm_unnois_avg_only$fpkm_average)

gxData_fpkm_DCM <- gxData_fpkm_unnois_avg[, sampleInfo$etiology == "DCM"]
gxData_fpkm_DCM <- gxData_fpkm_DCM %>%
  mutate(average = rowMeans(.))

gxData_fpkm_NF <- gxData_fpkm_unnois_avg[, sampleInfo$etiology == "NF"]
gxData_fpkm_NF <- gxData_fpkm_NF %>%
  mutate(average = rowMeans(.))

gxData_fpkm_DCM_avg_only <- gxData_fpkm_DCM[, ncol(gxData_fpkm_DCM), drop = FALSE]
gxData_fpkm_DCM_avg_only <- data.frame(gxData_fpkm_DCM_avg_only[match(
  row.names(gene_symbols_rows),
  row.names(gxData_fpkm_DCM_avg_only)
), ])
colnames(gxData_fpkm_DCM_avg_only)[1] <- "fpkm_DCM_average"

gxData_fpkm_NF_avg_only <- gxData_fpkm_NF[, ncol(gxData_fpkm_NF), drop = FALSE]
gxData_fpkm_NF_avg_only <- data.frame(gxData_fpkm_NF_avg_only[match(
  row.names(gene_symbols_rows),
  row.names(gxData_fpkm_NF_avg_only)
), ])

colnames(gxData_fpkm_NF_avg_only)[1] <- "fpkm_NF_average"


stat_data_updat <- data.frame(stat_data[match(
  row.names(gene_symbols_rows),
  row.names(stat_data)
), ])

gxData_fpkm_annotated <- gene_symbols_rows %>%
  mutate(
    fpkm_unnois_avg = gxData_fpkm_unnois_avg_only$fpkm_average,
    fpkm_DCM = gxData_fpkm_DCM_avg_only$fpkm_DCM_average,
    fpkm_NF = gxData_fpkm_NF_avg_only$fpkm_NF_average,
    stat_data_updat
  )

# Set the directory to the folder where you want save your txt file.

write.table(gxData_fpkm_annotated,
  file = "D:/Maastricht University/Period_2/Data Analysis/assignment/final_data.txt",
  sep = "\t"
)
