library(drug_net)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)

# get TCGA data
TCGA_data <- get_TCGA_RNAseq_data(project_name = "TCGA-BRCA", subtype = "Basal")

# read in GTEx data
GTEx_raw <- read.table("../BRCA_pipe/bulk-gex_v8_rna-seq_counts-by-tissue_gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2, header = T)

# combine data into single df
data <- TCGA_GTEx_combine(TCGA_data, GTEx_data = GTEx_raw)

# remove low activity genes
data_processed <- remove_low_activity_genes(data, top_x_samples = 0.50, min_samples = 0.10)