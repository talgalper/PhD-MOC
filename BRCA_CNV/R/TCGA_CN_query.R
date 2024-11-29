library(TCGAbiolinks)
library(SummarizedExperiment)
library(maftools)
library(tidyverse)
library(reshape2)
library(data.table)
library(DriverGenePathway)


getProjectSummary("TCGA-BRCA")

SNV_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

SNV_query_data <- getResults(SNV_query)

GDCdownload(SNV_query)
SNV_data <- GDCprepare(SNV_query)

save(SNV_data, SNV_query, SNV_query_data, file = "RData/BRCA_SNV_data.RData")

maf <- read.maf(SNV_data)

plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

DriverGenePathway(maf = maf, top = 10, removeNonMutated = TRUE)

#clincial <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")


CNV_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Copy Number Variation", 
  access = "open",
  sample.type = "Primary Tumor"
)

CNV_query_data <- getResults(CNV_query)

CNV_query_data_details <- lapply(CNV_query_data, table)

save(CNV_query, CNV_query_data, CNV_query_data_details, file = "RData/BRCA_CNV_query.RData")

GDCdownload(CNV_query)

CNV_data <- GDCprepare(CNV_query)

CNV_scores <- as.data.frame(assay(CNV_data, "copy_number"))
rownames(CNV_scores) <- gsub("\\.\\d+", "", rownames(CNV_scores))

# remove genes with NA values across all samples
NA_genes <- rownames(CNV_scores)[rowSums(is.na(CNV_scores)) == ncol(CNV_scores)]
CNV_scores <- CNV_scores[!rownames(CNV_scores) %in% NA_genes, ]

save(CNV_scores, file = "~/OneDrive - RMIT University/PhD/large_git_files/TCGA_CNV/geneLevel_CNV_NAfilt.RData")
load("~/OneDrive - RMIT University/PhD/large_git_files/TCGA_CNV/geneLevel_CNV_NAfilt.RData")

cnv_summary <- data.frame(
  gene = rownames(CNV_scores),  # Gene names
  num_gain = rowSums(CNV_scores > 2, na.rm = T),  # Count of gains for each gene
  num_loss = rowSums(CNV_scores < 2, na.rm = T),  # Count of losses for each gene
  num_none = rowSums(CNV_scores == 2, na.rm = T), # Count of neutral/no alterations for each gene
  NA_samples = rowSums(is.na(CNV_scores)) # Count the number of NA samples
)

cnv_summary$total_alterations <- cnv_summary$num_gain + cnv_summary$num_loss
cnv_summary$alterations_percentage <- (cnv_summary$total_alterations / (cnv_summary$total_alterations + cnv_summary$num_none)) * 100
cnv_summary$NA_percentage <- (cnv_summary$NA_samples / ncol(CNV_scores)) * 100

cnv_summary_filt <- cnv_summary[cnv_summary$NA_percentage <= 10, ]

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
                           filters = "ensembl_gene_id", 
                           values = cnv_summary_filt$gene, 
                           mart = ensembl)
ensembl_converted$description <- gsub("\\[.*?\\]", "", ensembl_converted$description)

cnv_summary_filt <- merge.data.table(ensembl_converted, cnv_summary_filt, by.x = "ensembl_gene_id", by.y = "gene", all.y = T)

save(cnv_summary_filt, file = "RData/BRCA_CNV_alterationsPerc.RData")





COSMIC <- fread("~/OneDrive - RMIT University/PhD/large_git_files/Cosmic_CompleteCNA_v100_GRCh38.tsv", sep = "\t")






library(RTCGAToolbox)

firehoseData <- getFirehoseData("BRCA", GISTIC = T, CNASNP = T, CNVSNP = T, CNASeq = T, Mutation = T)
CNV_GISTIC_data.allByGene <- getData(firehoseData, "GISTIC", platform = "AllByGene")
CNV_GISTIC_data.ThresholdedByGene <- getData(firehoseData, "GISTIC", platform = "ThresholdedByGene")
CNV_GISTIC_data.Peaks <- getData(firehoseData, "GISTIC", platform = "Peaks")

save(firehoseData, file = "RData/firehose_data.RData")






# download COSMIC data from: https://cancer.sanger.ac.uk/cosmic/download/cosmic/v100/completecna.
# This paper uses GISTIC2 data analysed by MutSig2CV which is probably how HHnet did it: https://www.mdpi.com/2072-6694/15/16/4090 
# Kylie paper filters for copy number polymorphisms which is masked copy number