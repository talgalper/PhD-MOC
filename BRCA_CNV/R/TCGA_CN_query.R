library(TCGAbiolinks)
library(SummarizedExperiment)
library(maftools)
library(tidyverse)
library(reshape2)


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

oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)

#clincial <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")



CNV_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Copy Number Variation", 
  access = "open",
  #data.type = "Gene Level Copy Number",
  sample.type = "Primary Tumor"
)

CNV_query_data <- getResults(CNV_query)

save(CNV_query, CNV_query_data, file = "RData/BRCA_CNV_geneLevelData.RData")

GDCdownload(CNV_query)

CNV_data <- GDCprepare(CNV_query)

CNV_scores <- as.data.frame(assay(CNV_data, "copy_number"))
rownames(CNV_scores) <- gsub("\\.\\d+", "", rownames(CNV_scores))

# remove genes with NA values across all samples
NA_genes <- rownames(CNV_scores)[rowSums(is.na(CNV_scores)) == ncol(CNV_scores)]
CNV_scores <- CNV_scores[!rownames(CNV_scores) %in% NA_genes, ]


CNV_plot_data <- CNV_scores
CNV_plot_data$gene <- rownames(CNV_plot_data)
CNV_plot_data <- melt(CNV_plot_data, id.vars = "gene")  # Convert to long format

# Categorize CNV values
CNV_plot_data$type <- ifelse(CNV_plot_data$value > 2, "Gain", 
                             ifelse(CNV_plot_data$value < 2, "Loss", "Neutral"))


ggplot(cnv_data_melt, aes(x = variable, y = gene, fill = type)) +
  geom_tile(color = "white") +  # Create a tile for each sample and gene
  scale_fill_manual(values = c("Gain" = "red", "Loss" = "green", "Neutral" = "grey")) +  # Color mapping
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate sample names
        axis.text.y = element_text(size = 7)) +  # Adjust text size for readability
  labs(x = "Samples", y = "Genes", fill = "Type")



library(RTCGAToolbox)

firehoseData <- getFirehoseData("BRCA", GISTIC = T, CNASNP = T, CNVSNP = T, CNASeq = T, Mutation = T)
CNV_GISTIC_data.allByGene <- getData(firehoseData, "GISTIC", platform = "AllByGene")
CNV_GISTIC_data.ThresholdedByGene <- getData(firehoseData, "GISTIC", platform = "ThresholdedByGene")
CNV_GISTIC_data.Peaks <- getData(firehoseData, "GISTIC", platform = "Peaks")

save(firehoseData, file = "RData/firehose_data.RData")






# download COSMIC data from: https://cancer.sanger.ac.uk/cosmic/download/cosmic/v100/completecna.
# This paper uses GISTIC2 data analysed by MutSig2CV which is probably how HHnet did it: https://www.mdpi.com/2072-6694/15/16/4090 
# Kylie paper filters for copy number polymorphisms which is masked copy number