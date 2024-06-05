### first run DE functions
library(biomaRt)
library(edgeR)
library(SummarizedExperiment)
library(tidyverse)
library(VennDiagram)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


### These results use BRCA.Normal samples as a control group ###

# Load data
load("RData/LumA/DE_data.RData")
load("RData/LumB/DE_data.RData")
load("RData/Her2/DE_data.RData")
load("RData/basal/DE_data.RData")
load("RData/TCGA_normal.RData")

# tidy env
counts_data <- list(LumA_unstranded = LumA_unstranded,
                    LumB_unstranded = LumB_unstranded,
                    Her2_unstranded = Her2_unstranded,
                    Basal_unstranded = Basal_unstranded,
                    Normal_unstranded = normal_unstranded)

rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded, normal_unstranded)


## luminal A results
lumA_QC <- filter_low_expr(counts_data$LumA_unstranded, counts_data$Normal_unstranded)
lumA_DE <- DE_analysis(counts_matrix = lumA_QC$counts_filt, sample_info = lumA_QC$sample_info)

## luminal B results
lumB_QC <- filter_low_expr(counts_data$LumB_unstranded, counts_data$Normal_unstranded)
lumB_DE <- DE_analysis(counts_matrix = lumB_QC$counts_filt, sample_info = lumB_QC$sample_info)

## Her2 results
Her2_QC <- filter_low_expr(counts_data$Her2_unstranded, counts_data$Normal_unstranded)
Her2_DE <- DE_analysis(counts_matrix = Her2_QC$counts_filt, sample_info = Her2_QC$sample_info)

## basal results
basal_QC <- filter_low_expr(counts_data$Basal_unstranded, counts_data$Normal_unstranded)
basal_DE <- DE_analysis(counts_matrix = basal_QC$counts_filt, sample_info = basal_QC$sample_info)


### These results use GTEx mammary tissue as a control group ###

# load and clear data
GTEx_data <- read.table("gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, -1]
rownames(GTEx_data) <- NULL

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames

# add GTEx data to data object and tidy env
counts_data["GTEx_data"] <- list(GTEx_ENS)
rm(GTEx_data, GTEx_ENS, rownames)


## luminal A results with GTEx as control
GTEx_lumA_QC <- filter_low_expr(counts_data$LumA_unstranded, counts_data$GTEx_data)
GTEx_lumA_DE <- DE_analysis(counts_matrix = GTEx_lumA_QC$counts_filt, sample_info = GTEx_lumA_QC$sample_info)

## luminal B results with GTEx as control
GTEx_lumB_QC <- filter_low_expr(counts_data$LumB_unstranded, counts_data$GTEx_data)
GTEx_lumB_DE <- DE_analysis(counts_matrix = GTEx_lumB_QC$counts_filt, sample_info = GTEx_lumB_QC$sample_info)

## Her2 results with GTEx as control
GTEx_Her2_QC <- filter_low_expr(counts_data$Her2_unstranded, counts_data$GTEx_data)
GTEx_Her2_DE <- DE_analysis(counts_matrix = GTEx_Her2_QC$counts_filt, sample_info = GTEx_Her2_QC$sample_info)

## basal results with GTEx as control
GTEx_basal_QC <- filter_low_expr(counts_data$Basal_unstranded, counts_data$GTEx_data)
GTEx_basal_DE <- DE_analysis(counts_matrix = GTEx_basal_QC$counts_filt, sample_info = GTEx_basal_QC$sample_info)



# compare logFC of markers
gene_id <- getBM(attributes = c("external_gene_name", "hgnc_symbol", "ensembl_gene_id"), 
                 filters = "external_gene_name", 
                 values = c("ESR1", "PGR", "ERBB2"), 
                 mart = ensembl)


LumA_subset <- lumA_DE$hits[rownames(lumA_DE$hits) %in% gene_id$ensembl_gene_id, ]
LumA_subset <- subset(LumA_subset, select = c("gene_id", "logFC"))
LumB_subset <- lumB_DE$hits[rownames(lumB_DE$hits) %in% gene_id$ensembl_gene_id, ]
LumB_subset <- subset(LumB_subset, select = c("gene_id", "logFC"))
Her2_subset <- Her2_DE$hits[rownames(Her2_DE$hits) %in% gene_id$ensembl_gene_id, ]
Her2_subset <- subset(Her2_subset, select = c("gene_id", "logFC"))
Basal_subset <- basal_DE$hits[rownames(basal_DE$hits) %in% gene_id$ensembl_gene_id, ]
Basal_subset <- subset(Basal_subset, select = c("gene_id", "logFC"))
GTEx_LumA_subset <- GTEx_lumA_DE$hits[rownames(GTEx_lumA_DE$hits) %in% gene_id$ensembl_gene_id, ]
GTEx_LumA_subset <- subset(GTEx_LumA_subset, select = c("gene_id", "logFC"))
GTEx_LumB_subset <- GTEx_lumB_DE$hits[rownames(GTEx_lumB_DE$hits) %in% gene_id$ensembl_gene_id, ]
GTEx_LumB_subset <- subset(GTEx_LumB_subset, select = c("gene_id", "logFC"))
GTEx_Her2_subset <- GTEx_Her2_DE$hits[rownames(GTEx_Her2_DE$hits) %in% gene_id$ensembl_gene_id, ]
GTEx_Her2_subset <- subset(GTEx_Her2_subset, select = c("gene_id", "logFC"))
GTEx_Basal_subset <- GTEx_basal_DE$hits[rownames(GTEx_basal_DE$hits) %in% gene_id$ensembl_gene_id, ]
GTEx_Basal_subset <- subset(GTEx_Basal_subset, select = c("gene_id", "logFC"))


marker_logFC <- merge(gene_id, LumA_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
colnames(marker_logFC)[3] <- "TCGA_lumA"
marker_logFC <- merge(marker_logFC, GTEx_LumA_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
colnames(marker_logFC)[4] <- "GTEx_lumA"

# Luminal B (LumB)
marker_logFC <- merge(marker_logFC, LumB_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
colnames(marker_logFC)[5] <- "TCGA_lumB"
marker_logFC <- merge(marker_logFC, GTEx_LumB_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
colnames(marker_logFC)[6] <- "GTEx_lumB"

# HER2-enriched (Her2)
marker_logFC <- merge(marker_logFC, Her2_subset, by.x = "ensembl_gene_id", by.y = "gene_id", all = T)
colnames(marker_logFC)[7] <- "TCGA_her2"
marker_logFC <- merge(marker_logFC, GTEx_Her2_subset, by.x = "ensembl_gene_id", by.y = "gene_id", all = T)
colnames(marker_logFC)[8] <- "GTEx_her2"

# Basal-like (Basal)
marker_logFC <- merge(marker_logFC, Basal_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
colnames(marker_logFC)[9] <- "TCGA_basal"
marker_logFC <- merge(marker_logFC, GTEx_Basal_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
colnames(marker_logFC)[10] <- "GTEx_basal"




venn.diagram(
  x = list(TCGA_control = lumA_DE$dif_exp_genes, GETx_control = rownames(GTEx_lumA_DE$dif_exp)),
  category.names = c("DE with TCGA Normal", "DE with GTEx Healthy"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("dodgerblue", "goldenrod1"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("dodgerblue", "goldenrod1"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "DE_comparison/lumA_TCGA_GTEx_comparison.png",
  disable.logging = TRUE
)

venn.diagram(
  x = list(TCGA_control = lumB_DE$dif_exp_genes, GETx_control = rownames(GTEx_lumB_DE$dif_exp)),
  category.names = c("DE with TCGA Normal", "DE with GTEx Healthy"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("dodgerblue", "goldenrod1"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("dodgerblue", "goldenrod1"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "DE_comparison/lumB_TCGA_GTEx_comparison.png",
  disable.logging = TRUE
)

venn.diagram(
  x = list(TCGA_control = Her2_DE$dif_exp_genes, GETx_control = rownames(GTEx_Her2_DE$dif_exp)),
  category.names = c("DE with TCGA Normal", "DE with GTEx Healthy"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("dodgerblue", "goldenrod1"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("dodgerblue", "goldenrod1"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "DE_comparison/Her2_TCGA_GTEx_comparison.png",
  disable.logging = TRUE
)

venn.diagram(
  x = list(TCGA_control = basal_DE$dif_exp_genes, GETx_control = rownames(GTEx_basal_DE$dif_exp)),
  category.names = c("DE with TCGA Normal", "DE with GTEx Healthy"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("dodgerblue", "goldenrod1"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("dodgerblue", "goldenrod1"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "DE_comparison/basal_TCGA_GTEx_comparison.png",
  disable.logging = TRUE
)
