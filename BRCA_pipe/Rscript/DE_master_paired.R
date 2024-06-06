### first run DE functions
library(biomaRt)
library(edgeR)
library(SummarizedExperiment)
library(tidyverse)
library(VennDiagram)
library(reshape2)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


### These results use BRCA.Normal samples as a control group ###

# load in paired sample info from Rscript/paired_analysis/paired_samples.R 
load("RData/paired/paired_subtypes.RData")

paired_samples <- list(lumA_paired = lumA_paired,
                       lumB_paired = lumB_paired,
                       Her2_paired = Her2_paired,
                       Basal_paired = Basal_paired)
rm(lumA_paired, lumB_paired, Her2_paired, Basal_paired, master)

# some samples had 2 healthy samples to one disease, needed to remove them
for (i in seq_along(paired_samples)) {
  subtype_info <- as.data.frame(paired_samples[i])
  subtype_info <- subtype_info[!duplicated(subtype_info[1]), ]
  colnames(subtype_info) <- sub(".*\\.", "", colnames(subtype_info))
  paired_samples[[i]] <- subtype_info
}
rm(i, subtype_info)


# read in counts data
load("RData/LumA/DE_data.RData")
load("RData/LumB/DE_data.RData")
load("RData/Her2/DE_data.RData")
load("RData/basal/DE_data.RData")
load("RData/TCGA_normal.RData")

counts_data <- list(LumA_unstranded = LumA_unstranded,
                    LumB_unstranded = LumB_unstranded,
                    Her2_unstranded = Her2_unstranded,
                    Basal_unstranded = Basal_unstranded,
                    Normal_unstranded = normal_unstranded)
# tidy env
rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded, normal_unstranded)


## luminal A results
lumA_QC <- paired_filter_low_expr(counts_data$LumA_unstranded, counts_data$Normal_unstranded, paired_samples$lumA_paired)
lumA_DE <- paired_DE_analysis(counts_matrix = lumA_QC$counts_filt, sample_info = lumA_QC$sample_info)

## luminal B results
lumB_QC <- paired_filter_low_expr(counts_data$LumB_unstranded, counts_data$Normal_unstranded, paired_samples$lumB_paired)
lumB_DE <- paired_DE_analysis(counts_matrix = lumB_QC$counts_filt, sample_info = lumB_QC$sample_info)

## Her2 results
Her2_QC <- paired_filter_low_expr(counts_data$Her2_unstranded, counts_data$Normal_unstranded, paired_samples$Her2_paired)
Her2_DE <- paired_DE_analysis(counts_matrix = Her2_QC$counts_filt, sample_info = Her2_QC$sample_info)

## basal results
basal_QC <- paired_filter_low_expr(counts_data$Basal_unstranded, counts_data$Normal_unstranded, paired_samples$Basal_paired)
basal_DE <- paired_DE_analysis(counts_matrix = basal_QC$counts_filt, sample_info = basal_QC$sample_info)


# save results
DE_results <- list(TCGA_lumA = list(hits = lumA_DE$hits, dif_exp = lumA_DE$dif_exp),
                   TCGA_lumB = list(hits = lumB_DE$hits, dif_exp = lumB_DE$dif_exp),
                   TCGA_Her2 = list(hits = Her2_DE$hits, dif_exp = Her2_DE$dif_exp),
                   TCGA_basal = list(hits = basal_DE$hits, dif_exp = basal_DE$dif_exp))

save(DE_results, file = "RData/DE_results_master_paired.RData")



# compare logFC of markers
gene_id <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
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


marker_logFC <- merge(gene_id, LumA_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
colnames(marker_logFC)[3] <- "TCGA_lumA"


# Luminal B (LumB)
marker_logFC <- merge(marker_logFC, LumB_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
colnames(marker_logFC)[5] <- "TCGA_lumB"


# HER2-enriched (Her2)
marker_logFC <- merge(marker_logFC, Her2_subset, by.x = "ensembl_gene_id", by.y = "gene_id", all = T)
colnames(marker_logFC)[7] <- "TCGA_her2"


# Basal-like (Basal)
marker_logFC <- merge(marker_logFC, Basal_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
colnames(marker_logFC)[9] <- "TCGA_basal"






# need to update with paired vs unpaired
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
