library(TCGAbiolinks)

# extract all TCGA-OV RNA-seq samples
TCGA_query <- GDCquery(project = "TCGA-OV",
                       access = "open", 
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq")

query_output <- getResults(TCGA_query)

# get clinical data
clinical_query <- GDCquery_clinic(project = "TCGA-OV",
                            type = "clinical")

clinical <- merge(query_output, clinical, by.x = "cases.submitter_id", by.y = "submitter_id")
clinical <- subset(clinical, select = c("cases", "cases.submitter_id", "figo_stage", 
                                                    "tissue_or_organ_of_origin", "sample_type"))

table(clinical_query$figo_stage)


# download data 
GDCdownload(TCGA_query, directory = "data/serous-OV/")

save(TCGA_query, file = "data/serous-OV/TCGA-OV_query.RData")
write.csv(clinical, "data/serous-OV/TCGA-OV_clinical.csv", row.names = F)

TCGA_OV_data <- GDCprepare(TCGA_query, summarizedExperiment = TRUE, directory = "data/TCGA-OV/")
TCGA_OV_data_unstranded <- assay(TCGA_OV_data, "unstranded")
rownames(TCGA_OV_data_unstranded) <- gsub("\\.\\d+", "", rownames(TCGA_OV_data_unstranded))
TCGA_OV_data_unstranded <- as.data.frame(TCGA_OV_data_unstranded)
save(TCGA_OV_data_unstranded, file = "data/serous-OV/TCGA-OV_unstranded.RData")
save(TCGA_OV_data, file = "~/OneDrive - RMIT University/PhD/large_git_files/TCGA-OV/TCGA-OV_sumExp.RData")

# read in GTEx ovarian healthy expression data
GTEx_raw <- read.table("data/serous-OV/gene_reads_v10_ovary.gct", skip = 2, header = T)
GTEx_data <- GTEx_raw[, -2]
GTEx_data <- column_to_rownames(GTEx_data, var = "Name")
rownames(GTEx_data) <- gsub("\\.\\d+", "", rownames(GTEx_data))
save(GTEx_data, file = "data/serous-OV/GTEx-OV_unstranded.RData")

## comapre MOC to TCGA-OV + GTEx-OV
library(SummarizedExperiment)
library(tidyverse)

# read in expression data
load("data/serous-OV/TCGA-OV_unstranded.RData")
load("data/serous-OV/GTEx-OV_unstranded.RData")
MOC_raw_counts <- read.csv("data/analysis_set_raw_counts.csv", row.names = 1)

# create sample info 
sample_info <- read.csv("data/All survival_CN_Aug18.csv")
sample_info <- sample_info[, c(1,2,4,5)]
sample_info$GAMUT_ID <- paste0("GAMuT_", sample_info$GAMUT_ID)
sample_info <- sample_info[sample_info$GAMUT_ID %in% colnames(MOC_raw_counts), ]
missing_samples <- setdiff(colnames(MOC_raw_counts), sample_info$GAMUT_ID)
sample_info <- rbind(sample_info, data.frame(GAMUT_ID = missing_samples,
                                             Classification = "UNK", 
                                             Grade = "UNK", 
                                             Stage = "UNK"))
sample_info <- sample_info[, c(1,2)]
colnames(sample_info) <- c("sample", "group")

sample_info <- rbind(
  sample_info,
  data.frame(sample = colnames(TCGA_OV_data_unstranded),
             group = rep("TCGA-OV", ncol(TCGA_OV_data_unstranded))),
  data.frame(sample = colnames(GTEx_data),
             group = rep("GTEx-OV", ncol(GTEx_data)))
)

all_expr_data <- merge(MOC_raw_counts, TCGA_OV_data_unstranded, by = "row.names")
all_expr_data <- column_to_rownames(all_expr_data, var = "Row.names")
all_expr_data <- merge(all_expr_data, GTEx_data, by = "row.names")
all_expr_data <- column_to_rownames(all_expr_data, var = "Row.names")

# filter low counts
counts_filt <- filterByExpr(all_expr_data, group = sample_info$group)
counts_filt <- all_expr_data[counts_filt, ]
low_exp_genes <- all_expr_data[!rownames(all_expr_data) %in% rownames(counts_filt), ]


# plot PCA + cluster circles
plot_PCA <- function(expr_data, sample_info, output_plot_data = T, circle_clust = F) {
  # convert expression data frame into CPM normalised + transposed matrix
  PCA_data <- cpm(as.matrix(expr_data), log = T)
  PCA_data <- t(PCA_data)
  
  pca <- prcomp(PCA_data, scale. = T, center = T)
  pca_data <- pca$x
  
  pca_var <- pca$sdev^2
  pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)
  
  pca_data <- as.data.frame(pca_data)
  
  # Merge sample mapping with PCA data
  pca_data <- merge(pca_data, sample_info, by.x = "row.names", by.y = "sample")
  
  # Create a custom colour palette for stages
  library(RColorBrewer)
  groups <- unique(sample_info[ ,2])
  num_colors <- length(groups)
  colours <- brewer.pal(n = num_colors, name = "Set1")
  names(colours) <- groups
  
  # Create the PCA plot with color mapping
  library(ggrepel)
  library(ggalt)
  
  if (isTRUE(circle_clust)) {
    PCA_plot <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
      geom_point(size = 4) +
      geom_encircle(aes(group = group), s_shape = 0, expand = 0.05, color = "black") +
      scale_color_manual(values = colours) +
      theme_bw() +
      labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
           y = paste0('PC2: ', pca_var_perc[2], ' %')) +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid = element_blank()
      )
  } else {
    PCA_plot <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
      geom_point(size = 4) +
      scale_color_manual(values = colours) +
      theme_bw() +
      labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
           y = paste0('PC2: ', pca_var_perc[2], ' %')) +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid = element_blank()
      )
  }

  print(PCA_plot)
  
  if (output_plot_data == T) {
    return(list(PCA_plot = PCA_plot, plot_data = pca_data, pca_var_perc = pca_var_perc))
  }
}

PCA_plot <- plot_PCA(expr_data = counts_filt, 
                     sample_info = sample_info, 
                     output_plot_data = T,
                     circle_clust = F)

save(PCA_plot, file = "RData/master-OV_PCA.RData")


# add labels to MOC samples
groups <- unique(sample_info[ ,2])
num_colors <- length(groups)
colours <- brewer.pal(n = num_colors, name = "Set1")
names(colours) <- groups

ggplot(PCA_plot$plot_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  geom_text_repel(
    data = subset(PCA_plot$plot_data, group == "MOC"),
    aes(label = Row.names),
    size = 5,
    show.legend = FALSE) +
  scale_color_manual(values = colours) +
  theme_bw() +
  labs(x = paste0('PC1: ', PCA_plot$pca_var_perc[1], ' %'),
       y = paste0('PC2: ', PCA_plot$pca_var_perc[2], ' %')) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    panel.grid = element_blank())








