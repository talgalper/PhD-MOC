### Differential expression analysis between TCGA-OV and GTEx-OV ###

library(edgeR)

load("data/serous-OV/GTEx-OV_unstranded.RData")
load("data/serous-OV/TCGA-OV_unstranded.RData")
MOC_raw_counts <- read.csv("data/analysis_set_raw_counts.csv", row.names = 1)
common_genes <- intersect(rownames(MOC_raw_counts), intersect(rownames(TCGA_OV_data_unstranded), rownames(GTEx_data)))

all_expr_data <- merge(TCGA_OV_data_unstranded, GTEx_data, by = "row.names")
all_expr_data <- tibble::column_to_rownames(all_expr_data, "Row.names")

all_expr_data <- all_expr_data[rownames(all_expr_data) %in% common_genes, ]

sample_info <- data.frame(sample = c(colnames(TCGA_OV_data_unstranded), colnames(GTEx_data)),
                          Classification = c(rep("TCGA-OV", ncol(TCGA_OV_data_unstranded)),
                                             rep("GTEx-OV", ncol(GTEx_data))))

# filter low counts
counts_filt <- filterByExpr(all_expr_data, group = sample_info$Classification)
counts_filt <- all_expr_data[counts_filt, ]
low_exp_genes <- all_expr_data[!rownames(all_expr_data) %in% rownames(counts_filt), ]

# plot PCA
PCA_plot <- plot_PCA(expr_data = counts_filt, 
                     sample_info = sample_info, 
                     output_plot_data = T,
                     circle_clust = F)

# perform DE analysis
DE_results <- DE_analysis(counts_matrix = counts_filt, 
                          sample_info = sample_info)

save(DE_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/MOC/DE_results_TCGA_vs_GTEx.RData")
load("~/OneDrive - RMIT University/PhD/large_git_files/MOC/DE_results_TCGA_vs_GTEx.RData")

# get DE counts summary
print(summary(decideTests(DE_results$qlf, p = 0.05, adjust = "fdr", lfc = 1)))


# create volcano plot
plot_data <- DE_results$toptags$table
plot_data$PValue[plot_data$PValue == 0] <- min(plot_data$PValue[plot_data$PValue!=0])
plot_data$logPValue <- -log10(plot_data$PValue)

library(EnhancedVolcano)

EnhancedVolcano(
  plot_data,
  lab = NA,
  x = "logFC",
  y = "PValue",
  labSize = 3,
  pCutoff = 1e-02,
  FCcutoff = 1,
  title = "MOC vs GTEx",
  legendPosition = "right",
  drawConnectors = T,
)


# convert to gene symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


DE_result_geneSymbol <- id_annot(ensembl, DE_results$dif_exp,
                                 input_type = "ensembl_gene_id",
                                 convert_to = c("external_gene_name", "description", "gene_biotype"))

DE_hits_geneSymbol <- id_annot(ensembl, DE_results$hits,
                               input_type = "ensembl_gene_id",
                               convert_to = c("external_gene_name", "description", "gene_biotype"))


write.csv(DE_result_geneSymbol, "DE/TCGA_vs_GTEx/TCGA_vs_GTEx_DE_results.csv", row.names = F)
write.csv(DE_hits_geneSymbol, "DE/TCGA_vs_GTEx/TCGA_vs_GTEx_DE_hits.csv", row.names = F)
