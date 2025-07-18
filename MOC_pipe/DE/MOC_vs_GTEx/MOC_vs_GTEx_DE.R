### Differential expression analysis between MOC and GTEx-OV ###

library(edgeR)

load("data/serous-OV/GTEx-OV_unstranded.RData")
MOC_raw_counts <- read.csv("data/analysis_set_raw_counts.csv", row.names = 1)
sample_info <- read.csv("data/All survival_CN_Aug18.csv")

colnames(MOC_raw_counts) <- sub("GAMuT_", "", colnames(MOC_raw_counts))

sample_info <- sample_info[sample_info$GAMUT_ID %in% colnames(MOC_raw_counts), ]
sample_info <- sample_info[, c(1,2,4,5)]
missing_samples <- setdiff(colnames(MOC_raw_counts), sample_info$GAMUT_ID)

sample_info <- rbind(sample_info, data.frame(GAMUT_ID = missing_samples,
                                             Classification = "UNK", 
                                             Grade = "UNK", 
                                             Stage = "UNK"))

sample_info <- sample_info[match(colnames(MOC_raw_counts), sample_info$GAMUT_ID), ]

# consolidate stages into single identifiers
sample_info$stage <- ifelse(sample_info$Stage %in% c("I", "IA", "IC"), "I",
                            ifelse(sample_info$Stage %in% c("II", "IIB"), "II",
                                   ifelse(sample_info$Stage %in% c("III", "IIIA", "IIIc", "IIIC"), "III",
                                          ifelse(sample_info$Stage %in% "IV", "IV",
                                                 ifelse(sample_info$Classification == "BEN", "BEN", "UNK")))))


colnames(sample_info)[1] <- "sample"
sample_info <- rbind(
  sample_info[,c(1,2,5)],
  data.frame(sample = colnames(GTEx_data),
             Classification = rep("GTEx-OV", ncol(GTEx_data)),
             stage = rep("Healthy", ncol(GTEx_data)))
)


all_expr_data <- merge(MOC_raw_counts, GTEx_data, by = "row.names")
all_expr_data <- tibble::column_to_rownames(all_expr_data, "Row.names")

keep_samples <- sample_info$sample[sample_info$Classification %in% c("MOC", "GTEx-OV")]
all_expr_data <- all_expr_data[, colnames(all_expr_data) %in% keep_samples]
sample_info <- sample_info[sample_info$sample %in% keep_samples, ]

# filter low counts
counts_filt <- filterByExpr(all_expr_data, group = sample_info$Classification)
counts_filt <- all_expr_data[counts_filt, ]
low_exp_genes <- all_expr_data[!rownames(all_expr_data) %in% rownames(counts_filt), ]


PCA_plot <- plot_PCA(expr_data = counts_filt, 
                     sample_info = sample_info, 
                     output_plot_data = T,
                     circle_clust = F,
                     shape = "stage", 
                     colour = "Classification")

# change colour of GTEx group to green
print(PCA_plot$PCA_plot + scale_color_manual(name = "Group", 
                                             values = c("#4DAF4A", "#E41A1C"), 
                                             labels = c("GTEx-OV", "MOC")))


# perform DE analysis
DE_results <- DE_analysis(counts_matrix = counts_filt, 
                          sample_info = sample_info)

save(DE_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/MOC/DE_results_MOC_vs_GTEx.RData")
load("~/OneDrive - RMIT University/PhD/large_git_files/MOC/DE_results_MOC_vs_GTEx.RData")

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


write.csv(DE_result_geneSymbol, "DE/MOC_vs_GTEx/MOC_vs_GTEx_DE_results.csv", row.names = F)
write.csv(DE_hits_geneSymbol, "DE/MOC_vs_GTEx/MOC_vs_GTEx_DE_hits.csv", row.names = F)
