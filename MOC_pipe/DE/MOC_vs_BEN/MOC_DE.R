### Perform pre-processing and DE analysis on MOC data ###

library(edgeR)

MOC_raw_counts <- read.csv("data/analysis_set_raw_counts.csv", row.names = 1)
sample_info <- read.csv("data/All survival_CN_Aug18.csv")

# change colnames to match sample info
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


# check distribution
hist(cpm(MOC_raw_counts, log = T))

old <- par(
  mar  = c(5, 4, 4, 8) + 0.1,
  xpd  = TRUE
)
stage_order <- unique(sample_info$stage[order(sample_info$stage)])
ord <- order(sample_info$stage)
colours <- RColorBrewer::brewer.pal(n = length(stage_order), name = "Set1")
names(colours) <- stage_order
par(font.lab = 2,
    mar = c(6,4,3,8))
boxplot(
  log2(MOC_raw_counts[, ord] + 1),
  las   = 2,
  col   = colours[as.character(sample_info$stage[ord])],
  #main  = "MOC Raw Counts (logCPM)",
  ylab  = "log2(counts + 1)",
  cex.lab = 1.2,
  cex.axis = 1.2)
legend(
  "topright",
  inset = c(-0.225, 0),
  legend = stage_order,
  fill = colours[stage_order],
  cex = 1.2,
  bty = "n",
  x.intersp = 0.2)
rm(colours, ord, stage_order, old)

# filter low counts
counts_filt <- filterByExpr(MOC_raw_counts, group = sample_info$Classification)
counts_filt <- MOC_raw_counts[counts_filt, ]
low_exp_genes <- MOC_raw_counts[!rownames(MOC_raw_counts) %in% rownames(counts_filt), ]

# plot PCA - read in from R/functions.R
PCA_plot <- plot_PCA(expr_data = counts_filt, 
                     sample_info = sample_info, 
                     output_plot_data = T,
                     circle_clust = F,
                     shape = "stage")

# remove UNK, BDL and EOM samples
keep_samples <- sample_info$GAMUT_ID[sample_info$Classification %in% c("MOC", "BEN")]
MOC_raw_counts_filt <- MOC_raw_counts[, colnames(MOC_raw_counts) %in% keep_samples]
sample_info_filt <- sample_info[sample_info$GAMUT_ID %in% keep_samples, ]

counts_filt <- filterByExpr(MOC_raw_counts_filt, group = sample_info_filt$Classification)
counts_filt <- MOC_raw_counts_filt[counts_filt, ]
low_exp_genes <- MOC_raw_counts_filt[!rownames(MOC_raw_counts_filt) %in% rownames(counts_filt), ]

hist(cpm(counts_filt, log = T))

old <- par(
  mar  = c(5, 4, 4, 8) + 0.1, 
  xpd  = TRUE    
)
stage_order <- unique(sample_info_filt$stage[order(sample_info_filt$stage)])
ord <- order(sample_info_filt$stage)
colours <- RColorBrewer::brewer.pal(n = length(stage_order), name = "Set1")
names(colours) <- stage_order
par(font.lab = 2,
    mar = c(6,4,3,8))
boxplot(
  log2(counts_filt[, ord] + 1),
  las = 2,
  col = colours[as.character(sample_info_filt$stage[ord])],
  #main = "MOC Raw Counts",
  ylab = "log2(counts + 1)",
  cex.lab = 1.2,
  cex.axis = 1.2
  )
legend(
  "topright",
  inset = c(-0.225, 0),
  legend = stage_order,
  fill = colours[stage_order],
  cex = 1.2,
  bty = "n",
  x.intersp = 0.2)
rm(colours, ord, stage_order, old)

PCA_plot <- plot_PCA(expr_data = counts_filt, 
                     sample_info = sample_info_filt, 
                     output_plot_data = T,
                     circle_clust = F,
                     shape = 'stage')


# perform DE analysis - read in from R/funcitons.R
DE_results <- DE_analysis(counts_matrix = counts_filt, 
                          sample_info = sample_info_filt)

save(DE_results, file = "DE/MOC_vs_BEN/DE_results.RData")
load("DE/MOC_vs_BEN/DE_results.RData")

# get DE counts summary
print(summary(decideTests(DE_results$qlf, p = 0.05, adjust = "fdr", lfc = 1)))
cat("Total DE:", nrow(DE_results$dif_exp))

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
  title = "MOC vs BEN",
  legendPosition = "right",
  drawConnectors = T,
)


# convert to gene symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# read in from R/functions.R
DE_result_geneSymbol <- id_annot_2(ensembl, DE_results$dif_exp,
                                 input_type = "ensembl_gene_id",
                                 convert_to = c("external_gene_name", "description", "gene_biotype"))

DE_hits_geneSymbol <- id_annot_2(ensembl, DE_results$hits,
                               input_type = "ensembl_gene_id",
                               convert_to = c("external_gene_name", "description", "gene_biotype"))


write.csv(DE_result_geneSymbol, "DE/MOC_vs_BEN/MOC_DE_results.csv", row.names = F)
write.csv(DE_hits_geneSymbol, "DE/MOC_vs_BEN/MOC_DE_hits.csv", row.names = F)


