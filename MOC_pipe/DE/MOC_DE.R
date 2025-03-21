### Perform pre-processing and DE analysis on MOC data ###

library(readxl)
library(edgeR)

MOC_raw_counts <- read.csv("data/analysis_set_raw_counts.csv", row.names = 1)
sample_info <- read.csv("data/All survival_CN_Aug18.csv")
#sample_info2 <- read_xlsx("data/sampleInformation.xlsx")

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

# filter low counts
counts_filt <- filterByExpr(MOC_raw_counts, group = sample_info$Classification)
counts_filt <- MOC_raw_counts[counts_filt, ]
low_exp_genes <- MOC_raw_counts[!rownames(MOC_raw_counts) %in% rownames(counts_filt), ]

# plot PCA + cluster circles
plot_PCA <- function(expr_data, sample_info, output_plot_data = T) {
  # convert expression data frame into CPM normalised + transposed matrix
  PCA_data <- cpm(as.matrix(expr_data), log = T)
  PCA_data <- t(PCA_data)
  
  pca <- prcomp(PCA_data, scale. = T, center = T)
  pca_data <- pca$x
  
  pca_var <- pca$sdev^2
  pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)
  
  pca_data <- as.data.frame(pca_data)
  
  # Merge sample mapping with PCA data
  pca_data <- merge(pca_data, sample_info, by.x = "row.names", by.y = "GAMUT_ID")
  
  # Create a custom colour palette for stages
  library(RColorBrewer)
  groups <- unique(sample_info[ ,2])
  num_colors <- length(groups)
  colours <- brewer.pal(n = num_colors, name = "Set1")
  names(colours) <- groups
  
  # Create the PCA plot with color mapping
  library(ggrepel)
  library(ggalt)
  
  PCA_plot <- ggplot(pca_data, aes(PC1, PC2, color = Classification)) +
    geom_point(size = 4) +
    geom_encircle(aes(group = Classification), s_shape = 0, expand = 0.05, color = "black") +
    scale_color_manual(values = colours) +
    theme_bw() +
    labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
         y = paste0('PC2: ', pca_var_perc[2], ' %')) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16)
    )
  
  print(PCA_plot)
  
  if (output_plot_data == T) {
    return(list(PCA_plot = PCA_plot, plot_data = pca_data, pca_var_perc = pca_var_perc))
  }
}

PCA_plot <- plot_PCA(expr_data = counts_filt, 
                     sample_info = sample_info, 
                     output_plot_data = T)


keep_samples <- sample_info$GAMUT_ID[sample_info$Classification %in% c("MOC", "BEN")]
MOC_raw_counts_filt <- MOC_raw_counts[, colnames(MOC_raw_counts) %in% keep_samples]
sample_info_filt <- sample_info[sample_info$GAMUT_ID %in% keep_samples, ]

counts_filt <- filterByExpr(MOC_raw_counts_filt, group = sample_info_filt$Classification)
counts_filt <- MOC_raw_counts_filt[counts_filt, ]
low_exp_genes <- MOC_raw_counts_filt[!rownames(MOC_raw_counts_filt) %in% rownames(counts_filt), ]

hist(cpm(counts_filt, log = T))

PCA_plot <- plot_PCA(expr_data = counts_filt, 
                     sample_info = sample_info, 
                     output_plot_data = T)



# unpaired DE analysis
DE_analysis <- function(counts_matrix, sample_info) {
  
  data <- DGEList(counts = counts_matrix, group = sample_info$Classification)
  
  cat("Normalising library size \n")
  
  data <- normLibSizes(data)
  
  design <- model.matrix(~Classification, data = sample_info)
  
  cat("Estimating dispersion \n")
  
  data <- estimateDisp(data)
  
  cat("Conduct quasi-likelihood F-test \n")
  
  fit <- glmQLFit(data, design)
  qlf <- glmQLFTest(fit)
  toptags <- topTags(qlf, n = Inf)
  
  dif_exp <- decideTests(qlf, p = 0.05, adjust = "fdr", lfc = 1)
  print(summary(dif_exp))
  dif_exp_genes <- rownames(data$counts)[as.logical(dif_exp)]
  hits <- toptags$table[toptags$table$FDR < 0.1, ]
  colnames <- colnames(hits)
  hits$gene_id <- rownames(hits)
  hits <- hits[,c("gene_id", colnames)]
  dif_exp <- hits[dif_exp_genes, ]
  cat("Total DE:", nrow(dif_exp))
  
  return(list(input_data = counts_matrix,
              sample_info = sample_info,
              hits = hits,
              dif_exp = dif_exp,
              toptags = toptags,
              DGEList = data,
              fit = fit,
              qlf = qlf,
              design = design))
}

DE_results <- DE_analysis(counts_matrix = MOC_raw_counts_filt, 
                          sample_info = sample_info_filt)

save(DE_results, file = "DE/DE_results.RData")
load("DE/DE_results.RData")

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
  title = "BRCA DE: TCGA vs GTEx",
  legendPosition = "right",
  drawConnectors = T,
)


# convert to gene symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_id_annot <- function(ensembl, data, col_id = 1, input_type, convert_to) {
  ensembl_annot <- getBM(attributes = c(input_type, convert_to), 
                             filters = input_type, 
                             values = data[, col_id], 
                             mart = ensembl)
  ensembl_annot$description <- gsub("\\[.*?\\]", "", ensembl_annot$description)
  
  data_annot <- merge(ensembl_annot, data, by.x = input_type, by.y = colnames(data)[col_id], all.y = T)
}

DE_result_geneSymbol <- ensembl_id_annot(ensembl, DE_results$dif_exp,
                                         input_type = "ensembl_gene_id",
                                         convert_to = c("external_gene_name", "description", "gene_biotype"))

DE_hits_geneSymbol <- ensembl_id_annot(ensembl, DE_results$hits,
                         input_type = "ensembl_gene_id",
                         convert_to = c("external_gene_name", "description", "gene_biotype"))


write.csv(DE_result_geneSymbol, "DE/MOC_DE_results.csv", row.names = F)
write.csv(DE_hits_geneSymbol, "DE/MOC_DE_hits.csv", row.names = F)



