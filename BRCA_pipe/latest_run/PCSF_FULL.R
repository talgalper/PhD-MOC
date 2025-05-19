library(PCSF)
library(tidyverse)
library(edgeR)

# Load data
load("RData/LumA/DE_data.RData")
load("RData/LumB/DE_data.RData")
load("RData/Her2/DE_data.RData")
load("RData/basal/DE_data.RData")

GTEx_data <- read.table("../BRCA_pipe/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
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
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})

control_info <- data.frame(sample = colnames(GTEx_ENS),
                           group = rep("control", ncol(GTEx_ENS)))
lumA_info <- data.frame(sample = colnames(LumA_unstranded),
                        group = rep("lumA", ncol(LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(LumB_unstranded),
                        group = rep("lumB", ncol(LumB_unstranded)))
her2_info <- data.frame(sample = colnames(Her2_unstranded),
                        group = rep("Her2", ncol(Her2_unstranded)))
basal_info <- data.frame(sample = colnames(Basal_unstranded),
                         group = rep("basal", ncol(Basal_unstranded)))
PCA_sample_info <- rbind(control_info, lumA_info, lumB_info, her2_info, basal_info)

# add sample type to sample_info
load("RData/TCGA_query.RData")
common <- common[, c(1,3)]
PCA_sample_info <- merge(PCA_sample_info, common, by.x = "sample", by.y = "cases", all.x = T)
PCA_sample_info$sample_type <- ifelse(is.na(PCA_sample_info$sample_type), "Healthy", PCA_sample_info$sample_type)

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(control_info, lumA_info, lumB_info, her2_info, basal_info)
rm(clinical, common, query_TCGA, subtypes)
gc()

hist(log(as.matrix(all_subtypes)), 
     xlab = "log(raw counts)"
)

# read in functions from "../BRCA_pipe/Rscript/DE_functions.R"
counts_filt <- filter_low_expr(disease_data = all_subtypes, 
                               control_data = GTEx_ENS)

hist(log(as.matrix(counts_filt$counts_filt)),
     xlab = "log(raw counts)"
)

## plot PCA
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
  pca_data <- merge(pca_data, sample_info, by.x = "row.names", by.y = "sample")
  
  # Create a custom colour palette for stages
  library(RColorBrewer)
  groups <- unique(PCA_sample_info[ ,2])
  num_colors <- length(groups)
  colours <- brewer.pal(n = num_colors, name = "Dark2")
  names(colours) <- groups
  
  # Create the PCA plot with color mapping
  library(ggrepel)
  library(ggalt)
  
  if (isTRUE(circle_clust)) {
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
        legend.text = element_text(size = 16),
        panel.grid = element_blank()
      )
  } else {
    PCA_plot <- ggplot(pca_data, aes(PC1, PC2, color = Classification)) +
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

PCA_plot <- plot_PCA(expr_data = counts_filt$counts_filt, 
                     sample_info = PCA_sample_info, 
                     output_plot_data = T,
                     circle_clust = F)


# remove solid tissue normal samples and re-plot PCA
STN_samples <- PCA_sample_info$sample[PCA_sample_info$sample_type == "Solid Tissue Normal"]
# remove STN samples from original data and re-filter low count genes
expr_data_filt <- filter_low_expr(disease_data = all_subtypes[, !colnames(all_subtypes) %in% STN_samples], 
                                  control_data = GTEx_ENS)

PCA_sample_info_filt <- PCA_sample_info[!PCA_sample_info$sample %in% STN_samples, ]

PCA_plot_filt <- plot_PCA(expr_data = expr_data_filt$counts_filt, 
                          sample_info = PCA_sample_info_filt, 
                          output_plot_data = T,
                          circle_clust = F)

save(PCA_plot, PCA_plot_filt, file = "~/OneDrive - RMIT University/PhD/large_git_files/DE_data/PCA_plot_data.RData")
load("~/OneDrive - RMIT University/PhD/large_git_files/DE_data/PCA_plot_data.RData")

## perform DE
sample_info_filt <- expr_data_filt$sample_info[!expr_data_filt$sample_info$sample %in% STN_samples, ]
DE_results <- DE_analysis(counts_matrix = expr_data_filt$counts_filt,
                          sample_info = sample_info_filt)

save(DE_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/DE_data/DE_results_STNfilt.RData")
load("~/OneDrive - RMIT University/PhD/large_git_files/DE_data/DE_results_STNfilt.RData")

# get DE counts summary
print(summary(decideTests(DE_results$qlf, p = 0.05, adjust = "fdr", lfc = 1)))

data <- DE_results$toptags$table
data$PValue[data$PValue == 0] <- min(data$PValue[data$PValue!=0])
logPValue <- -log10(data$PValue)
data$logPValue <- logPValue

library(EnhancedVolcano)
EnhancedVolcano(
  data,
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


# load in DE results from WGCNA where we already did this analysis for BRCA
#load("../WGCNA/BRCA/RData/DE_subset/dif_exp.RData")
dif_exp <- DE_results$dif_exp
DE_data <- subset(dif_exp, select = c("gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values

save(dif_exp, file = "latest_run/RData/STN_filt/dif_exp.RData")

# number of targets in DE data
targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2,4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]
temp <- DE_results$hits[DE_results$hits$gene_id %in% unique(targets$ensembl_gene_id), ]
temp <- merge(targets, temp, by.x = "ensembl_gene_id", by.y = "gene_id")



# convert to gene symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), 
                           filters = "ensembl_gene_id", 
                           values = DE_data$gene_id, 
                           mart = ensembl)
ensembl_converted$description <- gsub("\\[.*?\\]", "", ensembl_converted$description)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- DE_data[!DE_data$gene_id %in% ensembl_converted$ensembl_gene_id, ]

ensembl_converted <- ensembl_converted[ensembl_converted$external_gene_name != "", ]

novel_transcripts <- unmapped[grep("novel transcript", unmapped$description), ]
novel_proteins <- unmapped[grep("novel protein", unmapped$description), ]
pseudogene <- unmapped[grep("pseudogene", unmapped$description), ]

DE_data_geneSymbol <- merge(DE_data, ensembl_converted, by.x = "gene_id", by.y = "ensembl_gene_id")
DE_data_geneSymbol <- subset(DE_data_geneSymbol, select = c("external_gene_name", "logFC", "logFC_abs"))

save(DE_data_geneSymbol, unmapped, unrecognised, file = "latest_run/RData/DE_data_geneSymbol.RData")

# get interaction data from STRING
string_edge_data <- read.table("latest_run/intermediate/STRING network (physical) default edge - FULL.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("latest_run/intermediate/STRING network (physical) default node - FULL.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "display.name"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("display.name.x", "display.name.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

save(final_df, file = "latest_run/RData/STRING_PPI_FULL.RData")
load("latest_run/RData/STRING_PPI_FULL.RData")


# master human STRING network
#STRING <- read.table("latest_run/intermediate/9606.protein.physical.links.v12.0.txt", header = T)
#STRING_full <- read.table("latest_run/intermediate/9606.protein.physical.links.full.v12.0.txt", header = T)
#STRING_info <- read.table("latest_run/intermediate/9606.protein.info.v12.0.txt", header = T, sep = "\t")
#STRING_aliases <- read.table("../../../../OneDrive - RMIT University/PhD/large_git_files/9606.protein.aliases.v12.0.txt", header = T, "\t", fill = T)
#STRING_aliases <- STRING_aliases[STRING_aliases$source == "Ensembl_gene", ]
#
## Replace with preferred name
#STRING_updated <- STRING_full
#STRING_updated <- STRING_updated[STRING_updated$experiments != 0, ]
#STRING_updated <- STRING_updated[STRING_updated$combined_score >= 200, ]
#
#STRING_updated$protein1_gene <- STRING_info$preferred_name[match(STRING_updated$protein1, STRING_info$string_protein_id)]
#STRING_updated$protein2_gene <- STRING_info$preferred_name[match(STRING_updated$protein2, STRING_info$string_protein_id)]
#
#STRING_updated <- STRING_updated[, c("protein1_gene", "protein2_gene", "combined_score")]
#STRING_updated <- na.omit(STRING_updated)



# try PCSF with logFC scores from WGCNA results
#load("../WGCNA/BRCA/RData/all_default/venn_data.RData")
#common_genes <- Reduce(intersect, list(DE_genes, tumour_associated, top_kwithin, top_gene_membership))
#DE_data_subset <- DE_data[DE_data$gene_id %in% common_genes, ]



# set seed for reproducibility 
set.seed(1234)
# construct interactome
ppi <- construct_interactome(final_df)
# set terminals
terminals <- setNames(as.numeric(DE_data_geneSymbol$logFC_abs), DE_data_geneSymbol$external_gene_name)

# run PCSF with random noise
start_time <- Sys.time()
subnet <- PCSF_rand(ppi, terminals, n = 50, r = 0.1, b = 1, w = 2, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

plot.PCSF(subnet, node_label_cex = 15)

save(subnet, file = "latest_run/RData/PCSF_subnet_FULL.RData")
load("latest_run/RData/PCSF_subnet_FULL.RData")


# extract cluster data
clust <- components(subnet)
df <- data.frame(gene_id = names(clust$membership), cluster = factor(clust$membership))
betweenness <- betweenness(subnet) 
centrality <- degree(subnet) 
df$betweenness <- betweenness[as.character(df$gene_id)]
df$degree_centrality <- centrality[as.character(df$gene_id)]
df$betweenness <- as.integer(df$betweenness)
df$degree_centrality <- as.integer(df$degree_centrality)
df$prize <- V(subnet)$prize
df$type <- V(subnet)$type

rownames(df) <- 1:nrow(df)

df <- df[order(-df$degree_centrality), ]
rownames(df) <- NULL


# network enrichment
PCSF_enrich <- enrichment_analysis(subnet)

# Create a data frame with the enrichment results
enrichment_results <- PCSF_enrich$enrichment
enrichment_table <- data.frame(
  Cluster = enrichment_results$Cluster,
  Term = enrichment_results$Term,
  PValue = enrichment_results$P.value,
  Adjusted_Pvalue = enrichment_results$Adjusted.P.value,
  Genes = enrichment_results$Genes)



targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- unique(targets$ensembl_gene_id)

temp <- df[df$ensembl_gene_id %in% targets, ]

targets[!targets %in% df$gene_id]

