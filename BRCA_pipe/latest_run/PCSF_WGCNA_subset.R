library(PCSF)
library(tidyverse)
library(edgeR)

# Load data
load("RData/LumA/DE_data.RData")
load("RData/LumB/DE_data.RData")
load("RData/Her2/DE_data.RData")
load("RData/basal/DE_data.RData")
load("RData/TCGA_normal.RData")

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

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
collectGarbage()

# read in functions from "../BRCA_pipe/Rscript/DE_functions.R"
counts_filt <- filter_low_expr(disease_data = all_subtypes, 
                               control_data = GTEx_ENS)
hist(log(as.matrix(counts_filt$counts_filt)))

DE_results <- DE_analysis(counts_matrix = DE_counts_filt$counts_filt,
                          sample_info = DE_counts_filt$sample_info)

# load in DE results from WGCNA where we already did this analysis for BRCA
load("../WGCNA/BRCA/RData/DE_subset/dif_exp.RData")
DE_data <- subset(dif_exp, select = c("gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values


# try PCSF with logFC scores from WGCNA results
load("../WGCNA/BRCA/RData/all_default/signed/venn_data.RData")
common_genes <- Reduce(intersect, list(DE_genes, tumour_associated, top_kwithin, top_gene_membership))
DE_data <- DE_data[DE_data$gene_id %in% common_genes, ]


# convert to gene symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), 
                           filters = "ensembl_gene_id", 
                           values = DE_data$gene_id, 
                           mart = ensembl)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- DE_data[!DE_data$gene_id %in% ensembl_converted$ensembl_gene_id, ]

ensembl_converted <- ensembl_converted[ensembl_converted$external_gene_name != "", ]

novel_transcripts <- unmapped[grep("novel transcript", unmapped$description), ]
novel_proteins <- unmapped[grep("novel protein", unmapped$description), ]
pseudogene <- unmapped[grep("pseudogene", unmapped$description), ]

DE_data_geneSymbol <- merge(DE_data, ensembl_converted, by.x = "gene_id", by.y = "ensembl_gene_id")
DE_data_geneSymbol <- subset(DE_data_geneSymbol, select = c("external_gene_name", "logFC", "logFC_abs"))


load("latest_run/RData/STRING_PPI_FULL.RData")




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
plot.PCSFe(PCSF_enrich$subnet, edge_width = 8, node_size = 30, node_label_cex = 1)


# Create a data frame with the enrichment results
enrichment_results <- PCSF_enrich$enrichment
enrichment_table <- data.frame(
  Cluster = enrichment_results$Cluster,
  Term = enrichment_results$Term,
  PValue = enrichment_results$P.value,
  Adjusted_Pvalue = enrichment_results$Adjusted.P.value,
  Genes = enrichment_results$Genes)




split_data <- split(enrichment_table, enrichment_table$Cluster)

# Extract the top 3 rows for each cluster based on PValue
top3_rows <- lapply(split_data, function(df) {
  df[order(df$PValue), ][1:3, ]
})

# Combine the results back into a single data frame
top3_rows_df <- do.call(rbind, top3_rows)
