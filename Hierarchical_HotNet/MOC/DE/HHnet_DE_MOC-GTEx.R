library(PCSF)
library(tidyverse)
library(edgeR)

# read in data
MOC_raw_counts <- read.csv("../WGCNA/rna_seq_data/analysis_set_raw_counts_genenames.csv")
MOC_raw_counts_ENS <- read.csv("../WGCNA/rna_seq_data/analysis_set_raw_counts.csv")
ID_mapping <- data.frame(hgnc_symbol = MOC_raw_counts$X,
                         ensembl_id = MOC_raw_counts_ENS$X)
MOC_raw_counts_ENS <- column_to_rownames(MOC_raw_counts_ENS, "X")

sample_info <- read.csv("../WGCNA/rna_seq_data/All survival_CN_Aug18.csv")
sample_info$GAMUT_ID <- paste0("GAMuT_", sample_info$GAMUT_ID) # match IDs

# subset sample info to only those present in data matrix
MOC_samples <- data.frame(GAMUT_ID = colnames(MOC_raw_counts_ENS))
sample_info <- merge(MOC_samples, sample_info, by = "GAMUT_ID", all.x = T)
sample_info$Stage[sample_info$Grade == "BEN"] <- "BEN"

# remove EOM and BDL samples
sample_info_subset <- sample_info[sample_info$Classification != "EOM" & sample_info$Classification != "BDL", ]
sample_info_subset <- sample_info_subset[!is.na(sample_info_subset$GAMUT_ID), ]

# format for PCA
sample_info_subset <- subset(sample_info_subset, select = c("GAMUT_ID", "Classification"))
colnames(sample_info_subset) <- c("sample", "group")

# subset RNA-seq data
MOC_data <- MOC_raw_counts_ENS[, colnames(MOC_raw_counts_ENS) %in% sample_info_subset$sample[sample_info_subset$group == "MOC"]]
BEN_data <- MOC_raw_counts_ENS[, colnames(MOC_raw_counts_ENS) %in% sample_info_subset$sample[sample_info_subset$group == "BEN"]]

rm(sample_info, sample_info_subset)

# read in GTEx ovarian data
GTEx_data <- read.table("MOC/DE/data/gene_reads_ovary.gct", skip = 2)
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

# QC
counts_filt <- filter_low_expr(disease_data = MOC_data,
                               control_data = GTEx_ENS)

hist(log(as.matrix(counts_filt$counts_filt))) # looks good

# perform DE analysis on MOC vs GTEx ovairan
DE_results <- DE_analysis(counts_matrix = counts_filt$counts_filt,
                          sample_info = counts_filt$sample_info)

dif_exp <- DE_results$dif_exp
save(dif_exp, file = "MOC/DE/MOC_dif_exp.RData")
load("MOC/DE/MOC_dif_exp.RData")



library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), 
                           filters = "ensembl_gene_id", 
                           values = dif_exp$gene_id, 
                           mart = ensembl)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- dif_exp[!dif_exp$gene_id %in% ensembl_converted$ensembl_gene_id, ]

novel_transcripts <- unmapped[grep("novel transcript", unmapped$description), ]
novel_proteins <- unmapped[grep("novel protein", unmapped$description), ]
pseudogene <- unmapped[grep("pseudogene", unmapped$description), ]

ensembl_converted <- ensembl_converted[ensembl_converted$external_gene_name != "", ]

DE_data_geneSymbol <- merge(dif_exp, ensembl_converted, by.x = "gene_id", by.y = "ensembl_gene_id")
DE_data_geneSymbol <- subset(DE_data_geneSymbol, select = c("external_gene_name", "logFC"))

DE_data_geneSymbol$logFC_abs <- abs(DE_data_geneSymbol$logFC) # get absolute values


load("../BRCA_pipe/latest_run/RData/STRING_PPI_FULL.RData")
final_df <- final_df[, -3]

# create index and indexed edge list
gene_index <- as.data.frame(unique(c(final_df$node_1, final_df$node_2)))
gene_index$row_num <- seq.int(nrow(gene_index))
colnames(gene_index) <- c("ensembl_id", "index")
gene_index <- gene_index[c("index", "ensembl_id")]

edge_list_index <- data.frame(from = match(final_df$node_1, gene_index$ensembl_id),
                              to = match(final_df$node_2, gene_index$ensembl_id))

# create score file
DE_data <- subset(DE_data_geneSymbol, select = c("external_gene_name", "logFC"))
DE_data_abs <- subset(DE_data_geneSymbol, select = c("external_gene_name", "logFC_abs"))


write_tsv(gene_index, "MOC/DE/data/gene_index.tsv", col_names = F)
write_tsv(edge_list_index, "MOC/DE/data/edge_list_index.tsv", col_names = F)
write_tsv(final_df, "MOC/DE/data/STRING_network.tsv", col_names = F)
write_tsv(DE_data, "MOC/DE/data/logFC_scores.tsv", col_names = F)
write_tsv(DE_data_abs, "MOC/DE/data/logFC_scores_abs.tsv", col_names = F)




library(tidyverse)
library(igraph)



STRING_net <- read_graph("../../../../OneDrive - RMIT University/PhD/large_git_files/HHnet/STRING_hsa_physical_network.graphml", format = "graphml") # ubuntu
STRING_net <- as.undirected(STRING_net)

hh_results <- read_lines("MOC/DE/results/clusters_STRING_MOC_GTEx_logFC_abs.tsv", skip = 7)
hh_results <- str_split(hh_results, pattern = "\t")

# colour cluster members
V(STRING_net)$color <- ifelse(V(STRING_net)$`display name` %in% hh_results[[1]], "tomato", "white")

# subnetwork cluster
subnet <- induced_subgraph(graph = STRING_net, V(STRING_net)$`display name` %in% hh_results[[1]])
plot.igraph(subnet, asp = 0, vertex.size = 2, edge.arrow.size = 0.3, vertex.label.dist = 1, vertex.label = V(subnet)$`display name`)


write_graph(clust1_net, "MOC/DE/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
write_graph(subnet, "MOC/DE/results/hhnet_cluster1_net.graphml", format = "graphml")



df_subnet <- data.frame(display.name = V(subnet)$`display name`,
                        degree = degree(subnet),
                        betweenness = betweenness(subnet),
                        source = V(subnet)$color)
df_subnet$source <- ifelse(df_subnet$source == "tomato", "subnet", "STRING")
df_subnet <- df_subnet[order(-df_subnet$degree), ]


df_subnetNeighs <- data.frame(display.name = V(clust1_net)$`display name`,
                              degree = degree(clust1_net),
                              betweenness = betweenness(clust1_net),
                              source = V(clust1_net)$color)
df_subnetNeighs$source <- ifelse(df_subnetNeighs$source == "tomato", "subnet", "STRING")
df_subnetNeighs <- df_subnetNeighs[order(-df_subnetNeighs$degree), ]



