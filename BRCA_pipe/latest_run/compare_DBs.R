library(TCGAbiolinks)
library(ggVennDiagram)
library(biomaRt)
library(data.table)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


load("latest_run/RData/STN_filt/dif_exp.RData")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
                           filters = "ensembl_gene_id", 
                           values = dif_exp$gene_id, 
                           mart = ensembl)
ensembl_converted$description <- gsub("\\[.*?\\]", "", ensembl_converted$description)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- dif_exp[!dif_exp$gene_id %in% ensembl_converted$ensembl_gene_id, ]

dif_exp <- merge.data.table(ensembl_converted, dif_exp, by.x = "ensembl_gene_id", by.y = "gene_id", all.y = T)


GEPIA2_DE <- read.table("../../../../Downloads/table_degenes.txt", sep = "\t", header = T)

targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2,4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]



load("../../../../OneDrive - RMIT University/PhD/large_git_files/DE_data/DE_results_STNfilt.RData")
hits <- DE_results$hits
hits <- hits[hits$gene_id %in% targets$ensembl_gene_id, ]
hits <- merge(targets, hits, by.x = "ensembl_gene_id", by.y = "gene_id")

temp <- merge.data.table(targets, hits, by = "ensembl_gene_id", all.x = T)



OncoDB <- read.table("../../../../Downloads/expression_diff_BRCA.txt", header = T)


Xena <- read.csv("../../../../Downloads/DEG_results_Primary_Tumor vs. Normal_Tissue.csv", row.names = 1)
Xena <- Xena[Xena$logFC >= 1 | Xena$logFC <= -1, ]




venn_data <- list(GEPIA2 = GEPIA2_DE$Gene.Symbol,
                  customAnalysis = dif_exp$external_gene_name[dif_exp$external_gene_name != ""],
                  OncoDB = OncoDB$Gene,
                  Xena = rownames(Xena))

ggVennDiagram(venn_data)


BRCA_markers <- c("ERBB2", "TP53", "BRCA1", "BRCA2", "ESR1", "PGR", "PIK3CA", "PTEN",
                  "MYC", "GATA3", "CDH1", "AKT1", "CCND1", "FGFR1", "MUC1", "EGFR",
                  "NOTCH1", "KRAS", "BCL2", "MAPK1")

DE_datasets <- dif_exp
DE_datasets <- DE_datasets[DE_datasets$external_gene_name != "", ]

DE_datasets <- DE_datasets[ ,c(2,5)]
colnames(DE_datasets) <- c("gene_id", "MyData")
DE_datasets <- merge(DE_datasets, GEPIA2_DE, by.x = "gene_id", by.y = "Gene.Symbol", all = T)
DE_datasets <- DE_datasets[ ,-c(3,4,5,7)]
colnames(DE_datasets)[3] <- c("GEPIA2")
DE_datasets <- merge(DE_datasets, OncoDB, by.x = "gene_id", by.y = "Gene", all = T)
DE_datasets <- DE_datasets[ ,-c(4,5,6)]
colnames(DE_datasets)[4] <- c("OncoDB")
DE_datasets <- merge(DE_datasets, Xena, by.x = "gene_id", by.y = "row.names", all = T)
DE_datasets <- DE_datasets[ ,-c(6:10)]
colnames(DE_datasets)[5] <- c("Xena")

temp <- merge.data.table(targets, DE_datasets, by.x = "drugBank_target", by.y = "gene_id", all.x = T)
temp2 <- DE_datasets[DE_datasets$gene_id %in% BRCA_markers, ]

library(tidyverse)
temp <- DE_datasets
# remove rows with all NAs
temp <- temp[-c(which(rowSums(is.na(temp)) == ncol(temp))), ]
table(is.na(temp$gene_id))

temp[is.na(temp)] <- 0
temp <- temp[!duplicated(temp$gene_id), ]
rownames(temp) <- NULL
temp <- column_to_rownames(temp, "gene_id")
cor_matrix_pearson <- cor(temp, use = "pairwise.complete.obs", method = "pearson")
cor_matrix_spearman <- cor(temp, use = "pairwise.complete.obs", method = "spearman")



# GSEA analysis
DE_data_geneSymbol <- DE_data_geneSymbol[order(-DE_data_geneSymbol$logFC), ]
GSEA <- setNames(DE_data_geneSymbol$logFC, DE_data_geneSymbol$external_gene_name)

library(msigdbr)
library(fgsea)
msigdb_genesets <- msigdbr(species = "Homo sapiens", category = "H") # "H" for hallmark gene sets
pathways_list <- split(msigdb_genesets$gene_symbol, msigdb_genesets$gs_name)

set.seed(1234) 
fgsea_results <- fgsea(pathways = pathways_list,
                       stats = GSEA)

topPathways <- fgsea_results[order(fgsea_results$pval), ][1:10, ]

# Plot enrichment scores
ggplot(topPathways, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Top Enriched Pathways") +
  theme_minimal()

fgsea_results$leadingEdge_genes <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ", "))







samples <- DE_results$data$samples
samples <- rownames_to_column(samples)

data_tumour <- samples$rowname[samples$group == "disease"]
data_healthy <- samples$rowname[samples$group == "control"]

data <- DE_results$data$counts
data_tumour <- data[, colnames(data) %in% data_tumour]
data_healthy <- data[, colnames(data) %in% data_healthy]

data_tumour <- data_tumour[rownames(data_tumour) %in% "ENSG00000141510", ]
data_healthy <- data_healthy[rownames(data_healthy) %in% "ENSG00000141510", ]

hist(cpm(data_tumour, log = T), main = "tumour")
hist(cpm(data_healthy, log = T), main = "healthy")


summary(cpm(data_healthy, log = T))
summary(cpm(data_tumour, log = T))


hist(cpm(data_tumour, log = T))
hist(cpm(data_healthy, log = T))

