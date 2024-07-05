

# number of each stage in samples
stage_freq_df <- table(common$ajcc_pathologic_stage)
stage_freq_df <- as.data.frame(stage_freq_df)

ggplot(stage_freq_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Stage", y = "No. of Samples") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())



selected_barcodes <- common[common$Subtype_Selected == "BRCA.Her2", ]

stage_freq_df <- table(selected_barcodes$ajcc_pathologic_stage)
stage_freq_df <- as.data.frame(stage_freq_df)

ggplot(stage_freq_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Stage", y = "No. of Samples") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())





gene_expression <- query_output[query_output$type == "gene_expression", ]
mirna_expression <- query_output[query_output$type == "mirna_expression", ]

common_expression <- intersect(mirna_expression$id, gene_expression$id)

gene_expression_common <- gene_expression[gene_expression$id %in% common_expression, ]
mirna_expression_common <- mirna_expression[mirna_expression$id %in% common_expression, ]




# Load necessary libraries
library(VennDiagram)

# Create subsets based on different categories
gene_exp_ids <- unique(clinical_query$id[clinical_query$type == "gene_expression"])
mirna_exp_ids <- unique(clinical_query$id[clinical_query$type == "mirna_expression"])

gene_exp_data_ids <- unique(clinical_query$id[clinical_query$data_type == "Gene Expression Quantification"])
isoform_exp_data_ids <- unique(clinical_query$id[clinical_query$data_type == "Isoform Expression Quantification"])
mirna_exp_data_ids <- unique(clinical_query$id[clinical_query$data_type == "miRNA Expression Quantification"])

mirna_seq_ids <- unique(clinical_query$id[clinical_query$experimental_strategy == "miRNA-Seq"])
rna_seq_ids <- unique(clinical_query$id[clinical_query$experimental_strategy == "RNA-Seq"])

stage_I_ids <- unique(clinical_query$id[clinical_query$ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB")])
stage_II_ids <- unique(clinical_query$id[clinical_query$ajcc_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB")])
stage_III_ids <- unique(clinical_query$id[clinical_query$ajcc_pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")])
stage_IV_ids <- unique(clinical_query$id[clinical_query$ajcc_pathologic_stage == "Stage IV"])
stage_X_ids <- unique(clinical_query$id[clinical_query$ajcc_pathologic_stage == "Stage X"])

# Create Venn diagram
venn_list <- list(
  Gene_Expression_Data = gene_exp_data_ids,
  Isoform_Expression_Data = isoform_exp_data_ids,
  miRNA_Expression_Data = mirna_exp_data_ids,
  miRNA_Seq = mirna_seq_ids,
  RNA_Seq = rna_seq_ids
)

venn_list_stage <- list(
  Stage_I = stage_I_ids,
  Stage_II = stage_II_ids,
  Stage_III = stage_III_ids,
  Stage_IV = stage_IV_ids,
  Stage_X = stage_X_ids
)


venn.diagram(
  x = venn_list,
  category.names = names(venn_list),
  filename = "venn.png",
  output = TRUE,
  disable.logging = TRUE
)


gene_data <- gene_data[order(-gene_data$logFC), ]
upreg <- gene_data[1:20, ]
dnreg <- gene_data[(nrow(gene_data)-19):nrow(gene_data), ]
top_dif_exp <- rbind(upreg, dnreg)





# Luminal A Subtype
luminal_a <- data.frame(
  Gene_Target = c("ESR1", 
                  "PGR", 
                  "ERBB2", 
                  "CDK4",
                  "CDK6",
                  "MTOR",
                  "AKT1",
                  "ERK",
                  "SRC",
                  "FGFR1",
                  "FGFR2",
                  "FGFR3",
                  "FGFR4"),
  Example_Drugs = c("Tamoxifen, Fulvestrant", 
                    "Tamoxifen, Fulvestrant", 
                    "Trastuzumab, Pertuzumab", 
                    "Palbociclib, Ribociclib, Abemaciclib",
                    "Palbociclib, Ribociclib, Abemaciclib",
                    "Everolimus",
                    "Capivasertib, Ipatasertib",
                    "Trametinib, Cobimetinib",
                    "Dasatinib, Bosutinib",
                   "Erdafitinib, Infigratinib",
                   "Erdafitinib, Infigratinib",
                   "Erdafitinib, Infigratinib",
                   "Erdafitinib, Infigratinib")
)

LumA_genes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description", "uniprot_gn_id"), 
                       filters = "external_gene_name", 
                       values = luminal_a$Gene_Target, 
                       mart = ensembl)
LumA_genes$description <- gsub("\\s*\\[.*?\\]", "", LumA_genes$description)

luminal_a <- merge(LumA_genes, luminal_a, by.x = "external_gene_name", by.y = "Gene_Target")



# Luminal B Subtype
luminal_b <- data.frame(
  Gene_Target = c("ESR1", 
                  "PGR", 
                  "ERBB2", 
                  "CDK4",
                  "CDK6",
                  "MTOR",
                  "FGFR1",
                  "FGFR2",
                  "FGFR3",
                  "FGFR4",
                  "AKT1",
                  "ERK"),
  Example_Drugs = c("Tamoxifen, Fulvestrant", 
                    "Tamoxifen, Fulvestrant", 
                    "Trastuzumab, Pertuzumab", 
                    "Palbociclib, Ribociclib, Abemaciclib",
                    "Palbociclib, Ribociclib, Abemaciclib",
                    "Everolimus",
                    "Erdafitinib, Infigratinib",
                    "Erdafitinib, Infigratinib",
                    "Erdafitinib, Infigratinib",
                    "Erdafitinib, Infigratinib",
                    "Capivasertib, Ipatasertib")
)

LumB_genes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), 
                    filters = "external_gene_name", 
                    values = luminal_b$Gene_Target, 
                    mart = ensembl)
LumB_genes$description <- gsub("\\s*\\[.*?\\]", "", LumB_genes$description)

luminal_b <- merge(LumB_genes, luminal_b, by.x = "external_gene_name", by.y = "Gene_Target")


# HER2-Enriched Subtype
her2 <- data.frame(
  Gene_Target = c("ERBB2", 
                  "PI3K", 
                  "MTOR", 
                  "CDK4",
                  "CDK6"),
  Example_Drugs = c("Trastuzumab, Pertuzumab, Lapatinib", 
                    "Alpelisib, Buparlisib", 
                    "Everolimus", 
                    "Palbociclib, Ribociclib, Abemaciclib",
                    "Palbociclib, Ribociclib, Abemaciclib")
)

her2_genes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), 
                    filters = "external_gene_name", 
                    values = her2$Gene_Target, 
                    mart = ensembl)
her2_genes$description <- gsub("\\s*\\[.*?\\]", "", her2_genes$description)

her2 <- merge(her2_genes, her2, by.x = "external_gene_name", by.y = "Gene_Target")

# Basal-like Subtype (Triple-negative Breast Cancer)
basal_like <- data.frame(
  Gene_Target = c("PARP1", 
                  "EGFR", 
                  "CD274", 
                  "AR"),
  Example_Drugs = c("Olaparib, Talazoparib", 
                    "Cetuximab, Gefitinib, Erlotinib", 
                    "Atezolizumab, Pembrolizumab, Durvalumab", 
                    "Enzalutamide, Bicalutamide")
)

basal_genes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), 
                    filters = "external_gene_name", 
                    values = basal_like$Gene_Target, 
                    mart = ensembl)
basal_genes$description <- gsub("\\s*\\[.*?\\]", "", basal_genes$description)

basal_like <- merge(basal_genes, basal_like, by.x = "external_gene_name", by.y = "Gene_Target")

save(luminal_a, luminal_b, her2, basal_like, file = "RData/gene_targets.RData")



library(geneSynonym)

lumA_synon <- humanSyno(c("ESR1", "PGR", "ERBB2", "CDK4", "CDK6", 
                          "PI3K", "MTOR", "FGFR1", "FGFR2", "FGFR3", 
                          "FGFR4", "AKT", "ERK", "SRC"))

Basal_synon <- humanSyno(c("PARP", "EGFR", "PDL1", "AR", "PI3K", 
                           "FGFR1", "FGFR2", "FGFR3", "FGFR4", "AKT",
                           "SRC", "MEK", "MTOR"))


# find gene synonyms in ranked gene list
df <- final_gene_counts[0,]

for (i in seq_along(lumA_synon)) {
  i <- lumA_synon[i]
  syno <- unlist(i)
  syno_matched <- final_gene_counts[final_gene_counts$external_gene_name %in% syno, ]
  df <- rbind(df, syno_matched)
}


df <- gene_data[0,]

for (i in seq_along(lumA_synon)) {
  i <- lumA_synon[i]
  syno <- unlist(i)
  syno_matched <- gene_data[gene_data$external_gene_name %in% syno, ]
  df <- rbind(df, syno_matched)
}




synon <- humanSyno(c("ESR1", "PGR", "ERBB2", "CDK4", "CDK6", 
                     "PI3K", "MTOR", "FGFR1", "FGFR2", "FGFR3", 
                     "FGFR4", "AKT", "ERK", "SRC", "PARP", 
                     "PD-L1", "MEK", "ERBB3", "AR"))

synon <- melt(synon)
colnames(synon) <- c("synonym", "NCBI_gene", "input_term")




library(TCGAbiolinks)
library(SummarizedExperiment)

load("RData/TCGA_query.RData")

clinical <- GDCquery_clinic(project = "TCGA-BRCA",
                            type = "clinical")
clinical_query <- merge(query_output, clinical, by.x = "cases.submitter_id", by.y = "submitter_id")
clinical_query <- subset(clinical_query, select = c("cases", "cases.submitter_id", "ajcc_pathologic_stage", 
                                                    "tissue_or_organ_of_origin", "sample_type", "bcr_patient_barcode"))

subtypes <- PanCancerAtlas_subtypes()

master <- merge(clinical_query, subtypes, by.x = "cases", by.y = "pan.samplesID")
master <- subset(master, select = c("cases", "cases.submitter_id", "Subtype_Selected", "sample_type", 
                                    "ajcc_pathologic_stage", "tissue_or_organ_of_origin", "bcr_patient_barcode"))



table(duplicated(clinical_query$bcr_patient_barcode))
table(is.na(clinical_query$ajcc_pathologic_stage))
table(master$Subtype_Selected)

normals <- master[master$Subtype_Selected == "BRCA.Normal", ]
lumA <- master[master$Subtype_Selected == "BRCA.LumA", ]
lumB <- master[master$Subtype_Selected == "BRCA.LumB", ]
Basal <- master[master$Subtype_Selected == "BRCA.Basal", ]
Her2 <- master[master$Subtype_Selected == "BRCA.Her2", ]


paired_samples <- as.data.frame(subset(normals, select = c("bcr_patient_barcode")))
colnames(paired_samples)[1] <- "normal"

paired_samples <- merge(paired_samples, lumA, by.x = "normal", by.y = "bcr_patient_barcode", all = T)
paired_samples <- subset(paired_samples, select = c("normal", "ajcc_pathologic_stage"))
colnames(paired_samples)[2] <- "lumA"

paired_samples <- merge(paired_samples, lumB, by.x = "normal", by.y = "bcr_patient_barcode", all = T)
paired_samples <- subset(paired_samples, select = c("normal", "lumA", "ajcc_pathologic_stage"))
colnames(paired_samples)[3] <- "lumB"

paired_samples <- merge(paired_samples, Basal, by.x = "normal", by.y = "bcr_patient_barcode", all = T)
paired_samples <- subset(paired_samples, select = c("normal", "lumA", "lumB", "ajcc_pathologic_stage"))
colnames(paired_samples)[4] <- "basal"

paired_samples <- merge(paired_samples, Her2, by.x = "normal", by.y = "bcr_patient_barcode", all = T)
paired_samples <- subset(paired_samples, select = c("normal", "lumA", "lumB", "basal", "ajcc_pathologic_stage"))
colnames(paired_samples)[5] <- "Her2"

unparied <- paired_samples[is.na(paired_samples$lumA) & is.na(paired_samples$lumB) & is.na(paired_samples$basal) & is.na(paired_samples$Her2), ]
paired_samples <- paired_samples[!paired_samples$normal %in% unparied$normal, ]




library(edgeR)

load("RData/TCGA_query.RData")

table(is.na(common$Subtype_Selected))

master_query <- GDCquery(project = "TCGA-BRCA",
                         access = "open", 
                         data.category = "Transcriptome Profiling",
                         experimental.strategy = "RNA-Seq")
GDCdownload(master_query)

master_data <- GDCprepare(master_query, summarizedExperiment = T)

master_unstranded <- assay(master_data, "unstranded")  
rownames(master_unstranded) <- gsub("\\.\\d+", "", rownames(master_unstranded))
master_unstranded <- as.data.frame(master_unstranded)

subtype_subset <- master[master$cases %in% colnames(master_unstranded), ]
unstranded_subset <- master_unstranded[colnames(master_unstranded) %in% subtype_subset$cases]

subtype_subset <- subtype_subset[match(colnames(unstranded_subset), subtype_subset$cases), ]

group <- factor(subtype_subset$Subtype_Selected)

master_counts_filt <- filterByExpr(master_unstranded, group = group)
master_counts_filt <- master_unstranded[master_counts_filt, ]




# takes too long
plotMDS(master_counts_filt, top = 50)





# PCA
pca <- prcomp(t(master_counts_filt))
pca_data <- pca$x
pca_var <- pca$sdev^2

pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)


# Merge sample-stage mapping with PCA data
pca_data <- merge(pca_data, subtype_subset$, by.x = "row.names", by.y = "Sample")

# Create a custom color palette for stages
stage_colors <- c("normal" = "blue", "LumA" = "green", "LumB" = "red", "Her2" = "purple", "Basal" = "orange")

# Create the PCA plot with color mapping
ggplot(pca_data, aes(PC1, PC2, color = Stage)) +
  geom_point() +
  geom_text_repel(aes(label = row.names(pca_data)), size = 3) +  # Adjust the label size here
  scale_color_manual(values = stage_colors) + # Use the custom color palette
  theme_bw() +
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %'))



library(VennDiagram)

venn.diagram(
  x = list(paired_hits = paired_dif_exp$gene_id, unpaired_hits = unpaired_dif_exp$gene_id, wrong_paired_dif_exp$gene_id),
  category.names = c("paired genes", "unpaired genes", "paired F up"),
  filename = "paired_vs_unpaired.png",
  disable.logging = TRUE
)



hits <- read.table("intermediate/Basal/gene_list.txt")
unpaired_dif_exp <- hits

hits <- read.table("intermediate/paired/Basal/gene_list.txt")
paired_dif_exp <- hits

venn.diagram(
    x = list(unpaired_hits = unpaired_dif_exp$V1, paired_hits = paired_dif_exp$V1),
  category.names = c("Unpaired genes", "Paired genes"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("dodgerblue", "goldenrod1"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("dodgerblue", "goldenrod1"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1.5,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "paired_vs_unpaired.png",
  disable.logging = TRUE
)


venn.diagram(
  x = list(paired_hits = paired_dif_exp$gene_id, unpaired_hits = unpaired_dif_exp$gene_id, wrong_paired_dif_exp$gene_id),
  category.names = c("Paired genes", "Unpaired genes", "Paired F up"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("dodgerblue", "goldenrod1", "forestgreen"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("dodgerblue", "goldenrod1", "forestgreen"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1.5,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "paired_vs_unpaired.png",
  disable.logging = TRUE
)




final_gene_counts_paired <- final_gene_counts
final_gene_counts_unpaired <- read.csv("intermediate/LumA/final_gene_counts.csv")



venn.diagram(
  x = list(paired_hits = final_gene_counts$external_gene_name, unpaired_hits = final_gene_counts_unpaired$external_gene_name),
  category.names = c("paired", "unpaired"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("dodgerblue", "goldenrod1"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("dodgerblue", "goldenrod1"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1.5,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "paired_vs_unpaired_ranks.png",
  disable.logging = TRUE
)




library(ggplot2)
library(ggrepel)

# PCA

normal_filt <- filterByExpr(normal_unstranded)
normal_filt <- normal_unstranded[normal_filt, ]

cpm_normal <- cpm(normal_filt)


pca <- prcomp(t(cpm_normal), scale. = T)
pca_data <- pca$x
pca_var <- pca$sdev^2

pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)


# Create a custom color palette for stages
# Create the PCA plot with color mapping
ggplot(pca_data, aes(PC1, PC2)) +
  geom_point() +
  geom_text_repel(aes(label = row.names(pca_data)), size = 3) +  # Adjust the label size here
  theme_bw() +
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %'))



centroid <- colMeans(pca$x[,1:2])




GOI_hits <- hits[hits$gene_id %in% c("ENSG00000091831", "ENSG00000082175", "ENSG00000141736"), ]


edges <- as_edgelist(subnet)


rank_master <- merge(lumA_rank, lumB_rank, by = "external_gene_name", all = T)
rank_master <- merge(rank_master, Her2_rank, by = "external_gene_name", all = T)
rank_master <- merge(rank_master, basal_rank, by = "external_gene_name", all = T)
rownames(ranks) <- NULL

write.csv(rank_master, "../../../../Desktop/Chapter data/TCGA_norm_rank_master.csv")



rank_preview <- merge(lumA_rank[1:25, ], lumB_rank[1:25, ], by = "external_gene_name", all = T)
rank_preview <- merge(rank_preview, Her2_rank[1:25, ], by = "external_gene_name", all = T)
rank_preview <- merge(rank_preview, basal_rank[1:25, ], by = "external_gene_name", all = T)

rank_preview$lumA_rank <- as.integer(rank_preview$lumA_rank)
rank_preview$lumB_rank <- as.integer(rank_preview$lumB_rank)
rank_preview$Her2_rank <- as.integer(rank_preview$Her2_rank)
rank_preview$basal_rank <- as.integer(rank_preview$basal_rank)

rownames(rank_preview) <- NULL


library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(LuminalA = lumA_rank$external_gene_name, 
           LuminalB = lumB_rank$external_gene_name, 
           Her2 = Her2_rank$external_gene_name, 
           Basal = basal_rank$external_gene_name),
  category.names = c("Luminal A", "Luminal B", "Her2", "Basal"),
  col = "transparent",
  fill = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
  alpha = 0.5,
  cex = 2,
  main = "Gene Overlap Across Breast Cancer Subtypes",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-30, 30, -30, 30),
  cat.dist = c(0.05, 0.05, 0.05, 0.05),
  filename = "venn_rank.png",
  disable.logging = T
)

pdf("VennDiagram_ColorBlindFriendly.pdf", width = 10, height = 10)



library(plotly)
library(ggplot2)


plotData <- as.data.frame(table(OpenTargets_NCT_filtered$Subtype, useNA = "ifany"))

plotData <- as.data.frame(table(OpenTargets_unique$`Disease Name`))

pie(plotData$Freq, labels = plotData$Var1, main = "")

ggplot(plotData, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  theme_void() +
  labs(title = "Pie Chart with ggplot2")

plot_ly(plotData, labels = ~Var1, values = ~Freq, type = 'pie',
        textinfo = 'label+value') %>%
  layout(showlegend = FALSE)



ranks <- filtered_targets[!is.na(filtered_targets$lumA_rank) | !is.na(filtered_targets$lumB_rank) | !is.na(filtered_targets$Her2_rank) | !is.na(filtered_targets$basal_rank), ]
ranks <- subset(ranks, select = c("input_term", "synonym", "description", 
                                  "lumA_rank", "lumB_rank", "Her2_rank", "basal_rank"))



library(rDGIdb)

DGIdb <- queryDGIdb(filtered_targets$synonym)
results <- byGene(DGIdb)
detailed_results <- resultSummary(DGIdb)
x <- detailedResults(DGIdb)



plotData <- as.data.frame(table(ranks$Subtype, useNA = "ifany"))

plot_ly(plotData, labels = ~Var1, values = ~Freq, type = 'pie',
        textinfo = 'label+value') %>%
  layout(showlegend = FALSE)



#paired_backup <- unique_paired
#GTEx_backup <- unique_GTEx

unique_paired <- paired_backup
unique_GTEx <- GTEx_backup


unique_paired <- unique_paired[!duplicated(unique_paired$external_gene_name), ]
unique_paired[is.na(unique_paired)] <- "*"

table(grepl("\\*", unique_paired$lumA_logFC) & grepl("\\*", unique_paired$lumB_logFC) & grepl("\\*", unique_paired$Her2_logFC) & grepl("\\*", unique_paired$basal_logFC))
unique_paired <- unique_paired[!grepl("\\*", unique_paired$lumA_logFC) | !grepl("\\*", unique_paired$lumB_logFC) | !grepl("\\*", unique_paired$Her2_logFC) | !grepl("\\*", unique_paired$basal_logFC), ]

table(grepl("\\*", unique_paired$lumA_centrality) & grepl("\\*", unique_paired$lumB_centrality) & grepl("\\*", unique_paired$Her2_centrality) & grepl("\\*", unique_paired$basal_centrality))
unique_paired <- unique_paired[!grepl("\\*", unique_paired$lumA_centrality) | !grepl("\\*", unique_paired$lumB_centrality) | !grepl("\\*", unique_paired$Her2_centrality) | !grepl("\\*", unique_paired$basal_centrality), ]

table(grepl("\\*", unique_paired$lumA_rank) & grepl("\\*", unique_paired$lumB_rank) & grepl("\\*", unique_paired$Her2_rank) & grepl("\\*", unique_paired$basal_rank))
unique_paired <- unique_paired[!grepl("\\*", unique_paired$lumA_rank) | !grepl("\\*", unique_paired$lumB_rank) | !grepl("\\*", unique_paired$Her2_rank) | !grepl("\\*", unique_paired$basal_rank), ]


temp1 <- unique_paired[grepl("\\*", unique_paired$lumA_logFC) | grepl("\\*", unique_paired$lumB_logFC) | grepl("\\*", unique_paired$Her2_logFC) | grepl("\\*", unique_paired$basal_logFC), ]
temp1 <- unique_paired[grepl("\\*", unique_paired$lumA_centrality) | grepl("\\*", unique_paired$lumB_centrality) | grepl("\\*", unique_paired$Her2_centrality) | grepl("\\*", unique_paired$basal_centrality), ]



unique_GTEx <- unique_GTEx[!duplicated(unique_GTEx$external_gene_name), ]
unique_GTEx[is.na(unique_GTEx)] <- "*"

table(grepl("\\*", unique_GTEx$lumA_logFC) & grepl("\\*", unique_GTEx$lumB_logFC) & grepl("\\*", unique_GTEx$Her2_logFC) & grepl("\\*", unique_GTEx$basal_logFC))
unique_GTEx <- unique_GTEx[!grepl("\\*", unique_GTEx$lumA_logFC) | !grepl("\\*", unique_GTEx$lumB_logFC) | !grepl("\\*", unique_GTEx$Her2_logFC) | !grepl("\\*", unique_GTEx$basal_logFC), ]

table(grepl("\\*", unique_GTEx$lumA_centrality) & grepl("\\*", unique_GTEx$lumB_centrality) & grepl("\\*", unique_GTEx$Her2_centrality) & grepl("\\*", unique_GTEx$basal_centrality))
unique_GTEx <- unique_GTEx[!grepl("\\*", unique_GTEx$lumA_centrality) | !grepl("\\*", unique_GTEx$lumB_centrality) | !grepl("\\*", unique_GTEx$Her2_centrality) | !grepl("\\*", unique_GTEx$basal_centrality), ]

table(grepl("\\*", unique_GTEx$lumA_rank) & grepl("\\*", unique_GTEx$lumB_rank) & grepl("\\*", unique_GTEx$Her2_rank) & grepl("\\*", unique_GTEx$basal_rank))
unique_GTEx <- unique_GTEx[!grepl("\\*", unique_GTEx$lumA_rank) | !grepl("\\*", unique_GTEx$lumB_rank) | !grepl("\\*", unique_GTEx$Her2_rank) | !grepl("\\*", unique_GTEx$basal_rank), ]


temp2 <- unique_GTEx[grepl("\\*", unique_GTEx$lumA_logFC) | grepl("\\*", unique_GTEx$lumB_logFC) | grepl("\\*", unique_GTEx$Her2_logFC) | grepl("\\*", unique_GTEx$basal_logFC), ]
temp2 <- unique_GTEx[grepl("\\*", unique_GTEx$lumA_centrality) | grepl("\\*", unique_GTEx$lumB_centrality) | grepl("\\*", unique_GTEx$Her2_centrality) | grepl("\\*", unique_GTEx$basal_centrality), ]

temp2 <- temp1






table(is.na(filtered_targets$lumA_centrality) & is.na(filtered_targets$lumB_centrality) & is.na(filtered_targets$Her2_centrality) & is.na(filtered_targets$basal_centrality))


max_length <- max(length(ranks_unpaired$Drug.Name), length(ranks_unpaired$external_gene_name), length(ranks_paired$Drug.Name), length(ranks_paired$external_gene_name))



unpaired_drugs <- ranks_unpaired$Drug.Name
unparied_genes <- ranks_unpaired$external_gene_name
paired_drugs <- ranks_paired$Drug.Name
paired_genes <- ranks_paired$external_gene_name


length(unpaired_drugs) <- max_length
length(unparied_genes) <- max_length
length(paired_drugs) <- max_length
length(paired_genes) <- max_length

common <- data.frame(unpaired_drugs = unpaired_drugs,
                     unparied_genes = unparied_genes,
                     paired_drugs = paired_drugs,
                     paired_genes = paired_genes)


targets <- targets[!duplicated(targets$uniprot_gn_id), ]

targets_drug <- targets[!is.na(targets$druggability), ]

length(unique(targets$external_gene_name))
length(unique(targets_drug$external_gene_name))

targets$external_gene_name[!targets$external_gene_name %in% targets_drug$external_gene_name]


druggability <- subset(filtered_targets, select = c("external_gene_name", "pocket", "druggability", "max_hit", "struct_score"))

breast_cancer_targets <- c("PIK3CA", "ERBB2", "ESR1", "BRCA1", "BRCA2", "AR", 
                           "CDK4", "CDK6", "EGFR", "AKT1", "MTOR", "CHEK2", 
                           "ERBB3", "PTEN", "TP53", "MAOA", "KRAS")
druggability <- druggability[druggability$external_gene_name %in% breast_cancer_targets, ]


PCSF_lumA <- read.csv("intermediate/LumA/filterByExp/GTEx/PCSF_output.csv")
PCSF_lumB <- read.csv("intermediate/LumB/filterByExp/GTEx/PCSF_output.csv")
PCSF_Her2 <- read.csv("intermediate/Her2/filterByExp/GTEx/PCSF_output.csv")
PCSF_basal <- read.csv("intermediate/basal/filterByExp/GTEx/PCSF_output.csv")

PCSF_lumA <- as.data.frame(table(PCSF_lumA$cluster))
colnames(PCSF_lumA) <- c("cluster", "number_of_genes")
PCSF_lumB <- as.data.frame(table(PCSF_lumB$cluster))
colnames(PCSF_lumB) <- c("cluster", "number_of_genes")
PCSF_Her2 <- as.data.frame(table(PCSF_Her2$cluster))
colnames(PCSF_Her2) <- c("cluster", "number_of_genes")
PCSF_basal <- as.data.frame(table(PCSF_basal$cluster))
colnames(PCSF_basal) <- c("cluster", "number_of_genes")


PCSF_lumA <- read.csv("intermediate/paired/LumA/PCSF_output.csv")
PCSF_lumB <- read.csv("intermediate/paired/LumB/PCSF_output.csv")
PCSF_Her2 <- read.csv("intermediate/paired/Her2/PCSF_output.csv")
PCSF_basal <- read.csv("intermediate/paired/basal/PCSF_output.csv")

PCSF_lumA <- as.data.frame(table(PCSF_lumA$cluster))
colnames(PCSF_lumA) <- c("cluster", "number_of_genes")
PCSF_lumB <- as.data.frame(table(PCSF_lumB$cluster))
colnames(PCSF_lumB) <- c("cluster", "number_of_genes")
PCSF_Her2 <- as.data.frame(table(PCSF_Her2$cluster))
colnames(PCSF_Her2) <- c("cluster", "number_of_genes")
PCSF_basal <- as.data.frame(table(PCSF_basal$cluster))
colnames(PCSF_basal) <- c("cluster", "number_of_genes")





missing_genes_convert <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "uniprot_gn_id", "description"), 
                               filters = "ensembl_gene_id", 
                               values = missing_genes$gene_id, 
                               mart = ensembl)
missing_genes_convert <- merge(missing_genes_convert, missing_genes, by.x = "ensembl_gene_id", by.y = "gene_id", all = T)

# number of unrecognised terms
table(is.na(missing_genes_convert$uniprot_gn_id) & is.na(missing_genes_convert$description))

novel_transcripts <- missing_genes_convert[grep("novel transcript", missing_genes_convert$description), ]
novel_proteins <- missing_genes_convert[grep("novel protein", missing_genes_convert$description), ]
pseudogene <- missing_genes_convert[grep("pseudogene", missing_genes_convert$description), ]

temp <- filtered_targets[!is.na(filtered_targets$Her2_rank), ]
unique(temp$external_gene_name)



DE_results <- list(TCGA_lumA = lumA_DE$hits,
                   TCGA_lumB = lumB_DE$hits,
                   TCGA_Her2 = Her2_DE$hits,
                   TCGA_basal = basal_DE$hits,
                   GTEx_lumA = GTEx_lumA_DE$hits,
                   GTEx_lumB = GTEx_lumB_DE$hits,
                   GTEx_Her2 = GTEx_Her2_DE$hits,
                   GTEx_basal = GTEx_basal_DE$hits)

save(DE_results, file = "RData/DE_results_master.RData")


nrow(GTEx_lumA_QC$low_exp_genes)
temp <- merge(GTEx_ENS, LumA_unstranded, by = "row.names")



LumA_subset <- lumA_DE$hits[rownames(lumA_DE$hits) %in% gene_id$ensembl_gene_id, ]
LumB_subset <- lumB_DE$hits[rownames(lumB_DE$hits) %in% gene_id$ensembl_gene_id, ]
Her2_subset <- Her2_DE$hits[rownames(Her2_DE$hits) %in% gene_id$ensembl_gene_id, ]
Basal_subset <- basal_DE$hits[rownames(basal_DE$hits) %in% gene_id$ensembl_gene_id, ]
GTEx_LumA_subset <- GTEx_lumA_DE$hits[rownames(GTEx_lumA_DE$hits) %in% gene_id$ensembl_gene_id, ]
GTEx_LumB_subset <- GTEx_lumB_DE$hits[rownames(GTEx_lumB_DE$hits) %in% gene_id$ensembl_gene_id, ]
GTEx_Her2_subset <- GTEx_Her2_DE$hits[rownames(GTEx_Her2_DE$hits) %in% gene_id$ensembl_gene_id, ]
GTEx_Basal_subset <- GTEx_basal_DE$hits[rownames(GTEx_basal_DE$hits) %in% gene_id$ensembl_gene_id, ]


LumA_subset <- merge(gene_id, LumA_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
LumB_subset <- merge(gene_id, LumB_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
Her2_subset <- merge(gene_id, Her2_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
Basal_subset <- merge(gene_id, Basal_subset, by.x = "ensembl_gene_id", by.y = "gene_id")

GTEx_LumA_subset <- merge(gene_id, GTEx_LumA_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
GTEx_LumB_subset <- merge(gene_id, GTEx_LumB_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
GTEx_Her2_subset <- merge(gene_id, GTEx_Her2_subset, by.x = "ensembl_gene_id", by.y = "gene_id")
GTEx_Basal_subset <- merge(gene_id, GTEx_Basal_subset, by.x = "ensembl_gene_id", by.y = "gene_id")


temp <- lumA_DE$hits
temp <- temp[temp$logFC >= 1 | temp$logFC <= -1, ]




load("RData/DE_results_master.RData")

venn.diagram(
  x = list(DE_results$TCGA_basal$dif_exp$gene_id,
           DE_results$GTEx_basal$dif_exp$gene_id,
           basal_DE$dif_exp$gene_id),
  category.names = c("DE with TCGA Normal", "DE with GTEx Healthy", "Paired TCGA DE"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "DE_comparison/basal.png",
  disable.logging = TRUE
)



library(readxl)
study_data <- read_excel("../../../../Downloads/bsr-2021-2218_supp1/BSR-2021-2218_suppST3.xlsx")
study_data <- study_data[study_data$...1 %in% gene_id$external_gene_name, ]





missing_genes_convert <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "uniprot_gn_id", "description"), 
                               filters = "external_gene_name", 
                               values = missing_genes, 
                               mart = ensembl)
missing_genes_convert <- merge(missing_genes_convert, missing_genes, by.x = "external_gene_name", by.y = "x", all = T)

# number of unrecognised terms
table(is.na(missing_genes_convert$uniprot_gn_id) & is.na(missing_genes_convert$description))

novel_transcripts <- missing_genes_convert[grep("novel transcript", missing_genes_convert$description), ]
novel_proteins <- missing_genes_convert[grep("novel protein", missing_genes_convert$description), ]
pseudogene <- missing_genes_convert[grep("pseudogene", missing_genes_convert$description), ]


library(VennDiagram)
library(gridExtra)

### plot comparing results from the 3 pipelines
LumA_venn <- venn.diagram(
  x = list(DE_results$TCGA_lumA$dif_exp$gene_id,
           DE_results$GTEx_lumA$dif_exp$gene_id,
           DE_results_paired$TCGA_lumA$dif_exp$gene_id),
  main = "Luminal A",
  category.names = c("DE with TCGA Normal", "DE with GTEx Healthy", "Paired TCGA DE"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)

LumB_venn <- venn.diagram(
  x = list(DE_results$TCGA_lumB$dif_exp$gene_id,
           DE_results$GTEx_lumB$dif_exp$gene_id,
           DE_results_paired$TCGA_lumB$dif_exp$gene_id),
  main = "Luminal B",
  category.names = c("DE with TCGA Normal", "DE with GTEx Healthy", "Paired TCGA DE"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)

Her2_venn <- venn.diagram(
  x = list(DE_results$TCGA_Her2$dif_exp$gene_id,
           DE_results$GTEx_Her2$dif_exp$gene_id,
           DE_results_paired$TCGA_Her2$dif_exp$gene_id),
  main = "Her2",
  category.names = c("DE with TCGA Normal", "DE with GTEx Healthy", "Paired TCGA DE"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)

basal_venn <- venn.diagram(
  x = list(DE_results$TCGA_basal$dif_exp$gene_id,
           DE_results$GTEx_basal$dif_exp$gene_id,
           DE_results_paired$TCGA_basal$dif_exp$gene_id),
  main = "Basal",
  category.names = c("DE with TCGA Normal", "DE with GTEx Healthy", "Paired TCGA DE"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)

LumA_grob <- grobTree(LumA_venn)
LumB_grob <- grobTree(LumB_venn)
Her2_grob <- grobTree(Her2_venn)
basal_grob <- grobTree(basal_venn)

grid.arrange(
  grobs = list(
    LumA_grob,
    LumB_grob,
    Her2_grob,
    basal_grob
  ),
  ncol = 2)





LumA_TCGA_norm <- read.csv("intermediate/LumA/filterByExp/final_gene_counts.csv")
LumA_TCGA_norm <- rownames_to_column(LumA_TCGA_norm)
colnames(LumA_TCGA_norm)[1] <- "lumA_rank"
LumA_TCGA_norm <- subset(LumA_TCGA_norm, select = c("lumA_rank", "external_gene_name"))

LumB_TCGA_norm <- read.csv("intermediate/LumB/filterByExp/final_gene_counts.csv")
LumB_TCGA_norm <- rownames_to_column(LumB_TCGA_norm)
colnames(LumB_TCGA_norm)[1] <- "lumB_rank"
LumB_TCGA_norm <- subset(LumB_TCGA_norm, select = c("lumB_rank", "external_gene_name"))

Her2_TCGA_norm <- read.csv("intermediate/Her2/filterByExp/final_gene_counts.csv")
Her2_TCGA_norm <- rownames_to_column(Her2_TCGA_norm)
colnames(Her2_TCGA_norm)[1] <- "Her2_rank"
Her2_TCGA_norm <- subset(Her2_TCGA_norm, select = c("Her2_rank", "external_gene_name"))

basal_TCGA_norm <- read.csv("intermediate/basal/filterByExp/final_gene_counts.csv")
basal_TCGA_norm <- rownames_to_column(basal_TCGA_norm)
colnames(basal_TCGA_norm)[1] <- "basal_rank"
basal_TCGA_norm <- subset(basal_TCGA_norm, select = c("basal_rank", "external_gene_name"))



LumA_GTEx <- read.csv("intermediate/LumA/filterByExp/GTEx/final_gene_counts.csv")
LumA_GTEx <- rownames_to_column(LumA_GTEx)
colnames(LumA_GTEx)[1] <- "lumA_rank"
LumA_GTEx <- subset(LumA_GTEx, select = c("lumA_rank", "external_gene_name"))

LumB_GTEx <- read.csv("intermediate/LumB/filterByExp/GTEx/final_gene_counts.csv")
LumB_GTEx <- rownames_to_column(LumB_GTEx)
colnames(LumB_GTEx)[1] <- "lumB_rank"
LumB_GTEx <- subset(LumB_GTEx, select = c("lumB_rank", "external_gene_name"))

Her2_GTEx <- read.csv("intermediate/Her2/filterByExp/GTEx/final_gene_counts.csv")
Her2_GTEx <- rownames_to_column(Her2_GTEx)
colnames(Her2_GTEx)[1] <- "Her2_rank"
Her2_GTEx <- subset(Her2_GTEx, select = c("Her2_rank", "external_gene_name"))

basal_GTEx <- read.csv("intermediate/basal/filterByExp/GTEx/final_gene_counts.csv")
basal_GTEx <- rownames_to_column(basal_GTEx)
colnames(basal_GTEx)[1] <- "basal_rank"
basal_GTEx <- subset(basal_GTEx, select = c("basal_rank", "external_gene_name"))



LumA_paired <- read.csv("intermediate/paired/LumA/final_gene_counts.csv")
LumA_paired <- rownames_to_column(LumA_paired)
colnames(LumA_paired)[1] <- "lumA_rank"
LumA_paired <- subset(LumA_paired, select = c("lumA_rank", "external_gene_name"))

LumB_paired <- read.csv("intermediate/paired/LumB/final_gene_counts.csv")
LumB_paired <- rownames_to_column(LumB_paired)
colnames(LumB_paired)[1] <- "lumB_rank"
LumB_paired <- subset(LumB_paired, select = c("lumB_rank", "external_gene_name"))

Her2_paired <- read.csv("intermediate/paired/Her2/final_gene_counts.csv")
Her2_paired <- rownames_to_column(Her2_paired)
colnames(Her2_paired)[1] <- "Her2_rank"
Her2_paired <- subset(Her2_paired, select = c("Her2_rank", "external_gene_name"))

basal_paired <- read.csv("intermediate/paired/basal/final_gene_counts.csv")
basal_paired <- rownames_to_column(basal_paired)
colnames(basal_paired)[1] <- "basal_rank"
basal_paired <- subset(basal_paired, select = c("basal_rank", "external_gene_name"))


LumA_venn <- venn.diagram(
  x = list(LumA_TCGA_norm$external_gene_name,
           LumA_GTEx$external_gene_name,
           LumA_paired$external_gene_name),
  main = "Luminal A",
  category.names = c("TCGA Normal pipeline", "GTEx pipeline", "Paired pipeline"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)

LumB_venn <- venn.diagram(
  x = list(LumB_TCGA_norm$external_gene_name,
           LumB_GTEx$external_gene_name,
           LumB_paired$external_gene_name),
  main = "Luminal B",
  category.names = c("TCGA Normal pipeline", "GTEx pipeline", "Paired pipeline"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)

Her2_venn <- venn.diagram(
  x = list(Her2_TCGA_norm$external_gene_name,
           Her2_GTEx$external_gene_name,
           Her2_paired$external_gene_name),
  main = "Her2",
  category.names = c("TCGA Normal pipeline", "GTEx pipeline", "Paired pipeline"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)

basal_venn <- venn.diagram(
  x = list(basal_TCGA_norm$external_gene_name,
           basal_GTEx$external_gene_name,
           basal_paired$external_gene_name),
  main = "Basal",
  category.names = c("TCGA Normal pipeline", "GTEx pipeline", "Paired pipeline"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)

LumA_grob <- grobTree(LumA_venn)
LumB_grob <- grobTree(LumB_venn)
Her2_grob <- grobTree(Her2_venn)
basal_grob <- grobTree(basal_venn)

grid.arrange(
  grobs = list(
    LumA_grob,
    LumB_grob,
    Her2_grob,
    basal_grob
  ),
  ncol = 2)


lumA_common <- merge(LumA_TCGA_norm, LumA_GTEx, by = "external_gene_name", all = T)
lumA_common <- merge(lumA_common, LumA_paired, by = "external_gene_name", all = T)
colnames(lumA_common) <- c("external_gene_name", "TCGA_norm", "GTEx", "paired")
lumA_common[2:4] <- sapply(lumA_common[2:4], as.integer)

lumB_common <- merge(LumB_TCGA_norm, LumB_GTEx, by = "external_gene_name", all = T)
lumB_common <- merge(lumB_common, LumB_paired, by = "external_gene_name", all = T)
colnames(lumB_common) <- c("external_gene_name", "TCGA_norm", "GTEx", "paired")
lumB_common[2:4] <- sapply(lumB_common[2:4], as.integer)

Her2_common <- merge(Her2_TCGA_norm, Her2_GTEx, by = "external_gene_name", all = T)
Her2_common <- merge(Her2_common, Her2_paired, by = "external_gene_name", all = T)
colnames(Her2_common) <- c("external_gene_name", "TCGA_norm", "GTEx", "paired")
Her2_common[2:4] <- sapply(Her2_common[2:4], as.integer)

basal_common <- merge(basal_TCGA_norm, basal_GTEx, by = "external_gene_name", all = T)
basal_common <- merge(basal_common, basal_paired, by = "external_gene_name", all = T)
colnames(basal_common) <- c("external_gene_name", "TCGA_norm", "GTEx", "paired")
basal_common[2:4] <- sapply(basal_common[2:4], as.integer)

write.csv(lumA_common, "../../../../Desktop/lumA_common.csv")
write.csv(lumB_common, "../../../../Desktop/lumB_common.csv")
write.csv(Her2_common, "../../../../Desktop/Her2_common.csv")
write.csv(basal_common, "../../../../Desktop/basal_common.csv")

lumA_common <- merge(lumA_common, OpenTargets_NCT_filtered, by.x = "external_gene_name", by.y = "Target.Approved.Symbol")



TCGA_norm_venn <- venn.diagram(
  x = list(LumA_TCGA_norm$external_gene_name,
           LumB_TCGA_norm$external_gene_name,
           Her2_TCGA_norm$external_gene_name,
           basal_TCGA_norm$external_gene_name),
  main = "TCGA Normal Pipeline",
  category.names = c("luminal A", "Luminal B", "Her2", "Basal"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)


GTEx_venn <- venn.diagram(
  x = list(LumA_GTEx$external_gene_name,
           LumB_GTEx$external_gene_name,
           Her2_GTEx$external_gene_name,
           basal_GTEx$external_gene_name),
  main = "GTEx Pipeline",
  category.names = c("luminal A", "Luminal B", "Her2", "Basal"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)

paired_venn <- venn.diagram(
  x = list(LumA_paired$external_gene_name,
           LumB_paired$external_gene_name,
           Her2_paired$external_gene_name,
           basal_paired$external_gene_name),
  main = "GTEx Pipeline",
  category.names = c("luminal A", "Luminal B", "Her2", "Basal"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = NULL,
  disable.logging = TRUE
)


TCGA_norm_grob <- grobTree(TCGA_norm_venn)
GTEx_grob <- grobTree(GTEx_venn)
paired_grob <- grobTree(paired_venn)

grid.arrange(
  grobs = list(
    TCGA_norm_grob,
    GTEx_grob,
    paired_grob
    ),
  ncol = 2)


temp <- merge(pcsf_master, OpenTargets_NCT_filtered, by.x = "external_gene_name", by.y = "Target.Approved.Symbol")


temp <- filtered_targets[!duplicated(filtered_targets$external_gene_name), ]

# targets expressed in at least one subtype
removed <- temp[is.na(temp$lumA_logFC) & is.na(temp$lumB_logFC) & is.na(temp$Her2_logFC) & is.na(temp$basal_logFC), ]
temp <- temp[!(is.na(temp$lumA_logFC) & is.na(temp$lumB_logFC) & is.na(temp$Her2_logFC) & is.na(temp$basal_logFC)), ]

# targets in at least one subtype network
removed <- temp[is.na(temp$lumA_centrality) & is.na(temp$lumB_centrality) & is.na(temp$Her2_centrality) & is.na(temp$basal_centrality), ]
temp <- temp[!(is.na(temp$lumA_centrality) & is.na(temp$lumB_centrality) & is.na(temp$Her2_centrality) & is.na(temp$basal_centrality)), ]

# ranked at lesat once
removed <- temp[is.na(temp$lumA_rank) & is.na(temp$lumB_rank) & is.na(temp$Her2_rank) & is.na(temp$basal_rank), ]
temp <- temp[!(is.na(temp$lumA_rank) & is.na(temp$lumB_rank) & is.na(temp$Her2_rank) & is.na(temp$basal_rank)), ]

table(temp$druggability >= 0.5)



pcsf_master <- read.csv("intermediate/paired/LumA/PCSF_master_unique.csv", row.names = 1)
pcsf_master <- read.csv("intermediate/LumA/filterByExp/GTEx/PCSF_master_unique.csv", row.names = 1)




targets <- merge(IDs, OpenTargets_NCT_filtered, by.x = "external_gene_name", by.y = "Target.Approved.Symbol")



# add Fpocket druggability scores
af_drugability <- read.csv("../druggability_results/fpocket_druggability.csv")
targets <- merge(targets, af_drugability, by.x = "uniprot_gn_id", by.y = "uniprot_id", all.x = T)

# add PocketMiner scores
pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results_3.0.csv")
targets <- merge(targets, pocketminer_data, by.x = "uniprot_gn_id", by.y = "uniprot_id", all.x = T)

targets <- targets[order(-targets$druggability), ]
targets <- targets[!duplicated(targets$external_gene_name), ]

targets <- subset(targets, select = c("external_gene_name", "druggability", "max_hit", "struct_score"))
table(targets$druggability >= 0.5)




plot_data <- log(as.matrix(normal_unstranded))
plot(plot_data, xlab = "", ylab = "", main = "normal_unstranded")

plot_data <- log(as.matrix(GTEx_ENS))
plot(plot_data, xlab = "", ylab = "", main = "GTEx")

hist(log(as.matrix(GTEx_ENS)))



NIH_targets <- read.table("NIH_BRCA_approved_drugs.txt", sep = "\t")
colnames(NIH_targets)[1] <- "approved_drugs"
NIH_targets$approved_drugs <- toupper(NIH_targets$approved_drugs)

approved_openTargets <- merge(NIH_targets, OpenTargets, by.x = "approved_drugs", by.y = "Drug Name")
approved_openTargets <- approved_openTargets[!duplicated(approved_openTargets$approved_drugs), ]

NIH_targets$approved_drugs[!NIH_targets$approved_drugs %in% approved_openTargets$approved_drugs]


master_common <- merge(Her2_common, approved_openTargets, by.x = "external_gene_name", by.y = "Target Approved Name")


common_genes <- c(lumA_common$external_gene_name, lumB_common$external_gene_name, Her2_common$external_gene_name, basal_common$external_gene_name)
common_genes <- unique(common_genes)
table(unique(approved_openTargets$`Target Approved Symbol`) %in% common_genes)

