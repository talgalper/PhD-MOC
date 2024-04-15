

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
  x = list(paired_hits = paired_dif_exp$gene_id, unpaired_hits = unpaired_dif_exp$gene_id),
  category.names = c("paired genes", "unpaired genes"),
  filename = "paired_vs_unpaired.png",
  disable.logging = TRUE
)


venn.diagram(
  x = list(paired_hits = paired_dif_exp$gene_id, unpaired_hits = unpaired_dif_exp$gene_id),
  category.names = c("paired genes", "unpaired genes"),
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






