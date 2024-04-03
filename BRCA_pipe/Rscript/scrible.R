


stage_freq_df <- table(common$ajcc_pathologic_stage)
stage_freq_df <- as.data.frame(stage_freq_df)

ggplot(stage_freq_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Stage", y = "No. of Samples") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())


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
                  "CDK6"),
  Example_Drugs = c("Tamoxifen, Fulvestrant", 
                    "Tamoxifen, Fulvestrant", 
                    "Trastuzumab, Pertuzumab", 
                    "Palbociclib, Ribociclib, Abemaciclib",
                    "Palbociclib, Ribociclib, Abemaciclib")
)

LumA_genes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), 
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
                  "MTOR"),
  Example_Drugs = c("Tamoxifen, Fulvestrant", 
                    "Tamoxifen, Fulvestrant", 
                    "Trastuzumab, Pertuzumab", 
                    "Palbociclib, Ribociclib, Abemaciclib",
                    "Palbociclib, Ribociclib, Abemaciclib",
                    "Everolimus")
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



