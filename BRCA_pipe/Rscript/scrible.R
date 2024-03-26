


stage_freq_df <- table(common$ajcc_pathologic_stage)
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


