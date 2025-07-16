### series of exploratory analysis on MOC RS results ###


MOCvsBEN <- read.csv("results/HHnet_RS_ML_overlap.csv")
MOCvsGTEx <- read.csv("results/MOC_vs_GTEx/HHnet_RS_ML_overlap.csv")


# look at presence of pan-cancer genes in MOC data
pan_cancer_genes <- pan_cancer_genes <- c(
  "UBE2C", "KIF18B", "TROAP", "UHRF1", "KIF20A", "HJURP", "IQGAP3", "CDT1",
  "GTSE1", "SKA3", "TTK", "ASF1B", "NEK2", "CDC25C", "AURKB", "NCAPG",
  "PKMYT1", "CDC45", "TPX2", "BIRC5", "CCNB2", "FANCA", "ORC6", "MYBL2",
  "TOP2A", "POLQ", "BUB1B", "NCAPH", "MELK", "KIF4A", "MKI67", "DLGAP5",
  "CENPA", "MTFR2", "CEP55", "FOXM1", "SPC25", "ASPM", "BUB1", "MCM10",
  "CENPF", "FAM64A", "KIF14", "EZH2", "ERCC6L", "FAM111B", "E2F1", "MND1",
  "CLSPN", "GPC2", "MMP9", "C16orf59", "CDCA8", "KIF2C", "CDC20", "EME1"
)

temp <- MOCvsBEN[MOCvsBEN$external_gene_name %in% pan_cancer_genes, ]
temp <- MOCvsGTEx[MOCvsGTEx$external_gene_name %in% pan_cancer_genes, ]


# get a df with all key data
druggability <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")

temp <- merge(MOCvsBEN[,c(1:4,6,9)], druggability, by = c("external_gene_name", "uniprot_gn_id"), all.x = TRUE)
temp <- temp[order(temp$avg_rank), ]
rownames(temp) <- NULL
temp$CP_score <- round(temp$CP_score, digits = 3)
temp$`FP/CP` <- paste(temp$druggability, temp$CP_score, sep = "/")
temp2 <- temp[temp$counts < 5000, ]

temp <- merge(MOCvsGTEx[,c(1,2,3,4,6,9)], druggability, by = c("external_gene_name", "uniprot_gn_id"), all.x = TRUE)
temp <- temp[order(temp$avg_rank), ]
rownames(temp) <- NULL
temp$CP_score <- round(temp$CP_score, digits = 3)
temp$`FP/CP` <- paste(temp$druggability, temp$CP_score, sep = "/")
temp2 <- temp[temp$counts < 5000, ]

# test for normal distribution
shapiro.test(MOCvsBEN$counts)
shapiro.test(MOCvsBEN$Prediction_Score_rf)

# correlation between citation count and ML score
cor.test(MOCvsBEN$Prediction_Score_rf, MOCvsBEN$counts, method = "spearman")
cor.test(MOCvsGTEx$Prediction_Score_rf, MOCvsGTEx$counts, method = "spearman")


# correlation between all ML and Pubtator
feature_data_scores_appended <- id_annot(ensembl, feature_data_scores_appended[, c(1,105:108)], input_type = "uniprot_gn_id", convert_to = "external_gene_name")
temp <- merge(PubTator_counts, feature_data_scores_appended, by.x = "symbol", by.y = "external_gene_name")
temp <- temp[order(temp$Prediction_Score_rf, decreasing = TRUE), ]
cor.test(temp$Prediction_Score_rf, temp$counts, method = "spearman")




# gene ontology 
library(clusterProfiler)
library(org.Hs.eg.db)
library(progress)

MOCvsBEN <- read.csv("results/HHnet_RS_ML_overlap.csv")
# universe <- id_annot(data = rownames(counts_filt), input_type = "ensembl_gene_id", convert_to = "external_gene_name")
genes_ENS <- id_annot(ensembl, data = MOCvsBEN$external_gene_name, input_type = "external_gene_name", convert_to = "ensembl_gene_id")

GO <- enrichGO(genes_ENS$ensembl_gene_id, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP", universe = rownames(counts_filt))
result <- GO@result




# TCGA stage PCA
load("data/serous-OV/TCGA-OV_unstranded.RData")
load("data/serous-OV/GTEx-OV_unstranded.RData")
sample_info <- read.csv("data/serous-OV/TCGA-OV_clinical.csv")
MOC_raw_counts <- read.csv("data/analysis_set_raw_counts.csv", row.names = 1)

common_genes <- intersect(rownames(MOC_raw_counts), intersect(rownames(TCGA_OV_data_unstranded), rownames(GTEx_data)))

all_expr_data <- merge(TCGA_OV_data_unstranded, GTEx_data, by = "row.names")
all_expr_data <- tibble::column_to_rownames(all_expr_data, "Row.names")

all_expr_data_MOC_subset <- all_expr_data[rownames(all_expr_data) %in% common_genes, ]

sample_info$stage <- ifelse(sample_info$figo_stage %in% c("Stage I", "Stage IA", "Stage IC"), "I",
                            ifelse(sample_info$figo_stage %in% c("Stage II", "Stage IIB", "Stage IIC", "Stage IIA"), "II",
                                   ifelse(sample_info$figo_stage %in% c("Stage III", "Stage IIIA", "Stage IIIc", "Stage IIIC", "Stage IIIB"), "III",
                                          ifelse(sample_info$figo_stage %in% "Stage IV", "IV",
                                                 ifelse(is.na(sample_info$figo_stage), "UNK", NA)
                                                 ))))

temp <- data.frame(sample = c(colnames(TCGA_OV_data_unstranded), colnames(GTEx_data)),
                   Classification = c(rep("TCGA-OV", ncol(TCGA_OV_data_unstranded)), 
                                      rep("GTEx-OV", ncol(GTEx_data)))
                   )

sample_info <- merge(temp, sample_info, by.x = "sample", by.y = "cases", all = TRUE)
sample_info$stage <- ifelse(is.na(sample_info$stage), "GTEx", sample_info$stage)
rm(temp)


`TCGA-OV_PCA` <- plot_PCA(expr_data = all_expr_data, 
                          sample_info = sample_info, 
                          output_plot_data = TRUE, 
                          colour = "Classification", 
                          shape = "stage")
PCA_plot <- `TCGA-OV_PCA_MOC_subset`$PCA_plot


`TCGA-OV_PCA_MOC_subset` <- plot_PCA(expr_data = all_expr_data_MOC_subset,
                                     sample_info = sample_info,
                                     output_plot_data = TRUE,
                                     colour = "Classification",
                                     shape = "stage")
PCA_plot <- `TCGA-OV_PCA`$PCA_plot


# all data PCA
load("data/serous-OV/TCGA-OV_unstranded.RData")
load("data/serous-OV/GTEx-OV_unstranded.RData")
sample_info <- read.csv("data/serous-OV/TCGA-OV_clinical.csv")
MOC_raw_counts <- read.csv("data/analysis_set_raw_counts.csv", row.names = 1)
MOC_sample_info <- read.csv("data/All survival_CN_Aug18.csv")
colnames(MOC_raw_counts) <- sub("GAMuT_", "", colnames(MOC_raw_counts))
MOC_sample_info <- MOC_sample_info[MOC_sample_info$GAMUT_ID %in% colnames(MOC_raw_counts), ]
MOC_sample_info <- MOC_sample_info[, c(1,2,4,5)]

all_expr_data <- merge(TCGA_OV_data_unstranded, MOC_raw_counts, by = "row.names")
all_expr_data <- tibble::column_to_rownames(all_expr_data, "Row.names")
all_expr_data <- merge(all_expr_data, GTEx_data, by = "row.names")
all_expr_data <- tibble::column_to_rownames(all_expr_data, "Row.names")

sample_info$stage <- ifelse(sample_info$figo_stage %in% c("Stage I", "Stage IA", "Stage IC"), "I",
                            ifelse(sample_info$figo_stage %in% c("Stage II", "Stage IIB", "Stage IIC", "Stage IIA"), "II",
                                   ifelse(sample_info$figo_stage %in% c("Stage III", "Stage IIIA", "Stage IIIc", "Stage IIIC", "Stage IIIB"), "III",
                                          ifelse(sample_info$figo_stage %in% "Stage IV", "IV",
                                                 ifelse(is.na(sample_info$figo_stage), "UNK", NA)
                                          ))))

MOC_sample_info <- rbind(MOC_sample_info, data.frame(GAMUT_ID = setdiff(colnames(MOC_raw_counts), MOC_sample_info$GAMUT_ID),
                                                     Classification = "UNK", 
                                                     Grade = "UNK", 
                                                     Stage = "UNK"))
MOC_sample_info$stage <- ifelse(MOC_sample_info$Stage %in% c("I", "IA", "IC"), "I",
                            ifelse(MOC_sample_info$Stage %in% c("II", "IIB"), "II",
                                   ifelse(MOC_sample_info$Stage %in% c("III", "IIIA", "IIIc", "IIIC"), "III",
                                          ifelse(MOC_sample_info$Stage %in% "IV", "IV",
                                                 ifelse(MOC_sample_info$Classification == "BEN", "BEN", "UNK")))))


sample_info <- data.frame(sample = c(colnames(TCGA_OV_data_unstranded), colnames(MOC_raw_counts), colnames(GTEx_data)),
                          classification = c(rep("TCGA-OV", ncol(TCGA_OV_data_unstranded)), 
                                             rep("MOC", ncol(MOC_raw_counts)), 
                                             rep("GTEx-OV", ncol(GTEx_data))),
                          stage = c(sample_info$stage, 
                                    MOC_sample_info$stage, 
                                    rep("Healthy", ncol(GTEx_data)))
                          )

library(edgeR)
counts_filt <- filterByExpr(all_expr_data, group = sample_info$classification)
counts_filt <- all_expr_data[counts_filt, ]
low_exp_genes <- all_expr_data[!rownames(all_expr_data) %in% rownames(counts_filt), ]

all_OV_data_PCA <- plot_PCA(expr_data = all_expr_data,
                        sample_info = sample_info,
                        output_plot_data = TRUE,
                        colour = "classification",
                        shape = "stage",
                        shape_values = c(16,16,5,2,18,15,13))

print(all_OV_data_PCA$PCA_plot +
  scale_shape_manual(values = c(1,1,5,2,8,18,13))
)

# get structres that failed pre-checks
af_low_conf_struct <- read.csv("../Druggability_analysis/Fpocket/results_2024.05/af_low_conf_struct.csv")
af_low_conf_struct <- id_annot(ensembl,
                               data = af_low_conf_struct, 
                               col_id = 3, 
                               input_type = "uniprot_gn_id",
                               convert_to = "external_gene_name")
temp <- merge(temp, af_low_conf_struct, by = "external_gene_name", all.x = TRUE)


load("../Druggability_analysis/data_general/TTD_master.RData")
temp <- TTD_master[TTD_master$all_target_genes %in% "TP53" & TTD_master$DRUGTYPE %in% "Small molecular drug", ]



# venn diagram between GTEX and BEN DE
venn_data <- list(
  MOCvsBEN = MOCvsBEN$dif_exp$gene_id,
  MOCvsGTEx = MOCvsGTEx$dif_exp$gene_id
)

library(venn)
library(RColorBrewer)
venn(venn_data, 
     ellipse = T, 
     zcolor = brewer.pal(length(venn_data), name = "Dark2"),
     box = FALSE,
     ilabels = "counts",
     sncs = 2,
     ilcs = 2)





