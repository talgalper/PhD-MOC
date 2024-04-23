library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)

# load in TCGA query data and create master df with subtypes
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


# data check
table(duplicated(clinical_query$bcr_patient_barcode))
table(is.na(clinical_query$ajcc_pathologic_stage))
table(master$Subtype_Selected)
table(duplicated(master$bcr_patient_barcode))

# subset subtypes
normals <- master[master$Subtype_Selected == "BRCA.Normal", ]
lumA <- master[master$Subtype_Selected == "BRCA.LumA", ]
lumB <- master[master$Subtype_Selected == "BRCA.LumB", ]
Basal <- master[master$Subtype_Selected == "BRCA.Basal", ]
Her2 <- master[master$Subtype_Selected == "BRCA.Her2", ]

# get paried samples for each subtype
lumA_paired <- merge(normals, lumA, by = "bcr_patient_barcode", all = T)
lumA_paired <- lumA_paired[complete.cases(lumA_paired), ]
lumA_paired <- subset(lumA_paired, select = c("bcr_patient_barcode", "cases.x", "cases.y"))
colnames(lumA_paired) <- c("bcr_patient_barcode", "normal", "lumA")

lumB_paired <- merge(normals, lumB, by = "bcr_patient_barcode", all = T)
lumB_paired <- lumB_paired[complete.cases(lumB_paired), ]
lumB_paired <- subset(lumB_paired, select = c("bcr_patient_barcode", "cases.x", "cases.y"))
colnames(lumB_paired) <- c("bcr_patient_barcode", "normal", "lumB")

Her2_paired <- merge(normals, Her2, by = "bcr_patient_barcode", all = T)
Her2_paired <- Her2_paired[complete.cases(Her2_paired), ]
Her2_paired <- subset(Her2_paired, select = c("bcr_patient_barcode", "cases.x", "cases.y"))
colnames(Her2_paired) <- c("bcr_patient_barcode", "normal", "Her2")

Basal_paired <- merge(normals, Basal, by = "bcr_patient_barcode", all = T)
Basal_paired <- Basal_paired[complete.cases(Basal_paired), ]
Basal_paired <- subset(Basal_paired, select = c("bcr_patient_barcode", "cases.x", "cases.y"))
colnames(Basal_paired) <- c("bcr_patient_barcode", "normal", "Basal")


save(lumA_paired, lumB_paired, Her2_paired, Basal_paired, master, file = "RData/paired/paired_subtypes.RData")














#lumA <- as.data.frame(subset(normals, select = c("bcr_patient_barcode")))
#colnames(lumA)[1] <- "normal"
#rownames(paired_samples) <- paired_samples$normal
#
#paired_samples <- merge(paired_samples, lumA, by.x = "normal", by.y = "bcr_patient_barcode")
#paired_samples <- subset(paired_samples, select = c("normal", "ajcc_pathologic_stage"))
#colnames(paired_samples)[2] <- "lumA"
#
#paired_samples <- merge(paired_samples, lumB, by.x = "normal", by.y = "bcr_patient_barcode", all = T)
#paired_samples <- subset(paired_samples, select = c("normal", "lumA", "ajcc_pathologic_stage"))
#colnames(paired_samples)[3] <- "lumB"
#
#paired_samples <- merge(paired_samples, Basal, by.x = "normal", by.y = "bcr_patient_barcode", all = T)
#paired_samples <- subset(paired_samples, select = c("normal", "lumA", "lumB", "ajcc_pathologic_stage"))
#colnames(paired_samples)[4] <- "basal"
#
#paired_samples <- merge(paired_samples, Her2, by.x = "normal", by.y = "bcr_patient_barcode", all = T)
#paired_samples <- subset(paired_samples, select = c("normal", "lumA", "lumB", "basal", "ajcc_pathologic_stage"))
#colnames(paired_samples)[5] <- "Her2"
#
## get unpaired samples
#unparied <- paired_samples[is.na(paired_samples$lumA) & is.na(paired_samples$lumB) & is.na(paired_samples$basal) & is.na(paired_samples$Her2), ]
#
## remove unpaired samples from main df
#paired_samples <- paired_samples[!paired_samples$normal %in% unparied$normal, ]
#
## get samples that have >2 pairings + remove them
#multi_paired <- paired_samples[rowSums(!is.na(paired_samples)) == 3, ]
#paired_samples <- paired_samples[!paired_samples$normal %in% multi_paired$normal, ]
#
#
## subset clinical data with paired samples
#master_subet <- master[master$bcr_patient_barcode %in% paired_samples$normal, ]





