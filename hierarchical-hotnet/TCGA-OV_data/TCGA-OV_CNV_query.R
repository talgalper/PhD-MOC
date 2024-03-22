library(TCGAbiolinks)
library(SummarizedExperiment)
library(RobustRankAggreg)
library(dplyr)
library(tidyverse)
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


#### Get TCGA-OV data from TCGA ####
query <- GDCquery(project = "TCGA-OV",
                  data.category = "Copy Number Variation")

query <- getResults(query)



# extract specifically gene level copy number/copy number estimate
ASCAT3_query <- GDCquery(project = 'TCGA-OV',
                         data.category = 'Copy Number Variation',
                         experimental.strategy = 'Genotyping Array',
                         data.type = 'Gene Level Copy Number',
                         workflow.type = 'ASCAT3',
                         access = 'open')

ASCAT3_query_result <- getResults(ASCAT3_query)
GDCdownload(ASCAT3_query)

ASCAT3_query_result <- subset(ASCAT3_query_result, select = c("submitter_id", "cases"))

# get clinical data
clinical_data <- GDCquery_clinic("TCGA-OV", type = "clinical")
clinical_data <- subset(clinical_data, select = c("project", "submitter_id", "figo_stage", "tissue_or_organ_of_origin"))

# get biospecimen data
biospecimen_data <- GDCquery_clinic("TCGA-OV", type = "Biospecimen")
biospecimen_data <- subset(biospecimen_data, select = c("sample_id", "sample_type", "submitter_id", "tissue_type"))
biospecimen_data$merger <- substr(biospecimen_data$submitter_id, 1, 12) # create column to merge by

# merge by submitter ID
TCGA_OV_master <- merge(clinical_data, biospecimen_data, by.x = "submitter_id", by.y = "merger")
TCGA_OV_master <- subset(TCGA_OV_master, select = c("submitter_id", "submitter_id.y", "figo_stage", 
                                                    "sample_type", "tissue_or_organ_of_origin"))

# add barcodes to data
barcode_subset <- subset(ASCAT3_query_result, select = c("cases", "submitter_id"))
barcode_subset$clinical_id <- substr(barcode_subset$cases, 1, 12)
TCGA_OV_master <- merge(TCGA_OV_master, barcode_subset[, c(1, 3)], by.x = "submitter_id", by.y = "clinical_id")

# subset incomplete entries
incomplete_entries <- TCGA_OV_master[!complete.cases(TCGA_OV_master), ]

# subset stage IV barcodes
stage_IV_barcodes <- TCGA_OV_master[TCGA_OV_master$figo_stage == "Stage IV", ]
stage_IV_barcodes <- na.omit(stage_IV_barcodes)
stage_IV_barcodes <- unique(stage_IV_barcodes$cases)

# rerun GDC query but only for stage IV barcodes
stage_IV_query <- GDCquery(project = 'TCGA-OV',
                           data.category = 'Copy Number Variation',
                           experimental.strategy = 'Genotyping Array',
                           data.type = 'Gene Level Copy Number',
                           workflow.type = 'ASCAT3',
                           access = 'open',
                           barcode = stage_IV_barcodes)

# create matrix of CNV scores for selected barcodes
stage_IV_assays <- GDCprepare(stage_IV_query, summarizedExperiment = T)
stage_IV_data <- assay(stage_IV_assays, "copy_number")
stage_IV_data <- as.data.frame(stage_IV_data)

# subset and remove any genes with NA values
NAgenes <- stage_IV_data[!complete.cases(stage_IV_data), ]
stage_IV_data <- stage_IV_data[complete.cases(stage_IV_data), ]

rownames(stage_IV_data) <- sub("\\..*", "", rownames(stage_IV_data))

# convert gene ensembl to gene symbol
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = rownames(stage_IV_data),
                      mart = ensembl)

# get terms that could not be converted
unmatched_terms <- gene_symbols$ensembl_gene_id[gene_symbols$external_gene_name == ""]


hist(apply(ASCAT3_data_subset, 1, mean))
hist(apply(ASCAT3_data_subset, 1, median))
table(apply(ASCAT3_data_subset, 1, median))

hist(apply(ASCAT3_data_subset, 1, sd))


