library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)

# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-OV')

# building a query
query_TCGA <- GDCquery(project = 'TCGA-OV',
                       data.category = 'Copy Number Variation')
output_query_TCGA <- getResults(query_TCGA)


# build a query to retrieve gene expression data 
CNV_TCGA_disease <- GDCquery(project = 'TCGA-OV',
                       data.category = 'Copy Number Variation',
                       data.type = 'Copy Number Segment',
                       experimental.strategy = 'Genotyping Array',
                       workflow.type = 'DNAcopy',
                       access = 'open',
                       sample.type = 'Primary Tumor')

output_CNV_TCGA_disease <- getResults(CNV_TCGA_disease)

# repeat for 'normal' tissue (collected next to the tumor)
CNV_TCGA_normal <- GDCquery(project = 'TCGA-OV',
                               data.category = 'Copy Number Variation',
                               data.type = 'Copy Number Segment',
                               experimental.strategy = 'Genotyping Array',
                               workflow.type = 'DNAcopy',
                               access = 'open',
                               sample.type = 'Solid Tissue Normal')

output_CNV_TCGA_normal <- getResults(CNV_TCGA_normal)

# download data 
GDCdownload(CNV_TCGA_disease)
GDCdownload(CNV_TCGA_normal)

# prepare data
tcga_ov_data_disease <- GDCprepare(CNV_TCGA_disease, summarizedExperiment = TRUE)
tcga_ov_data_normal <- GDCprepare(CNV_TCGA_normal, summarizedExperiment = TRUE)


# get clinical data
clinical_data <- GDCquery_clinic("TCGA-OV", type = "clinical")
clinical_data <- subset(clinical_data, select = c("project", "submitter_id", "figo_stage", "tissue_or_organ_of_origin"))

biospecimen_data <- GDCquery_clinic("TCGA-OV", type = "Biospecimen")
biospecimen_data <- subset(biospecimen_data, select = c("sample_id", "sample_type", "submitter_id"))

normal_master <- merge(tcga_ov_data_normal, biospecimen_data, by.x = "GDC_Aliquot", by.y = "sample_id")








