library(TCGAbiolinks)
library(SummarizedExperiment)
library(RobustRankAggreg)
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

ASCAT3_data <- GDCprepare(ASCAT3_query, summarizedExperiment = T)
ASCAT3_data_raw <- assay(ASCAT3_data)
ASCAT3_data_raw <- as.data.frame(ASCAT3_data_raw)


# subset and remove any genes with NA values
NAgenes <- ASCAT3_data_raw[!complete.cases(ASCAT3_data_raw), ]
ASCAT3_data_subset <- ASCAT3_data_raw[complete.cases(ASCAT3_data_raw), ]

rownames(ASCAT3_data_subset) <- sub("\\..*", "", rownames(ASCAT3_data_subset))

# convert gene ensembl to gene symbol
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = rownames(ASCAT3_data_subset),
                      mart = ensembl)

