library(TCGAbiolinks)
library(tidyverse)


getProjectSummary('TCGA-BRCA')

# extract all BRCA RNA-seq samples
query_TCGA <- GDCquery(project = "TCGA-BRCA",
                       access = "open", 
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq")

query_output <- getResults(query_TCGA)

# get clinical data
clinical <- GDCquery_clinic(project = "TCGA-BRCA",
                            type = "clinical")

clinical_query <- clinical[complete.cases(clinical$ajcc_pathologic_stage), ]
clinical_query <- merge(query_output, clinical_query, by.x = "cases.submitter_id", by.y = "submitter_id")
clinical_query <- subset(clinical_query, select = c("cases", "cases.submitter_id", "ajcc_pathologic_stage", 
                                                    "tissue_or_organ_of_origin", "sample_type"))

table(clinical_query$ajcc_pathologic_stage)

# get subtype information
subtypes <- PanCancerAtlas_subtypes()

# samples that are mapped with subtype data
common <- merge(clinical_query, subtypes, by.x = "cases", by.y = "pan.samplesID")
common <- subset(common, select = c("cases", "Subtype_Selected", "sample_type", "ajcc_pathologic_stage"))

# download all data (optional)
data_download <- GDCquery(project = "TCGA-BRCA",
                          access = "open", 
                          data.category = "Transcriptome Profiling",
                          experimental.strategy = "RNA-Seq",
                          barcode = common$cases)
GDCdownload(data_download)

# save all objects for later if required
save(query_TCGA, clinical, subtypes, common, file = "RData/TCGA_query.RData")



