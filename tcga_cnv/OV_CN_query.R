library(TCGAbiolinks)
library(SummarizedExperiment)
library(RobustRankAggreg)
library(tidyverse)



query <- GDCquery(project = "TCGA-OV",
                  data.category = "Copy Number Variation")

query <- getResults(query)


## ASCAT3 TCGA CNV data
ASCAT3_query <- GDCquery(project = 'TCGA-OV',
                         data.category = 'Copy Number Variation',
                         experimental.strategy = 'Genotyping Array',
                         data.type = 'Gene Level Copy Number',
                         workflow.type = 'ASCAT3',
                         access = 'open')

ASCAT3 <- getResults(ASCAT3_query)

GDCdownload(ASCAT3_query)

ASCAT3_data <- GDCprepare(ASCAT3_query, summarizedExperiment = T)

ASCAT3_data_raw <- assay(ASCAT3_data, "copy_number")
ASCAT3_data_raw <- as.data.frame(ASCAT3_data_raw)


rownames(ASCAT3_data_raw) <- gsub("\\.\\d+", "", rownames(ASCAT3_data_raw))


## ASCAT2 TCGA CNV data
ASCAT2_query <- GDCquery(project = 'TCGA-OV',
                         data.category = 'Copy Number Variation',
                         experimental.strategy = 'Genotyping Array',
                         data.type = 'Gene Level Copy Number',
                         workflow.type = 'ASCAT2',
                         access = 'open')

ASCAT2 <- getResults(ASCAT2_query)

GDCdownload(ASCAT2_query)

ASCAT2_data <- GDCprepare(ASCAT2_query, summarizedExperiment = T)

ASCAT2_data_raw <- assay(ASCAT2_data, "copy_number")
ASCAT2_data_raw <- na.omit(ASCAT2_data_raw)
ASCAT2_data_raw <- as.data.frame(ASCAT2_data_raw)

rownames(ASCAT2_data_raw) <- gsub("\\.\\d+", "", rownames(ASCAT2_data_raw))


## ABSOLUTE LiftOver TCGA CNV data
ABSOLUTELiftOver_query <- GDCquery(project = 'TCGA-OV',
                         data.category = 'Copy Number Variation',
                         experimental.strategy = 'Genotyping Array',
                         data.type = 'Gene Level Copy Number',
                         workflow.type = 'ABSOLUTE LiftOver',
                         access = 'open')

ABSOLUTELiftOver <- getResults(ABSOLUTELiftOver_query)

GDCdownload(ABSOLUTELiftOver_query)

ABSOLUTELiftOver_data <- GDCprepare(ABSOLUTELiftOver_query, summarizedExperiment = T)

ABSOLUTELiftOver_data_raw <- assay(ABSOLUTELiftOver_data, "copy_number")
ABSOLUTELiftOver_data_raw <- na.omit(ABSOLUTELiftOver_data_raw)
ABSOLUTELiftOver_data_raw <- as.data.frame(ABSOLUTELiftOver_data_raw)

rownames(ASCAT2_data_raw) <- gsub("\\.\\d+", "", rownames(ASCAT2_data_raw))




# get genes with all NA values
rows_all_na <- apply(ASCAT3_data_raw, 1, function(row) all(is.na(row)))
rows_all_na <- rownames(ASCAT3_data_raw)[rows_all_na]

# check how many NAs each gene has
table(rowSums(is.na(ASCAT3_data_raw)))


# calculate the mean of each row, excluding the NA values
ASCAT3_mean_CN <- ASCAT3_data_raw[rowSums(is.na(ASCAT3_data_raw)) <= ceiling(ncol(ASCAT3_data_raw)*0.1), ] # remove genes with NAs in more than 10% of the samples
ASCAT3_mean_CN <- rowMeans(ASCAT3_mean_CN, na.rm = T)
ASCAT3_mean_CN <- as.data.frame(ASCAT3_mean_CN)
ASCAT3_mean_CN <- rownames_to_column(ASCAT3_mean_CN, "gene_id")
colnames(ASCAT3_mean_CN)[2] <- "mean_cn"

barplot(ASCAT3_mean_CN$mean_cn)

#selected_columns <- sapply(ASCAT3_data_raw, function(x) quantile(x, 0.75) > 2)
#ASCAT3_CN_scores <- ASCAT3_data_raw[, selected_columns, drop = FALSE]

# subset copy numbers in Q4
ASCAT3_high_mean_CN <- ASCAT3_mean_CN[ASCAT3_mean_CN$mean_cn >= quantile(ASCAT3_mean_CN$mean_cn, 0.75), ]

# 0:1 normalisation to differentiate the data a bit
ASCAT3_high_mean_CN$normalised <- (ASCAT3_high_mean_CN$mean_cn - min(ASCAT3_high_mean_CN$mean_cn)) / (max(ASCAT3_high_mean_CN$mean_cn) - min(ASCAT3_high_mean_CN$mean_cn)) 


ASCAT3_high_mean_CN$gene_id <- gsub("\\.\\d+", "", ASCAT3_high_mean_CN$gene_id)

write.csv(ASCAT3_high_mean_CN, "ASCAT3_mean_CN.csv", row.names = F)




