library(edgeR)
library(tidyverse)

tumour_ec <- read.table("RNAseqDB/data/brca-rsem-count-tcga-t.txt", sep = "\t", header = T)
tumour_ec <- column_to_rownames(tumour_ec, "Hugo_Symbol")
tumour_ec <- tumour_ec[,-1]
normal_ec <- read.table("RNAseqDB/data/brca-rsem-count-tcga.txt", sep = "\t", header = T)
normal_ec <- column_to_rownames(normal_ec, "Hugo_Symbol")
normal_ec <- normal_ec[,-1]
GTEx_ec <- read.table("RNAseqDB/data/breast-rsem-count-gtex.txt", sep = "\t", header = T)
GTEx_ec <- column_to_rownames(GTEx_ec, "Hugo_Symbol")
GTEx_ec <- GTEx_ec[,-1]

data_filt <- filter_low_expr(tumour_ec, GTEx_ec)
DE_results <- DE_analysis(counts_matrix = data_filt$counts_filt, sample_info = data_filt$sample_info)
print_summary(DE_results)

dif_exp <- DE_results$dif_exp


# separate RNAseqDB data based on subtype
load("RData/LumA/DE_data.RData")
load("RData/LumB/DE_data.RData")
load("RData/Her2/DE_data.RData")
load("RData/basal/DE_data.RData")
rm(normal_unstranded)

colnames(tumour_ec) <- gsub("\\.", "-", colnames(tumour_ec))

lumA_data <- tumour_ec[, colnames(tumour_ec) %in% colnames(LumA_unstranded)]
lumB_data <- tumour_ec[, colnames(tumour_ec) %in% colnames(LumB_unstranded)]
Her2_data <- tumour_ec[, colnames(tumour_ec) %in% colnames(Her2_unstranded)]
basal_data <- tumour_ec[, colnames(tumour_ec) %in% colnames(Basal_unstranded)]

subtype_data <- list(lumA_data = lumA_data,
                     lumB_data = lumB_data,
                     Her2_data = Her2_data,
                     basal_data = basal_data)

rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded,
   lumA_data, lumB_data, Her2_data, basal_data)

subtype_filt_data <- list()
subtype_DE_results <- list()
for (i in seq_along(subtype_data)) {
  subtype <- names(subtype_data)[i]
  subtype <- gsub("_data", "", subtype)
  cat("\n\n#### Working on", subtype, "subtype ####\n")
  
  dat.filt <- filter_low_expr(subtype_data[[i]], GTEx_ec)
  DE.results <- DE_analysis(counts_matrix = dat.filt$counts_filt, sample_info = dat.filt$sample_info)
  
  subtype_filt_data[[subtype]] <- dat.filt
  subtype_DE_results[[subtype]] <- DE.results
  
  rm(dat.filt, DE.results, subtype, i)
}

subtype_DE_results$lumA$hits[rownames(subtype_DE_results$lumA$hits) %in% c("ESR1", "PGR", "ERBB2"), ]
subtype_DE_results$lumB$hits[rownames(subtype_DE_results$lumB$hits) %in% c("ESR1", "PGR", "ERBB2"), ]
subtype_DE_results$Her2$hits[rownames(subtype_DE_results$Her2$hits) %in% c("ESR1", "PGR", "ERBB2"), ]
subtype_DE_results$basal$hits[rownames(subtype_DE_results$basal$hits) %in% c("ESR1", "PGR", "ERBB2"), ]

print_summary(subtype_DE_results$basal)

# compare to original data
load("RData/DE_results_master.RData")

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes <- DE_results$GTEx_lumA$dif_exp
gene_id <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                 filters = "ensembl_gene_id", 
                 values = genes$gene_id, 
                 mart = ensembl)

gene_id <- merge(gene_id, genes, all.y = T, by.x = "ensembl_gene_id", by.y = "gene_id")

table(subtype_DE_results$lumA$dif_exp$gene_id %in% gene_id$external_gene_name)

OpenTargets_NCT <- read.csv("OpenTargets_data/OpenTargets_NCT.csv", row.names = 1)
OpenTargets_NCT_filtered <- OpenTargets_NCT[OpenTargets_NCT$Impact.On.Cacner.Progression != "no", ]
OpenTargets_NCT_filtered <- OpenTargets_NCT_filtered[!duplicated(OpenTargets_NCT_filtered$Target.Approved.Symbol), ]

common <- subtype_DE_results$lumA$dif_exp[subtype_DE_results$lumA$dif_exp$gene_id %in% OpenTargets_NCT_filtered$Target.Approved.Symbol, ]
common <- gene_id[gene_id$external_gene_name %in% OpenTargets_NCT_filtered$Target.Approved.Symbol, ]



# number of common starting genes
lumA_filt <- filter_low_expr(LumA_unstranded, GTEx_ENS)

gene_id <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                 filters = "ensembl_gene_id", 
                 values = rownames(lumA_filt$counts_filt), 
                 mart = ensembl)
table(rownames(tumour_ec) %in% gene_id$external_gene_name)



