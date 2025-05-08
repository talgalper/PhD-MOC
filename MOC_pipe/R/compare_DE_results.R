### Compare and analyse DE results ###
library(edgeR)

load("DE/MOC_vs_BEN/DE_results.RData")
MOC_vs_BEN <- DE_results
load("~/OneDrive - RMIT University/PhD/large_git_files/MOC/DE_results_MOC_vs_GTEx.RData")
MOC_vs_GTEx <- DE_results
load("~/OneDrive - RMIT University/PhD/large_git_files/MOC/DE_results_TCGA_vs_GTEx.RData")
TCGA_vs_GTEx <- DE_results
load("~/OneDrive - RMIT University/PhD/large_git_files/MOC/DE_results_TCGA_vs_GTEx(full).RData")
TCGA_vs_GTEx_full <- DE_results
rm(DE_results)

print(summary(decideTests(MOC_vs_BEN$qlf, p = 0.05, adjust = "fdr", lfc = 1)))
print(summary(decideTests(MOC_vs_GTEx$qlf, p = 0.05, adjust = "fdr", lfc = 1)))
print(summary(decideTests(TCGA_vs_GTEx$qlf, p = 0.05, adjust = "fdr", lfc = 1)))
print(summary(decideTests(TCGA_vs_GTEx_full$qlf, p = 0.05, adjust = "fdr", lfc = 1)))


MOC_vs_BEN_difExp <- read.csv("DE/MOC_vs_BEN/MOC_DE_results.csv")
MOC_vs_GTEx_difExp <- read.csv("DE/MOC_vs_GTEx/MOC_vs_GTEx_DE_results.csv")
TCGA_vs_GTEx_difExp <- read.csv("DE/TCGA_vs_GTEx/TCGA_vs_GTEx_DE_results.csv")
TCGA_vs_GTEx_difExp_full <- read.csv("DE/TCGA_vs_GTEx/full/TCGA_vs_GTEx_DE_results.csv")

venn_data <- list(
  MOC_vs_BEN = MOC_vs_BEN_difExp$ensembl_gene_id,
  MOC_vs_GTEx = MOC_vs_GTEx_difExp$ensembl_gene_id,
  TCGA_vs_GTEx = TCGA_vs_GTEx_difExp$ensembl_gene_id)

library(venn)
library(RColorBrewer)
venn(venn_data, 
     ellipse = F, 
     zcolor = brewer.pal(length(venn_data), name = "Dark2"),
     box = FALSE,
     ilabels = "counts",
     sncs = 2,
     ilcs = 2)


common <- merge(MOC_vs_BEN_difExp[,c(1:5)], MOC_vs_GTEx_difExp[,c(1,5)], by = "ensembl_gene_id")
common <- merge(common, TCGA_vs_GTEx_difExp[,c(1,5)], by = "ensembl_gene_id")
colnames(common)[5:7] <- c("MOC_vs_BEN_logFC", "MOC_vs_GTEx_logFC", "TCGA_vs_GTEx_logFC")


# get genes unique to MOC vs BEN
MOC_unique <- setdiff(MOC_vs_BEN_difExp$ensembl_gene_id, union(MOC_vs_GTEx_difExp$ensembl_gene_id, TCGA_vs_GTEx_difExp$ensembl_gene_id))
MOC_unique <- MOC_vs_BEN_difExp[MOC_vs_BEN_difExp$ensembl_gene_id %in% MOC_unique, ]

# get unique to MOC vs BEN & MOC vs GTEx
MOC_unique <- setdiff(union(MOC_vs_BEN_difExp$ensembl_gene_id, MOC_vs_GTEx_difExp$ensembl_gene_id), TCGA_vs_GTEx_difExp$ensembl_gene_id)


MOC_markers <- c("MUC5AC", "MUC2", "CDX2", "KRT20", "HER2", "CLDN18", "WT1", "PAX8", "CA125", "THBS2", "TAGLN")

MOC_vs_BEN_markers <- MOC_vs_BEN_difExp[MOC_vs_BEN_difExp$external_gene_name %in% MOC_markers, ]
MOC_vs_GTEx_markers <- MOC_vs_GTEx_difExp[MOC_vs_GTEx_difExp$external_gene_name %in% MOC_markers, ]
TCGA_vs_GTEx_markers <- TCGA_vs_GTEx_difExp[TCGA_vs_GTEx_difExp$external_gene_name %in% MOC_markers, ]



temp <- merge(MOC_vs_BEN_markers[,c(2,5)], MOC_vs_GTEx_markers[,c(2,5)], by = "external_gene_name", all = T)
colnames(temp)[c(2,3)] <- c("MOCvsBEN_logFC", "MOCvsGTEx_logFC")
temp <- merge(temp, TCGA_vs_GTEx_markers[,c(2,5)], by = "external_gene_name", all = T)
colnames(temp)[4] <- "TCGAvsGTEx_logFC"
temp <- merge(data.frame(external_gene_name = MOC_markers), temp, all = T)


venn_data <- list(
  MOC_vs_BEN = MOC_vs_BEN_markers$ensembl_gene_id,
  MOC_vs_GTEx = MOC_vs_GTEx_markers$ensembl_gene_id,
  TCGA_vs_GTEx = TCGA_vs_GTEx_markers$ensembl_gene_id)

library(venn)
library(RColorBrewer)
venn(venn_data, 
     ellipse = T, 
     zcolor = brewer.pal(length(venn_data), name = "Dark2"),
     box = FALSE,
     ilabels = "counts",
     sncs = 2,
     ilcs = 2)



     
     
     
     
     







