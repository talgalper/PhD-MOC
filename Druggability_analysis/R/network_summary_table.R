

pocket_scores <- merge(Fpocket_scores, PocketMiner_scores, by = c("ID", "uniprot_id"))
pocket_scores <- pocket_scores[order(-pocket_scores$druggability)]

# TP53, STAT3, MYC, TGFB1, TGFB2 and MET
non_druggable <- data.frame(uniprot_id = c("P04637", "P40763", "P01106", "P01137", "P61812", "P08581"), 
                            gene_id = c("TP53", "STAT3", "MYC", "TGFB1", "TGFB2", "MET"))
non_druggable <- merge(non_druggable, pocket_scores, by = "uniprot_id", all.x = T)
non_druggable <- subset(non_druggable, select = c("gene_id", "druggability", "cryptic_pocket"))



### create master DF with all network metrics for FDA targets ###


all_db_targets <- read.csv("data_general/target_all_dbs.csv")
FDA_drug_targets <- unique(all_db_targets$drugBank_target)


library(igraph)
library(tidyverse)
library(data.table)
library(biomaRt)


# load in and convert DE data to gene symbols
load("../BRCA_pipe/latest_run/RData/STN_filt/dif_exp.RData")

load("../BRCA_pipe/latest_run/RData/STN_filt/PCSF_results.RData")

HHnet_result <- fread("../Hierarchical_HotNet/BRCA/STN_filt/results/df_subnetNeighs.csv")
HHnetsubnet_result <- fread("../Hierarchical_HotNet/BRCA/STN_filt/results/df_subnet.csv")


# load in and format WGCNA data
load("../WGCNA/BRCA/RData/STN_filt/all_kwithin.RData")
load("../WGCNA/BRCA/RData/STN_filt/all_bwnet.RData")

load("../../../../OneDrive - RMIT University/PhD/large_git_files/WGCNA/TCGA_GTEx_filt_norm.RData") # Mac
load("../../../../Desktop/WGCNA_BRCA_large_files/TCGA_GTEx_filt_norm.RData") # Ubuntu

load("../WGCNA/BRCA/RData/STN_filt/venn_data.RData")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

WGCNAsigned_modules <- as.data.frame(bwnet$colors)
WGCNAsigned_modules <- rownames_to_column(WGCNAsigned_modules)
ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
                           filters = "ensembl_gene_id", 
                           values = names(bwnet$colors), 
                           mart = ensembl)
ensembl_converted$description <- gsub("\\s*\\[.*?\\]", "", ensembl_converted$description)


WGCNAsigned_modules <- merge(WGCNAsigned_modules, ensembl_converted, by.x = "rowname", by.y = "ensembl_gene_id", all.x = T)
colnames(WGCNAsigned_modules) <- c("ensembl_id", "module", "external_gene_name", "description", "gene_biotype")

tumour_associated <- as.data.frame(tumour_associated)
tumour_associated <- merge(tumour_associated, WGCNAsigned_modules, by.x = "tumour_associated", by.y = "ensembl_id", all.x = T)

load("../WGCNA/BRCA/RData/STN_filt/venn_data.RData")
common_genes <- Reduce(intersect, list(DE_genes, tumour_associated, top_kwithin, top_gene_membership))

library(ggVennDiagram)
venn_data <- list(PCSF = df$ensembl_gene_id,
                  WGCNA = common_genes,
                  `HHnet Neighs` = HHnet_result$ensembl_gene_id,
                  HHnet = HHnetsubnet_result$ensembl_gene_id)

ggVennDiagram(venn_data,
              set_size = 8,
              label_size = 8) + 
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 17))


# combine all the data for the FDA drug targets
targets <- read.csv("data_general/target_all_dbs.csv")
targets <- unique(targets$drugBank_target)

master <- data.frame(target = targets)
master <- merge(master, WGCNAsigned_modules, by.x = "target", by.y = "external_gene_name", all.x = T)
master <- subset(master, select = c("target", "ensembl_id", "description", "module"))
master <- merge(master, HHnet_result, by.x = "target", by.y = "external_gene_name", all.x = T)
master <- master[, -c(5,6,7,9:13)]
colnames(master)[5] <- "HHnet_degree" 
master <- merge(master, df, by.x = "target", by.y = "external_gene_name", all.x = T)
master <- master[, -c(6:10,12,13)]
colnames(master)[6] <- "PCSF_degree" 
master <- merge(master, dif_exp, by.x = "ensembl_id", by.y = "gene_id", all.x = T)
master <- master[, -c(8:11)]
master <- merge(master, kWithin, by.x = "ensembl_id", by.y = "row.names", all.x = T)
master <- master[, -c(8,10,11)]
master$tumour_associated <- ifelse(master$target %in% tumour_associated$external_gene_name, "yes", "no")

write.csv(master, "~/Desktop/temp.csv")
