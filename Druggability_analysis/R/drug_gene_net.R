library(igraph)
library(tidyverse)
library(data.table)

BioSNAP_edgelist <- fread("data_general/ChG-Miner_miner-chem-gene.tsv", sep = "\t")
drug_mapping <- fread("data_general/drug-mappings.tsv", sep = "\t")
drug_mapping <- drug_mapping[, c(1,2)]
drug_mapping <- drug_mapping[drug_mapping$drugbankId != "null" & drug_mapping$name != "null", ]
drug_mapping <- unique(drug_mapping)

temp <- drug_mapping[duplicated(drug_mapping$drugbankId), ] 
temp2 <- drug_mapping[duplicated(drug_mapping$drugbankId, fromLast = T), ]
temp <- merge(temp, temp2, all = T)
rm(temp2)

BioSNAP_edgelist <- merge(BioSNAP_edgelist, drug_mapping, by.x = "#Drug", by.y = "drugbankId", all.x = T)
BioSNAP_edgelist <- BioSNAP_edgelist[-1, ]


library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "asia")

uniprot_converted <- getBM(attributes = c("uniprot_gn_id", "external_gene_name"), 
                           filters = "uniprot_gn_id", 
                           values = unique(BioSNAP_edgelist$Gene), 
                           mart = ensembl)

unmapped <- uniprot_converted[uniprot_converted$external_gene_name == "", ]
unrecognised <- BioSNAP_edgelist[!BioSNAP_edgelist$Gene %in% uniprot_converted$uniprot_gn_id, ]

BioSNAP_edgelist <- merge(BioSNAP_edgelist, uniprot_converted, by.x = "Gene", by.y = "uniprot_gn_id", all.x = T)
BioSNAP_edgelist$external_gene_name[is.na(BioSNAP_edgelist$external_gene_name)] <- BioSNAP_edgelist$Gene[is.na(BioSNAP_edgelist$external_gene_name)]

BioSNAP_edgelist <- BioSNAP_edgelist[, c(3,4)]
BioSNAP_edgelist$name <- toupper(BioSNAP_edgelist$name)


DrugBank <- read.csv("DrugBank/DrugBank_targets_ENS.csv")
DrugBank <- DrugBank[!is.na(DrugBank$DrugBank_ID), ]
DrugBank$drug <- toupper(DrugBank$drug)

BioSNAP_edgelist_filt <- BioSNAP_edgelist[BioSNAP_edgelist$name %in% DrugBank$drug, ]
not_in_BioSNAP <- DrugBank[!DrugBank$drugBank_target %in% BioSNAP_edgelist$external_gene_name, ]


BioSNAP_BRCA_graph <- graph_from_edgelist(as.matrix(BioSNAP_edgelist_filt))
BioSNAP_BRCA_graph <- as_undirected(BioSNAP_BRCA_graph)

plot.igraph(BioSNAP_BRCA_graph, asp = 0, edge.arrow.size = 0.3, vertex.size = 3, vertex.label.dist = 1)








BioSNAP_edgelist <- read_tsv("data_general/ChG-Miner_miner-chem-gene.tsv")
DrugBank <- read.csv("DrugBank/DrugBank_targets_ENS.csv")
DrugBank <- DrugBank[!is.na(DrugBank$DrugBank_ID), ]

table(unique(DrugBank$DrugBank_ID) %in% unique(BioSNAP_edgelist$`#Drug`))

not_in_BioSNAP <- DrugBank[!DrugBank$DrugBank_ID %in% BioSNAP_edgelist$`#Drug`, ]


BioSNAP_edgelist_filt <- BioSNAP_edgelist[BioSNAP_edgelist$`#Drug` %in% unique(DrugBank$DrugBank_ID), ]

durgbank_dictionary <- subset(DrugBank, select = c("drug", "DrugBank_ID"))
durgbank_dictionary <- unique(durgbank_dictionary)

BioSNAP_edgelist_filt <- merge(BioSNAP_edgelist_filt, durgbank_dictionary, by.x = "#Drug", by.y = "DrugBank_ID")

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

uniprot_converted <- getBM(attributes = c("uniprot_gn_id", "external_gene_name"), 
                           filters = "uniprot_gn_id", 
                           values = BioSNAP_edgelist_filt$Gene, 
                           mart = ensembl)

unmapped <- uniprot_converted[uniprot_converted$external_gene_name == "", ]
unrecognised <- BioSNAP_edgelist_filt[!BioSNAP_edgelist_filt$Gene %in% uniprot_converted$uniprot_gn_id, ]

BioSNAP_edgelist_filt <- merge(BioSNAP_edgelist_filt, uniprot_converted, by.x = "Gene", by.y = "uniprot_gn_id", all.x = T)

BioSNAP_edgelist_filt$external_gene_name <- ifelse(is.na(BioSNAP_edgelist_filt$external_gene_name), BioSNAP_edgelist_filt$Gene, BioSNAP_edgelist_filt$external_gene_name)
BioSNAP_edgelist_filt <- BioSNAP_edgelist_filt[, 3:4]
BioSNAP_edgelist_filt$drug <- toupper(BioSNAP_edgelist_filt$drug)

BioSNAP_graph <- graph_from_edgelist(as.matrix(BioSNAP_edgelist_filt))
BioSNAP_graph <- as.undirected(BioSNAP_graph)


plot.igraph(BioSNAP_graph, asp = 0, edge.arrow.size = 0.3, vertex.size = 3, vertex.label.dist = 1)
