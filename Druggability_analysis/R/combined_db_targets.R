#### Combine data from all drug-gene interaction databases to get most reliable interactions ###
library(tidyverse)


NCI_targeted <- read.table("NCI/NCI_targetted.txt", sep = "\t")
NCI_targeted$V1 <- toupper(NCI_targeted$V1)
NCI_targeted$V1 <- trimws(NCI_targeted$V1)


# read in data from drug-gene interaciton databases
DGIdb <- read_tsv("DGIdb/web_dwnld_2024.07/interactions.tsv")
DGIdb$drug_claim_name <- toupper(DGIdb$drug_claim_name)

#OpenTargets <- read_tsv("OpenTargets/old/OpenTargets_data.tsv") # uses Python GraphQL query
OpenTargets <- read_tsv("OpenTargets/OpenTargets_web_dwnld_2024.07.tsv") # web dwnld from safari

DrugBank <- read.csv("DrugBank/DrugBank_targets_ENS.csv")
DrugBank$drug <- toupper(DrugBank$drug)
DrugBank$drug <- trimws(DrugBank$drug)

# combine drugbank data with targeted therapy drugs. includes replacement terms for ambigous drug names
DrugBank_matched <- merge(NCI_targeted, DrugBank, by.x = "V1", by.y = "drug", all.x = T, all.y = T)
DrugBank_matched$Targeted <- ifelse(DrugBank_matched$V1 %in% NCI_targeted$V1, TRUE, FALSE)
colnames(DrugBank_matched)[1] <- "drugName"

temp <- DrugBank_matched[is.na(DrugBank_matched$drugBank_target), ] # unmatched terms
targets_annot <- DrugBank_matched[, -5] # remove Targeted column

# merge with OpenTargets data
targets_annot <- merge(targets_annot, OpenTargets, by = "drugName", all.x = T)

# subset drug terms based unique gene targets
#targets_annot <- targets_annot[!duplicated(targets_annot$drugName) | !duplicated(targets_annot$drugBank_target) | is.na(targets_annot$symbol), ]
colnames(targets_annot)[12] <- "OpenTarget_target" # rename column "symbol"

# merge with DGIdb data
targets_annot <- merge(targets_annot, DGIdb, by.x = "drugName", by.y = "drug_name", all.x = T)
colnames(targets_annot)[19] <- "DGIdb_target" # rename column "gene_name"
targets_annot <- targets_annot[, -17] # remove gene_claim_name

# subset drug terms based unique gene targets across all db's
rows_with_na <- apply(targets_annot[, c("drugName", "drugBank_target", "OpenTarget_target", "DGIdb_target")], 1, function(x) any(is.na(x)))
unique_rows <- !duplicated(targets_annot[!rows_with_na, c("drugName", "drugBank_target", "OpenTarget_target", "DGIdb_target")])
targets_annot_unique <- rbind(targets_annot[rows_with_na, ], targets_annot[!rows_with_na, ][unique_rows, ])

# remove non anti-neoplastic and small molecule drugs
targets_annot_unique <- targets_annot_unique[
  (targets_annot_unique$anti_neoplastic != FALSE | is.na(targets_annot_unique$anti_neoplastic)) & 
    (!(targets_annot_unique$type %in% c("Antibody", "Antibody drug conjugate")) | is.na(targets_annot_unique$type)), ]

# susbet with only gene names for easy viewing
targets_annot_summary <- subset(targets_annot_unique, select = c("drugName", "drugBank_target", "OpenTarget_target", "DGIdb_target",
                                                                 "actionType", "interaction_type", "type", "anti_neoplastic", "interaction_score"))


# get drugs where target is identified in all databases
target_all_dbs <- targets_annot_unique[
  !is.na(targets_annot_unique$drugBank_target) &
    !is.na(targets_annot_unique$OpenTarget_target) &
    !is.na(targets_annot_unique$DGIdb_target) &
    targets_annot_unique$drugBank_target == targets_annot_unique$OpenTarget_target &
    targets_annot_unique$OpenTarget_target == targets_annot_unique$DGIdb_target,
]


write.csv(target_all_dbs, "data_general/target_all_dbs.csv", row.names = F)


