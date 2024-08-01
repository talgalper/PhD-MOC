#### Combine and format NCI known BRCA drug targets ####

# read in NCI data
NCI_approved <- read.table("NCI/NCI_approved.txt", sep = "\t")
NCI_approved$V1 <- toupper(NCI_approved$V1)

NCI_targeted <- read.table("NCI/NCI_targetted.txt", sep = "\t")
NCI_targeted$V1 <- toupper(NCI_targeted$V1)


# combine targets
NCI_master <- data.frame(Drugs = unique(c(NCI_approved$V1, NCI_targeted$V1)))
NCI_master$Targeted <- NCI_master$Drugs %in% NCI_targeted$V1
NCI_master$Approved <- NCI_master$Drugs %in% NCI_approved$V1
NCI_master$Common <- NCI_master$Targeted == NCI_master$Approved

write.csv(NCI_master, "NCI/NCI_approved_BRCA_drugs.csv", row.names = F)


