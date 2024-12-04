library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

drug_scores <- read.csv("data_general/druggability_scores.csv")


uniprot_gn_id <- getBM(attributes = c( "uniprot_gn_id", "external_gene_name", "description"), 
                       filters = "uniprot_gn_id", 
                       values = drug_scores$uniprot_id, 
                       mart = ensembl)

uniprotswissprot <- getBM(attributes = c( "uniprotswissprot", "external_gene_name", "description"),
                          filters = "uniprotswissprot",
                          values = drug_scores$uniprot_id,
                          mart = ensembl)


IDs_converted <- merge(uniprot_gn_id, uniprotswissprot, by.x = "uniprot_gn_id", by.y = "uniprotswissprot", all = T)

IDs_converted$external_gene_name <- ifelse(IDs_converted$external_gene_name.x == "" | is.na(IDs_converted$external_gene_name.x),
                                           IDs_converted$external_gene_name.y, IDs_converted$external_gene_name.x)
IDs_converted <- subset(IDs_converted, select = c("uniprot_gn_id", "external_gene_name", "description.x"))
colnames(IDs_converted)[3] <- "description"
IDs_converted$description <- gsub("\\s*\\[.*?\\]", "", IDs_converted$description)

unmapped <- IDs_converted[IDs_converted$external_gene_name == "", ]
unmapped <- na.omit(unmapped)
unrecognised <- drug_scores[!drug_scores$uniprot_id %in% IDs_converted$uniprot_gn_id, ]

IDs_converted <- IDs_converted[IDs_converted$external_gene_name != "", ]
IDs_converted <- IDs_converted[!is.na(IDs_converted$external_gene_name), ]

novel_transcripts <- unmapped[grep("novel transcript", unmapped$description), ]
novel_proteins <- unmapped[grep("novel protein", unmapped$description), ]
pseudogene <- unmapped[grep("pseudogene", unmapped$description), ]

drug_scores <- merge(IDs_converted, drug_scores, by.x = "uniprot_gn_id", by.y = "uniprot_id", all.y = T)
drug_scores <- unique(drug_scores)

write.csv(drug_scores, "data_general/druggability_scores_annot.csv", row.names = F)


drug_scores <- read.csv("data_general/druggability_scores_annot.csv")
AF_Fpocket <- read.csv("Fpocket/results_2024.05/fpocket_druggability_full.csv")
SM_Fpocket <- read.csv("Fpocket/SWISSMODEL/fpocket_druggability_full.csv")

AKT <- c(
"P31749",
"P31751",
"Q9Y243"
)

AKT_AF <- AF_Fpocket[AF_Fpocket$uniprot_id %in% AKT, ]
AKT_SM <- SM_Fpocket[SM_Fpocket$uniprot_id %in% AKT, ]
