library(ggplot2)

Fpocket_AF <- read.csv("Fpocket/results_2024.05/fpocket_druggability.csv")
# remove duplicate uniprot IDs, keep highest druggability
Fpocket_AF <- Fpocket_AF[order(-Fpocket_AF$druggability), ] 
Fpocket_AF <- Fpocket_AF[!duplicated(Fpocket_AF$uniprot_id), ]

Fpocket_SM <- read.csv("Fpocket/SWISSMODEL/fpocket_druggability.csv")

Fpocket_AF$method <- rep("AF", nrow(Fpocket_AF))
Fpocket_SM$method <- rep("SM", nrow(Fpocket_SM))

Fpocket_AF <- subset(Fpocket_AF, select = c("uniprot_id", "druggability", "method"))
Fpocket_SM <- subset(Fpocket_SM, select = c("uniprot_id", "druggability", "method"))


Fpocket_data <- rbind(Fpocket_AF, Fpocket_SM)

# plot distribution of druggability scores
ggplot(Fpocket_data, aes(x=method, y=druggability, fill=method)) +
  geom_boxplot() +
  geom_jitter(width=0.2, alpha=0.5) +
  labs(title="Comparison of Druggability Scores between SWISSMODEL and AlphaFold",
       x="Method", y="Druggability Score") +
  theme_minimal()

ggplot(Fpocket_data, aes(x=druggability, fill=method)) +
  geom_density(alpha=0.5) +
  labs(title="Density Plot of Druggability Scores",
       x="Druggability Score", y="Density") +
  theme_minimal()



# common uniprot IDs between two databases
common <- intersect(Fpocket_AF$uniprot_id, Fpocket_SM$uniprot_id)
Fpocket_data_common <- Fpocket_data[Fpocket_data$uniprot_id %in% common, ]

ggplot(Fpocket_data_common, aes(x=method, y=druggability, fill=method)) +
  geom_boxplot() +
  geom_jitter(width=0.2, alpha=0.5) +
  labs(title="Comparison of Druggability Scores between SWISSMODEL and AlphaFold",
       x="Method", y="Druggability Score") +
  theme_minimal()

wilcox.test(druggability ~ method, data = Fpocket_data_common)



SM_data <- Fpocket_data_common[Fpocket_data_common$method == "SM", ]
AF_data <- Fpocket_data_common[Fpocket_data_common$method == "AF", ]
merged_data <- merge(SM_data, AF_data, by = "uniprot_id", suffixes = c("_SM", "_AF"))
merged_data$difference <- abs(merged_data$druggability_SM - merged_data$druggability_AF)
merged_data <- merged_data[order(-merged_data$difference), ]

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

uniprot_converted <- getBM(attributes = c("uniprot_gn_id", "external_gene_name", "description"), 
                           filters = "uniprot_gn_id", 
                           values = merged_data$uniprot_id, 
                           mart = ensembl)
uniprot_converted$description <- gsub("\\s*\\[.*?\\]", "", uniprot_converted$description)

merged_data <- merge(uniprot_converted, merged_data, by.x = "uniprot_gn_id", by.y = "uniprot_id", all.y = T)
merged_data[merged_data == ""] <- NA
merged_data <- merged_data[order(-merged_data$difference), ]


targets <- c(
  "AKT1",
  "AKT2",
  "AKT3",
  "CDK4",
  "CDK6",
  "CYP19A1",
  "EGFR",
  "ERBB2",
  "ESR1",
  "PARP1",
  "PARP2",
  "PARP3",
  "PIK3CA")

compare_targets <- merged_data[merged_data$external_gene_name %in% targets, ]
compare_targets <- compare_targets[, -c(5,7)]
compare_targets <- compare_targets[order(compare_targets$external_gene_name), ]
rownames(compare_targets) <- NULL
write.csv(compare_targets, "../../../../Desktop/temp.csv")




