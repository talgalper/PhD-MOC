library(ggplot2)

PocketMiner_AF <- read.csv("PocketMiner/results/pocketminer_results_4.0.csv")
# remove duplicate uniprot IDs, keep highest druggability
PocketMiner_AF <- PocketMiner_AF[order(-PocketMiner_AF$max_hit), ] 
PocketMiner_AF <- PocketMiner_AF[!duplicated(PocketMiner_AF$uniprot_id), ]

PocketMiner_SM <- read.csv("PocketMiner/results/SWISSMODEL/pocketminer_results.csv")

PocketMiner_AF$method <- rep("AF", nrow(PocketMiner_AF))
PocketMiner_SM$method <- rep("SM", nrow(PocketMiner_SM))

PocketMiner_AF <- PocketMiner_AF[, -1]
PocketMiner_SM <- PocketMiner_SM[, -1]


PocketMiner_data <- rbind(PocketMiner_AF, PocketMiner_SM)

# plot distribution of druggability scores
ggplot(PocketMiner_data, aes(x=method, y=max_hit, fill=method)) +
  geom_boxplot() +
  geom_jitter(width=0.2, alpha=0.5) +
  labs(title="Comparison of Cryptic Pocket Scores between SWISSMODEL and AlphaFold",
       x="Method", y="Druggability Score") +
  theme_minimal()

ggplot(PocketMiner_data, aes(x=druggability, fill=method)) +
  geom_density(alpha=0.5) +
  labs(title="Density Plot of Druggability Scores",
       x="Druggability Score", y="Density") +
  theme_minimal()



# common uniprot IDs between two databases
common <- intersect(PocketMiner_AF$uniprot_id, PocketMiner_SM$uniprot_id)
PocketMiner_data_common <- PocketMiner_data[PocketMiner_data$uniprot_id %in% common, ]
PocketMiner_data_common <- PocketMiner_data_common[, -3]

ggplot(PocketMiner_data_common, aes(x=method, y=max_hit, fill=method)) +
  geom_boxplot() +
  geom_jitter(width=0.2, alpha=0.5) +
  labs(title="Comparison of Cryptic Pocket Scores between SWISSMODEL and AlphaFold",
       x="Method", y="Druggability Score") +
  theme_minimal()

wilcox.test(max_hit ~ method, data = PocketMiner_data_common)



SM_data <- PocketMiner_data_common[PocketMiner_data_common$method == "SM", ]
AF_data <- PocketMiner_data_common[PocketMiner_data_common$method == "AF", ]
merged_data <- merge(SM_data, AF_data, by = "uniprot_id", suffixes = c("_SM", "_AF"))
merged_data$difference <- abs(merged_data$max_hit_SM - merged_data$max_hit_AF)
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
compare_targets <- compare_targets[order(-compare_targets$difference), ]
rownames(compare_targets) <- NULL
write.csv(compare_targets, "../../../../Desktop/temp.csv")




