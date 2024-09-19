library(ggplot2)

Fpocket_AF <- read.csv("Fpocket/results_2024.05/fpocket_druggability.csv")
Fpocket_AF <- Fpocket_AF[order(-Fpocket_AF$druggability), ] # remove duplicate uniprot IDs, keep highest druggability
Fpocket_AF <- Fpocket_AF[!duplicated(Fpocket_AF$uniprot_id), ]
Fpocket_SM <- read.csv("Fpocket/SWISSMODEL/fpocket_druggability.csv")
Fpocket_AF$method <- rep("AF", nrow(Fpocket_AF))
Fpocket_SM$method <- rep("SM", nrow(Fpocket_SM))
Fpocket_SM$struct_score <- Fpocket_SM$struct_score * 100 # convert to percentage
Fpocket_AF <- subset(Fpocket_AF, select = c("uniprot_id", "druggability", "struct_score", "method"))
Fpocket_SM <- subset(Fpocket_SM, select = c("uniprot_id", "druggability", "struct_score", "method"))
Fpocket_data <- rbind(Fpocket_AF, Fpocket_SM)

PocketMiner_AF <- read.csv("PocketMiner/results/pocketminer_results_4.0.csv")
PocketMiner_AF <- PocketMiner_AF[order(-PocketMiner_AF$max_hit), ] # remove duplicate uniprot IDs, keep highest druggability
PocketMiner_AF <- PocketMiner_AF[!duplicated(PocketMiner_AF$uniprot_id), ]
PocketMiner_SM <- read.csv("PocketMiner/results/SWISSMODEL/pocketminer_results.csv")
PocketMiner_AF$method <- rep("AF", nrow(PocketMiner_AF))
PocketMiner_SM$method <- rep("SM", nrow(PocketMiner_SM))
PocketMiner_AF <- PocketMiner_AF[, -1]
PocketMiner_SM <- PocketMiner_SM[, -1]
PocketMiner_data <- rbind(PocketMiner_AF, PocketMiner_SM)


# plot distribution of druggble vs undruggable proteins
temp <- merge(Fpocket_data, PocketMiner_data, by = "uniprot_id", all = T)
colnames(temp)[5] <- "CP_score"
temp$druggability[is.na(temp$druggability)] <- 0
temp$CP_score[is.na(temp$CP_score)] <- 0
temp$highest_score <- pmax(temp$druggability, temp$CP_score)
temp <- temp[order(-temp$highest_score), ]
temp <- temp[!duplicated(temp$uniprot_id), ]
colnames(temp)[4] <- "Fpocket_method"
colnames(temp)[7] <- "PocketMiner_method"
temp <- temp[, -6]

# Create the source column
temp$source <- ifelse(temp$highest_score == temp$druggability, 
                    temp$Fpocket_method, 
                    temp$PocketMiner_method)

# density plot
ggplot(temp, aes(x=highest_score, fill=source)) +
  geom_density(alpha=0.5) +
  labs(title="Density Plot of Druggable Proteins",
       x="Highest Score", y="Density") +
  geom_vline(xintercept = 0.7, color = "green", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
  theme_minimal()

# frequency plot
ggplot(temp, aes(x=highest_score, color=source)) +
  geom_freqpoly(binwidth = 0.05, size = 1) +
  labs(title="Frequency Plot of Druggable Proteins",
       x="Highest Score", y="Frequency") +
  geom_vline(xintercept = 0.7, color = "green", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = seq(min(temp$highest_score), max(temp$highest_score), by = 0.10)) +
  theme_minimal()


library(VennDiagram)

AF_data <- merge(Fpocket_AF, PocketMiner_AF, by = "uniprot_id")
AF_data <- subset(AF_data, select = c("uniprot_id", "druggability", "max_hit"))
AF_data$highest_score <- pmax(AF_data$druggability, AF_data$max_hit)
AF_data$source <- ifelse(AF_data$highest_score == AF_data$druggability, "Fpocket", "PocketMiner")

SM_data <- merge(Fpocket_SM, PocketMiner_SM, by = "uniprot_id")
SM_data <- subset(SM_data, select = c("uniprot_id", "druggability", "max_hit"))
SM_data$highest_score <- pmax(SM_data$druggability, SM_data$max_hit)
SM_data$source <- ifelse(SM_data$highest_score == SM_data$druggability, "Fpocket", "PocketMiner")

temp2 <- merge(AF_data, SM_data, by = "uniprot_id")

D_ND <- temp2$uniprot_id[temp2$highest_score.x >= 0.5 & temp2$highest_score.y < 0.5]
ND_D <- temp2$uniprot_id[temp2$highest_score.x < 0.5 & temp2$highest_score.y >= 0.5]

unchanged <- temp2$uniprot_id[!temp2$uniprot_id %in% c(D_ND, ND_D)]

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

name_converted <- getBM(attributes = c("external_gene_name", "uniprot_gn_id", "description"), 
                           filters = "external_gene_name", 
                           values = FDA_drug_targets, 
                           mart = ensembl)



Fpocket_D_ND <- temp2$uniprot_id[temp2$druggability.x >= 0.5 & temp2$druggability.y < 0.5]
Fpocket_ND_D <- temp2$uniprot_id[temp2$druggability.x < 0.5 & temp2$druggability.y >= 0.5]

Cryptic_D_ND <- temp2$uniprot_id[temp2$max_hit.x >= 0.5 & temp2$max_hit.y < 0.5]
Cryptic_ND_D <- temp2$uniprot_id[temp2$max_hit.x < 0.5 & temp2$max_hit.y >= 0.5]

unchanged <- temp2$uniprot_id[!temp2$uniprot_id %in% c(Fpocket_D_ND, Fpocket_ND_D, Cryptic_D_ND, Cryptic_ND_D)]
unchanged <- temp2$uniprot_id

# higher score given to already druggable structure
table(temp2$highest_score.x > 0.5 & temp2$highest_score.y > 0.5 & temp2$highest_score.y > temp2$highest_score.x)


venn.diagram(
  x = list(D_ND, ND_D, unchanged),
  category.names = c("drug -> undrug", "undrug -> drug", "unchanged"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("dodgerblue", "lightcoral", "mediumseagreen"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("dodgerblue", "lightcoral", "mediumseagreen"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1.5,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "temp.png",
  disable.logging = TRUE
)



