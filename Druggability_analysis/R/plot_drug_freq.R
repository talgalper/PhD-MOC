library(ggplot2)

Fpocket_AF <- read.csv("Fpocket/results_2024.05/fpocket_druggability.csv")
Fpocket_AF <- Fpocket_AF[order(-Fpocket_AF$druggability), ] # remove duplicate uniprot IDs, keep highest druggability
Fpocket_AF <- Fpocket_AF[!duplicated(Fpocket_AF$uniprot_id), ]
Fpocket_SM <- read.csv("Fpocket/SWISSMODEL/fpocket_druggability.csv")
Fpocket_AF$method <- rep("AF", nrow(Fpocket_AF))
Fpocket_SM$method <- rep("SM", nrow(Fpocket_SM))
Fpocket_AF <- subset(Fpocket_AF, select = c("uniprot_id", "druggability", "method"))
Fpocket_SM <- subset(Fpocket_SM, select = c("uniprot_id", "druggability", "method"))
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


# plot distribution of druggble vs undruggable p\roteins
temp <- merge(Fpocket_data, PocketMiner_data, by = "uniprot_id", all = T)
colnames(temp)[4] <- "CP_score"
temp$druggability[is.na(temp$druggability)] <- 0
temp$CP_score[is.na(temp$CP_score)] <- 0
temp$highest_score <- pmax(temp$druggability, temp$CP_score)
temp <- temp[order(-temp$highest_score), ]
temp <- temp[!duplicated(temp$uniprot_id), ]
colnames(temp)[3] <- "Fpocket_method"
colnames(temp)[6] <- "PocketMiner_method"
temp <- temp[, -5]

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

temp <- merge(AF_data, SM_data, by = "uniprot_id")



table(duplicated(temp$uniprot_id))




