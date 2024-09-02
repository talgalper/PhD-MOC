
# plot distribution of druggble vs undruggable proteins

temp <- merge(Fpocket_data, PocketMiner_data, by = "uniprot_id", all = T)
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
