library(tidyverse)
library(gridExtra)
library(reshape2)


### citaiton data analysis ###


citation_scores <- read.csv("intermediate/citation_scores_3.0.csv")

description <- getBM(attributes = c("external_gene_name", "description"), 
                     filters = "external_gene_name", 
                     values = citation_scores$gene_id, 
                     mart = ensembl)

description$description <- gsub("\\s*\\[.*?\\]", "", description$description)

citation_scores <- merge(description, citation_scores, by.x = "external_gene_name", by.y = "gene_id")


# Melt the data frame to long format for easy plotting
citation_scores <- column_to_rownames(citation_scores, "gene_id")
melted_data <- melt(citation_scores)
melted_data$value <- lapply(melted_data$value, function(col) log(col + 1))
melted_data$value <- as.numeric(melted_data$value)

options(scipen = 999)

ggplot(melted_data, aes(x = variable, y = value)) +
  labs(x = "", y = "Citation counts") +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 13, colour = "black"),  
        axis.title.y = element_text(size = 13, colour = "black", vjust = 6),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  ylim(0, 400000)

ggplot(melted_data, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(breaks = round(seq(min(melted_data$value), max(melted_data$value), by = 2), 1)) +
  labs(x = "Value (log + 1)", y = "Density") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 13, colour = "black"),  
        axis.title.y = element_text(size = 13, colour = "black", vjust = 6),
        axis.title.x = element_text(size = 13, color = "black", vjust = -6),
        plot.margin = unit(c(1,1,1,1), "cm"))


# Function to check for unique values between columns
check_unique_values <- function(df) {
  num_cols <- ncol(df)
  unique_pairs <- combn(num_cols, 2)  # Get all combinations of column pairs
  
  unique_values <- list()
  for (i in 1:ncol(unique_pairs)) {
    col1 <- unique_pairs[1, i]
    col2 <- unique_pairs[2, i]
    
    unique_values[[paste(colnames(df)[col1], colnames(df)[col2], sep = "_vs_")]] <- 
      setdiff(df[[col1]], df[[col2]])  # Find unique values between columns
  }
  
  return(unique_values)
}

# Check for unique values between columns in the dataframe
unique_vals_between_cols <- check_unique_values(master)



# Function to find genes unique to a single column
find_unique_genes <- function(df) {
  unique_genes <- list()
  for (col in colnames(df)) {
    other_cols <- setdiff(colnames(df), col)
    unique_genes[[col]] <- rownames(df)[which(rowSums(df[col] == 1) == sum(df[other_cols] == 0))]
  }
  return(unique_genes)
}

# Find genes unique to each column
unique_genes_per_column <- find_unique_genes(master)

for (col in names(unique_genes_per_column)) {
  cat("Genes unique to", col, ":", ifelse(length(unique_genes_per_column[[col]]) > 0, unique_genes_per_column[[col]], "None"), "\n")
}







