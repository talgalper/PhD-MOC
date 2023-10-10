library(pheatmap)
library(sva)

sample_info <- read.csv("~/Desktop/final_copy/pcsf_kylie/raw_data/All survival_CN_Aug18.csv")
sample_info <- subset(sample_info, select = c("GAMUT_ID", "Grade", "Stage"))

## Check for batch effects
be_data <- data_log_norm

# Convert GAMUT_ID to character in sample_info to handle non-numeric values
sample_info$GAMUT_ID <- as.character(sample_info$GAMUT_ID)

# Create a mapping between GAMUT_ID and Stage
id_to_stage <- sample_info[, c("GAMUT_ID", "Stage")]

# Initialize a vector to store new column names
new_colnames <- character(length(colnames(be_data)))

# Loop through each column in data_log_norm
for (i in seq_along(colnames(be_data))) {
  col_name <- colnames(be_data)[i]
  # Extract GAMUT_ID from the column name
  gamut_id <- sub("GAMuT_", "", col_name)
  
  # Find the corresponding stage from the mapping
  stage <- id_to_stage$Stage[id_to_stage$GAMUT_ID == gamut_id]
  
  # If a corresponding stage is found, create the new column name
  if (!is.na(stage[1])) {
    new_col_name <- paste(col_name, stage[1], sep = "_")
    new_colnames[i] <- new_col_name
  } else {
    new_colnames[i] <- col_name
  }
}

# Assign the new column names to data_log_norm
colnames(be_data) <- new_colnames


# Reorder the columns of the data frame
colnames(be_data) <- gsub("[AaBbCc]", "", colnames(be_data))
stage_labels <- sub(".*_", "", colnames(be_data))
ordering_index <- order(stage_labels)
be_data <- be_data[, ordering_index]


pheatmap(cor(data_log_norm))






sample_data <- data.frame(str_split_fixed(colnames(be_data), "_", 3))

# Rejoin the first two columns
rejoined_columns <- apply(sample_data[, c("X1", "X2")], 1, function(row) {
  paste(row, collapse = "_")
})

# Convert the result to a data frame
rejoined_df <- data.frame(Rejoined = rejoined_columns)

sample_data <- cbind(rejoined_df, sample_data$X3)
colnames(sample_data) <- c("sample", "stage")

sample_data_pca <- data.frame(sample_data, prcomp(be_data, scale. = T, center = T)$rotation)

ggplot(sample_data_pca, aes(
  x = PC2,
  y = PC3,
  colour = stage)) + geom_point(size = 7)



combat_data <- ComBat(be_data, sample_data_pca$stage, par.prior = T)


combat_data_pca <- data.frame(sample_data, prcomp(combat_data, scale. = T, center = T)$rotation)
ggplot(combat_data_pca, aes(
  x = PC2,
  y = PC3,
  colour = stage)) + geom_point(size = 7)

