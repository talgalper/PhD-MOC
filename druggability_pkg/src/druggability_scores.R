### converts .txt output file from fpocket into organised df ###
if (!require("readr")) {
  install.packages("readr")
}

library(readr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]


fpocket_format <- function(txt_file){
  # read in .txt file
  df <- read.delim(txt_file, sep = ",", header = F)
  
  # Create a list to store the new columns
  new_columns <- list()
  
  # Loop through every 20th row in column 1 and move it to a new column
  for (i in seq(20, nrow(df), 20)) {
    new_col <- df[(i-19):i, 1]
    new_columns[[length(new_columns)+1]] <- new_col
  }
  
  # Bind the new columns together
  new_df <- do.call(cbind, new_columns)
  new_df <- as.data.frame(new_df)
  
  # rename columns
  new_df <- new_df[-1,]
  new_df <- as.data.frame(new_df) # make sure its a data frame
  colnames(new_df) <- paste0("pocket_", 1:ncol(new_df))
  
  # rename rows
  new_df <- cbind(new_df, do.call("rbind", strsplit(as.character(new_df$`pocket_1`), ":")))
  new_df <- new_df[, -ncol(new_df)]
  row_names <- new_df[, ncol(new_df)]
  new_df <- new_df[, -ncol(new_df)]
  new_df <- as.data.frame(new_df) # make sure its a data frame
  colnames(new_df) <- paste0("pocket_", 1:ncol(new_df))
  rownames(new_df) <- row_names
  
  # loop through every cell in the data frame and remove text before ":"
  for (i in 1:nrow(new_df)) {
    for (j in 1:ncol(new_df)) {
      new_df[i, j] <- gsub(".*:", "", new_df[i, j])
      new_df[i, j] <- gsub("\t", "", new_df[i, j])
    }
  }
  
  # convert all cells from character to numeric
  for (col in names(new_df)) {
    new_df[[col]] <- as.numeric(new_df[[col]])
  }
  
  return(new_df)
}


# List all files in directory "scores"
files <- list.files(input_dir)

# Create a data frame to store the results
results <- data.frame(filename = character(),
                      pocket = logical(),
                      druggability = logical(),
                      stringsAsFactors = FALSE)

# Loop through each file
for (i in seq_along(files)) {
  file <- files[i]
  # Extract the UniProt ID from the filename
  uniprot_id <- sub("^[^-]+-([^-]+)-.*", "\\1", file)
  
  print(paste0("Formatting file ", i, " of ", length(files), ": ", file))
  
  # Read in the data using the function
  data <- fpocket_format(paste0("results/scores/", file))
  
  # Find the highest value in the first row
  pocket_max <- max(data[1,])
  
  # Find the highest value in the second row
  druggability_max <- max(data[2,])
  
  # Count the number of values above 0.4 in row 2
  num_drug_pockets <- sum(data[2,] >= 0.4)
  
  # Add the results to the data frame
  results <- rbind(results, data.frame(filename = file,
                                       uniprot_id = uniprot_id,
                                       pocket = pocket_max,
                                       druggability = druggability_max,
                                       num_drug_pockets = num_drug_pockets))
}


# Order the table by highest druggability scores
results <- results[order(-results$druggability),]

write_csv(results, "results/fpocket_druggability.csv")

