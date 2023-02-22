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
  colnames(new_df) <- new_df[1,]
  new_df <- new_df[-1,]
  
  # rename rows
  new_df <- cbind(new_df, do.call("rbind", strsplit(as.character(new_df$`Pocket 1 :`), ":")))
  new_df <- new_df[, -ncol(new_df)]
  rownames(new_df) <- new_df[, ncol(new_df)]
  new_df <- new_df[, -ncol(new_df)]
  
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
                      both_greater_than_0.4 = logical(),
                      stringsAsFactors = FALSE)

# Loop through each file
for (file in files) {
  # Extract the UniProt ID from the filename
  uniprot_id <- gsub("^.*_(\\w{6})\\.txt$", "\\1", file)
  
  # Read in the data using the function
  data <- fpocket_format(paste0("results/scores/", file))
  
  # Find the highest value in the first row
  pocket_max <- max(data[1,])
  
  # Find the highest value in the second row
  druggability_max <- max(data[2,])
  
  # Check if both column 2 and 3 are greater than or equal to 0.4
  both_greater_than_0.4 <- pocket_max >= 0.4 & druggability_max >= 0.4
  
  # Count the number of values above 0.4 in row 2
  num_drug_pockets <- sum(data[2,] >= 0.4)
  
  # Add the results to the data frame
  results <- rbind(results, data.frame(filename = file,
                                       uniprot_id = uniprot_id,
                                       pocket = pocket_max,
                                       druggability = druggability_max,
                                       both_greater_than_0.4 = both_greater_than_0.4,
                                       num_drug_pockets = num_drug_pockets))
}


write_csv(results, "results/fpocket_druggability.csv")

