### converts .txt output file from fpocket into organised df ###
if (!require("readr")) {
  install.packages("readr")
}

library(readr)

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
files <- list.files("results/scores/")

# Create a data frame to store the results
results <- data.frame(file_id = character(),
                      uniprot_id = character(),
                      ID = character(),
                      pocket = logical(),
                      druggability = logical(),
                      stringsAsFactors = FALSE)

# Loop through each file
for (i in seq_along(files)) {
  file <- files[i]

  if (file.size(paste0("results/scores/", file)) == 0) {
    print(paste0("File ", i, " of ", length(files), " has no pockets", ": ", file))
    next
  }
  
  print(paste0("Formatting file ", i, " of ", length(files), ": ", file))
  
  split_name <- strsplit(file, split = "-")
  split_name <- unlist(split_name)
  id <- paste0(split_name[2], "-", split_name[3])
  
  # Extract the UniProt ID from the filename
  uniprot_id <- sub("^[^-]+-([^-]+)-.*", "\\1", file)
  
  # file id identifier
  file_id <- gsub("_info.txt", "", file)
  
  # Read in the data using the function
  data <- fpocket_format(paste0("results/scores/", file))
  
  # Find the highest value in the first row
  pocket_max <- max(data[1,])
  
  # Find the highest value in the second row
  druggability_max <- max(data[2,])
  
  # Count the number of values above 0.4 in row 2
  num_drug_pockets <- sum(data[2,] >= 0.4)
  
  # Add the results to the data frame
  results <- rbind(results, data.frame(file_id = file_id,
                                       uniprot_id = uniprot_id,
                                       ID = id,
                                       pocket = pocket_max,
                                       druggability = druggability_max,
                                       num_drug_pockets = num_drug_pockets))
}

# read in results from af_struct_conf.R
af_struct_conf <- read.csv("results/af_struct_score.csv")

# add the structure confidence score to the final df
results_master <- merge(af_struct_conf, results, by.x = "uniprot_id",  by.y = "file_id")
results_master <- results_master[, -c(2,3,4,6,7)]

# Order the table by highest druggability scores
results_master <- results_master[order(-results_master$druggability), ]
rownames(results_master) <- NULL

write_csv(results_master, "results/fpocket_druggability.csv")




