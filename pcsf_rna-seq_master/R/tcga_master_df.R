### if run sequentially this script will provide a single data frame for all your tcga rna-seq data ###
# repeat for each stage

library(tidyr)
library(dplyr)
library(readr)

# set name of target dir
stage <- "stage_IV"


# initialise file names, directories and dataframes
filepath <- file.path("tcga/data", stage)
num_of_files <- list.files(filepath)
num_of_files <- as.character(print(length(num_of_files)))
output_filename <- paste0("tcga/", stage, "_master_df.csv")

col_counter <- 1

files <- list.files(file.path("tcga/data", stage))

temp_df <- read.delim(file.path("tcga/data", stage, files[1]), header = T, skip = 1)
temp_df <- temp_df[-c(1:4), ]

master_df <- data.frame(gene_id = temp_df$gene_id,
                        gene_name = temp_df$gene_name)

## create a master df from all the tcga files
for (i in seq_along(files)) {
  file <- files[i]
  print(paste0("Processing file ", i, " of ", length(files)))
  
  df <- read.delim(file.path("tcga/data", stage, file), header = T, skip = 1)
  # delete first 4 rows 
  df <- df[-c(1:4), ]
  # subset desired columns
  df <- subset(df, select = c("gene_id", "gene_name", "unstranded", "tpm_unstranded", "fpkm_unstranded"))
  
  col_name <- paste0("unstranded_", col_counter)
  master_df[[col_name]] <- df$unstranded
  col_counter <- col_counter + 1
}


write.csv(master_df, output_filename)
