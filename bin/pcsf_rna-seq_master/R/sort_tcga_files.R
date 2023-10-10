### This script requires "clinical" and "gdc_sample_sheet" downloaded from TCGA cart

library(tidyr)
library(dplyr)

# read in "clinical" file from TCGA and subset desired columns
df <- read.delim("tcga/data/sample/clinical.cart.2023-05-11/clinical.tsv", sep = "\t", header = T)
df <- subset(df, select = c("case_id", "case_submitter_id", "figo_stage"))
stage_data <- df[!duplicated(df), ] # remove duplicate rows
stage_data <- stage_data[!grepl("'--", stage_data$figo_stage),] # remove all rows containing '--

# read in "gdc_sample_sheet" from TCGA and subset desired columns
file_id <- read.delim("tcga/data/sample/gdc_sample_sheet.2023-05-11.tsv", sep = "\t", header = T)
file_id <- subset(file_id, select = c("File.Name", "Case.ID"))
colnames(file_id) <- c("file_name", "case_submitter_id")
file_id$case_submitter_id <- gsub("(.*),.*", "\\1", file_id$case_submitter_id) # remove duplicate case submitter id

# create df with col = filename and stage
stage_files <- left_join(file_id, stage_data, by = "case_submitter_id")
stage_files <- subset(stage_files, select = c("file_name", "figo_stage"))

# view the number of each stage
table(stage_files$figo_stage)

# split data based on stage
stage_list <- split(stage_files, stage_files$figo_stage)


# Assign each data frame to its own variable
for (i in seq_along(stage_list)) {
  assign(paste0(names(stage_list[i])), stage_list[[i]])
}

# consolidate stage sub categories
stage_I <- `Stage IC`
stage_II <- rbind(`Stage IIA`, `Stage IIB`, `Stage IIC`)
stage_III <- rbind(`Stage IIIA`, `Stage IIIB`, `Stage IIIC`)
stage_IV <- `Stage IV`

stage_list_2 <- list(stage_I, stage_II, stage_III, stage_IV)

# clean up global env
rm(list = c("Stage IC", "Stage IIA", "Stage IIB", "Stage IIC", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV"))

# source and destination directories
src_dir <- "tcga/data/rna_seq"
dest_dirs <- c("tcga/data/stage_I", "tcga/data/stage_II", "tcga/data/stage_III", "tcga/data/stage_IV")

# Loop over each data frame and copy files to corresponding directories
for (i in seq_along(stage_list_2)) {
  files_to_move <- stage_list_2[[i]]$file_name
  dest_dir <- dest_dirs[i]
  file.copy(from = file.path(src_dir, files_to_move), to = dest_dir)
}


# use "tcga_master_df.R" to create single df for stage II & III + I & IV



