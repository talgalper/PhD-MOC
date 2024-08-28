library(tidyverse)
library(bio3d)

# read in metadata from SWISS-MODEL download
meta <- read.table("../../../../Desktop/SWISS-MODEL_Repository/INDEX", sep = "\t", header = T)

# subset structure with highest qmeandisco_global score for each uniprot ID
meta_subset <- meta
meta_subset$uniprot_id <- sub("-\\d+$", "", meta_subset$UniProtKB_ac)
#meta_subset <- meta_subset[order(-meta_subset$qmeandisco_global), ]
#meta_subset <- meta_subset[!duplicated(meta_subset$uniprot_id), ]


# function to split uniprot ID into paired chars
split_into_pairs <- function(input_string) {
  # Split the string into individual characters
  chars <- strsplit(input_string, NULL)[[1]]
  # Combine the characters into groups of two
  pairs <- sapply(seq(1, length(chars), by = 2), function(i) {
    paste0(chars[i], ifelse(i + 1 <= length(chars), chars[i + 1], ""))
  })
  return(pairs)
}

# create location ID to navigate SWISSMODEL directory
library(progress)
pb <- progress_bar$new(
  format = "  Splitting [:bar] :percent eta: :eta",
  total = length(meta_subset$uniprot_id), clear = FALSE)

loc_id <- list()
for (string in meta_subset$uniprot_id) {
  result <- split_into_pairs(string)
  
  loc_id[[string]] <- result
  
  rm(string, result)
  pb$tick()
}


library(progress)
pb <- progress_bar$new(
  format = "  Copying files [:bar] :percent eta: :eta",
  total = length(unique(meta_subset$uniprot_id)), clear = FALSE)
# move selected structures to new directory, with fancy new progress bar
for (i in seq_along(unique(meta_subset$uniprot_id))) {
  uniprot_id <- unique(meta_subset$uniprot_id)[i]
  dir <- file.path("../../../../Desktop/SWISS-MODEL_Repository", 
                   loc_id[[uniprot_id]][1], 
                   loc_id[[uniprot_id]][2], 
                   loc_id[[uniprot_id]][3], 
                   "swissmodel")
  files <- list.files(dir, full.names = TRUE)
  
  # Create a directory for each Uniprot ID in the output directory
  output_dir <- file.path("../../../../Desktop/SWISSMODEL_Fpocket/structures", uniprot_id)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Loop over each file in the directory and move it to the corresponding Uniprot ID folder
  for (file in files) {
    file.copy(file, file.path(output_dir, basename(file)))
  }
  
  rm(uniprot_id, dir, files, i, file, output_dir)
  pb$tick()
}
rm(pb)





## This code selects for structures with highest qmeandisco_global score indicated earlier
library(progress)
pb <- progress_bar$new(
  format = "  Copying files [:bar] :percent eta: :eta",
  total = length(meta_subset$uniprot_id), clear = FALSE)
# loop over all uniprot IDs and identify dir where files are
for (i in seq_along(meta_subset$uniprot_id)) {
  coordinate_id <- meta_subset$coordinate_id[i]
  uniprot_id <- meta_subset$uniprot_id[i]
  dir <- file.path("../../../../Desktop/SWISS-MODEL_Repository", loc_id[[uniprot_id]][1], loc_id[[uniprot_id]][2], loc_id[[uniprot_id]][3], "swissmodel")
  files <- list.files(dir, full.names = T)
  # loop over each PDB file in dir and get the one that matches to the highest qmeandisco_global score
  for (file in files) {
    if (grepl(coordinate_id, basename(file))) {
      file.copy(file, file.path("../../../../Desktop/SWISSMODEL_Fpocket/structures", paste0(uniprot_id, ".pdb"))) # choose output dir here
    }
  }
  rm(coordinate_id, uniprot_id, dir, files, i, file)
  pb$tick()
}
rm(pb)
