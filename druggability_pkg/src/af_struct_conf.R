### removes AlphaFold structures with confidence scores below 50% ###

if (!require("bio3d")) {
  install.packages("bio3d")
}

if (!require("RCurl")) {
  install.packages("RCurl")
}

library(bio3d)
library(RCurl)

args <- commandArgs(trailingOnly = TRUE)
pdb_dir <- args[1] # path to structres dir

# Get a list of all the PDB files in the directory
pdb_files <- list.files(pdb_dir, pattern = "\\.pdb$", full.names = TRUE)

# Create an empty data frame to store the results
result_df <- data.frame(filename = character(), atom_residue = logical(), struct_score = numeric())
removed_df <- data.frame(filename = character(), atom_residue = logical(), struct_score = numeric())
  
# Loop over the pdb files
for (i in seq_along(pdb_files)) {
  
  pdb_file <- pdb_files[i]
  
  print(paste0("Formatting file ", i, " of ", length(pdb_files), ": ", pdb_file))
  
  # Load the pdb file into R using the `read.pdb` function
  pdb <- read.pdb(pdb_file)
  
  # Convert the parsed pdb data into a data frame
  pdb_df <- as.data.frame(pdb$atom) 
  
  # Check if all the "b" column values are the same for each value in the "resno" column
  resno_b_df <- subset(pdb_df, select = c("resno", "b"))
  unique_resno <- unique(resno_b_df$resno)
  b_values <- sapply(unique_resno, function(resno) unique(resno_b_df[resno_b_df$resno == resno, "b"]))
  atom_residue <- all(sapply(b_values, function(b) length(b) == 1))
  
  # Calculate the percentage of "b" column values that are greater than or equal to 50
  struct_score <- mean(pdb_df$b)
  
  # removes file prefix
  pdb_file <- gsub(".*/", "", pdb_file)
  
  # Add the result to the appropriate data frame
  if (struct_score >= 50) {
    result_df <- rbind(result_df, data.frame(filename = pdb_file, atom_residue = atom_residue, struct_score = struct_score))
  } else {
    removed_df <- rbind(removed_df, data.frame(filename = pdb_file, atom_residue = atom_residue, struct_score = struct_score))
  }
}

write.csv(result_df, "results/af_struct_score.csv", row.names = F)
write.csv(removed_df, "results/af_low_conf_struct.csv", row.names = F)


# for shell script: echo 
# "Warning!: check af_struct_score.csv for files with varying atomic confidence scores per amino acid residue"
# if any values in the atom_residue column = FALSE

