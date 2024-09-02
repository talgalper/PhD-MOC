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
pdb_dir <- args[1] # path to structures dir
pdb_dir <- "structures"

# Get a list of all the PDB files in the directory
pdb_files <- list.files(pdb_dir, pattern = "\\.pdb$", full.names = TRUE, recursive = TRUE) # adjusted for new dir format
pdb_names <- list.files(pdb_dir, full.names = F)

# Create an empty data frame to store the results
result_df <- data.frame(filename = character(), file_id = character(), uniprot_id = character(), atom_residue = logical(), struct_score = numeric())
removed_df <- data.frame(filename = character(),file_id = character(),  uniprot_id = character(), atom_residue = logical(), struct_score = numeric())

library(progress)
pb <- progress_bar$new(
  format = "  Formatting PDB's & structure score [:bar] :percent eta: :eta",
  total = length(pdb_files), clear = FALSE)

# Loop over the pdb files
for (i in seq_along(pdb_files)) {
  
  pdb_file <- pdb_files[i]
  pdb_name <- pdb_names[i]
  
  # Load the pdb file into R using the `read.pdb` function
  pdb <- read.pdb(pdb_file)
  
  # Convert the parsed pdb data into a data frame
  pdb_df <- as.data.frame(pdb$atom)
  pdb_df <- pdb_df[pdb_df$type == "ATOM", ] # remove non standard residues
  pdb_df <- pdb_df[is.na(pdb_df$insert), ]  # remove inserts
  
  # if structure has multiple chains, subset only residues from chain A
  if (length(unique(pdb_df$chain)) > 1) {
    first_chain <- pdb_df$chain[1]
    pdb_df <- pdb_df[pdb_df$chain == first_chain, ]
  }
  
  # Check if all the "b" column values are the same for each value in the "resno" column
  resno_b_df <- subset(pdb_df, select = c("resno", "b"))
  unique_resno <- unique(resno_b_df$resno)
  b_values <- sapply(unique_resno, function(resno) unique(resno_b_df[resno_b_df$resno == resno, "b"]))
  atom_residue <- all(sapply(b_values, function(b) length(b) == 1))
  
  # Calculate the percentage of "b" column values that are greater than or equal to 50
  struct_score <- mean(pdb_df$b)
  
  # isolate uniprot id from file name
  uniprot_id <- sub("*.pdb", "", pdb_name)
  
  # Add the result to the appropriate data frame
  if (struct_score >= 0.5) {
    result_df <- rbind(result_df, data.frame(filename = pdb_name, file_id = pdb_file, uniprot_id = uniprot_id, atom_residue = atom_residue, struct_score = struct_score))
  } else {
    removed_df <- rbind(removed_df, data.frame(filename = pdb_name, file_id = pdb_file, uniprot_id = uniprot_id, atom_residue = atom_residue, struct_score = struct_score))
  }
  
  # Create a formatted PDB file and write it out
  pdb$atom <- pdb_df
  xyz <- as.numeric(c(t(pdb_df[, c("x", "y", "z")])))
  dim(xyz) <- c(1, length(xyz))
  pdb$xyz <- xyz
  pdb$calpha <- pdb_df$elety == "CA"
  write.pdb(pdb, file = paste0("structures_formatted/", uniprot_id, ".pdb"))
  
  rm(pdb, pdb_df, i, atom_residue, uniprot_id, b_values, unique_resno, resno_b_df, 
     struct_score, pdb_file, pdb_name, xyz)
  pb$tick()
}
rm(pb)

# check if there are any FALSE values in the atom_residue column
if (any(result_df$atom_residue == FALSE)) {
  message("Warning!: Some atoms within a single residue have varying confidence scores.
          Advised to check af_struct_score.csv for these residues (look for the FLASE values).
          Likely multiple chains in PDB.")
}

write.csv(result_df, "results/af_struct_score.csv", row.names = F)
write.csv(removed_df, "results/af_low_conf_struct.csv", row.names = F)



# move low struct score files to separate dir
result_df <- read.csv("results/af_struct_score.csv") # read in if file already exists

for (file in removed_df$file_id) {
  file_id <- sub(".*\\/([A-Z0-9]+)", "\\1", file)
  file.rename(file, paste0("low_conf_struct/", file_id))
  
  rm(file, file_id)
}


pdb <- read.pdb("structures/P20592.pdb")
pdb <- pdb$atom








