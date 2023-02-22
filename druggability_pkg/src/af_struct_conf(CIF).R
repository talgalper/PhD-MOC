### removes AlphaFold structures with confidence scores below 50% ###

if (!require("bio3d")) {
  install.packages("bio3d")
}

if (!require("RCurl")) {
  install.packages("RCurl")
}

library(bio3d)
library(RCurl)

cif_dir <- "input_structures/"

# Get a list of all the CIF files in the directory
cif_files <- list.files(cif_dir, pattern = "\\.cif$", full.names = TRUE)

# Create an empty data frame to store the results
result_df <- data.frame(file_name = character(), atom_residue = logical(), struct_conf = numeric())
removed_df <- data.frame(file_name = character(), atom_residue = logical(), struct_conf = numeric())

# Loop over the CIF files
for (cif_file in cif_files) {
  # Load the CIF file into R using the `read.cif` function
  cif <- read.cif(cif_file)
  
  # Convert the parsed CIF data into a data frame
  cif_df <- as.data.frame(cif$atom) 
  
  # Check if all the "b" column values are the same for each value in the "resno" column
  resno_b_df <- subset(cif_df, select = c("resno", "b"))
  unique_resno <- unique(resno_b_df$resno)
  b_values <- sapply(unique_resno, function(resno) unique(resno_b_df[resno_b_df$resno == resno, "b"]))
  atom_residue <- all(sapply(b_values, function(b) length(b) == 1))
  
  # Calculate the percentage of "b" column values that are greater than or equal to 50
  struct_conf <- mean(cif_df$b)
  
  # Add the result to the appropriate data frame
  if (struct_conf >= 50) {
    result_df <- rbind(result_df, data.frame(file_name = cif_file, atom_residue = atom_residue, struct_conf = struct_conf))
  } else {
    removed_df <- rbind(removed_df, data.frame(file_name = cif_file, atom_residue = atom_residue, struct_conf = struct_conf))
  }
}

write.csv(result_df, "results/af_struct_conf.csv")
write.csv(removed_df, "results/af_low_conf_struct.csv")



# read in final csv in and check the atoms_residies column for FALSE cells.
# if found, print "Warning!: check file $FALSE_cif_filenames, varying atomic confidence scores per amino acid residue




cif <- read.cif("Desktop/example_cif/AF-A0A0A0MS04-F1-model_v4.cif")

# Convert the parsed CIF data into a data frame
cif_df <- as.data.frame(cif$atom) 



