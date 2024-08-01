#!/bin/bash

# Define source and destination directories
source_dir="Desktop/AlphaFold_data/"
dest_dir="Desktop/druggability_pkg/structures/"

# Find all .pdb.gz files in the source directory and move them to the destination directory
find "$source_dir" -name "*.pdb.gz" -exec mv {} "$dest_dir" \;

# Change to the destination directory
cd "$dest_dir"

# Unzip all the .pdb.gz files in the destination directory
gunzip *.pdb.gz

# Optional: Remove the .pdb.gz files after unzipping
# rm *.pdb.gz

echo "Files moved and unzipped successfully."
