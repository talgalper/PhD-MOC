#!/bin/bash


### runs fpocket over all .pdb files in a directory ###


# Check if directory path is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <path_to_data_dir>"
    exit 1
fi

# Change to the specified directory
cd $1

# Create output directory
output_dir="../results/structural"
if [ ! -d $output_dir ]; then
    mkdir $output_dir
fi

# Loop through all PDB files in the directory
for file in *.pdb; do
    echo "Analysing file $file"
    # Run fpocket on each file
    fpocket -f "$file"
    
    # Move fpocket output directory to output directory
    fpocket_structures="${file%.*}_out"
    if [ -d $fpocket_structures ]; then
        mv $fpocket_structures $output_dir
        
        # Move txt files to fpocket_scores directory
        scores_dir="../results/scores"
        if [ ! -d $scores_dir ]; then
            mkdir $scores_dir
        fi
        txt_files="$output_dir/$fpocket_structures/*.txt"
        if [ -f $txt_files ]; then
            mv $txt_files $scores_dir
        fi
    fi
done

