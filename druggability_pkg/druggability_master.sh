#!/bin/bash


# example command: ./druggability_master.sh path/to/structures/
#                                   or
#                  sh druggability_master.sh path/to/structures/



### master script for druggability scoring###
set -e

mkdir -p results
mkdir -p structures
mkdir -p low_conf_struct

## fpocket setup
sh src/fpocket_setup.sh

# check to see if there are pdb files in the structures folder
structures="structures/"
if [ -z "$(find $structures -maxdepth 1 -type f -name '*.pdb' -print -quit 2>/dev/null)" ]; then
    echo "structure directory does not contain any .pdb files, please add at least 1 to begin"
    exit 1
fi

## remove structures with less than 50% confidence
# creates list of files with confidence >= 50%
Rscript src/af_struct_conf.R structures

# separates all files with confidence <50% and removes quotation marks
filename_list=$(cut -d ',' -f 1 results/af_low_conf_struct.csv | tail -n +2 | sed 's/"//g')

# loop through the filenames and move the corresponding files to the new directory
for filename in $filename_list; do
  echo "Moving file $filename to low_conf_struct/"
  mv structures/"$filename" low_conf_struct/
done

## start fpocket on all structures in structure directory
echo "Starting fpocket on structures"
sh src/fpocket_run.sh structures/

## format scores
Rscript src/druggability_scores.R results/scores

echo "Results saved as fpocket_druggability.csv in results dir"


## potential additions
# find best paraemter(number) for small molecule binding interface
# can have a green zone 70-100% & 50-70%
# highlight pymol druggability stuctures greater than 40 with green
# search literature for druggability threshold

