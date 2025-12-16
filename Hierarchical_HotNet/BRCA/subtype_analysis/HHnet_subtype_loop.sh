#!/usr/bin/env bash

# requires existing matrix and beta file in lumA/intermediate/STRING

for subtype in lumA lumB Her2 basal; do

echo "===== Starting analysis for: $subtype ====="

# IMPORTANT: Make sure to double check and update paths and filenames
# Set the paths for your data, intermediate, and results directories
data=$PWD/BRCA/subtype_analysis/$subtype/data
intermediate=$PWD/BRCA/subtype_analysis/$subtype/intermediate
results=$PWD/BRCA/subtype_analysis/$subtype/results

# Create directories
mkdir -p $data
mkdir -p $intermediate
mkdir -p $results

# Set the name of your network and score file
network=STRING
score=${subtype}_logFC

# Copy required STRING files from lumA into the current subtype directory
mkdir -p $intermediate/"$network"
cp $PWD/BRCA/subtype_analysis/lumA/intermediate/"$network"/beta.txt $intermediate/"$network"/
cp $PWD/BRCA/subtype_analysis/lumA/intermediate/"$network"/similarity_matrix.h5 $intermediate/"$network"/

# Filenames for your input files
index_gene_file=gene_index.tsv
edge_list_file=edge_list_index.tsv
score_file=logFC_scores_abs.tsv

num_permutations=100

# This script uses gnu parallel to parallelize some of the scripts. If gnu
# parallel is not already installed on your system, then you can install it from
# https://www.gnu.org/software/parallel/. You can change the num_cores variable
# to specify the number of cores for your system.

##### MAKE SURE THIS IS CORRECT #####
num_cores=48

mkdir -p $intermediate/"$network"
mkdir -p $intermediate/"$network"_"$score"

# Compile Fortran module (if needed)
#cd ../src
#f2py -c fortran_module.f95 -m fortran_module > /dev/null
#cd ..

# Permute scores
echo "Permuting scores for $subtype..."
cp $data/"$score_file" $intermediate/"$network"_"$score"/scores_0.tsv

python src/find_permutation_bins.py \
    -gsf $intermediate/"$network"_"$score"/scores_0.tsv \
    -igf $data/"$index_gene_file" \
    -elf $data/"$edge_list_file" \
    -ms  1000 \
    -o   $intermediate/"$network"_"$score"/score_bins.tsv

parallel -u -j $num_cores --bar \
    python src/permute_scores.py \
        -i  $intermediate/"$network"_"$score"/scores_0.tsv \
        -bf $intermediate/"$network"_"$score"/score_bins.tsv \
        -s  {} \
        -o  $intermediate/"$network"_"$score"/scores_{}.tsv \
    ::: `seq $num_permutations`

# Construct hierarchies
echo "Constructing hierarchies for $subtype..."
parallel -u -j $num_cores --bar \
    python src/construct_hierarchy.py \
        -smf  $intermediate/"$network"/similarity_matrix.h5 \
        -igf  $data/"$index_gene_file" \
        -gsf  $intermediate/"$network"_"$score"/scores_{}.tsv \
        -helf $intermediate/"$network"_"$score"/hierarchy_edge_list_{}.tsv \
        -higf $intermediate/"$network"_"$score"/hierarchy_index_gene_{}.tsv \
    ::: `seq 0 $num_permutations`

# Process hierarchies
echo "Processing hierarchies for $subtype..."
python src/process_hierarchies.py \
    -oelf $intermediate/"$network"_"$score"/hierarchy_edge_list_0.tsv \
    -oigf $intermediate/"$network"_"$score"/hierarchy_index_gene_0.tsv \
    -pelf $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/hierarchy_edge_list_"$i".tsv "; done) \
    -pigf $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/hierarchy_index_gene_"$i".tsv "; done) \
    -lsb  10 \
    -cf   $results/clusters_"$network"_"$score".tsv \
    -pl   $network $score \
    -pf   $results/sizes_"$network"_"$score".pdf \
    -nc   $num_cores

done

# Clean up copied STRING files (keep original files in lumA)
echo "Cleaning up temporary STRING files..."
for subtype in "lumB" "Her2" "basal"; do
    rm -f $PWD/BRCA/subtype_analysis/$subtype/intermediate/STRING/beta.txt
    rm -f $PWD/BRCA/subtype_analysis/$subtype/intermediate/STRING/similarity_matrix.h5
done
