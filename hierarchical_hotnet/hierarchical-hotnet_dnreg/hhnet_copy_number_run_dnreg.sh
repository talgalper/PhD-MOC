#!/usr/bin/env bash

data=$PWD/data
intermediate=$PWD/intermediate
results=$PWD/results

num_permutations=100

# This script uses gnu parallel to parallizes some of the scripts.  If gnu
# parallel is not already installed on your system, then you can install it from
# https://www.gnu.org/software/parallel/.  You can change the num_cores variable
# to specify the number of cores for your system.

num_cores=8

# Compile Fortran module.
#cd ../src
#f2py -c fortran_module.f95 -m fortran_module > /dev/null
#cd ..

################################################################################
#
#   Prepare data.
#
################################################################################

# Create data, intermediate data and results, and results directories.
mkdir -p $data
mkdir -p $intermediate
mkdir -p $results

for network in network_1
do
    mkdir -p $intermediate/"$network"
done

for network in network_1
do
    for score in scores_1
    do
        mkdir -p $intermediate/"$network"_"$score"
    done
done

################################################################################
#
#   Construct similarity matrices.
#
################################################################################

echo "Construct similarity matrices..."

for network in network_1
do
    python src/construct_similarity_matrix.py \
        -i   $data/dnreg_edge_list_index.tsv \
        -o   $intermediate/similarity_matrix.h5 \
        -bof $intermediate/beta.txt
done

################################################################################
#
#   Permute data.
#
################################################################################

echo "Permuting scores..."

for network in network_1
do
    for score in scores_1
    do
        cp $data/dnreg_protein_to_score.tsv $intermediate/scores/scores_0.tsv

        python src/find_permutation_bins.py \
            -gsf $intermediate/scores/scores_0.tsv \
            -igf $data/dnreg_index_to_protein.tsv \
            -elf $data/dnreg_edge_list_index.tsv \
            -ms  1000 \
            -o   $intermediate/score_bins.tsv

        parallel -u -j $num_cores --bar \
            python src/permute_scores.py \
                -i  $intermediate/scores/scores_0.tsv \
                -bf $intermediate/score_bins.tsv \
                -s  "$i" \
                -o  $intermediate/scores/scores_"$i".tsv
            ::: `seq $num_permutations`
    done
done

################################################################################
#
#   Construct hierarchies.
#
################################################################################

echo "Constructing hierarchies..."

for network in network_1
do
    for score in scores_1
    do
        parallel -u -j $num_cores --bar \
            python src/construct_hierarchy.py \
                -smf  $intermediate/similarity_matrix.h5 \
                -igf  $data/dnreg_index_to_protein.tsv \
                -gsf  $intermediate/scores/scores_"$i".tsv \
                -helf $intermediate/scores/hierarchy_dnreg_edge_list_index_"$i".tsv \
                -higf $intermediate/scores/hierarchy_dnreg_index_to_protein_"$i".tsv
            ::: `seq 0 $num_permutations`
    done
done

################################################################################
#
#   Process hierarchies.
#
################################################################################

echo "Processing hierarchies..."

# This example uses -lsb/--lower_size_bound 1 because it is a small toy example
# with 25 vertices.  Use larger value (default is 10) for larger graphs.
for network in network_1
do
    for score in scores_1
    do
        python src/process_hierarchies.py \
            -oelf $intermediate/scores/hierarchy_dnreg_edge_list_index_0.tsv \
            -oigf $intermediate/scores/hierarchy_dnreg_index_to_protein_0.tsv \
            -pelf $(for i in `seq $num_permutations`; do echo " $intermediate/$scores/hierarchy_dnreg_edge_list_index_"$i".tsv "; done) \
            -pigf $(for i in `seq $num_permutations`; do echo " $intermediate/$scores/hierarchy_dnreg_index_to_protein_"$i".tsv "; done) \
            -lsb  1 \
            -cf   $results/clusters.tsv \
            -pl   "Network" "Score" \
            -pf   $results/sizes.pdf \
            -nc   $num_cores
    done
done

################################################################################
#
#   Perform consensus.
#
################################################################################

echo "Performing consensus..."

python src/perform_consensus.py \
    -cf  $results/clusters.tsv \
    -igf $data/dnreg_index_to_protein.tsv \
    -elf $data/dnreg_edge_list_index.tsv \
    -n   network \
    -s   scores \
    -t   2 \
    -cnf $results/consensus_nodes.tsv \
    -cef $results/consensus_edges.tsv
