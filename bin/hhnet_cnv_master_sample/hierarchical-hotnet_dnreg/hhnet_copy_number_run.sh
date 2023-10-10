#!/usr/bin/env bash

data=$PWD/data
intermediate=$PWD/intermediate
results=$PWD/results

num_permutations=100

# This script uses gnu parallel to parallizes some of the scripts.  If gnu
# parallel is not already installed on your system, then you can install it from
# https://www.gnu.org/software/parallel/.  You can change the num_cores variable
# to specify the number of cores for your system.

num_cores=10

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
        -i   $data/"$network"_edge_list.tsv \
        -o   $intermediate/"$network"/similarity_matrix.h5 \
        -bof $intermediate/"$network"/beta.txt
done

################################################################################
#
#   Permute data.
#
################################################################################

# This example does not use permuted networks, but these commands show how to
# generate them.
echo "Permuting scores..."

for network in network_1
do
    for score in scores_1
    do
        cp $data/"$score".tsv $intermediate/"$network"_"$score"/scores_0.tsv

        python src/find_permutation_bins.py \
            -gsf $intermediate/"$network"_"$score"/scores_0.tsv \
            -igf $data/"$network"_index_gene.tsv \
            -elf $data/"$network"_edge_list.tsv \
            -ms  1000 \
            -o   $intermediate/"$network"_"$score"/score_bins.tsv

        parallel -u -j $num_cores --bar \
            python src/permute_scores.py \
                -i  $intermediate/"$network"_"$score"/scores_0.tsv \
                -bf $intermediate/"$network"_"$score"/score_bins.tsv \
                -s  {} \
                -o  $intermediate/"$network"_"$score"/scores_{}.tsv \
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
                -smf  $intermediate/"$network"/similarity_matrix.h5 \
                -igf  $data/"$network"_index_gene.tsv \
                -gsf  $intermediate/"$network"_"$score"/scores_{}.tsv \
                -helf $intermediate/"$network"_"$score"/hierarchy_edge_list_{}.tsv \
                -higf $intermediate/"$network"_"$score"/hierarchy_index_gene_{}.tsv \
            ::: `seq 0 $num_permutations`
    done
done

################################################################################
#
#   Process hierarchies.
#
################################################################################

echo "Processing hierarchies..."

for network in network_1
do
    for score in scores_1
    do
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
done

################################################################################
#
#   Perform consensus.
#
################################################################################

echo "Performing consensus..."

python src/perform_consensus.py \
    -cf  $results/clusters_network_1_scores_1.tsv \
    -igf $data/network_1_index_gene.tsv \
    -elf $data/network_1_edge_list.tsv \
    -n   network_1 \
    -s   scores_1 \
    -t   2 \
    -cnf $results/consensus_nodes.tsv \
    -cef $results/consensus_edges.tsv
