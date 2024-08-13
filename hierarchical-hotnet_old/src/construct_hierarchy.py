#!/usr/bin/python

# Load modules.
import math, numpy as np
import sys, argparse

from hierarchical_clustering import tarjan_HD, strongly_connected_components
from common import combined_similarity_matrix
from hhio import load_matrix, load_index_gene, load_gene_score, save_edge_list, save_index_gene

# Parse arguments.
def get_parser():
    description = 'Construct the hierarchical decomposition of the SCCs of the Hierarchical HotNet similarity matrix.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-smf', '--similarity_matrix_file', type=str, required=True, help='HH similarity matrix filename')
    parser.add_argument('-smn', '--similarity_matrix_name', type=str, required=False, default='PPR', help='HH similarity matrix name')
    parser.add_argument('-igf', '--index_gene_file', type=str, required=True, help='Index-gene filename')
    parser.add_argument('-gsf', '--gene_score_file', type=str, required=False, help='Gene-score filename')
    parser.add_argument('-lt', '--log_transform', action='store_true', help='Log transform scores')
    parser.add_argument('-st', '--score_threshold', type=float, required=False, default=float('nan'), help='Score threshold')
    parser.add_argument('-helf', '--hierarchy_edge_list_file', type=str, required=True, help='Hierarchy edge list filename')
    parser.add_argument('-higf', '--hierarchy_index_gene_file', type=str, required=True, help='Hierarchy index-gene filename')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    return parser

# Run script.
def run(args):
    # Load data.
    if args.verbose:
        print('Loading data...')

    index_to_gene, gene_to_index = load_index_gene(args.index_gene_file)
    P = load_matrix(args.similarity_matrix_file, args.similarity_matrix_name)

    if args.gene_score_file:
        gene_to_score = load_gene_score(args.gene_score_file)
        if not math.isnan(args.score_threshold):
            if args.log_transform:
                gene_to_score = dict((gene, score) for gene, score in gene_to_score.items() if score<=args.score_threshold)
            else:
                gene_to_score = dict((gene, score) for gene, score in gene_to_score.items() if score>=args.score_threshold)
        if args.log_transform:
            gene_to_score = dict((gene, -math.log10(score) if score!=1.0 else 0.0) for gene, score in gene_to_score.items())
    else:
        gene_to_score = dict((gene, 1.0) for gene in gene_to_index)

    # Process data.
    if args.verbose:
        print('Processing data...')

    S, common_index_to_gene, common_gene_to_index = combined_similarity_matrix(P, gene_to_index, gene_to_score)

    # If the digraph associated with S is not strongly connected, then restrict to a largest strongly connected component.
    components = strongly_connected_components(S)
    if len(components)>1:
        component = sorted(max(components, key=len))
        S = S[np.ix_(component, component)]
        common_index_to_gene = dict((i+1, common_index_to_gene[j+1]) for i, j in enumerate(component))
        common_gene_to_index = dict((gene, i) for i, gene in common_index_to_gene.items())

    # Construct hierarchical decomposition.
    if args.verbose:
        print('Constructing hierarchical decomposition...')

    T = tarjan_HD(np.asarray(S, dtype=np.float32), reverse=True, verbose=args.verbose)

    # Save results.
    if args.verbose:
        print('Saving results...')

    save_edge_list(args.hierarchy_edge_list_file, T)
    save_index_gene(args.hierarchy_index_gene_file, common_index_to_gene)

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
