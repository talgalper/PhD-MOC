import tensorflow as tf
from models import MQAModel
import numpy as np
from glob import glob
import mdtraj as md
import os
from tqdm import tqdm

from validate_performance_on_xtals import process_strucs, predict_on_xtals

def make_predictions(pdb_paths, model, nn_path, debug=False, output_basename=None):
    '''
        pdb_paths : list of pdb paths
        model : MQAModel corresponding to network in nn_path
        nn_path : path to checkpoint files
    '''
    strucs = [md.load(s) for s in pdb_paths]
    X, S, mask = process_strucs(strucs)
    if debug:
        np.save(f'{output_basename}_X.npy', X)
        np.save(f'{output_basename}_S.npy', S)
        np.save(f'{output_basename}_mask.npy', mask)
    predictions = predict_on_xtals(model, nn_path, X, S, mask)
    return predictions


# get list of files to be run from structures directory
structure_list = os.listdir("../structures/")


# loop pocketminer over all files in the structures directory. Uses tqdm to track progress
# output saved after each iteration to results directory
for structure in tqdm(structure_list):
    # main method
    if __name__ == '__main__':
        # TO DO - provide input pdb(s), output name, and output folder
        strucs = [
            f'../structures/{structure}']
        
        name_split = structure.split("-")
        uniprot_id = name_split[1]
        subunit_id = name_split[2]

        output_name = uniprot_id + "-" + subunit_id
        
        output_folder_txt = '../results/txt_results/'
        output_folder_npy = '../results/npy_results/'

        # debugging mode can be turned on to output protein features and sequence
        debug = False

        # Load MQA Model used for selected NN network
        nn_path = "../models/pocketminer.index"
        DROPOUT_RATE = 0.1
        NUM_LAYERS = 4
        HIDDEN_DIM = 100
        model = MQAModel(node_features=(8, 50), edge_features=(1, 32),
                         hidden_dim=(16, HIDDEN_DIM),
                         num_layers=NUM_LAYERS, dropout=DROPOUT_RATE)


        if debug:
            output_basename = f'{output_folder}/{output_name}'
            predictions = make_predictions(strucs, model, nn_path, debug=True, output_basename=output_basename)
        else:
            predictions = make_predictions(strucs, model, nn_path)

        # output filename can be modified here
        np.save(f'{output_folder_npy}/{output_name}-preds.npy', predictions)
        np.savetxt(os.path.join(output_folder_txt,f'{output_name}-predictions.txt'), predictions, fmt='%.4g', delimiter='\n')
