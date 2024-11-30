#!/usr/bin/env python3

import sys

# List of required modules
required_modules = [
    'requests',
    'pandas',
    'numpy',
    'tqdm',
    'nltk',
    'argparse',
    'pickle'
]

# Function to check for required modules
def check_modules(modules):
    missing_modules = []
    for module in modules:
        try:
            __import__(module)
        except ImportError:
            missing_modules.append(module)
    if missing_modules:
        print("The following required modules are missing:")
        for mod in missing_modules:
            print(f"- {mod}")
            
        response = input("Do you want to install the missing modules now? (y/n): ")
        if response.lower() == 'y':
            import subprocess
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', *missing_modules])
        else:
            sys.exit(1)

# Check for required modules
check_modules(required_modules)

import requests
import pandas as pd
import numpy as np
from tqdm import tqdm
import nltk
from nltk.corpus import words
import argparse
import urllib.parse
import pickle
import os

# Ensure the NLTK 'words' corpus is downloaded
try:
    nltk.data.find('corpora/words')
except LookupError:
    print("Downloading NLTK 'words' corpus...")
    nltk.download('words')
print(f'Download location: {nltk.data.find("corpora/words")}')

def get_citation_count(search_query):
    """
    Query PubMed for citation count
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        'db': 'pubmed',
        'term': search_query,
        'retmax': 1,
        'usehistory': 'y'
    }

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        count = int(response.text.split("<Count>")[1].split("</Count>")[0])
        return count
    except requests.exceptions.RequestException as e:
        print(f"Network error while retrieving citation count for query '{search_query}': {e}")
        return None
    except IndexError:
        print(f"Error: Unable to retrieve citation count for query '{search_query}'")
        print("Response:", response.text)
        return None

def get_gene_citations(input_file, cosmic_file, out_filename):
    """
    Process gene citations from input file
    """
    try:
        # Read the input file
        print("Reading input file...")
        df = pd.read_csv(input_file)
        
        if 'external_gene_name' not in df.columns:
            raise ValueError("Input file must contain 'external_gene_name' column")
        
        if 'description' not in df.columns:
            raise ValueError("Input file must contain 'description' column")

        # Print the number of genes to process
        print(f"Number of genes to process: {len(df)}")
        
        # Create PubMed_search term and swap small gene names with description    
        df['PubMed_search'] = np.where(
            df['external_gene_name'].str.len() > 2,
            df['external_gene_name'],
            df['description']
        )

        # Drop NA rows and get list of genes to search
        df = df.dropna(subset=['external_gene_name'])

        # Remove commas and everything that follows in 'PubMed_search' column
        df['PubMed_search'] = df['PubMed_search'].str.split(",").str[0]

        # Define MeSH query component
        mesh_query = '(Neoplasms[MH] AND Humans[MH])'
        
        # Initialise results DataFrame
        results = pd.DataFrame(columns=['external_gene_name', 'description', 'search_id', 'MeSH_count'])
        
        # Check if checkpoint file exists
        checkpoint_file = 'checkpoint.pkl'
        if os.path.exists(checkpoint_file):
            print("\nCheckpoint file found. Resuming from last saved point...")
            with open(checkpoint_file, 'rb') as f:
                results = pickle.load(f)
            # Get the list of genes already processed
            processed_genes = set(results['external_gene_name'])
            # Filter out already processed genes
            df = df[~df['external_gene_name'].isin(processed_genes)]
            print(f"Number of genes remaining to process: {len(df)}")
        else:
            # Initialize results DataFrame
            results = pd.DataFrame(columns=['external_gene_name', 'description', 'search_id', 'MeSH_count'])
            processed_genes = set()

        # Get citations
        print("\nQuerying PubMed citations...")
        checkpoint_interval = 100  # Save checkpoint every 100 iterations
        iteration = 0

        for index, row in tqdm(df.iterrows(), total=len(df)):
            gene = row['PubMed_search']
            external_gene_name = row['external_gene_name']
            description = row['description']
            MeSH_query = f'{gene}[TIAB] AND "{mesh_query}"'
            MeSH_count = get_citation_count(MeSH_query)
            results = pd.concat([results, pd.DataFrame({
                'external_gene_name': [external_gene_name],
                'description': [description],
                'search_id': [gene], 
                'MeSH_count': [MeSH_count]
            })], ignore_index=True)
            
            iteration += 1
            # Save checkpoint at specified intervals
            if iteration % checkpoint_interval == 0:
                with open(checkpoint_file, 'wb') as f:
                    pickle.dump(results, f)

        # Save final checkpoint after loop completion
        with open(checkpoint_file, 'wb') as f:
            pickle.dump(results, f)
        print("Checkpoint file saved.")
        
        # Check for missing entries and update score
        def update_nan_values(count_results):
            nan_count = count_results['MeSH_count'].isna().sum()
            if nan_count > 0:
                print(f"\nUpdating {nan_count} NaN values...")
            else:
                print("\nNo NaN values to update.")
                return count_results

            for index, row in tqdm(count_results[count_results['MeSH_count'].isna()].iterrows(), total=nan_count):
                gene = row['search_id']
                MeSH_query = f'"{gene}"[TIAB] AND {mesh_query}'
                new_value = get_citation_count(MeSH_query)

                # Update the DataFrame
                count_results.at[index, 'MeSH_count'] = new_value

            print(f"Number of NaN values after update: {count_results['MeSH_count'].isna().sum()}")
            return count_results

        # Use the function to update NaN values
        updated_results = update_nan_values(results)
        updated_results['MeSH_count'] = updated_results['MeSH_count'].astype(int, errors='ignore')
        
        # Save checkpoint
        results_TIAB = updated_results

        # Function to check if searched gene is an English word
        english_words = set(word.upper() for word in words.words())

        def is_english_word(word):
            if isinstance(word, str):
                return word.upper() in english_words
            else:
                return False

        # Number of English words in set
        num_english_words = len([gene for gene in results_TIAB['search_id'] if gene.upper() in english_words])
        print(f"\nNumber of English words in gene set: {num_english_words}")

        # Add 'andGene_count' column
        results_TIAB['andGene_count'] = np.nan
        # Check for English words and rerun with "gene" following gene name
        print("\nUpdating counts for English words with 'gene' appended...")
        for idx, row in tqdm(results_TIAB.iterrows(), total=len(results_TIAB)):
            gene = row['search_id']

            # Only process if it's an English word
            if is_english_word(gene):
                try:
                    # Construct new query with gene[TIAB]
                    new_query = f'"{gene} gene"[TIAB] AND {mesh_query}'
                    count = get_citation_count(new_query)
                    results_TIAB.loc[idx, 'andGene_count'] = count
                except Exception as e:
                    print(f"Error processing gene {gene}: {str(e)}")
                    results_TIAB.loc[idx, 'andGene_count'] = np.nan

        results_TIAB['andGene_count'] = results_TIAB['andGene_count'].astype('Int64')
        
        # Add 'description_count' column
        results_TIAB['description_count'] = np.nan
        # Use descriptions for English words 
        print("\nUpdating counts using descriptions for English words...")
        for idx, row in tqdm(results_TIAB.iterrows(), total=len(results_TIAB)):
            gene = row['search_id']

            if is_english_word(gene):
                try:
                    description = row['description']
                    # Construct query with description
                    desc_query = f'"{description}"[TIAB] AND {mesh_query}'
                    count = get_citation_count(desc_query)
                    results_TIAB.loc[idx, 'description_count'] = count

                except Exception as e:
                    print(f"Error processing description for gene {gene}: {str(e)}")
                    results_TIAB.loc[idx, 'description_count'] = np.nan

        results_TIAB['description_count'] = results_TIAB['description_count'].astype('Int64')
        
        # Read COSMIC hallmark genes from the provided file
        COSMIC_hallmark_genes_df = pd.read_csv(cosmic_file, sep='\t')
        COSMIC_hallmark_genes = COSMIC_hallmark_genes_df['GENE_SYMBOL'].tolist()
        COSMIC_hallmark_genes = list(set(COSMIC_hallmark_genes))

        # Get "best score" 
        def get_best_score(row):
            # Check if the gene is a cancer hallmark gene
            if row['external_gene_name'] in COSMIC_hallmark_genes:
                return row['MeSH_count']
            else:
                # If both 'andGene_count' and 'description_count' are NA, use 'MeSH_count'
                if pd.isna(row['andGene_count']) and pd.isna(row['description_count']):
                    return row['MeSH_count']
                # Otherwise, return the highest non-NA value between 'andGene_count' and 'description_count'
                return max(
                    row['andGene_count'] if not pd.isna(row['andGene_count']) else float('-inf'),
                    row['description_count'] if not pd.isna(row['description_count']) else float('-inf')
                )

        # Add the new column to results_TIAB
        results_TIAB['best_score'] = results_TIAB.apply(get_best_score, axis=1)
        
        # Save results
        output_file = out_filename
        results_TIAB.to_csv(output_file, index=False)
        print(f"\nResults saved to: {output_file}")
        
        # Remove the checkpoint file after successful completion
        if os.path.exists(checkpoint_file):
            os.remove(checkpoint_file)
            print("Checkpoint file removed.")
            
        return results_TIAB
        
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Process gene citations from CSV file')
    parser.add_argument('input_file', help='Path to input CSV file containing gene names (external_gene_name) and description')
    parser.add_argument('cosmic_file', help='Path to COSMIC hallmark genes CSV file')
    parser.add_argument('output_file', help='Output file name to save results') 
    args = parser.parse_args()
    
    get_gene_citations(args.input_file, args.cosmic_file, args.output_file)

if __name__ == "__main__":
    main()
