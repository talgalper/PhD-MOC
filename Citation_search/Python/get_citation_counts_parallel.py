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
    'pickle',
    'concurrent.futures',
    'threading',
    'time'
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
import concurrent.futures  # For parallelization
import threading  # For thread-safe progress bar and rate limiting
import time  # For rate limiting

# Ensure the NLTK 'words' corpus is downloaded
try:
    nltk.data.find('corpora/words')
except LookupError:
    print("Downloading NLTK 'words' corpus...")
    nltk.download('words')
print(f'Download location: {nltk.data.find("corpora/words")}')

# Rate limiting semaphore
RATE_LIMIT = 3  # Maximum 3 requests per second
semaphore = threading.Semaphore(RATE_LIMIT)
lock = threading.Lock()  # For thread-safe operations

def get_citation_count(search_query):
    """
    Query PubMed for citation count in parallel with rate limiting
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        'db': 'pubmed',
        'term': search_query,
        'retmax': 1,
        'usehistory': 'y'
    }

    # Acquire semaphore before making a request
    with semaphore:
        try:
            # Enforce rate limit by sleeping for the appropriate duration
            time.sleep(1 / RATE_LIMIT)
            response = requests.get(url, params=params, timeout=10)
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
        
        # Initialize results DataFrame
        results_columns = ['external_gene_name', 'description', 'search_id', 'MeSH_count']
        results = pd.DataFrame(columns=results_columns)
        
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
            processed_genes = set()

        # Function to process each gene (for parallel execution)
        def process_gene(row):
            gene = row['PubMed_search']
            external_gene_name = row['external_gene_name']
            description = row['description']
            MeSH_query = f'{gene}[TIAB] AND "{mesh_query}"'
            MeSH_count = get_citation_count(MeSH_query)
            return {
                'external_gene_name': external_gene_name,
                'description': description,
                'search_id': gene,
                'MeSH_count': MeSH_count
            }

        # Get citations in parallel
        print("\nQuerying PubMed citations in parallel with rate limiting (3req/sec)...")
        max_workers = 3  # Since we can make at most 3 requests per second
        progress_bar = tqdm(total=len(df))

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for index, row in df.iterrows():
                futures.append(executor.submit(process_gene, row))

            results_list = []
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                results_list.append(result)
                with lock:
                    progress_bar.update(1)
            progress_bar.close()

        # Convert results to DataFrame
        new_results = pd.DataFrame(results_list)
        # Combine with any previously saved results
        results = pd.concat([results, new_results], ignore_index=True)
        # Save checkpoint
        with open(checkpoint_file, 'wb') as f:
            pickle.dump(results, f)
        print("Checkpoint file saved.")
        
        # Proceed with the rest of the processing using 'results'

        # Check for missing entries and update score
        def update_nan_values(count_results):
            nan_count = count_results['MeSH_count'].isna().sum()
            if nan_count > 0:
                print(f"\nUpdating {nan_count} NaN values...")
            else:
                print("\nNo NaN values to update.")
                return count_results

            # Prepare the DataFrame with NaN MeSH_count
            nan_df = count_results[count_results['MeSH_count'].isna()].copy()

            # Function to update MeSH_count for NaN entries
            def update_gene(row):
                gene = row['search_id']
                MeSH_query = f'"{gene}"[TIAB] AND {mesh_query}'
                new_value = get_citation_count(MeSH_query)
                return row.name, new_value

            # Update in parallel
            print("\nUpdating NaN MeSH_counts with rate limiting...")
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {executor.submit(update_gene, row): idx for idx, row in nan_df.iterrows()}
                for future in concurrent.futures.as_completed(futures):
                    idx, new_value = future.result()
                    count_results.at[idx, 'MeSH_count'] = new_value

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
        # Prepare data for parallel processing
        english_word_indices = results_TIAB[results_TIAB['search_id'].apply(is_english_word)].index

        # Function to update 'andGene_count'
        def update_andGene_count(idx_row):
            idx, row = idx_row
            gene = row['search_id']
            if pd.isna(row['andGene_count']):
                try:
                    new_query = f'"{gene} gene"[TIAB] AND {mesh_query}'
                    count = get_citation_count(new_query)
                    return idx, count
                except Exception as e:
                    print(f"Error processing gene {gene}: {str(e)}")
                    return idx, np.nan
            else:
                return idx, row['andGene_count']

        # Update 'andGene_count' with rate limiting
        print("\nUpdating counts for English words with 'gene' appended with rate limiting...")
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(update_andGene_count, (idx, results_TIAB.loc[idx])): idx for idx in english_word_indices}
            for future in concurrent.futures.as_completed(futures):
                idx, count = future.result()
                results_TIAB.at[idx, 'andGene_count'] = count

        results_TIAB['andGene_count'] = results_TIAB['andGene_count'].astype('Int64')
        
        # Add 'description_count' column
        results_TIAB['description_count'] = np.nan

        # Function to update 'description_count'
        def update_description_count(idx_row):
            idx, row = idx_row
            gene = row['search_id']
            if pd.isna(row['description_count']):
                try:
                    description = row['description']
                    desc_query = f'"{description}"[TIAB] AND {mesh_query}'
                    count = get_citation_count(desc_query)
                    return idx, count
                except Exception as e:
                    print(f"Error processing description for gene {gene}: {str(e)}")
                    return idx, np.nan
            else:
                return idx, row['description_count']

        # Update 'description_count' with rate limiting
        print("\nUpdating counts using descriptions for English words with rate limiting...")
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(update_description_count, (idx, results_TIAB.loc[idx])): idx for idx in english_word_indices}
            for future in concurrent.futures.as_completed(futures):
                idx, count = future.result()
                results_TIAB.at[idx, 'description_count'] = count

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
