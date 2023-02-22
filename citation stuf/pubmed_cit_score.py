import requests
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt


def get_citation_count(protein_id):
    # Build the API request URL
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={protein_id}&retmax=1&usehistory=y"

    # Send the request to the API
    response = requests.get(url)

    # Parse the response to retrieve the count of citations
    count = int(response.text.split("<Count>")[1].split("</Count>")[0])

    return count

# Create an argument parser
parser = argparse.ArgumentParser()

parser.add_argument("protein_list", help="CSV file containing the list of proteins")
parser.add_argument("--out_dir", help="Directory where the results will be saved", default=os.getcwd())
parser.add_argument("--plot", help="Plot the distribution of the citation scores", action='store_true')

args = parser.parse_args()

# Read the protein list from the CSV file
protein_df = pd.read_csv(args.protein_list)
protein_list = protein_df['protein_id'].tolist()

# Create an empty dataframe to store the results
results = pd.DataFrame(columns=['protein_id', 'citation_score'])

# Iterate over the list of proteins and retrieve the citation count for each
for protein in protein_list:
    citation_count = get_citation_count(protein)
    results = results.append({'protein_id': protein, 'citation_score': citation_count}, ignore_index=True)

# Calculate the minimum and maximum citation scores
min_citation_score = results['citation_score'].min()
max_citation_score = results['citation_score'].max()

# Calculate the range of the citation scores
citation_score_range = max_citation_score - min_citation_score

# Apply the min-max scaling method to the citation scores
results['normalized_citation_score'] = (results['citation_score'] - min_citation_score) / citation_score_range

# Create a folder in the output directory to store output data
if not os.path.exists(f"{args.out_dir}/citation_results"):
    os.makedirs(f"{args.out_dir}/citation_results")

if args.plot:
    plt.scatter(results['protein_id'], results['citation_score'])
    plt.title("Distribution of Citation Scores")
    plt.xlabel("Protein ID")
    plt.ylabel("Citation Score")
    plt.savefig(f"{args.out_dir}/citation_results/citation_distribution_plot.png")

# Save the results dataframe to a CSV file
results.to_csv(f"{args.out_dir}/citation_results/citation_scores.csv", index=False)
