import requests
from xml.etree import ElementTree
import re
import argparse
import pandas as pd
import os


parser = argparse.ArgumentParser(description='Creates list of citaiton scores for list of proteins')
parser.add_argument('protein_list', type=str, help='CSV file with containing list of proteins, with header "protein_id"')
parser.add_argument('--output_dir', help='output file path', default=os.getcwd())
args = parser.parse_args()

protein_df = pd.read_csv(args.protein_list, header=0)
protein_list = protein_df['protein_id'].tolist()


def protein_citation_score(protein_name):
    # Use the PubMed API to search for articles that mention the protein
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={protein_name}'
    response = requests.get(url)
    root = ElementTree.fromstring(response.content)
    ids = [elem.text for elem in root.iterfind("./IdList/Id")]
    if not ids:
        return 0
    id_str = ",".join(ids)

    # Use the PubMed API to retrieve the articles and their citation count
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={id_str}&retmode=xml'
    response = requests.get(url)
    root = ElementTree.fromstring(response.content)

    # Initialize a variable to store the citation score
    citation_score = 0

    # Iterate through the articles and add the number of citations to the score
    for article in root.iterfind("./PubmedArticle"):
        citation = article.find("./PubmedData/ReferenceList/Reference/Citation")
        if citation is not None:
            citation_str = citation.text
            match = re.search(r'(\d+)', citation_str)
            if match:
                citation_score += int(match.group(1))

    return citation_score


# Create an empty dataframe to store the results
results = pd.DataFrame(columns=['protein', 'citation_score'])

for protein in protein_list:
        citation_count = protein_citation_score(protein)
        results = results.append({'protein': protein, 'citation_score': citation_count}, ignore_index=True)

# normalised scores
min_citation_score = results['citation_score'].min()
max_citation_score = results['citation_score'].max()

citation_score_range = max_citation_score - min_citation_score

results['normalized_citation_score'] = (results['citation_score'] - min_citation_score) / citation_score_range

results.to_csv(f"{args.output_dir}/citation_scores.csv", index=False)

