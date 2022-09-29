# setwd to PCSF folder

# PCSF one time setup
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
devtools::install_github("IOR-Bioinformatics/PCSF", repos=BiocManager::repositories(),
                         dependencies=TRUE, type="source", force=TRUE)


# start

library(tidyverse)
library(PCSF)

# read and format data.
# deframe() sets protein names to the frame of the df.
terminals <- read.table("data/snp_max_consequence_moc_noBen_scores.tsv") %>%
  tibble::deframe()

# Need to provide a data frame composed of three columns, where each row corresponds 
# to an edge in which the first element is a head node, the second element is a tail node, and the last element represents the cost of the edge.
df <- read.table("data/wgcna-10th-power-MOC-edgelist-cost-filterDot01.txt", header = T)

# construct interactome and solve PCSF
ppi <- construct_interactome(df)

subnet <- PCSF(ppi, terminals, w = 1, b = 1, mu = 0.0005)

# plot 
plot.PCSF(subnet)


# enrichment analysis not working yet
library(topGO)

subnet <- PCSF_rand(ppi, snpNOBEN, n = 10, r = 0.1, w = 2, b = 1, mu = 5e-04)

gene_universe <- V(ppi)$name
res <- enrichment_analysis(subnet, mode=1, gene_universe)

res <- enrichment_analysis(subnet)

plot(res$subent, edge_width = 8, node_size = 30, node_label_cex = 1)


### error ###
'''> res <- enrichment_analysis(subent, mode=1, gene_universe)
  Performing enrichment analysis...

  Enrichment is being performed by topGO package ...

Building most specific GOs .....
	( 0 GO terms found. )

Build GO DAG topology ..........
	( 0 GO terms and 0 relations. )
Nothing to do:
Error in split.default(names(sort(nl)), f.index) : 
  first argument must be a vector
In addition: Warning messages:
1: In cluster_edge_betweenness(subnet) :
  At core/community/edge_betweenness.c:486 : Membership vector will be selected based on the highest modularity score.
2: In cluster_edge_betweenness(subnet) :
  At core/community/edge_betweenness.c:493 : Modularity calculation with weighted edge betweenness community detection might not make sense -- 
  modularity treats edge weights as similarities while edge betwenness treats them as distances.'''


