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
gene_universe <- V(ppi)$name
res <- enrichment_analysis(subnet, mode=1, gene_universe)




