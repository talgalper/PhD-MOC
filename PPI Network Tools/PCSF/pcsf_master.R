# setwd to PCSF folder

# PCSF one time setup
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
devtools::install_github("IOR-Bioinformatics/PCSF", repos=BiocManager::repositories(),
                         dependencies=TRUE, type="source", force=TRUE)

## PCSF parameters

#' @param ppi An interaction network, an \pkg{igraph} object.
#' @param terminals  A list of terminal genes with prizes to be analyzed in the PCSF context.
#' A named \code{numeric} vector, where terminal genes are named same as in the interaction network
#' and numeric values correspond to the importance of the gene within the study.
#' @param w A \code{numeric} value for tuning the number of trees in the output. A default value is 2.
#' @param b A \code{numeric} value for tuning the node prizes. A default value is 1.
#' @param mu A \code{numeric} value for a hub penalization. A default value is 0.0005.
#' @param dummies A list of nodes that are to connected to the root of the tree. If missing the root will be connected to all terminals.

library(tidyverse)
library(PCSF)

# read and format data.
# deframe() sets protein names to the frame of the df.

terminals <- read.table("data/snp_max_consequence_moc_noBen_scores.tsv") %>%
  tibble::deframe()

# Need to provide a data frame composed of three columns, where each row corresponds 
# to an edge in which the first element is a head node, the second element 
# is a tail node, and the last element represents the cost of the edge.

df <- read.table("data/wgcna-10th-power-MOC-edgelist-cost-filterDot01.txt", header = T)

# construct interactome and solve PCSF

ppi <- construct_interactome(df)

subnet <- PCSF(ppi, terminals, w = 1, b = 1, mu = 0.0005)


## edit plot features using following parameters. 
#' @param style A \code{boolean} value to determine the visualization style of the network,
#' where \code{0} plots the \code{static} network and \code{1} plots the \code{dynamic}
#' network. The default valu is 0.
#' @param edge_width A \code{numeric} value to emphasize a maximum edge width. A default value is 5.
#' This value must be greater than 1.
#' @param node_size A \code{numeric} value to emphasize a maximum node size. A default value is 40.
#' This value must be greater than 10.
#' @param node_label_cex A \code{numeric} value to set a node label size. A default value is 30.
#' @param Steiner_node_color A \code{string} to set a color for \code{Steiner} nodes.
#' A default value is "lightblue".
#' @param Terminal_node_color A \code{string} to set a color for \code{terminal} nodes.
#' @param Steiner_node_legend A \code{string} to set a legend for \code{Steiner} nodes.
#' A default legend is "Steiner".
#' @param Terminal_node_legend A \code{string} to set a legend for \code{terminal} nodes.


plot.PCSF(subnet)


# not working yet
library(topGO)
gene_universe <- V(ppi)$name
res <- enrichment_analysis(subnet, mode=1, gene_universe)




