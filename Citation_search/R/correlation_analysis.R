library(data.table)
library(igraph)

# PubTator3 counts and get total counts per gene
PubTator <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/citaiton_counts_ognsmAnnot.csv") # mac
PubTator$tax_id <- as.character(PubTator$tax_id)
PubTator$entrezgene_id <- as.character(PubTator$entrezgene_id)
PubTator[, combined := ifelse(symbol == "", entrezgene_id, symbol)]
PubTator_counts <- PubTator[
  , .(counts = sum(count)),
  by = combined
]
PubTator_counts <- PubTator_counts[order(-PubTator_counts$counts), ]
colnames(PubTator_counts)[1] <- "symbol"

# PCSF
load("~/Documents/GitHub/PhD-MOC/BRCA_pipe/latest_run/RData/STN_filt/PCSF_results.RData")
PCSF <- df
rm(df)

# HHnet
HHnet <- read.csv("../Hierarchical_HotNet/BRCA/STN_filt/results/df_subnet.csv")
HHnet_enrich <- read.csv("../Hierarchical_HotNet/BRCA/STN_filt/results/df_subnetNeighs.csv")

# STRING
STRING_net <- fread("../Hierarchical_HotNet/STRING_data/STRING_physical_geneSymbol.csv")
STRING_net <- STRING_net[!duplicated(t(apply(STRING_net, 1, sort))), ] # get rid of doubled up edges i.e. (a,b) (b,a) same edge weight
STRING_net <- graph_from_data_frame(STRING_net, directed = F)
STRING <- data.frame(symbol = V(STRING_net)$name,
                     degree = degree(STRING_net),
                     betweenness = betweenness(STRING_net),
                     closeness = closeness(STRING_net),
                     eigen_centrality = eigen_centrality(STRING_net)$vector)
STRING <- STRING[order(-STRING$degree), ]
rm(STRING_net)
gc()

summary(PubTator_counts$counts)
summary(PCSF$degree_centrality)
summary(HHnet$degree)

# normalise data
PCSF$degree_norm <- (PCSF$degree_centrality - min(PCSF$degree_centrality)) / (max(PCSF$degree_centrality) - min(PCSF$degree_centrality))
HHnet$degree_norm <- (HHnet$degree - min(HHnet$degree)) / (max(HHnet$degree) - min(HHnet$degree))
HHnet_enrich$degree_norm <- (HHnet_enrich$degree - min(HHnet_enrich$degree)) / (max(HHnet_enrich$degree) - min(HHnet_enrich$degree))
STRING$degree_norm <- (STRING$degree - min(STRING$degree)) / (max(STRING$degree) - min(STRING$degree))
PubTator_counts$counts_norm <- log10(PubTator_counts$counts)

hist(PubTator_counts$counts_norm)
hist(PCSF$degree_norm)
hist(HHnet$degree_norm)
hist(HHnet_enrich$degree_norm)


PCSF <- merge.data.table(as.data.table(PCSF), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
PCSF <- PCSF[order(-PCSF$degree_centrality), ]
HHnet <- merge.data.table(as.data.table(HHnet), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
HHnet <- HHnet[order(-HHnet$degree), ]
HHnet_enrich <- merge.data.table(as.data.table(HHnet_enrich), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
HHnet_enrich <- HHnet_enrich[order(-HHnet_enrich$degree), ]
STRING <- merge.data.table(as.data.table(STRING), PubTator_counts, by = "symbol", all.x = T)
STRING <- STRING[order(-STRING$degree), ]

# Spearman correlation - spearman for non-normal data
temp <- cor.test(
  PCSF$degree_norm, 
  PCSF$counts_norm,
  method = "spearman",
  use = "complete.obs"
)
temp$estimate^2
# repeat using kendall method - apparently better at handling ties
cor.test(
  PCSF$degree_norm, 
  PCSF$counts_norm,
  method = "kendall"
)

# repeat for HHnet
temp <- cor.test(
  HHnet$degree_norm, 
  HHnet$counts_norm,
  method = "spearman",
  use = "complete.obs"
)
temp$estimate^2
cor.test(
  HHnet$degree_norm, 
  HHnet$counts_norm,
  method = "kendall"
)

temp <- cor.test(
  HHnet_enrich$degree_norm, 
  HHnet_enrich$counts_norm,
  method = "spearman",
  use = "complete.obs"
)
temp$estimate^2
cor.test(
  HHnet_enrich$degree_norm, 
  HHnet_enrich$counts_norm,
  method = "kendall"
)

# repeat for STRING
temp <- cor.test(
  STRING$degree_norm, 
  STRING$counts_norm,
  method = "spearman",
  use = "complete.obs"
)
temp$estimate^2
cor.test(
  STRING$degree_norm, 
  STRING$counts_norm,
  method = "kendall"
)


# visualise
plot(PCSF$degree_centrality, PCSF$counts,
     xlab = "Degree", ylab = "Citation Count",
     main = "Correlation between Degree and Citation Count",
     pch = 19)
abline(lm(counts ~ degree_centrality, data = PCSF), col = "red")

plot(HHnet$degree, HHnet$counts,
     xlab = "Degree", ylab = "Citation Count",
     main = "Correlation between Degree and Citation Count",
     pch = 19)
abline(lm(counts ~ degree, data = HHnet), col = "red")

plot(HHnet_enrich$degree, HHnet_enrich$counts,
     xlab = "Degree", ylab = "Citation Count",
     main = "Correlation between Degree and Citation Count",
     pch = 19)
abline(lm(counts ~ degree, data = HHnet_enrich), col = "red")

plot(STRING$degree, STRING$counts,
     xlab = "Degree", ylab = "Citation Count",
     main = "Correlation between Degree and Citation Count",
     pch = 19)
abline(lm(counts ~ degree, data = STRING), col = "red")


