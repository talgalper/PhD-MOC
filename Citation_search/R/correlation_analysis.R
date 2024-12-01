

PCSF_counts <- read.csv("results/PCSF_citation_counts[TIAB].csv")

PCSF_degree <- read.csv("data/PCSF_results.csv")
PCSF_degree <- PCSF_degree[,c(2,7)]

PCSF_data <- merge(PCSF_degree, PCSF_counts, by = "external_gene_name", all = T)
PCSF_data <- PCSF_data[!is.na(PCSF_data$external_gene_name) & PCSF_data$external_gene_name != "", ]
PCSF_data <- PCSF_data[order(-PCSF_data$degree_centrality), ]

plot(PCSF_data$degree_centrality, PCSF_data$best_score,
     xlab = "Degree Centrality",
     ylab = "Citation Counts (best score)",
     main = "Scatter Plot of Degree Centrality vs. Citation Counts")


hist(log(PCSF_data$degree_centrality), main = "Histogram of Degree Centrality")
hist(log(PCSF_data$best_score), main = "Histogram of Citation Counts")

ad.test(PCSF_data$degree_centrality)
ad.test(PCSF_data$best_score)

cor.test(PCSF_data$degree_centrality, PCSF_data$MeSH_count, method = "spearman")





HHnet_counts <- read.csv("results/HHnet_DEneighs_citation_counts.csv")

HHnet_degree <- read.csv("data/df_subnetNeighs.csv")
HHnet_degree <- HHnet_degree[,c(2,5)]

HHnet_data <- merge(HHnet_degree, HHnet_counts, by = "external_gene_name", all = T)
HHnet_data <- HHnet_data[!is.na(HHnet_data$external_gene_name) & HHnet_data$external_gene_name != "", ]
HHnet_data <- HHnet_data[order(-HHnet_data$degree), ]


plot(HHnet_data$degree, HHnet_data$best_score,
     xlab = "Degree Centrality",
     ylab = "Citation Counts (best score)",
     main = "Scatter Plot of Degree Centrality vs. Citation Counts")


hist(log(HHnet_data$degree), main = "Histogram of Degree Centrality")
hist(log(HHnet_data$best_score), main = "Histogram of Citation Counts")

ad.test(HHnet_data$degree)
ad.test(HHnet_data$best_score)

cor.test(HHnet_data$degree, HHnet_data$MeSH_count, method = "spearman")



targets <- unique(BRCA_drugTargets_citation_counts$external_gene_name)

temp <- HHnet_data[HHnet_data$external_gene_name %in% targets, ]
temp2 <- PCSF_data[PCSF_data$external_gene_name %in% targets, ]


