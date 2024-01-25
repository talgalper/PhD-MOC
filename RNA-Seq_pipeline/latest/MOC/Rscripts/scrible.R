hits <- read.csv("intermediate/edgeR/edgeR_hits.csv", row.names = 1)

barplot(hits$logFC, 
        ylim = c(-6, 8),
        main = "Distribution of logFC values")

barplot(hits$logCPM, 
        ylim = c(0, 14),
        main = "Distribution of logCPM values")


barplot(rowSums(df),
        main = "Dsitribution of raw counts")

CPMdf <- cpm(df)

barplot(colSums(CPMdf),
        ylim = c(0, 1200000),
        main = "Distribution of sample CPM")
        

options(scipen = 999)


PCSF_results <- read.csv("results/PCSF_results.csv")


results <- read.csv("results/MOC_PCSF_master_unique.csv", row.names = 1)

description <- getBM(attributes = c("external_gene_name","description"),
                  filters = "external_gene_name", 
                  values = results$external_gene_name, 
                  mart = ensembl)



results <- subset(results, select = c("external_gene_name", "druggability", "num_drug_pockets"))

results <- merge(description, results, by = "external_gene_name")


results <- subset(kylie_pcsf_master_ranked, select = c("external_gene_name", "betweenness", "degree_centrality", "druggability", "citation_score", "combined_score"))
  
  

description <- getBM(attributes = c("external_gene_name","description"),
                     filters = "external_gene_name", 
                     values = results$external_gene_name, 
                     mart = ensembl)


results <- merge(description, results)
results$description <- gsub("\\s*\\[.*?\\]", "", results$description)



num_drug_poc <- read.csv("../../../../../../Desktop/lick my puss.csv", row.names = 1)
num_drug_poc <- num_drug_poc[1:25, ]

matched <- num_drug_poc[num_drug_poc$external_gene_name %in% final_gene_counts$external_gene_name, ]


ggplot(kylie_pcsf_master, aes(x = druggability, y = num_drug_pockets)) +
  geom_point() +
  theme(panel.background = element_blank()) +
  xlab("Largest druggable pocket") +
  ylab("Number of druggable pockets")


matched <- aggregate_ranks[1:25, ][!aggregate_ranks$external_gene_name[1:25] %in% gene_counts$gene[1:25], ]





all_run <- read.csv("../../../../../../Desktop/lick my puss.csv", row.names = 1)

matched <- all_run[1:25, ][!all_run$external_gene_name[1:25] %in% gene_counts$gene[1:25], ]




matched <- aggregate_ranks[1:100, ][aggregate_ranks$external_gene_name[1:100] %in% gene_counts$gene, ]


plot_data <- subset(kylie_pcsf_master, select = c("max_hit", "num_hits"))
plot_data <- plot_data[plot_data$max_hit >= 0.7, ]



ggplot(plot_data, aes(x = max_hit, y = num_hits)) +
  geom_point() +
  theme(panel.background = element_blank(),
        axis.title.x = element_text(size = 12, colour = "black", vjust = -3),
        axis.title.y = element_text(size = 12, colour = "black", vjust = 3),
        axis.text = element_text(size = 10, colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Largest cryptic pocket") +
  ylab("Number of cryptic pockets")


ggplot(kylie_pcsf_master, aes(x = druggability, y = max_hit)) +
  geom_point() +
  theme(panel.background = element_blank(),
        axis.title.x = element_text(size = 12, colour = "black", vjust = -3),
        axis.title.y = element_text(size = 12, colour = "black", vjust = 3),
        axis.text = element_text(size = 10, colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Largest druggable pocket") +
  ylab("Largest cryptic pocket")

model <- lm(kylie_pcsf_master$max_hit ~ kylie_pcsf_master$druggability)
summary(model)$r.squared


x <- read.csv("results/PCSF_results.csv")
description <- getBM(attributes = c("external_gene_name", "description"), 
                     filters = "external_gene_name", 
                     values = x$gene_id, 
                     mart = ensembl)

description$description <- gsub("\\s*\\[.*?\\]", "", description$description)

x <- merge(description, x, by.x = "external_gene_name", by.y = "gene_id")


data <- data.frame(names = c("a", "b", "c", NA, "e", "f"),
                   values = c(1, 2, NA, 4, 5, 6),
                   values2 = c(1, 2, 3, NA, 5, 6))





master <- data.frame(basic = gene_counts_basic$external_gene_name[1:25], 
                     NDP = final_gene_counts$external_gene_name[1:25])

master <- cbind(master, gene_counts$gene[1:25])
colnames(master)[3] <- "all"

master <- cbind(master, aggregate_ranks$external_gene_name[1:25])
colnames(master)[6] <- "RRA_all"



ranks <- data.frame(logFC = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$logFC)],
                    betweenness = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$betweenness)],
                    degree_centrality = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$degree_centrality)],
                    druggability = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$druggability)],
                    num_drug_pockets = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$num_drug_pockets)],
                    up_reg = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$logFC)], 
                    cryptic_pocket = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$max_hit)],
                    num_cryp_pockets = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$num_hits)],
                    citation = kylie_pcsf_master$external_gene_name[order(kylie_pcsf_master$citation_score)])


DE_scores <- read.csv("intermediate/PCSF_data/pcsf_kylie_score.csv")



ggplot(checkpoint_4, aes(x = druggability, y = cryp_pocket)) +
  geom_point()


enrichment_results <- read.csv("results/enrichR_results.csv", row.names = 1)
enrichment_results <- enrichment_results[order(enrichment_results$Cluster, enrichment_results$PValue), ]
enrichment_results <- enrichment_results[!duplicated(enrichment_results$Cluster) | duplicated(enrichment_results$Cluster, fromLast = TRUE), ]


enrichment_results <- enrichment_results %>%
  arrange(Cluster, PValue) %>%
  group_by(Cluster) %>%
  slice_head(n = 2) %>%
  ungroup()



library(RCy3)
createNetworkFromIgraph(kylie_subnet)
exportNetwork("../../../../../../Desktop/network", "graphML")





af_druggability <- subset(kylie_pcsf_master, select = c("external_gene_name", "druggability", "num_drug_pockets", "struct_score")) 

description <- getBM(attributes = c("external_gene_name", "description"), 
                     filters = "external_gene_name", 
                     values = af_drugability$external_gene_name, 
                     mart = ensembl)
description$description <- gsub("\\s*\\[.*?\\]", "", description$description)

af_druggability <- merge(description, af_druggability, by = "external_gene_name")
af_druggability <- af_druggability[order(-af_druggability$druggability), ]


af_druggability <- af_druggability[!duplicated(af_druggability$external_gene_name), ]


library(rDGIdb)

DGIdb <- queryDGIdb(kylie_pcsf_master_ranked$external_gene_name[1:25])
results <- byGene(DGIdb)
detailed_results <- resultSummary(DGIdb)




data <- read.csv("intermediate/citation_scores.csv")
logData <- log(data$citation_score)


hist(logData,
     xlab = "Citaion score",
     xlim = c(0, 12),
     main = "")

barplot(data$citation_score)



data <- read.csv("results/MOC_PCSF_master_unique.csv", row.names = 1)

hist(data$betweenness)
hist(log(data$betweenness))



hist(colSums(kylie_data[2:ncol(kylie_data)]))



attributes <- listAttributes(ensembl)
attributes_subset <- attributes[attributes$page == "feature_page", ]



library(progress)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attributes <- listAttributes(ensembl)

results <- data.frame(original_term = checkpoint_1$external_gene_name)

pb <- progress_bar$new(total = length(attributes$name))
for (attribute in attributes$name) {
  x <- getBM(attributes = c("external_gene_name", attribute), 
               filters = "external_gene_name", 
               values = "ANPEP", 
               mart = ensembl)
  
  results <- merge(results, x, by.x = "original_gene_name", by.y = "external_gene_name")
  
  pb$tick()
}


for (attribute in attributes$name) {
  for (gene in checkpoint_1$external_gene_name) {
    
    x <- getBM(attributes = c("external_gene_name", attribute), 
               filters = "external_gene_name", 
               values = gene, 
               mart = ensembl)
  }
}


library(convertid)

z <- checkpoint_1[19, ]
z <- z[, 1:2]

pb <- progress_bar$new(total = length(attributes$name))
for (attribute in attributes$name) {
  x <- convert.bm(dat = z, id = "external_gene_name", 
                  biom.data.set = "human", 
                  biom.attributes = c("external_gene_name", attribute), 
                  biom.filter = "external_gene_name")
  
  pb$tick()
}



x <- convert.bm(dat = z, id = "external_gene_name", 
                biom.data.set = "human")

x <- likely_symbol(syms = "ANPEP", orgnsm = "human", output = c("likely", "symbols", "all"), prev_sym = F)



x <- alias2Symbol("CD13", species = "Hs", exp)


x <- getBM(attributes = c("external_gene_name", "entrezgene_id"), 
           filters = "external_gene_name", 
           values = "ANPEP", 
           mart = ensembl)


kylie_pcsf_master <- subset(kylie_pcsf_master, select = c("external_gene_name", "max_hit", "num_hits"))

gene_description <- getBM(attributes = c("external_gene_name", "description"), 
                          filters = "external_gene_name", 
                          values = kylie_pcsf_master$external_gene_name, 
                          mart = ensembl)
gene_description$description <- gsub("\\s*\\[.*?\\]", "", gene_description$description)


kylie_pcsf_master <- merge(gene_description, kylie_pcsf_master, by = "external_gene_name")
colnames(kylie_pcsf_master) <- c("Gene", "Description", "Cryptic pocket score", "Number of cryptic pockets")
write.csv(kylie_pcsf_master, "~/Desktop/pocketminer_scores.csv")


missing_genes_converted <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                                 filters = "ensembl_gene_id", 
                                 values = missing_genes$gene_id, 
                                 mart = ensembl)

gene_description$description <- gsub("\\s*\\[.*?\\]", "", gene_description$description)



RS <- read.csv("~/Desktop/rank_sensitivity.csv")
RRA <- read.csv("~/Desktop/RRA.csv")

x <- RS25[RS25$external_gene_name %in% RRA25$external_gene_name, ]
y <- RRA25[RRA25$external_gene_name %in% RS25$external_gene_name, ]

`!RRA25` <- RS25[!(RS25$external_gene_name %in% RRA25$external_gene_name), ]
`!RS25` <- RRA25[!(RRA25$external_gene_name %in% RS25$external_gene_name), ]

RS25 <- RS[1:25, ]
RRA25 <- RRA[1:25, ]

common25 <- merge(RRA25, RS25, by = c("external_gene_name", "description"))
uncommon25 <- anti_join(RS25, RRA25)



data <- read.csv("~/Desktop/all_variable_sens.csv", row.names = 1)

DGIdb <- queryDGIdb(RRA25$external_gene_name)
results <- byGene(DGIdb)
detailed_results <- resultSummary(DGIdb)
x <- detailedResults(DGIdb)





