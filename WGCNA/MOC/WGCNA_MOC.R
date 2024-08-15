library(tidyverse)
library(WGCNA)
library(matrixStats)
library(reshape2)
library(edgeR)
library(DESeq2)

# read in data
MOC_raw_counts <- read.csv("rna_seq_data/analysis_set_raw_counts_genenames.csv")
MOC_raw_counts_ENS <- read.csv("rna_seq_data/analysis_set_raw_counts.csv")
ID_mapping <- data.frame(hgnc_symbol = MOC_raw_counts$X,
                         ensembl_id = MOC_raw_counts_ENS$X)
MOC_raw_counts_ENS <- column_to_rownames(MOC_raw_counts_ENS, "X")

# shows duplicated terms when using hgnc symbols
# will need to use Ensembles for proceeding analysis
temp <- MOC_raw_counts[duplicated(MOC_raw_counts$X), ]
temp2 <- MOC_raw_counts[duplicated(MOC_raw_counts$X, fromLast = T), ]
temp <- rbind(temp, temp2)

sample_info <- read.csv("rna_seq_data/All survival_CN_Aug18.csv")
sample_info$GAMUT_ID <- paste0("GAMuT_", sample_info$GAMUT_ID) # match IDs

# subset sample info to only those present in data matrix
MOC_samples <- data.frame(GAMUT_ID = colnames(MOC_raw_counts_ENS))
sample_info <- merge(MOC_samples, sample_info, by = "GAMUT_ID", all.x = T)
sample_info$Stage[sample_info$Grade == "BEN"] <- "BEN"

# remove EOM and BDL samples
sample_info_subset <- sample_info[sample_info$Classification != "EOM" & sample_info$Classification != "BDL", ]
sample_info_subset <- sample_info_subset[!is.na(sample_info_subset$GAMUT_ID), ]

# format for PCA
sample_info_subset <- subset(sample_info_subset, select = c("GAMUT_ID", "Classification"))
colnames(sample_info_subset) <- c("sample", "group")

# subset RNA-seq data
MOC_data <- MOC_raw_counts_ENS[, colnames(MOC_raw_counts_ENS) %in% sample_info_subset$sample[sample_info_subset$group == "MOC"]]
BEN_data <- MOC_raw_counts_ENS[, colnames(MOC_raw_counts_ENS) %in% sample_info_subset$sample[sample_info_subset$group == "BEN"]]

# QC
counts_filt <- filter_low_expr(tumour_matrix = MOC_data,
                               control_matrix = BEN_data)

# perform DE if havent already
DE_results <- DE_analysis(counts_matrix = counts_filt,
                          sample_info = sample_info_subset)
dif_exp <- DE_results$dif_exp

save(dif_exp, file = "MOC/RData/MOC_dif_exp.RData")

# normalisation (transposes matrix)
MOC_data_norm <- vst_norm(counts_filt)

# plot PCA
MOC_PCA <- plot_PCA(expr_data = MOC_data_norm,
                    sample_info = sample_info_subset,
                    plot_tree = T,
                    output_plot_data = T)


# choose soft thresholding power
sft <- pickSoftThreshold(MOC_data_norm,
                         blockSize = 45000,
                         verbose = 2)

sft <- sft$fitIndices

library(gridExtra)
library(grid)
a1 <- ggplot(sft, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit') +
  theme_classic()

a2 <- ggplot(sft, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
rm(a1, a2)

save(sft, file = "MOC/RData/combined_sft.RData")
load("MOC/RData/combined_sft.RData")



# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
bwnet <- blockwiseModules(MOC_data_norm,
                          maxBlockSize = 45000,
                          power = 6,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3,
                          saveTOMs = FALSE)
elapsed_time <- Sys.time() - start_time
print(paste0("Elapsed time: ", elapsed_time))


# Plot the dendrogram
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


save(bwnet, sample_info, file = "MOC/RData/bwnet.RData")
load("MOC/RData/bwnet.RData")




## trait correlation
library(tidyverse)
library(CorLevelPlot)

sample_info <- read.csv("rna_seq_data/All survival_CN_Aug18.csv")
sample_info$GAMUT_ID <- paste0("GAMuT_", sample_info$GAMUT_ID) # match IDs

# subset sample info to only those present in data matrix
MOC_samples <- data.frame(GAMUT_ID = colnames(MOC_raw_counts_ENS))
sample_info <- merge(MOC_samples, sample_info, by = "GAMUT_ID", all.x = T)
sample_info$Stage[sample_info$Grade == "BEN"] <- "BEN"

# remove EOM and BDL samples
sample_info_subset <- sample_info[sample_info$Classification != "EOM" & sample_info$Classification != "BDL", ]
sample_info_subset <- sample_info_subset[!is.na(sample_info_subset$GAMUT_ID), ]

sample_info_subset <- subset(sample_info_subset, select = c("GAMUT_ID", "Classification", "Stage"))

stage_I <- c("I","IA","IB","IC")
stage_II <- c("II", "IIA", "IIB", "IIC")
stage_III <- c("III", "IIIA", "IIIB", "IIIC", "IIIc")
stage_IV <- c("IV")


# Create a function to consolidate stages
consolidate_stage <- function(stage) {
  if (stage %in% stage_I) {
    return("I")
  } else if (stage %in% stage_II) {
    return("II")
  } else if (stage %in% stage_III) {
    return("III")
  } else if (stage %in% stage_IV) {
    return("IV")
  } else {
    return(stage)
  }
}

# Apply the function to the Stage column
sample_info_subset$Consolidated_Stage <- sapply(sample_info_subset$Stage, consolidate_stage)

traits.state <- binarizeCategoricalColumns.forPlots(sample_info_subset$Classification)
traits.subtype <- binarizeCategoricalColumns.forPlots(sample_info_subset$Consolidated_Stage)
traits <- cbind(traits.state, traits.subtype)
rownames(traits) <- c(colnames(MOC_data), colnames(BEN_data))

moduleTrait_cor <- cor(bwnet$MEs, traits, use = "p")
moduleTrait_cor_pvals <- corPvalueStudent(moduleTrait_cor, nSamples = nrow(MOC_data_norm))

heatmap_data <- merge(bwnet$MEs, traits, by = "row.names")

heatmap_data <- column_to_rownames(heatmap_data, "Row.names")
heatmap_data <- heatmap_data[, -48] # remove duplucate BEN column

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[46:51],
             y = names(heatmap_data)[2:45],
             col = c("blue1", "skyblue", "white", "pink", "red"))



# intramodular connectivity
colours <- labels2colors(bwnet$colors)
kWithin <- intramodularConnectivity.fromExpr(MOC_data_norm, colours, power = 6)
rownames(kWithin) <- colnames(MOC_data_norm)
kWithin <- kWithin[order(-kWithin$kWithin), ]
save(kWithin, file = "MOC/RData/all_kwithin.RData")
load("MOC/RData/all_kwithin.RData")

# Does the same thing as above (intra modular connectivity)
module.membership.measure <- cor(bwnet$MEs, MOC_data_norm, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nrow(MOC_data_norm))
module.membership.measure.pvals <- as.data.frame(t(module.membership.measure.pvals))

# gene significance
gene.signf.corr <- cor(MOC_data_norm, traits$data.MOC, use = 'p')
gene.signf.corr.pvals <- as.data.frame(corPvalueStudent(gene.signf.corr, nrow(MOC_data_norm)))
gene.signf.corr.pvals <- rownames_to_column(gene.signf.corr.pvals)
gene.signf.corr.pvals <- gene.signf.corr.pvals[order(gene.signf.corr.pvals$V1), ]
rownames(gene.signf.corr.pvals) <- NULL
colnames(gene.signf.corr.pvals) <- c("gene_id", "pvalue")

# get top 10 genes for connectivity for each module
top_connectivity_genes = list()
for (module in unique(bwnet$colors)) {
  moduleGenes = names(bwnet$colors)[bwnet$colors == module]
  moduleKWithin = kWithin[moduleGenes, ]
  topNumGenes <- ceiling(length(moduleGenes) * 0.10)
  topModuleGenes = head(order(moduleKWithin$kWithin, decreasing = TRUE), topNumGenes)
  top_connectivity_genes[[module]] = moduleGenes[topModuleGenes]
  
  rm(moduleGenes, module, topModuleGenes, topNumGenes)
}


library(reshape2)
top_connectivity_genes <- melt(top_connectivity_genes)
colnames(top_connectivity_genes) <- c("ensembl_id", "module")
#top_connectivity_genes$DE <- top_connectivity_genes$ensembl_id %in% dif_exp$gene_id


# perform GO and pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(progress)
pb <- progress_bar$new(total = length(unique(bwnet$colors)))

# run GO enrichment
# turn this into a function at some point
all_GO <- list()
for (module in unique(bwnet$colors)) {
  genes <- names(bwnet$colors)[bwnet$colors %in% module]
  GO <- enrichGO(genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
  
  all_GO[[module]] <- GO
  
  rm(genes, GO, module)
  pb$tick()
}

save(all_GO, file = "../../../../OneDrive - RMIT University/PhD/large_git_files/WGCNA/MOC/WGCNA_GO_BP.RData")


GO_formatted <- data.frame()
for (i in seq_along(all_GO)) {
  module <- all_GO[[i]]
  module_name <- names(all_GO)[i]
  
  result <- module@result
  result_top <- head(result, 5)
  result_top$module <- rep(module_name, nrow(result_top))
  result_top$`-log(p.adjust)` <- -log(result_top$p.adjust)
  
  GO_formatted <- rbind(GO_formatted, result_top)
  
  rm(module, i, result, result_top)
}

# Function to convert GeneRatio to numeric
convert_gene_ratio <- function(gene_ratio) {
  # Split the string by "/"
  parts <- strsplit(gene_ratio, "/")[[1]]
  # Convert the parts to numeric and calculate the ratio
  ratio <- as.numeric(parts[1]) / as.numeric(parts[2])
  return(ratio)
}

GO_formatted$GeneRatio.num <- sapply(GO_formatted$GeneRatio, convert_gene_ratio)
save(GO_formatted, file = "BRCA/RData/all_default/signed/all_WGCNA_GO_BP.RData")

ggplot(data = GO_formatted, aes(x = module, y = Description, 
                                color = `p.adjust`, size = GeneRatio.num)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  labs(size = "Gene Ratio") +
  ggtitle("GO enrichment analysis (BP)")



## KEGG pwathway
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
pb <- progress_bar$new(total = length(unique(bwnet$colors)))

all_KEGG <- list()
for (module in unique(bwnet$colors)) {
  genes <- names(bwnet$colors)[bwnet$colors %in% module]
  genes_converted <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                           filters = "ensembl_gene_id", 
                           values = genes, 
                           mart = ensembl,
                           verbose = F)
  
  KEGG <- enrichKEGG(genes_converted$entrezgene_id, organism = "hsa", keyType = "ncbi-geneid")
  
  all_KEGG[[module]] <- KEGG
  
  rm(genes, genes_converted, KEGG, module)
  pb$tick()
}

KEGG_formatted <- data.frame()
for (i in seq_along(all_KEGG)) {
  module <- all_KEGG[[i]]
  module_name <- names(all_KEGG)[i]
  
  result <- module@result
  result_top <- head(result, 5)
  result_top$module <- rep(module_name, nrow(result_top))
  result_top$`-log(p.adjust)` <- -log(result_top$p.adjust)
  
  KEGG_formatted <- rbind(KEGG_formatted, result_top)
  
  rm(module, i, result, result_top)
}

# uses convert_gene_ratio function
KEGG_formatted$GeneRatio.num <- sapply(KEGG_formatted$GeneRatio, convert_gene_ratio)
save(KEGG_formatted, file = "BRCA/RData/all_default/signed/all_WGCNA_KEGG.RData")

ggplot(data = KEGG_formatted, aes(x = module, y = Description, 
                                  color = `p.adjust`, size = GeneRatio.num)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  labs(size = "Gene Ratio") +
  ggtitle("KEGG enrichment analysis")


## cross section of tumour modules and DE genes
load("MOC/RData/MOC_dif_exp.RData")

DE_genes_bwnet <- bwnet$colors[names(bwnet$colors) %in% dif_exp$gene_id]
DE_genes_bwnet <- as.data.frame(table(DE_genes_bwnet))
colnames(DE_genes_bwnet) <- c("module", "DE_genes")

temp <- as.data.frame(table(bwnet$colors))
colnames(temp) <- c("module", "total_size")
DE_genes_bwnet <- merge(temp, DE_genes_bwnet, by = "module", all = T)
DE_genes_bwnet[is.na(DE_genes_bwnet)] <- 0
DE_genes_bwnet$`proportion(%)` <- DE_genes_bwnet$DE_genes/DE_genes_bwnet$total_size * 100
DE_genes_bwnet <- DE_genes_bwnet[order(-DE_genes_bwnet$`proportion(%)`), ]


## venn diagram for cross 
DE_genes <- dif_exp$gene_id
tumour_associated <- names(bwnet$colors)[!bwnet$colors %in% c("tan", "salmon", "turquoise", "magenta", "pink")] # modules here are sig associated to control group
top_kwithin <- top_connectivity_genes$ensembl_id
top_gene_membership <- gene.signf.corr.pvals$gene_id[1:(length(gene.signf.corr.pvals$gene_id) * 0.1)]
save(DE_genes, tumour_associated, top_kwithin, top_gene_membership, file = "MOC/RData/venn_data.RData")

kWithin[rownames(kWithin) %in% "ENSG00000141510", ] #TP53

load("BRCA/RData/all_default/venn_data.RData")
library(VennDiagram)

venn.diagram(
  x = list(DE_genes = DE_genes, 
           tumour_associated = tumour_associated, 
           kwithin = top_kwithin,
           top_tumour_membership = top_gene_membership),
  category.names = c("DE genes", "Tumour associated", "Top10% Kwithin", "Top10% MM"),
  col = "transparent",  # set the color of the intersections to transparent
  fill = c("dodgerblue", "goldenrod1", "lightcoral", "mediumseagreen"),  # set colors for each category
  alpha = 0.5,  # set the transparency level of the circles
  cat.col = c("dodgerblue", "goldenrod1", "lightcoral", "mediumseagreen"),  # set colors for category labels
  cat.fontfamily = "Arial",  # set the font family for category labels
  cat.fontface = "bold",  # set the font face for category labels
  cat.fontsize = 10,  # set the font size for category labels
  cex = 1.5,  # increase the size of the circles
  margin = 0.1,  # set the margin size (proportion of the plot)
  filename = "MOC/RData/consensus_genes.png",
  disable.logging = TRUE
)


common_genes <- Reduce(intersect, list(DE_genes, tumour_associated, top_kwithin, top_gene_membership))

genes_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                         filters = "ensembl_gene_id", 
                         values = common_genes, 
                         mart = ensembl,
                         verbose = F)



GO <- enrichGO(common_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_formatted <- GO@result


###############################################################################
# Will now perform WGCNA on MOC and BEN groups separately for module     
# preservation analysis                                                    
###############################################################################
counts_filt <- filter_low_expr(tumour_matrix = MOC_data,
                               control_matrix = BEN_data,
                               sep = T)

MOC_data_norm <- vst_norm(counts_filt$tumour)
BEN_data_norm <- vst_norm(counts_filt$control)

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

network_modules <- function(WGCNA_data, Power) {
  start_time <- Sys.time()
  bwnet <- blockwiseModules(WGCNA_data,
                            maxBlockSize = 45000,
                            power = Power,
                            mergeCutHeight = 0.25,
                            numericLabels = FALSE,
                            randomSeed = 1234,
                            verbose = 3,
                            saveTOMs = FALSE)
  elapsed_time <- Sys.time() - start_time
  cat("Elapsed time: ")
  print(elapsed_time)
  
  # Plot the dendrogram
  plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                      c("unmerged", "merged"),
                      dendroLabels = FALSE,
                      addGuide = TRUE,
                      hang= 0.03,
                      guideHang = 0.05)
  
  return(bwnet)
}

MOC_bwnet <- network_modules(MOC_data_norm, Power = 6)
BEN_bwnet <- network_modules(BEN_data_norm, Power = 6)

save(tumour_bwnet, file = "MOC/RData/MOC_bwnet.RData")
save(control_bwnet, file = "MOC/RData/BEN_bwnet.RData")

load("MOC/RData/MOC_bwnet.RData")
load("MOC/RData/BEN_bwnet.RData")


# module preservation analysis
multidata <- multiData(Control = control_data, 
                       Tumour = tumour_data)
multicolour <- list(Control = control_bwnet$colors,
                    Tumour = tumour_bwnet$colors)

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
preserved_modules <- modulePreservation(multiData = multidata,
                                        multiColor = multicolour,
                                        dataIsExpr = T,
                                        quickCor = 1,
                                        randomSeed = 1234,
                                        verbose = 3,
                                        nPermutations = 100,
                                        maxModuleSize = max(max(table(tumour_bwnet$colors)), 
                                                            max(table(control_bwnet$colors))),                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time

# plot results
modulePreservation_plt <- plot_preserved_modules(preserved_modules)

save(preserved_modules, modulePreservation_plt, file = "BRCA/RData/all_default/modulePreservation(n=100).RData")














### old script for MOC ###
# combine stages
ben <- read.csv("rna_seq_data/ben_master_df.csv")
stage_I <- read.csv("rna_seq_data/stage_I_master_df.csv")
stage_II <- read.csv("rna_seq_data/stage_II_master_df.csv")
stage_III <- read.csv("rna_seq_data/stage_III_master_df.csv")
stage_IV <- read.csv("rna_seq_data/stage_IV_master_df.csv")

kylie_data <- merge(ben, stage_I, by = "X")
kylie_data <- merge(kylie_data, stage_II, by = "X")
kylie_data <- merge(kylie_data, stage_III, by = "X")
kylie_data <- merge(kylie_data, stage_IV, by = "X")

colnames(kylie_data)[1] <- "gene_id"
rownames(kylie_data) <- kylie_data$gene_id
kylie_data <- kylie_data[, -1] # remove gene ids for calcs


# create stage ID for each sample
stage_info <- data.frame(
  Sample = colnames(kylie_data),
  Stage = character(length(colnames(kylie_data)))
)

stages <- c("ben", "stage_I", "stage_II", "stage_III", "stage_IV")

for (stage in stages) {
  stage_data <- get(stage)
  stage_info$Stage[colnames(kylie_data) %in% colnames(stage_data)] <- stage
}


## Outlier detection

# Identify outliers. Searches for 0 variance genes and missing entries
gsg <- goodSamplesGenes(t(kylie_data))
summary(gsg)
gsg$allOK

# This method visualises outliers as cluster tree
htree <- hclust(dist(t(kylie_data)))
plot(htree)


# PCA
pca <- prcomp(t(kylie_data))
pca_data <- pca$x
pca_var <- pca$sdev^2

pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)


# Merge sample-stage mapping with PCA data
pca_data <- merge(pca_data, stage_info, by.x = "row.names", by.y = "Sample")

# Create a custom color palette for stages
stage_colors <- c("ben" = "blue", "stage_I" = "green", "stage_II" = "red", "stage_III" = "purple", "stage_IV" = "orange")

# Create the PCA plot with color mapping
ggplot(pca_data, aes(PC1, PC2, color = Stage)) +
  geom_point() +
  geom_text_repel(aes(label = row.names(pca_data)), size = 3) +  # Adjust the label size here
  scale_color_manual(values = stage_colors) + # Use the custom color palette
  theme_bw() +
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %'))


outliers <- c("GAMuT_23091", "GAMuT_41828",  # BEN
              "GAMuT_IC257") # Stage I


kylie_data_subset <- kylie_data[, !colnames(kylie_data) %in% outliers]
kylie_data_subset <- t(kylie_data_subset)
## Normalisation
# log transformation normalisation

# does this need to be done?
# if so should it go before or after removal of low activity genes?
data_log_norm <- log(kylie_data_subset + 1) 


hist(data_log_norm)



par(mfrow = c(1, 2))

hist(colSums(kylie_data_subset),
     main = "raw data hist")

hist(colSums(data_log_norm),
     main = "log transformed hist")

par(mfrow = c(1, 2))

barplot(colSums(kylie_data_subset),
        main = "raw data library size")

barplot(colSums(data_log_norm),
        main = "log transformed library size")

# top x% of samples
num_genes <- nrow(data_log_norm)
top_threshold <- ceiling(num_genes * 0.50)

# min required samples
num_samples <- ncol(data_log_norm)
min_required_samples <- ceiling(num_samples * 0.10) # plot sensitivity



# Extract the gene names and data columns
gene_names <- rownames(data_log_norm)
data_columns <- data_log_norm

# Initialise an empty data frame
results_df <- data.frame(gene_id = gene_names)

# Iterate through each data column and perform the test
for (col_id in seq_along(data_columns)) {
  # Get the column name and the data
  col_name <- colnames(data_columns)[col_id]
  col_data <- data_columns[, col_id]
  
  # Sort the data and mark genes within the top_threshold rows as TRUE
  sorted_indices <- order(col_data, decreasing = TRUE)
  top_indices <- sorted_indices[1:top_threshold]
  
  # Create a logical vector for gene presence in top_threshold
  gene_presence <- rep(FALSE, nrow(data_log_norm))
  gene_presence[top_indices] <- TRUE
  
  # Add to the results data frame
  results_df <- cbind(results_df, gene_presence)
}

# Rename columns 
colnames(results_df) <- c("gene_id", colnames(data_columns))


# Calculate the number of TRUE values for each gene across columns
gene_counts <- rowSums(results_df[, -1]) # Exclude the first column ("gene_id")
failed_genes <- gene_counts < min_required_samples
print(summary(failed_genes))

# Remove failed genes from data
wgcna_data <- subset(data_log_norm, !(rownames(data_log_norm) %in% gene_names[failed_genes]))

wgcna_data <- t(wgcna_data)




# Choose a set of soft-threshold powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(wgcna_data,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft_data <- sft$fitIndices


# Visualise to pick power
a1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# Convert matrix to numeric
wgcna_data[] <- sapply(wgcna_data, as.numeric)

soft_power <- 12 # selected based on number of samples
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(wgcna_data,
                          maxBlockSize = 15000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor

module_eigengenes <- bwnet$MEs

table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)



# clean up sample_info and filter for samples in wgcna_data
sample_info <- read.csv("rna_seq_data/All survival_CN_Aug18.csv")
sample_info <- subset(sample_info, select = c("GAMUT_ID", "Grade", "Stage"))

sample_info$GAMUT_ID <- paste0("GAMuT_", sample_info$GAMUT_ID)

wgcna_samples <- rownames(wgcna_data)

sample_info_filtered <- sample_info %>% 
  filter(GAMUT_ID %in% wgcna_samples)




# create list of stage ids 
stage_ids <- c("I","IA","IB","IC", 
               "II", "IIA", "IIB", "IIC", 
               "III", "IIIA", "IIIB", "IIIC", "IIIc", 
               "IV")

# remove unwanted samples
sample_info_filtered <- sample_info_filtered %>%
  filter(is.na(Stage) | Stage %in% stage_ids) %>%
  mutate(Stage = gsub("[AaBbCc]", "", Stage, ignore.case = TRUE)) 

# don't need grade column any more
sample_info_filtered <- sample_info_filtered[, -2]
rownames(sample_info_filtered) <- sample_info_filtered$GAMUT_ID


traits <- sample_info_filtered %>% 
  mutate(stage.I.vs.all = ifelse(Stage == 'I', 1, 0)) %>% 
  select(3)



sample_info_filtered$Stage <- factor(sample_info_filtered$Stage, levels = c("I", "II", "III", "IV"))

sample_traits <- binarizeCategoricalColumns(sample_info_filtered$Stage,
                                            includePairwise = FALSE,
                                            includeLevelVsAll = TRUE,
                                            minCount = 1)

traits <- cbind(traits, sample_traits)


# Define numbers of genes and samples
nSamples <- nrow(wgcna_data)
nGenes <- ncol(wgcna_data)


module_trait_corr <- cor(module_eigengenes, traits, use = 'p')
module_trait_corr_pvals <- corPvalueStudent(module_trait_corr, nSamples)


# visualize module-trait association as a heatmap
heatmap_data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap_data)

heatmap_data <- heatmap_data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[11:14],
             y = names(heatmap_data)[1:10],
             col = c("blue1", "skyblue", "white", "pink", "red"))



module_gene_mapping <- as.data.frame(bwnet$colors)
red_module <- module_gene_mapping %>% 
  filter(`bwnet$colors` == 'red') %>% 
  rownames()


# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module_membership_measure <- cor(module_eigengenes, wgcna_data, use = 'p')
module_membership_measure_pvals <- corPvalueStudent(module_membership_measure, nSamples)
module_membership_measure_pvals <- as.data.frame(t(module_membership_measure_pvals))


# Calculate the gene significance and associated p-values 
gene_signf_corr <- cor(wgcna_data, traits$data.IV.vs.all, use = 'p')
gene_signf_corr_pvals <- corPvalueStudent(gene_signf_corr, nSamples)
gene_signf_corr_pvals <- as.data.frame(gene_signf_corr_pvals)
gene_signf_corr_pvals <- rownames_to_column(gene_signf_corr_pvals)
gene_signf_corr_pvals <- gene_signf_corr_pvals[order(gene_signf_corr_pvals$V1), ]
rownames(gene_signf_corr_pvals) <- NULL


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

stage_IV_hub_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                            filters = "ensembl_gene_id", 
                            values = gene_signf_corr_pvals$rowname, 
                            mart = ensembl)

stage_IV_hub_genes <- merge(stage_IV_hub_genes, gene_signf_corr_pvals, by.x = "ensembl_gene_id", by.y = "rowname")

stage_IV_hub_genes <- stage_IV_hub_genes[order(stage_IV_hub_genes$V1), ]




