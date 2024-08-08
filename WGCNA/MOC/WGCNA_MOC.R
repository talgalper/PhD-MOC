library(tidyverse)
library(WGCNA)
library(matrixStats)
library(reshape2)

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




# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
bwnet <- blockwiseModules(all_wgcna_data,
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
rownames(traits) <- c(colnames(all_subtypes), colnames(GTEx_ENS))
traits <- subset(traits, select = c("data.control", "data.tumour", "data.lumA", "data.lumB", "data.Her2","data.basal", "data.GTEx")) # reorder columns

moduleTrait_cor <- cor(bwnet$MEs, traits, use = "p")
moduleTrait_cor_pvals <- corPvalueStudent(moduleTrait_cor, nSamples = nrow(all_wgcna_data))

heatmap_data <- merge(bwnet$MEs, traits, by = "row.names")

heatmap_data <- column_to_rownames(heatmap_data, "Row.names")

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[18:23],
             y = names(heatmap_data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))








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




