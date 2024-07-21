
# Load multiWGCNA
library(multiWGCNA)
library(dplyr)

# Load data
wgcna_data_subset = read.csv('old _data/wgcna_data_subset.csv')
sampleTable = read.csv('old _data/sampleTable.csv')

# Clean up sampleTable
goodSampleTable = sampleTable[,c('cases', 'Subtype_Selected', 'ajcc_pathologic_stage')] %>% 
  setNames(c('Sample', 'Subtype', 'Stage'))
goodSampleTable$Sample = gsub('-', '.', goodSampleTable$Sample)
goodSampleTable$Stage = gsub('Stage ', '', goodSampleTable$Stage)

# Set up conditions
conditions1 = sort(unique(goodSampleTable$Subtype))
conditions2 = sort(unique(goodSampleTable$Stage))

# these four stages lead to error likely due to having too few samples
conditions2_fixed = conditions2[!conditions2 %in% c('IB', 'II', 'III', 'X')] 

# Clean up datExpr
rownames(wgcna_data_subset) = wgcna_data_subset$X
wgcna_data_subset = wgcna_data_subset[,-1]

# Check input
wgcna_data_subset[1:5,1:5]
head(goodSampleTable)

# Check that all column names match between sampleTable and datExpr
any(is.na(match(colnames(wgcna_data_subset), goodSampleTable$Sample)))

# Construct networks
galper_networks = constructNetworks(wgcna_data_subset, 
                                    goodSampleTable, 
                                    conditions1, 
                                    conditions2_fixed,
                                    networkType = "unsigned", 
                                    power = 10,
                                    minModuleSize = 40, 
                                    maxBlockSize = 45000,
                                    reassignThreshold = 0, 
                                    minKMEtoStay = 0.7,
                                    mergeCutHeight = 0.10, 
                                    numericLabels = TRUE,
                                    pamRespectsDendro = FALSE, 
                                    verbose=3,
                                    saveTOMs = FALSE)
