

wgcna_data_subset <- wgcna_data[sample(rownames(wgcna_data), 1500),]



selectedBarcodes <- c(colnames(LumA_unstranded), colnames(normal_unstranded))

sampleTable <- common[ ,-3]
sampleTable <- sampleTable[sampleTable$cases %in% selectedBarcodes, ]
rownames(sampleTable) <- sampleTable[ , 1]
sampleTable$ajcc_pathologic_stage <- gsub( "Stage ", "", sampleTable$ajcc_pathologic_stage)

# define variables for `constructNetworks` function
se <- SummarizedExperiment(assay = list(LumA_unstranded = wgcna_data_subset), 
                           colData = list(Sample = sampleTable$cases,
                                          Status = sampleTable$Subtype_Selected,
                                          Stage = sampleTable$ajcc_pathologic_stage),
                           rowData = rownames(wgcna_data_subset))

sampleTable <- DataFrame(sampleTable)

conditions1 = unique(sampleTable[,2])
conditions2 = unique(sampleTable[,3])


# Construct the combined networks and all the sub-networks
LumA_networks <- constructNetworks(se, sampleTable, conditions1, conditions2,
                                    networkType = "unsigned", power = 10,
                                    minModuleSize = 40, maxBlockSize = 25000,
                                    reassignThreshold = 0, minKMEtoStay = 0.7,
                                    mergeCutHeight = 0.10, numericLabels = TRUE,
                                    pamRespectsDendro = FALSE, verbose=3,
                                    saveTOMs = FALSE)











