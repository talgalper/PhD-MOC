
selected_barcodes <- c(lumA_paired$normal, lumA_paired$lumA)

# drop duplciate samples
selected_barcodes <- selected_barcodes[!duplicated(selected_barcodes) & !duplicated(selected_barcodes, fromLast = TRUE)]


LumA_query <- GDCquery(project = "TCGA-BRCA",
                       access = "open", 
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq",
                       barcode = selected_barcodes)

GDCdownload(LumA_query)

LumA_data <- GDCprepare(LumA_query, summarizedExperiment = T)


LumA_unstranded <- assay(LumA_data, "unstranded")  
rownames(LumA_unstranded) <- gsub("\\.\\d+", "", rownames(LumA_unstranded))
LumA_unstranded <- as.data.frame(LumA_unstranded)



subtype_subset <- master[master$cases %in% colnames(LumA_unstranded), ]
subtype_subset <- subtype_subset[match(colnames(LumA_unstranded), subtype_subset$cases), ]

group <- factor(subtype_subset$Subtype_Selected)

counts_filt <- filterByExpr(LumA_unstranded, group = group)
counts_filt <- LumA_unstranded[counts_filt, ]



#### Differential expression analysis ####

data <- DGEList(counts = counts_filt, group = group)

design <- model.matrix(~group)

# Estimate a common negative binomial dispersion parameter for a DGE dataset with a general experimental design
common <- estimateGLMCommonDisp(data, design, verbose = T)

# Estimate the abundance-dispersion trend by Cox-Reid approximate profile likelihood.
trend <- estimateGLMTrendedDisp(common, design)

# Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each tag, 
# with expression levels specified by a log-linear model.
tagwise <- estimateGLMTagwiseDisp(trend, design)

# Fit a negative binomial generalized log-linear model to the read counts for each gene. 
# Conduct genewise statistical tests for a given coefficient or coefficient contrast.
fit <- glmFit(tagwise, design)

# Conduct genewise statistical tests for a given coefficient or coefficient contrast.
lrt <- glmLRT(fit, coef = 2)

# Extract the most differentially expressed genes (or sequence tags) from a test object, 
# ranked either by p-value or by absolute log-fold-change.
toptags <- topTags(lrt, n = Inf)

# Identify which genes are significantly differentially expressed from 
# an edgeR fit object containing p-values and test statistics.
dif_exp <- decideTestsDGE(lrt, p = 0.05, adjust = "fdr", lfc = 1)
print(summary(dif_exp))

dif_exp_genes <- rownames(tagwise)[as.logical(dif_exp)]

# create a results df
hits <- toptags$table[toptags$table$FDR < 0.1, ]
hits <- rownames_to_column(hits)
rownames(hits) <- hits$rowname
colnames(hits)[1] <- "gene_id"

dif_exp <- hits[dif_exp_genes, ]

# Plot the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
plotBCV(tagwise, main = "biological coefficient of variation")

# Make a mean-difference plot of two libraries of count data with smearing of points with very low counts, 
# especially those that are zero for one of the columns.
plotSmear(lrt, de.tags = dif_exp_genes, main = "Mean-difference plot")

# plot Pvalues of different logFC scores
ggplot(hits, aes(x=logFC, y=-log(FDR))) + geom_point() + labs(title = "Adjusted logFC")

paired_hits <- toptags$table
paired_hits_subset <- paired_hits[paired_hits$FDR < 0.1, ]



