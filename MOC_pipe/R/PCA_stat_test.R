### Stat test to compare PCA groups (difference/similarity) ###


library(MVN)

pca <- as.matrix(cpm(counts_filt, log = TRUE))
pca <- prcomp(t(counts_filt), center = TRUE, scale. = TRUE)
pca_scores <- pca$x[, 1:3]

# Performing MANOVA using the PCs as response variables.
manova_model <- manova(pca_scores ~ sample_info$Classification)
summary(manova_model, test = "Wilks")

# Check normality of residuals univariately
resid_manova <- residuals(manova_model)
apply(resid_manova, 2, function(x) shapiro.test(x)$p.value)

# Check multivariate normality
library(MVN)
# Test multivariate normality using Mardia's test
mvn_result <- mvn(data = resid_manova, mvnTest = "mardia")
print(mvn_result$multivariateNormality)

library(biotools)
# test homogeneity of covariances.
boxm_test <- boxM(pca_scores, grouping = sample_info$Classification)
print(boxm_test)





library(vegan)
# Use the distance matrix (e.g., Euclidean) from the top PCs or the original expression data:
distance_matrix <- dist(pca_scores)
permanova_result <- adonis2(distance_matrix ~ sample_info$Classification, data = as.data.frame(pca_scores))
print(permanova_result)



# pull out eigenvectors for MOC samples
pca <- prcomp(t(cpm(counts_filt, log = TRUE)), scale. = TRUE, center = TRUE)
loadings_PC1 <- pca$rotation[, "PC1"]
top_genes_PC1 <- sort(abs(loadings_PC1), decreasing = TRUE)[1:25]
print(top_genes_PC1)

top_genes_PC1 <- id_annot(data = as.data.frame(top_genes_PC1), 
                          col_id = 0,
                          input_type = "ensembl_gene_id", 
                          convert_to = c("external_gene_name", "description"))

top_genes_PC1 <- top_genes_PC1[order(top_genes_PC1$top_genes_PC1, decreasing = TRUE), ]
rownames(top_genes_PC1) <- NULL


# check to see if any PCs correlate with technical variables
library(edgeR)
sample_info <- read.csv("data/All survival_CN_Aug18.csv")
sample_info <- sample_info[sample_info$GAMUT_ID %in% colnames(counts_filt), ]

pca <- prcomp(t(cpm(counts_filt, log = TRUE)), scale. = TRUE, center = TRUE)
pcs <- pca$x[,1:5]
cor.mat <- sapply(sample_info[,c(colnames(sample_info))],
                  function(x) cor(as.numeric(factor(x)), pcs[,1]))
print(sort(abs(cor.mat), decreasing = TRUE))

# surrogate variables
library(sva)
mod <- model.matrix(~ Classification, data=sample_info)
svobj <- sva(cpm(counts_filt, log = TRUE), mod)
length(svobj$sv)




