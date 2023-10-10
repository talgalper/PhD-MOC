library(DESeq2)


gtex_data <- read.table("gtex/gene_reads_2017-06-05_v8_ovary.gct", skip = 2)

# set col names of gtex data
colnames(gtex_data) <- gtex_data[1, ]
gtex_data <- gtex_data[-1, ]

# rename cols
colnames(gtex_data)[2] <- "gene_id"
colnames(gtex_data)[3] <- "gene_name"
gtex_data <- gtex_data[, -1]


tcga_stage_II <- read.csv("tcga/stage_II_master_df.csv", row.names = 1)

# remove number after dot point in ensembl id in both files
gtex_data$gene_id <- gsub("\\.[0-9]*$", "", gtex_data$gene_id)
tcga_stage_II$gene_id <- gsub("\\.[0-9]*$", "", tcga_stage_II$gene_id)


df <- inner_join(tcga_stage_II, gtex_data, by = c("gene_id", "gene_name"))
rownames(df) <- df$gene_id
df <- df[, -c(1:2)]
df[] <- apply(df, 2, as.numeric)

countdata <- as.matrix(df)

# Select the columns that start with "fpkm"
fpkm_cols <- grep("^fpkm", colnames(df))

# Assign the fpkm columns to the cancer group
cancer <- colnames(df)[fpkm_cols]

# Assign the rest of the columns to the healty group
healthy <- setdiff(colnames(df), cancer)

# Create the group variable
condition <- factor(c(rep("cancer", length(cancer)), rep("healthy", length(healthy))))

coldata <- data.frame(row.names=colnames(countdata), condition)

countdata[] <- apply(countdata, 2, as.integer)

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~condition)

keep <- rowSums(counts(dds)) > 1 
dds <- dds[keep,]

dds <- DESeq(dds)

rld <- rlogTransformation(dds) 
head(assay(rld)) 
hist(assay(rld))








