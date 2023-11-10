
library(biomaRt)
library(tidyverse)
library(GenomicRanges)
library(gridExtra)

cnv_data <- read.csv("TCGA_disease_data.csv", row.names = 1)


# load ensemble database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


## weighted seg mean function
regions <- paste0(cnv_data$Chromosome, ":", cnv_data$Start, "-", cnv_data$End)

annotated_regions <- getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
                           filters = "chromosomal_region",
                           values = regions,
                           mart = ensembl)

load("annotated_regions.RData")


# delete rows with no gene names
annotated_regions <- annotated_regions[annotated_regions$external_gene_name !="", ]
# calculate gene sizes
annotated_regions <- annotated_regions %>%
  mutate(gene_size = end_position - start_position)

# sort by gene size and keep only the largest gene for each gene name
annotated_regions <- annotated_regions %>%
  arrange(desc(gene_size)) %>%
  group_by(external_gene_name) %>%
  filter(row_number() == 1) %>%
  ungroup()

# remove the gene size column
annotated_regions <- annotated_regions[, -5]

# create loc_id
grouped_data <- group_by(cnv_data, Chromosome)
cnv_data <- mutate(grouped_data, segment_number = row_number())
cnv_data <- mutate(cnv_data, loc_id = paste(Chromosome, segment_number, sep = "."))
cnv_data <- cnv_data[ ,-7]



# create GRanges object for cnv_data
cnv_gr <- GRanges(
  seqnames = cnv_data$Chromosome,
  IRanges(start = cnv_data$Start, end = cnv_data$End),
  loc_id = cnv_data$loc_id,
  seg_mean = cnv_data$Segment_Mean
)

# create GRanges object for annotated data
annot_gr <- GRanges(
  seqnames = annotated_regions$chromosome_name,
  IRanges(start = annotated_regions$start_position, end = annotated_regions$end_position),
    gene_id = annotated_regions$external_gene_name
)


# find overlapping genes
overlaps <- findOverlaps(cnv_gr, annot_gr)
overlap_genes <- as.character(annot_gr$gene_id[subjectHits(overlaps)])
overlap_df <- annotated_regions[annotated_regions$external_gene_name %in% overlap_genes, ]

# create data frame with overlapping loc_ids for each gene
overlap_df <- data.frame(gene_id = overlap_genes, loc_id = cnv_gr[queryHits(overlaps)]$loc_id)
overlap_df <- aggregate(loc_id ~ gene_id, overlap_df, paste, collapse = ";")


# create frequency table of the number of loc_id's per gene
num_loc_ids <- table(sapply(strsplit(as.character(overlap_df$loc_id), ";"), length))
print(num_loc_ids)

# create df of genes that do not overlap with any segment
non_overlap_genes <- annotated_regions$external_gene_name[!annotated_regions$external_gene_name %in% overlap_genes]
non_overlap_df <- annotated_regions[annotated_regions$external_gene_name %in% non_overlap_genes, ]


# separate genes based on how many segments they overlap and assign names to the dfs
overlap_df$loc_id_count <- sapply(strsplit(overlap_df$loc_id, ";"), length)
dfs <- list()
for (i in unique(overlap_df$loc_id_count)) {
  df_name <- paste0("df_", i)  # generate unique name for each df
  dfs[[df_name]] <- subset(overlap_df, loc_id_count == i)
}


perform_weighting <- function(df, num_overlaps) {
  
  # Select the overlaps dataframe with the specified number of overlaps
  overlaps <- df[[as.character(num_overlaps)]]
  
  if (!is.null(overlaps)) {
    
    # Split loc_id into overlap_a, overlap_b, etc. columns
    overlap_cols <- paste0("overlap_", rep(1:num_overlaps))
    overlaps <- overlaps %>% separate(loc_id, overlap_cols, sep = ";")
    
    # Add columns with the length of each overlap
    length_cols <- paste0("length_", rep(1:num_overlaps))
    for (i in seq_along(overlap_cols)) {
      overlaps <- overlaps %>%
        mutate(!!length_cols[i] := cnv_data$End[match(!!sym(overlap_cols[i]), cnv_data$loc_id)] - 
                 cnv_data$Start[match(!!sym(overlap_cols[i]), cnv_data$loc_id)] + 1)
    }
    
    # Calculate weighted segment mean for each overlap
    weighted_mean_cols <- paste0("weighted_seg_mean_", rep(1:num_overlaps))
    for (i in seq_along(overlap_cols)) {
      overlaps <- overlaps %>%
        mutate(!!weighted_mean_cols[i] := (!!sym(length_cols[i]) / rowSums(select(overlaps, starts_with("length")))) * 
                 cnv_data$Segment_Mean[match(!!sym(overlap_cols[i]), cnv_data$loc_id)])
    }
    
    # Calculate adjusted segment mean 
    overlaps <- overlaps %>%
      mutate(adjusted_seg_mean = rowMeans(select(overlaps, starts_with("weighted_seg_mean"))))
    
    # Remove unnecessary columns
    overlaps <- subset(overlaps, select = c("gene_id", "adjusted_seg_mean"))
    
    # Rename columns to be consistent
    colnames(overlaps) <- c("gene_id", "seg_mean")
    
    return(overlaps)
    
  } else {
    return(NULL)
  }
  
}


names(dfs) <- unique(as.character(overlap_df$loc_id_count))

for (i in names(dfs)) {
  dfs[[i]] <- subset(overlap_df, loc_id_count == as.integer(i))
}

results <- list()

for (i in seq_along(dfs)) {
  print(paste0("Progress: ", i, "/", length(dfs)))
  df <- dfs[i]
  results[[i]] <- perform_weighting(df, names(dfs[i]))
}

final_df <- do.call(rbind, results)

return(final_df)

# apparently how you get copy number scores from segment means
copy_number <- c()
for (i in final_df$seg_mean) {
  x <- 2*(2^i)
  copy_number <- append(copy_number, x)
}










## unweighted seg mean function
unweighted_seg_mean <- function(cnv_data) {
  
  regions <- paste0(cnv_data$Chromosome, ":", cnv_data$Start, "-", cnv_data$End)
  
  annotated_regions <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                             filters = "chromosomal_region",
                             values = regions,
                             mart = ensembl)
  
  # delete rows with no gene names
  annotated_regions <- annotated_regions[annotated_regions$hgnc_symbol !="", ]
  # calculate gene sizes
  annotated_regions <- annotated_regions %>%
    mutate(gene_size = end_position - start_position)
  
  # sort by gene size and keep only the largest gene for each gene name
  annotated_regions <- annotated_regions %>%
    arrange(desc(gene_size)) %>%
    group_by(hgnc_symbol) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  # remove the gene size column
  annotated_regions <- annotated_regions %>%
    select(-gene_size)
  
  # create loc_id
  grouped_data <- group_by(cnv_data, Chromosome)
  cnv_data <- mutate(grouped_data, segment_number = row_number())
  cnv_data <- mutate(cnv_data, loc_id = paste(Chromosome, segment_number, sep = "."))
  cnv_data <- cnv_data[ ,-7]
  
  
  
  # create GRanges object for cnv_data
  cnv_gr <- GRanges(
    seqnames = cnv_data$Chromosome,
    IRanges(start = cnv_data$Start, end = cnv_data$End),
    loc_id = cnv_data$loc_id,
    seg_mean = cnv_data$Segment_Mean
  )
  
  # create GRanges object for annotated data
  annot_gr <- GRanges(
    seqnames = annotated_regions$chromosome_name,
    IRanges(start = annotated_regions$start_position, end = annotated_regions$end_position),
    gene_id = annotated_regions$hgnc_symbol
  )
  
  
  # find overlapping genes
  overlaps <- findOverlaps(cnv_gr, annot_gr)
  overlap_genes <- as.character(annot_gr$gene_id[subjectHits(overlaps)])
  overlap_df <- annotated_regions[annotated_regions$hgnc_symbol %in% overlap_genes, ]
  
  # create data frame with overlapping loc_ids for each gene
  overlap_df <- data.frame(gene_id = overlap_genes, loc_id = cnv_gr[queryHits(overlaps)]$loc_id)
  overlap_df <- aggregate(loc_id ~ gene_id, overlap_df, paste, collapse = ";")
  
  
  # create frequency table of the number of loc_id's per gene
  num_loc_ids <- table(sapply(strsplit(as.character(overlap_df$loc_id), ";"), length))
  print(num_loc_ids)
  
  
  overlap_df$loc_id_count <- sapply(strsplit(overlap_df$loc_id, ";"), length)
  #overlap_df <- subset(overlap_df, loc_id_count != 1) # removes genes that only overlap 1 seg
  overlap_df <- subset(overlap_df, select = c("gene_id", "loc_id"))
  
  # remove all loc_ids except first
  overlap_df$loc_id <- gsub(";.*", "", overlap_df$loc_id)
  
  overlaps <- merge(overlap_df, cnv_data[c("loc_id", "Segment_Mean")], by = "loc_id", all.x = TRUE)
  
  overlaps <- subset(overlaps, select = c("gene_id", "Segment_Mean"))
  
  colnames(overlaps) <- c("gene_id", "seg_mean")
  
  return(overlaps)
}




