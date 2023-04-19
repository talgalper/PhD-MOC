### This Rscript will annotate, weight and combine copy number segment mean data ###

# check BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# required packages
required_packages <- c("ggplot2", "biomaRt", "dplyr", "tidyr", "GenomicRanges", "gridExtra", "knitr")

# Check if the required packages are installed, if not then install them
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
    library(package, character.only = TRUE)
  } else {
    library(package, character.only = TRUE)
  }
}


args <- commandArgs(trailingOnly = TRUE)
cnv_dir <- args[1]


# load ensemble database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


## weighted seg mean function
weight_seg_mean <- function(cnv_data) {
  
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
  
  # create df of genes that do not overlap with any segment
  non_overlap_genes <- annotated_regions$hgnc_symbol[!annotated_regions$hgnc_symbol %in% overlap_genes]
  non_overlap_df <- annotated_regions[annotated_regions$hgnc_symbol %in% non_overlap_genes, ]
  
  
  
  # seperate genes based on how many segments they overlap
  overlap_df$loc_id_count <- sapply(strsplit(overlap_df$loc_id, ";"), length)
  dfs <- list()
  for (i in unique(overlap_df$loc_id_count)) {
    dfs[[i]] <- subset(overlap_df, loc_id_count == i)
  }
  
  
  
  # This code checks to see if this data frame is available within dfs[]
  if (length(dfs) >= 1 && !is.null(dfs[[1]])) {
    # assign segment means to each gene based on overlaps
    overlaps_1 <- dfs[[1]]
    # Merge no_overlaps_1 with cnv_data to assign segment_mean scores to loc_ids
    overlaps_1 <- merge(overlaps_1, cnv_data[c("loc_id", "Segment_Mean")], by = "loc_id", all.x = TRUE)
    
    # remove unnecessary columns
    overlaps_1 <- subset(overlaps_1, select = c("gene_id", "Segment_Mean"))
    
    # rename columns to that all are the same
    colnames(overlaps_1) <- c("gene_id", "seg_mean")
  }
  
  
  
  if (length(dfs) >= 2 && !is.null(dfs[[2]])) {
    
    # Select the second overlaps dataframe
    overlaps_2 <- dfs[[2]]
    
    # Split loc_id into overlap_a and overlap_b columns
    overlaps_2 <- overlaps_2 %>% separate(loc_id, c("overlap_a", "overlap_b"), sep = ";")
    
    # Add columns with the length of each overlap
    overlaps_2 <- overlaps_2 %>%
      mutate(length_a = cnv_data$End[match(overlap_a, cnv_data$loc_id)] - cnv_data$Start[match(overlap_a, cnv_data$loc_id)] + 1,
             length_b = cnv_data$End[match(overlap_b, cnv_data$loc_id)] - cnv_data$Start[match(overlap_b, cnv_data$loc_id)] + 1)
    
    # Calculate weighted segment mean for overlap_a
    overlaps_2 <- overlaps_2 %>%
      mutate(weighted_seg_mean_a = (length_a/(length_a + length_b)) * cnv_data$Segment_Mean[match(overlap_a, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_b
    overlaps_2 <- overlaps_2 %>%
      mutate(weighted_seg_mean_b = (length_b/(length_a + length_b)) * cnv_data$Segment_Mean[match(overlap_b, cnv_data$loc_id)])
    
    # Calculate adjusted segment mean 
    overlaps_2 <- overlaps_2 %>%
      mutate(adjusted_seg_mean = (weighted_seg_mean_a + weighted_seg_mean_b) / 2)
    
    # remove unnecessary columns
    overlaps_2 <- subset(overlaps_2, select = c("gene_id", "adjusted_seg_mean"))
    
    # rename columns to that all are the same
    colnames(overlaps_2) <- c("gene_id", "seg_mean")
  }
  
  
  
  if (length(dfs) >= 3 && !is.null(dfs[[3]])) {
    overlaps_3 <- dfs[[3]]
    
    # Split loc_id into overlap_a and overlap_b columns
    overlaps_3 <- overlaps_3 %>% separate(loc_id, c("overlap_a", "overlap_b", "overlap_c"), sep = ";")
    
    # Add columns with the length of each overlap
    overlaps_3 <- overlaps_3 %>%
      mutate(length_a = cnv_data$End[match(overlap_a, cnv_data$loc_id)] - cnv_data$Start[match(overlap_a, cnv_data$loc_id)] + 1,
             length_b = cnv_data$End[match(overlap_b, cnv_data$loc_id)] - cnv_data$Start[match(overlap_b, cnv_data$loc_id)] + 1,
             length_c = cnv_data$End[match(overlap_c, cnv_data$loc_id)] - cnv_data$Start[match(overlap_c, cnv_data$loc_id)] + 1)
    
    # Calculate weighted segment mean for overlap_a
    overlaps_3 <- overlaps_3 %>%
      mutate(weighted_seg_mean_a = (length_a/(length_a + length_b + length_c)) * cnv_data$Segment_Mean[match(overlap_a, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_b
    overlaps_3 <- overlaps_3 %>%
      mutate(weighted_seg_mean_b = (length_b/(length_a + length_b + length_c)) * cnv_data$Segment_Mean[match(overlap_b, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_b
    overlaps_3 <- overlaps_3 %>%
      mutate(weighted_seg_mean_c = (length_c/(length_a + length_b + length_c)) * cnv_data$Segment_Mean[match(overlap_c, cnv_data$loc_id)])
    
    
    # Calculate adjusted segment mean 
    overlaps_3 <- overlaps_3 %>%
      mutate(adjusted_seg_mean = (weighted_seg_mean_a + weighted_seg_mean_b + weighted_seg_mean_c) / 3)
    
    
    # remove unnecessary columns
    overlaps_3 <- subset(overlaps_3, select = c("gene_id", "adjusted_seg_mean"))
    
    # rename columns to that all are the same
    colnames(overlaps_3) <- c("gene_id", "seg_mean")
  }
  
  
  
  if (length(dfs) >= 4 && !is.null(dfs[[4]])) {
    overlaps_4 <- dfs[[4]]
    
    # Split loc_id into overlap_a and overlap_b columns
    overlaps_4 <- overlaps_4 %>% separate(loc_id, c("overlap_a", "overlap_b", "overlap_c", "overlap_d"), sep = ";")
    
    # Add columns with the length of each overlap
    overlaps_4 <- overlaps_4 %>%
      mutate(length_a = cnv_data$End[match(overlap_a, cnv_data$loc_id)] - cnv_data$Start[match(overlap_a, cnv_data$loc_id)] + 1,
             length_b = cnv_data$End[match(overlap_b, cnv_data$loc_id)] - cnv_data$Start[match(overlap_b, cnv_data$loc_id)] + 1,
             length_c = cnv_data$End[match(overlap_c, cnv_data$loc_id)] - cnv_data$Start[match(overlap_c, cnv_data$loc_id)] + 1,
             length_d = cnv_data$End[match(overlap_d, cnv_data$loc_id)] - cnv_data$Start[match(overlap_d, cnv_data$loc_id)] + 1)
    
    # Calculate weighted segment mean for overlap_a
    overlaps_4 <- overlaps_4 %>%
      mutate(weighted_seg_mean_a = (length_a/(length_a + length_b + length_c + length_d)) * cnv_data$Segment_Mean[match(overlap_a, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_b
    overlaps_4 <- overlaps_4 %>%
      mutate(weighted_seg_mean_b = (length_b/(length_a + length_b + length_c + length_d)) * cnv_data$Segment_Mean[match(overlap_b, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_c
    overlaps_4 <- overlaps_4 %>%
      mutate(weighted_seg_mean_c = (length_c/(length_a + length_b + length_c + length_d)) * cnv_data$Segment_Mean[match(overlap_c, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_d
    overlaps_4 <- overlaps_4 %>%
      mutate(weighted_seg_mean_d = (length_d/(length_a + length_b + length_c + length_d)) * cnv_data$Segment_Mean[match(overlap_d, cnv_data$loc_id)])
    
    # Calculate adjusted segment mean 
    overlaps_4 <- overlaps_4 %>%
      mutate(adjusted_seg_mean = (weighted_seg_mean_a + weighted_seg_mean_b + weighted_seg_mean_c + weighted_seg_mean_d) / 4)
    
    
    # remove unnecessary columns
    overlaps_4 <- subset(overlaps_4, select = c("gene_id", "adjusted_seg_mean"))
    
    # rename columns to that all are the same
    colnames(overlaps_4) <- c("gene_id", "seg_mean")
  }
  
  
  
  if (length(dfs) >= 5 && !is.null(dfs[[5]])) {
    overlaps_5 <- dfs[[5]]
    
    # Split loc_id into overlap_a and overlap_b columns
    overlaps_5 <- overlaps_5 %>% separate(loc_id, c("overlap_a", "overlap_b", "overlap_c", "overlap_d", "overlap_e"), sep = ";")
    
    # Add columns with the length of each overlap
    overlaps_5 <- overlaps_5 %>%
      mutate(length_a = cnv_data$End[match(overlap_a, cnv_data$loc_id)] - cnv_data$Start[match(overlap_a, cnv_data$loc_id)] + 1,
             length_b = cnv_data$End[match(overlap_b, cnv_data$loc_id)] - cnv_data$Start[match(overlap_b, cnv_data$loc_id)] + 1,
             length_c = cnv_data$End[match(overlap_c, cnv_data$loc_id)] - cnv_data$Start[match(overlap_c, cnv_data$loc_id)] + 1,
             length_d = cnv_data$End[match(overlap_d, cnv_data$loc_id)] - cnv_data$Start[match(overlap_d, cnv_data$loc_id)] + 1,
             length_e = cnv_data$End[match(overlap_e, cnv_data$loc_id)] - cnv_data$Start[match(overlap_e, cnv_data$loc_id)] + 1)
    
    # Calculate weighted segment mean for overlap_a
    overlaps_5 <- overlaps_5 %>%
      mutate(weighted_seg_mean_a = (length_a/(length_a + length_b + length_c + length_d + length_e)) * cnv_data$Segment_Mean[match(overlap_a, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_b
    overlaps_5 <- overlaps_5 %>%
      mutate(weighted_seg_mean_b = (length_b/(length_a + length_b + length_c + length_d + length_e)) * cnv_data$Segment_Mean[match(overlap_b, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_c
    overlaps_5 <- overlaps_5 %>%
      mutate(weighted_seg_mean_c = (length_c/(length_a + length_b + length_c + length_d + length_e)) * cnv_data$Segment_Mean[match(overlap_c, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_d
    overlaps_5 <- overlaps_5 %>%
      mutate(weighted_seg_mean_d = (length_d/(length_a + length_b + length_c + length_d + length_e)) * cnv_data$Segment_Mean[match(overlap_d, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_e
    overlaps_5 <- overlaps_5 %>%
      mutate(weighted_seg_mean_e = (length_e/(length_a + length_b + length_c + length_d + length_e)) * cnv_data$Segment_Mean[match(overlap_e, cnv_data$loc_id)])
    
    
    # Calculate adjusted segment mean 
    overlaps_5 <- overlaps_5 %>%
      mutate(adjusted_seg_mean = (weighted_seg_mean_a + weighted_seg_mean_b + weighted_seg_mean_c + weighted_seg_mean_d + weighted_seg_mean_e) / 5)
    
    
    # remove unnecessary columns
    overlaps_5 <- subset(overlaps_5, select = c("gene_id", "adjusted_seg_mean"))
    
    # rename columns to that all are the same
    colnames(overlaps_5) <- c("gene_id", "seg_mean")
  }
  
  
  
  if (length(dfs) >= 6 && !is.null(dfs[[6]])) {
    overlaps_6 <- dfs[[6]]
    
    
    # Split loc_id into overlap_a and overlap_b columns
    overlaps_6 <- overlaps_6 %>% separate(loc_id, c("overlap_a", "overlap_b", "overlap_c", "overlap_d", "overlap_e", "overlap_f"), sep = ";")
    
    # Add columns with the length of each overlap
    overlaps_6 <- overlaps_6 %>%
      mutate(length_a = cnv_data$End[match(overlap_a, cnv_data$loc_id)] - cnv_data$Start[match(overlap_a, cnv_data$loc_id)] + 1,
             length_b = cnv_data$End[match(overlap_b, cnv_data$loc_id)] - cnv_data$Start[match(overlap_b, cnv_data$loc_id)] + 1,
             length_c = cnv_data$End[match(overlap_c, cnv_data$loc_id)] - cnv_data$Start[match(overlap_c, cnv_data$loc_id)] + 1,
             length_d = cnv_data$End[match(overlap_d, cnv_data$loc_id)] - cnv_data$Start[match(overlap_d, cnv_data$loc_id)] + 1,
             length_e = cnv_data$End[match(overlap_e, cnv_data$loc_id)] - cnv_data$Start[match(overlap_e, cnv_data$loc_id)] + 1,
             length_f = cnv_data$End[match(overlap_f, cnv_data$loc_id)] - cnv_data$Start[match(overlap_f, cnv_data$loc_id)] + 1)
    
    # Calculate weighted segment mean for overlap_a
    overlaps_6 <- overlaps_6 %>%
      mutate(weighted_seg_mean_a = (length_a/(length_a + length_b + length_c + length_d + length_e + length_f)) * cnv_data$Segment_Mean[match(overlap_a, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_b
    overlaps_6 <- overlaps_6 %>%
      mutate(weighted_seg_mean_b = (length_b/(length_a + length_b + length_c + length_d + length_e + length_f)) * cnv_data$Segment_Mean[match(overlap_b, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_c
    overlaps_6 <- overlaps_6 %>%
      mutate(weighted_seg_mean_c = (length_c/(length_a + length_b + length_c + length_d + length_e + length_f)) * cnv_data$Segment_Mean[match(overlap_c, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_d
    overlaps_6 <- overlaps_6 %>%
      mutate(weighted_seg_mean_d = (length_d/(length_a + length_b + length_c + length_d + length_e + length_f)) * cnv_data$Segment_Mean[match(overlap_d, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_e
    overlaps_6 <- overlaps_6 %>%
      mutate(weighted_seg_mean_e = (length_e/(length_a + length_b + length_c + length_d + length_e + length_f)) * cnv_data$Segment_Mean[match(overlap_e, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_f
    overlaps_6 <- overlaps_6 %>%
      mutate(weighted_seg_mean_f = (length_f/(length_a + length_b + length_c + length_d + length_e + length_f)) * cnv_data$Segment_Mean[match(overlap_f, cnv_data$loc_id)])
    
    
    # Calculate adjusted segment mean 
    overlaps_6 <- overlaps_6 %>%
      mutate(adjusted_seg_mean = (weighted_seg_mean_a + weighted_seg_mean_b + weighted_seg_mean_c + weighted_seg_mean_d + weighted_seg_mean_e + weighted_seg_mean_f) / 6)
    
    
    # remove unnecessary columns
    overlaps_6 <- subset(overlaps_6, select = c("gene_id", "adjusted_seg_mean"))
    
    # rename columns to that all are the same
    colnames(overlaps_6) <- c("gene_id", "seg_mean")
  }
  
  
  
  if (length(dfs) >= 7 && !is.null(dfs[[7]])) {
    overlaps_7 <- dfs[[7]]
    
    
    # Split loc_id into overlap_a and overlap_b columns
    overlaps_7 <- overlaps_7 %>% separate(loc_id, c("overlap_a", "overlap_b", "overlap_c", "overlap_d", "overlap_e", "overlap_f", "overlap_g"), sep = ";")
    
    # Add columns with the length of each overlap
    overlaps_7 <- overlaps_7 %>%
      mutate(length_a = cnv_data$End[match(overlap_a, cnv_data$loc_id)] - cnv_data$Start[match(overlap_a, cnv_data$loc_id)] + 1,
             length_b = cnv_data$End[match(overlap_b, cnv_data$loc_id)] - cnv_data$Start[match(overlap_b, cnv_data$loc_id)] + 1,
             length_c = cnv_data$End[match(overlap_c, cnv_data$loc_id)] - cnv_data$Start[match(overlap_c, cnv_data$loc_id)] + 1,
             length_d = cnv_data$End[match(overlap_d, cnv_data$loc_id)] - cnv_data$Start[match(overlap_d, cnv_data$loc_id)] + 1,
             length_e = cnv_data$End[match(overlap_e, cnv_data$loc_id)] - cnv_data$Start[match(overlap_e, cnv_data$loc_id)] + 1,
             length_f = cnv_data$End[match(overlap_f, cnv_data$loc_id)] - cnv_data$Start[match(overlap_f, cnv_data$loc_id)] + 1,
             length_g = cnv_data$End[match(overlap_g, cnv_data$loc_id)] - cnv_data$Start[match(overlap_g, cnv_data$loc_id)] + 1)
    
    # Calculate weighted segment mean for overlap_a
    overlaps_7 <- overlaps_7 %>%
      mutate(weighted_seg_mean_a = (length_a/(length_a + length_b + length_c + length_d + length_e + length_f + length_g)) * cnv_data$Segment_Mean[match(overlap_a, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_b
    overlaps_7 <- overlaps_7 %>%
      mutate(weighted_seg_mean_b = (length_b/(length_a + length_b + length_c + length_d + length_e + length_f + length_g)) * cnv_data$Segment_Mean[match(overlap_b, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_c
    overlaps_7 <- overlaps_7 %>%
      mutate(weighted_seg_mean_c = (length_c/(length_a + length_b + length_c + length_d + length_e + length_f + length_g)) * cnv_data$Segment_Mean[match(overlap_c, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_d
    overlaps_7 <- overlaps_7 %>%
      mutate(weighted_seg_mean_d = (length_d/(length_a + length_b + length_c + length_d + length_e + length_f + length_g)) * cnv_data$Segment_Mean[match(overlap_d, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_e
    overlaps_7 <- overlaps_7 %>%
      mutate(weighted_seg_mean_e = (length_e/(length_a + length_b + length_c + length_d + length_e + length_f + length_g)) * cnv_data$Segment_Mean[match(overlap_e, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_f
    overlaps_7 <- overlaps_7 %>%
      mutate(weighted_seg_mean_f = (length_f/(length_a + length_b + length_c + length_d + length_e + length_f + length_g)) * cnv_data$Segment_Mean[match(overlap_f, cnv_data$loc_id)])
    
    # Calculate weighted segment mean for overlap_g
    overlaps_7 <- overlaps_7 %>%
      mutate(weighted_seg_mean_g = (length_g/(length_a + length_b + length_c + length_d + length_e + length_f + length_g)) * cnv_data$Segment_Mean[match(overlap_g, cnv_data$loc_id)])
    
    
    # Calculate adjusted segment mean 
    overlaps_7 <- overlaps_7 %>%
      mutate(adjusted_seg_mean = (weighted_seg_mean_a + weighted_seg_mean_b + weighted_seg_mean_c + weighted_seg_mean_d + weighted_seg_mean_e + weighted_seg_mean_f + weighted_seg_mean_g) / 7)
    
    
    # remove unnecessary columns
    overlaps_7 <- subset(overlaps_7, select = c("gene_id", "adjusted_seg_mean"))
    
    # rename columns to that all are the same
    colnames(overlaps_7) <- c("gene_id", "seg_mean")
  }
  
  
  # List of overlaps data frames to check and rbind
  overlaps_list <- list("overlaps_1", "overlaps_2", "overlaps_3", "overlaps_4", "overlaps_5", "overlaps_6", "overlaps_7")
  
  # Initialize empty list to store overlapping data frames
  overlaps_to_rbind <- list()
  
  # Loop through each overlaps data frame and check if it exists in the workspace
  for (i in seq_along(overlaps_list)) {
    overlaps_name <- paste0("overlaps_", i)
    if (exists(overlaps_name)) {
      overlaps_to_rbind[[length(overlaps_to_rbind) + 1]] <- get(overlaps_name)
    }
  }
  
  # rbind overlapping data frames
  overlaps_master <- do.call(rbind, overlaps_to_rbind)
  
  
  return(overlaps_master)
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


## loop fucntions over cnv files and combine them into single dfs
copy_number_files <- list.files(cnv_dir)

weighted_master_df <- data.frame(gene_id = character())
unweighted_master_df <- data.frame(gene_id = character())


# for weighted segment means
for (i in seq_along(copy_number_files)) {
  file <- copy_number_files[i]
  cnv_data <- read.table(paste0(cnv_dir, file), header=TRUE, stringsAsFactors=FALSE, sep="\t")
  print(paste0("processing file ", i, " of ", length(copy_number_files), ": ", file))
  
  # Get the weighted seg means
  weighted_seg_means <- weight_seg_mean(cnv_data)
  
  # Rename the seg_mean column to include an identifier
  seg_mean_col_name <- paste0("seg_mean_", i)
  names(weighted_seg_means)[2] <- seg_mean_col_name
  
  # Add the weighted seg means to the combined dataframe
  weighted_master_df <- merge(weighted_master_df, weighted_seg_means, by = "gene_id", all = TRUE)
}

# for unweighted segment means
for (i in seq_along(copy_number_files)) {
  file <- copy_number_files[i]
  cnv_data <- read.table(paste0(cnv_dir, file), header=TRUE, stringsAsFactors=FALSE, sep="\t")
  print(paste0("processing file ", i, " of ", length(copy_number_files), ": ", file))
  
  # Get the weighted seg means
  unweighted_seg_means <- unweighted_seg_mean(cnv_data)
  
  # Rename the seg_mean column to include an identifier
  seg_mean_col_name <- paste0("seg_mean_", i)
  names(unweighted_seg_means)[2] <- seg_mean_col_name
  
  # Add the weighted seg means to the combined dataframe
  unweighted_master_df <- merge(unweighted_master_df, unweighted_seg_means, by = "gene_id", all = TRUE)
}


weighted_master_df <- read.csv("results/weighted_seg_means_master.csv", row.names = 1)
unweighted_master_df <- read.csv("results/unweighted_seg_means_master.csv", row.names = 1)

## Average out the segment means per gene across all samples
# set gene_id to row names
rownames(weighted_master_df) <- weighted_master_df$gene_id
weighted_master_df <- weighted_master_df[,-1]

rownames(unweighted_master_df) <- unweighted_master_df$gene_id
unweighted_master_df <- unweighted_master_df[,-1]


avg_weighted_master_df <- data.frame(row.names = rownames(weighted_master_df), avg_seg_mean = apply(weighted_master_df, 1, mean))
avg_unweighted_master_df <- data.frame(row.names = rownames(unweighted_master_df), avg_seg_mean = apply(unweighted_master_df, 1, mean))


write.csv(avg_weighted_master_df, "results/avg_weighted_seg_mean.csv")
write.csv(avg_unweighted_master_df, "results/avg_unweighted_seg_mean.csv")


