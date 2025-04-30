### MOC pipeline Functions ###

plot_PCA <- function(expr_data, sample_info, output_plot_data = TRUE, circle_clust = FALSE, 
                     label_group = NULL, colour = "Classification", shape = NULL) {
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(ggalt)
  library(RColorBrewer)
  library(rlang)
  
  # Convert expression data frame into CPM-normalized + transposed matrix
  PCA_data <- cpm(as.matrix(expr_data), log = TRUE)
  PCA_data <- t(PCA_data)
  
  pca <- prcomp(PCA_data, scale. = TRUE, center = TRUE)
  pca_data <- as.data.frame(pca$x)
  
  pca_var <- pca$sdev^2
  pca_var_perc <- round(pca_var/sum(pca_var) * 100, digits = 2)
  
  # Merge sample mapping with PCA data
  pca_data <- merge(pca_data, sample_info, by.x = "row.names", by.y = 1)
  
  # Create a custom colour palette for stages
  groups <- unique(sample_info[ ,colour])
  num_colors <- length(groups)
  colours <- brewer.pal(n = num_colors, name = "Set1")
  names(colours) <- groups
  
  # Set up the base aesthetic mapping depending on whether a shape variable is provided
  base_aes <- if (!is.null(shape)) {
    aes(PC1, PC2, color = !!sym(colour), shape = !!sym(shape))
  } else {
    aes(PC1, PC2, color = !!sym(colour))
  }
  
  # Build the ggplot object based on whether we want cluster encircling or not
  if (isTRUE(circle_clust)) {
    PCA_plot <- ggplot(pca_data, base_aes) +
      geom_point(size = 4) +
      geom_encircle(aes(group = !!sym(colour)), s_shape = 0, expand = 0.05, color = "black") +
      scale_color_manual(values = colours) +
      theme_bw() +
      labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
           y = paste0('PC2: ', pca_var_perc[2], ' %')) +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18, colour = "black"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid = element_blank()
      )
  } else {
    PCA_plot <- ggplot(pca_data, base_aes) +
      geom_point(size = 4) +
      scale_color_manual(values = colours) +
      theme_bw() +
      labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
           y = paste0('PC2: ', pca_var_perc[2], ' %')) +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18, colour = "black"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid = element_blank()
      )
  }
  
  # Optionally add labels for a specific group if label_group is provided
  if (!is.null(label_group)) {
    PCA_plot <- PCA_plot + 
      geom_text_repel(
        data = subset(pca_data, Classification == label_group),
        aes(label = Row.names),
        size = 5,
        show.legend = FALSE
      )
  }
  
  # Print and optionally return the plot and plot data
  print(PCA_plot)
  
  if (output_plot_data) {
    return(list(PCA_plot = PCA_plot, plot_data = pca_data, pca_var_perc = pca_var_perc))
  }
}



# unpaired DE analysis
DE_analysis <- function(counts_matrix, sample_info) {
  
  data <- DGEList(counts = counts_matrix, group = sample_info$Classification)
  
  cat("Normalising library size \n")
  
  data <- normLibSizes(data)
  
  design <- model.matrix(~Classification, data = sample_info)
  
  cat("Estimating dispersion \n")
  
  data <- estimateDisp(data)
  
  cat("Conduct quasi-likelihood F-test \n")
  
  fit <- glmQLFit(data, design)
  qlf <- glmQLFTest(fit)
  toptags <- topTags(qlf, n = Inf)
  
  dif_exp <- decideTests(qlf, p = 0.05, adjust = "fdr", lfc = 1)
  print(summary(dif_exp))
  dif_exp_genes <- rownames(data$counts)[as.logical(dif_exp)]
  hits <- toptags$table[toptags$table$FDR < 0.1, ]
  colnames <- colnames(hits)
  hits$gene_id <- rownames(hits)
  hits <- hits[,c("gene_id", colnames)]
  dif_exp <- hits[dif_exp_genes, ]
  cat("Total DE:", nrow(dif_exp))
  
  return(list(input_data = counts_matrix,
              sample_info = sample_info,
              hits = hits,
              dif_exp = dif_exp,
              toptags = toptags,
              DGEList = data,
              fit = fit,
              qlf = qlf,
              design = design))
}


id_annot <- function(ensembl, data, col_id = 1, input_type, convert_to) {
  library(tibble)
  
  if (is.character(data)) {
    ensembl_annot <- getBM(
      attributes = c(input_type, convert_to), 
      filters = input_type, 
      values = data, 
      mart = ensembl
    )
    
    if ("description" %in% convert_to) {
      ensembl_annot$description <- gsub("\\[.*?\\]", "", ensembl_annot$description)
    }
    
    # Create a matrix of NA's with rows equal to the length of 'data' and
    missing_df <- data.frame(
      matrix(NA, nrow = length(data), ncol = length(c(input_type, convert_to))),
      stringsAsFactors = FALSE
    )
    # Set the column names to match the ones used in getBM:
    colnames(missing_df) <- c(input_type, convert_to)
    # Fill in the input identifier column with the values in 'data'
    missing_df[, input_type] <- data
    
    # Combine the returned annotations with the missing rows data frame
    data_annot <- rbind(ensembl_annot, missing_df)
    data_annot <- data_annot[!duplicated(data_annot[,input_type]), ]
    
  } else {
    if (col_id == 0) {
      data <- rownames_to_column(data)
      col_id <- 1
    }
    
    ensembl_annot <- getBM(
      attributes = c(input_type, convert_to), 
      filters = input_type, 
      values = data[, col_id], 
      mart = ensembl
    )
    
    if (convert_to %in% "description") {
      ensembl_annot$description <- gsub("\\[.*?\\]", "", ensembl_annot$description)
    }    
    
    data_annot <- merge(
      ensembl_annot, 
      data, 
      by.x = input_type, 
      by.y = colnames(data)[col_id], 
      all.y = TRUE
    )
  }
  
  return(data_annot)
}





id_annot_2 <- function(ensembl, data, col_id = 1,
                       input_type, convert_to) {
  # deps: biomaRt::getBM, tibble::rownames_to_column()
  
  strip_desc <- function(df) {
    if ("description" %in% convert_to && "description" %in% names(df)) {
      df$description <- gsub("\\[.*?\\]", "", df$description)
    }
    df
  }
  
  # helper to build the "missing" rows data.frame
  make_missing <- function(keys) {
    miss <- data.frame(
      matrix(NA, 
             nrow = length(keys),
             ncol = 1 + length(convert_to)),
      stringsAsFactors = FALSE
    )
    colnames(miss) <- c(input_type, convert_to)
    miss[[input_type]] <- keys
    miss
  }
  
  # FLAG: converting FROM uniprot_gn_id into anything
  special_in  <- (input_type == "uniprot_gn_id")
  # FLAG: converting TO uniprot_gn_id (gene→UniProt)
  special_out <- ("uniprot_gn_id" %in% convert_to && 
                    input_type != "uniprot_gn_id")
  
  # —— WORK WITH character VECTOR ——  
  if (is.character(data)) {
    keys <- data
    
    # 1) INPUT_TYPE = uniprot_gn_id → union of Swiss‐Prot + GN lookups
    if (special_in) {
      # Swiss‐Prot first
      bm1 <- biomaRt::getBM(
        attributes = c("uniprotswissprot", convert_to),
        filters    = "uniprotswissprot",
        values     = keys,
        mart       = ensembl
      )
      bm1 <- bm1[bm1$uniprotswissprot != "" & 
                   !is.na(bm1$uniprotswissprot), , drop=FALSE]
      names(bm1)[1] <- input_type
      bm1 <- strip_desc(bm1)
      
      # then other UniProt GN
      bm2 <- biomaRt::getBM(
        attributes = c("uniprot_gn_id", convert_to),
        filters    = "uniprot_gn_id",
        values     = keys,
        mart       = ensembl
      )
      bm2 <- bm2[bm2$uniprot_gn_id != "" & 
                   !is.na(bm2$uniprot_gn_id), , drop=FALSE]
      names(bm2)[1] <- input_type
      bm2 <- strip_desc(bm2)
      
      # union, preferring bm1 rows
      used <- bm1[[input_type]]
      annot <- rbind(bm1, bm2[! bm2[[input_type]] %in% used, , drop=FALSE])
    }
    
    # 2) GENE→uniprot_gn_id: collapse Swiss vs GN, plus any extra attrs
    else if (special_out) {
      other_attrs <- setdiff(convert_to, "uniprot_gn_id")
      attrs <- c(input_type,
                 "uniprotswissprot",
                 "uniprot_gn_id",
                 other_attrs)
      bm <- biomaRt::getBM(
        attributes = attrs,
        filters    = input_type,
        values     = keys,
        mart       = ensembl
      )
      bm <- strip_desc(bm)
      
      # split & collapse per key
      by_key <- split(bm, bm[[input_type]])
      out_list <- lapply(by_key, function(df_sub) {
        # pick uniprot_gn_id
        swiss <- unique(na.omit(
          ifelse(df_sub$uniprotswissprot=="", NA, df_sub$uniprotswissprot)
        ))
        if (length(swiss)>0) {
          upid <- paste(swiss, collapse=";")
        } else {
          gids <- unique(na.omit(
            ifelse(df_sub$uniprot_gn_id=="", NA, df_sub$uniprot_gn_id)
          ))
          upid <- paste(gids, collapse=";")
        }
        # build one row
        row <- setNames(
          as.list(rep(NA, 1+length(convert_to))),
          c(input_type, convert_to)
        )
        row[[input_type]]      <- df_sub[[input_type]][1]
        row[["uniprot_gn_id"]] <- upid
        
        # collapse any other attrs too
        for (attr in other_attrs) {
          vals <- unique(na.omit(df_sub[[attr]]))
          row[[attr]] <- if (length(vals)>0) 
            paste(vals, collapse=";") 
          else 
            NA
        }
        as.data.frame(row, stringsAsFactors=FALSE)
      })
      annot <- do.call(rbind, out_list)
    }
    
    # 3) DEFAULT: single-query lookup
    else {
      annot <- biomaRt::getBM(
        attributes = c(input_type, convert_to),
        filters    = input_type,
        values     = keys,
        mart       = ensembl
      )
      annot <- strip_desc(annot)
    }
    
    # 4) bind in missing rows & dedupe
    miss    <- make_missing(keys)
    all_out <- rbind(annot, miss)
    data_annot <- all_out[!duplicated(all_out[[input_type]]), , drop=FALSE]
  }
  
  # —— WORK WITH data.frame ——  
  else {
    df <- data
    if (col_id == 0) {
      df <- tibble::rownames_to_column(df)
      col_id <- 1
    }
    keys <- df[[ colnames(df)[col_id] ]]
    
    # same 3‐branch logic, but merging at the end
    if (special_in) {
      bm1 <- biomaRt::getBM(
        attributes = c("uniprotswissprot", convert_to),
        filters    = "uniprotswissprot",
        values     = keys,
        mart       = ensembl
      )
      bm1 <- bm1[bm1$uniprotswissprot!="" & 
                   !is.na(bm1$uniprotswissprot), , drop=FALSE]
      names(bm1)[1] <- input_type
      bm1 <- strip_desc(bm1)
      
      bm2 <- biomaRt::getBM(
        attributes = c("uniprot_gn_id", convert_to),
        filters    = "uniprot_gn_id",
        values     = keys,
        mart       = ensembl
      )
      bm2 <- bm2[bm2$uniprot_gn_id!="" & 
                   !is.na(bm2$uniprot_gn_id), , drop=FALSE]
      names(bm2)[1] <- input_type
      bm2 <- strip_desc(bm2)
      
      used  <- bm1[[input_type]]
      annot <- rbind(bm1, bm2[! bm2[[input_type]] %in% used, , drop=FALSE])
    }
    else if (special_out) {
      other_attrs <- setdiff(convert_to, "uniprot_gn_id")
      attrs <- c(input_type,
                 "uniprotswissprot",
                 "uniprot_gn_id",
                 other_attrs)
      bm <- biomaRt::getBM(
        attributes = attrs,
        filters    = input_type,
        values     = keys,
        mart       = ensembl
      )
      bm <- strip_desc(bm)
      
      by_key <- split(bm, bm[[input_type]])
      out_list <- lapply(by_key, function(df_sub) {
        swiss <- unique(na.omit(
          ifelse(df_sub$uniprotswissprot=="", NA, df_sub$uniprotswissprot)
        ))
        if (length(swiss)>0) {
          upid <- paste(swiss, collapse=";")
        } else {
          gids <- unique(na.omit(
            ifelse(df_sub$uniprot_gn_id=="", NA, df_sub$uniprot_gn_id)
          ))
          upid <- paste(gids, collapse=";")
        }
        row <- setNames(
          as.list(rep(NA, 1+length(convert_to))),
          c(input_type, convert_to)
        )
        row[[input_type]]      <- df_sub[[input_type]][1]
        row[["uniprot_gn_id"]] <- upid
        for (attr in other_attrs) {
          vals <- unique(na.omit(df_sub[[attr]]))
          row[[attr]] <- if (length(vals)>0)
            paste(vals, collapse=";")
          else
            NA
        }
        as.data.frame(row, stringsAsFactors=FALSE)
      })
      annot <- do.call(rbind, out_list)
    }
    else {
      annot <- biomaRt::getBM(
        attributes = c(input_type, convert_to),
        filters    = input_type,
        values     = keys,
        mart       = ensembl
      )
      annot <- strip_desc(annot)
    }
    
    # merge back to original df
    data_annot <- merge(
      annot, df,
      by.x  = input_type,
      by.y  = colnames(df)[col_id],
      all.y = TRUE,
      sort  = FALSE
    )
  }
  
  return(data_annot)
}
















