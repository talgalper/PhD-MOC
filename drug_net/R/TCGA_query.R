
## Note: Only extracts "open" source samples. Downloads may be large depending on project size.
## Data will be downloaded into a new directory called GDCdata

#projects <- getGDCprojects() # view projects, use project name from "id" column



#' Retreive transcriptomic data from TCGA
#' 
#' @param project_name The TCGA project name
#' @param subtype The cancer subtype
#' @param view_clincal Boolean value to view clinical figo stage of samples
#' @param manifest_file Opt to keep manifest file that comes with downloaded GDC data
#' @return Data frame with unstranded counts
#' @examples
#' data <- get_TCGA_RNAseq_data(TCGA-BRCA, LumB)
#' @export
get_TCGA_RNAseq_data <- function(project_name, subtype, view_clinical = FALSE, manifest_file = FALSE) {
  
  query_TCGA <- GDCquery(project = project_name,
                         access = "open", 
                         data.category = "Transcriptome Profiling")
  
  query_output <- getResults(query_TCGA) # make initial query
  
  clinical <- GDCquery_clinic(project = project_name, # get the clinical data
                              type = "clinical")
  
  # print out the stage data (optional)
  if (view_clinical != FALSE) {
    clinical <- merge(query_output, clinical, by.x = "cases.submitter_id", by.y = "submitter_id")
    clinical <- subset(clinical, select = c("cases.submitter_id", "ajcc_pathologic_stage", "tissue_or_organ_of_origin"))
    table(clinical$ajcc_pathologic_stage)
  }
  
  # downlaod and prepare data. optional subtype parameter included
  if (missing(subtype)){
    if (tolower(readline("No subtype input, get all project data? (y/n): ") == "y")) {
      GDCdownload(query_TCGA)
      
      if (manifest_file == FALSE) {
        file.remove("MANIFEST.txt")
      }
      
      data <- GDCprepare(query_TCGA, summarizedExperiment = T)
      
      unstranded <- assay(data, "unstranded")  
      rownames(unstranded) <- gsub("\\.\\d+", "", rownames(unstranded))
      unstranded <- as.data.frame(unstranded)
      
      return(unstranded)
    } else {
      cat("For subtype info: PanCancerAtlas_subtypes()")
    }
    
  } else {
    subtypes <- PanCancerAtlas_subtypes()
    if (subtype %in% subtypes$Subtype_mRNA) {
      common <- query_output[query_output$cases %in% subtypes$pan.samplesID, ]
      common <- merge(query_output, subtypes, by.x = "cases", by.y = "pan.samplesID")
      selected_barcodes <- subset(common, select = c("cases", "Subtype_mRNA", "sample_type"))
      selected_barcodes <- selected_barcodes[selected_barcodes$Subtype_mRNA == subtype, ]
      
      subtype_query <- GDCquery(project = "TCGA-BRCA",
                                access = "open",
                                data.category = "Transcriptome Profiling",
                                barcode = selected_barcodes$cases)
      
      GDCdownload(subtype_query)
      
      if (manifest_file == FALSE) {
        file.remove("MANIFEST.txt")
      }
      
      data <- GDCprepare(subtype_query, summarizedExperiment = T)
      
      unstranded <- assay(data, "unstranded")  
      rownames(unstranded) <- gsub("\\.\\d+", "", rownames(unstranded))
      unstranded <- as.data.frame(unstranded)
      
      return(unstranded)
    } else {
      cat("Subtype not found, please double check notation or select a different one.", "\n")
      cat("Example: BRCA.Her2")
    }
  }
}


#' Show subtypes for a TCGA project
#' 
#' @param project Name of TCGA project
#' @return List of subtypes associated with TCGA project
#' @examples
#' subtypes <- get_subtypes(BRCA)
#' @export
get_subtypes <- function(project) {
  subtypes <- PanCancerAtlas_subtypes()
  subtypes <- subset(subtypes, select = c("cancer.type", "Subtype_mRNA"))
  project_subtypes <- subtypes[subtypes$cancer.type == project, ]
  
  if (nrow(project_subtypes) == 0){
    cat("No subtypes found for project name: ", project)
  } else {
    cat(project, " subtypes include: ", paste(unique(project_subtypes$Subtype_mRNA), collapse = ", "))
  }
}
