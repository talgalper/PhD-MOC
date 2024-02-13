library(TCGAbiolinks)


### Retreive transcriptoimic data from TCGA ###
## project name: TCGA project name i.e. "TCGA-BRCA"
## subtype: cancer subtype input using correct TCGA notaiton i.e. "BRCA.Her2"
## get projects: boolean input to check current project details
## view_clinical: boolean input to check figo stage of selected samples

## note: only extracts "open" source samples



get_TCGA_RNAseq_data <- function(project_name, subtype, get_projects = FALSE, view_clinical = FALSE) {
  
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
    if (subtype %in% subtypes$Subtype_Selected) {
      common <- query_output[query_output$cases %in% subtypes$pan.samplesID, ]
      common <- merge(query_output, subtypes, by.x = "cases", by.y = "pan.samplesID")
      selected_barcodes <- subset(common, select = c("cases", "Subtype_Selected", "sample_type"))
      selected_barcodes <- selected_barcodes[selected_barcodes$Subtype_Selected == subtype, ]
      
      subtype_query <- GDCquery(project = "TCGA-BRCA",
                                access = "open",
                                data.category = "Transcriptome Profiling",
                                barcode = selected_barcodes$cases)
      
      GDCdownload(subtype_query)
      
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

get_subtypes <- function(project) {
  subtypes <- PanCancerAtlas_subtypes()
  subtypes <- subset(subtypes, select = c("cancer.type", "Subtype_Selected"))
  project_subtypes <- 
}



get_TCGA_RNAseq_data(project_name = "TCGA-BRCA", subtype = "BRCA.")
