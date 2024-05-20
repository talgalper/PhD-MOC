# Load necessary packages
library(shiny)
library(TCGAbiolinks)
library(edgeR)
library(dplyr)
library(tibble)
library(SummarizedExperiment)

# Define the UI
ui <- fluidPage(
  titlePanel("TCGA Data Query and Analysis for Breast Cancer"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Query TCGA Data"),
      hr(),
      checkboxInput("save_query", "Save query data", value = TRUE),
      tags$p("Save the TCGA query object as a local file to save time for future runs.",
             style = "font-size: 12px; margin-top: -10px;"),
      hr(),
      checkboxInput("search_local", "Search for local query"),
      tags$p("Search for the TCGA query object in the local directory to load it if available.",
             style = "font-size: 12px; margin-top: -10px;"),
      hr(),
      selectInput("project", "Select Project", 
                  choices = c("TCGA-BRCA")),
      selectInput("subtype", "Select Subtype",
                  choices = c("BRCA.LumA", "BRCA.LumB", "BRCA.Her2", "BRCA.Basal")),
      selectInput("access", "Select Access Type", 
                  choices = c("open", "controlled")),
      selectInput("data_category", "Select Data Category", 
                  choices = c("Transcriptome Profiling")),
      selectInput("exp_strategy", "Select Experimental Strategy", 
                  choices = c("RNA-Seq")),
      checkboxInput("download_data", "Download GDC Data"),
      tags$p("TCGA files must be in GDCdata directory. Download only when required",
             style = "font-size: 12px; margin-top: -10px;"),
      actionButton("query_tcga", "Query TCGA"),
      checkboxInput("save_output", "Save output", value = TRUE),
      hr(),
      h3("Filter Genes"),
      numericInput("cpm_value", "CPM Threshold", value = 1, min = 0, step = 0.1),
      sliderInput("gene_proportion", "Proportion of Samples", min = 0, max = 1, value = 0.5, step = 0.05),
      actionButton("filter_genes", "Filter Genes")
    ),
    
    mainPanel(
      fluidRow(
        column(12,
               tabsetPanel(
                 tabPanel("Query Summary",
                          plotOutput("library_size_hist"),
                          verbatimTextOutput("data_summary")),
                 
                 tabPanel("Filtered Genes",
                          tableOutput("filtered_genes")),
                 
                 tabPanel("About",
                          h4("Description"),
                          p("This application allows you to query TCGA data for breast cancer projects, 
                            filter genes based on CPM thresholds, and analyze the data interactively."),
                          p("Use the controls on the left to query data and filter genes. The results will be displayed
                            in the respective tabs."),
                          p("More sections and functionalities can be added as needed.")
                 )
               )
        )
      )
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # Reactive values to store data
  rv <- reactiveValues(
    subtype_unstranded = NULL,
    normal_unstranded = NULL,
    merged_df = NULL,
    filtered_genes = NULL,
    error_message = NULL,
    query_TCGA = NULL,
    query_output = NULL,
    clinical_query = NULL,
    common = NULL
  )
  
  observeEvent(input$save_query, {
    if (input$save_query) {
      updateCheckboxInput(session, "search_local", value = FALSE)
    }
  })
  
  observeEvent(input$search_local, {
    if (input$search_local) {
      updateCheckboxInput(session, "save_query", value = FALSE)
    }
  })
  
  observeEvent(input$query_tcga, {
    tryCatch({
      withProgress(message = 'Querying TCGA Data...', value = 0, {
        if (input$search_local && file.exists("app_data/TCGA_query.RData")) {
          incProgress(0.1, detail = "Loading TCGA query data")
          load("app_data/TCGA_query.RData")
          rv$query_TCGA <- query_TCGA
          rv$query_output <- query_output
          rv$clinical_query <- clinical_query
          rv$common <- common
          
          selected_barcodes <- rv$common[rv$common$Subtype_Selected == input$subtype, ]
          
          incProgress(0.3, detail = "Querying subtype")
          subtype_query <- GDCquery(project = input$project,
                                    access = input$access, 
                                    data.category = input$data_category,
                                    experimental.strategy = input$exp_strategy,
                                    barcode = selected_barcodes$cases)
          
          if (input$download_data) {
            incProgress(0.5, detail = "Downloading data: May take a few minutes...")
            GDCdownload(subtype_query)
          }
          
          incProgress(0.7, detail = "Extracting results")
          subtype_data <- GDCprepare(subtype_query, summarizedExperiment = TRUE)
          rv$subtype_unstranded <- as.data.frame(assay(subtype_data, "unstranded"))
          
        } else {
          incProgress(0.1, detail = "Initializing query")
          
          query_TCGA <- GDCquery(project = input$project,
                                 access = input$access, 
                                 data.category = input$data_category,
                                 experimental.strategy = input$exp_strategy)
          
          query_output <- getResults(query_TCGA)
          
          clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
          
          clinical_query <- clinical[complete.cases(clinical$ajcc_pathologic_stage), ]
          clinical_query <- merge(query_output, clinical_query, by.x = "cases.submitter_id", by.y = "submitter_id")
          clinical_query <- subset(clinical_query, select = c("cases", "cases.submitter_id", "ajcc_pathologic_stage", 
                                                              "tissue_or_organ_of_origin", "sample_type"))
          
          table(clinical_query$ajcc_pathologic_stage)
          
          subtypes <- PanCancerAtlas_subtypes()
          
          common <- merge(clinical_query, subtypes, by.x = "cases", by.y = "pan.samplesID")
          common <- subset(common, select = c("cases", "Subtype_Selected", "sample_type", "ajcc_pathologic_stage"))
          
          rv$query_TCGA <- query_TCGA
          rv$query_output <- query_output
          rv$clinical_query <- clinical_query
          rv$common <- common
          
          incProgress(0.2, detail = "Querying subtype")
          selected_barcodes <- rv$common[rv$common$Subtype_Selected == input$subtype, ]
          
          subtype_query <- GDCquery(project = input$project,
                                    access = input$access, 
                                    data.category = input$data_category,
                                    experimental.strategy = input$exp_strategy,
                                    barcode = selected_barcodes$cases)
          
          if (input$download_data) {
            incProgress(0.4, detail = "Downloading data: May take a few minutes...")
            GDCdownload(subtype_query)
          }
          
          incProgress(0.6, detail = "Extracting results")
          subtype_data <- GDCprepare(subtype_query, summarizedExperiment = TRUE)
          rv$subtype_unstranded <- as.data.frame(assay(subtype_data, "unstranded"))
          
          if (input$save_query) {
            if (!dir.exists("app_data")) {
              dir.create("app_data")
            }
            save(query_TCGA, query_output, clinical_query, common,
                 file = "app_data/TCGA_query.RData")
          }
        }
        
        if (input$save_output) {
          save(subtype_unstranded, file = paste0("app_data/", input$subtype, "_unstranded.RData"))
        }
        
        incProgress(1, detail = "Finished")
        Sys.sleep(3)
      })
    }, error = function(e) {
      rv$error_message <- paste("Error: ", e$message)
      showModal(modalDialog(
        title = "Error",
        rv$error_message,
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  observeEvent(input$filter_genes, {
    tryCatch({
      withProgress(message = 'Filtering Genes...', value = 0, {
        req(rv$subtype_unstranded)
        subtype_unstranded <- rv$subtype_unstranded
        
        normal_unstranded <- subtype_unstranded # Placeholder, replace with actual normal data
        
        colnames(subtype_unstranded) <- paste("subtype", 1:ncol(subtype_unstranded), sep = "_")
        colnames(normal_unstranded) <- paste("normal", 1:ncol(normal_unstranded), sep = "_")
        
        merged_df <- merge(subtype_unstranded, normal_unstranded, by = "row.names")
        merged_df <- column_to_rownames(merged_df, "Row.names")
        
        gene_id <- rownames(merged_df)
        
        incProgress(0.5, detail = "Filtering genes")
        keepTheseGenes <- (rowSums(cpm(merged_df) > input$cpm_value) >= ncol(merged_df) * input$gene_proportion)
        print(summary(keepTheseGenes))
        
        merged_df <- cbind(gene_id, merged_df)
        
        removedGenes <- merged_df$gene_id[!keepTheseGenes]
        removedGenes <- as.data.frame(removedGenes)
        colnames(removedGenes)[1] <- "gene_id"
        
        merged_df <- merged_df[keepTheseGenes, ]
        merged_df <- merged_df[, -1]
        
        rv$merged_df <- merged_df
        rv$filtered_genes <- merged_df
        
        incProgress(1, detail = "Finished")
        Sys.sleep(2)
      })
    }, error = function(e) {
      rv$error_message <- paste("Error: ", e$message)
      showModal(modalDialog(
        title = "Error",
        rv$error_message,
        easyClose = TRUE,
        footer = NULL
      ))
    })
  })
  
  # Display the filtered genes
  output$filtered_genes <- renderTable({
    req(rv$filtered_genes)
    head(rv$filtered_genes)
  })
  
  # Display the library size histogram
  output$library_size_hist <- renderPlot({
    req(rv$subtype_unstranded)
    hist(colSums(rv$subtype_unstranded), 
         main = "Library Size Histogram",
         xlab = "Library Size",
         ylab = "Frequency",
         col = "skyblue",
         border = "white")
  })
  
  # Display the data summary
  output$data_summary <- renderPrint({
    req(rv$subtype_unstranded)
    data <- rv$subtype_unstranded
    list(
      rows = nrow(data),
      columns = ncol(data),
      head = head(data)
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
