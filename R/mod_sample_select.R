# === Module: mod_sample_select ===

#' Sample Selection Module
#'
#' This module enables selection of samples from uploaded or demo data.
#' Users can view the design table and choose samples for downstream analyses.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param loaded_data_rv list storing raw_counts samples normailized counts and species
#' @param filtered_data_rv a reactive data with filtered tables.
#' @param filtered_dds_rv a reactive filtered dds object
#' @param dds_rv reactive value full dds object
#' @return Updates reactive values: `filtered_data_rv` (filtered data) and `filtered_dds_rv` (filtered DESeq2 object)
#' @importFrom DT renderDT datatable
#' @importFrom shinythemes shinytheme 
#' @export
mod_sample_select <- function(input, output, session, dds_rv, loaded_data_rv, filtered_data_rv, filtered_dds_rv) {
  observe({
    req(loaded_data_rv()$samples)
    samples <- loaded_data_rv()$samples
    print("Sample select inputs being populated...")  # Confirm execution
    updateSelectInput(session, "filter_sample_type", choices = unique(samples$sample_type))
    updateSelectInput(session, "filter_source", choices = unique(samples$source))
    updateSelectInput(session, "filter_EED_status", choices = unique(samples$EED_status))
    updateSelectInput(session, "filter_data_type", choices = unique(samples$data_type))
    updateSelectInput(session, "filter_NF1_status", choices = unique(samples$NF1_status))
    updateSelectInput(session, "filter_SUZ12_status", choices = unique(samples$SUZ12_status))
    updateSelectInput(session, "filter_PRC2_status", choices = unique(samples$PRC2_status))
    updateSelectInput(session, "filter_TP53_status", choices = unique(samples$TP53_status))
    updateSelectInput(session, "sample_select", choices = rownames(samples), selected = NULL)
    tryCatch({
      cat("Loaded data reactive values:\n")
      str(loaded_data_rv())
      print(colnames(samples))
    }, error = function(e) {
      message("Debug print failed: ", e$message)
    })
  })

  

  observeEvent(input$run_filter, {
    req(loaded_data_rv(), dds_rv())
    filtering_data <- loaded_data_rv()
    
    filters <- list(
      data_type = input$filter_data_type,
      source = input$filter_source,
      EED = input$filter_EED_status,
      sample_type = input$filter_sample_type,
      SUZ12 = input$filter_SUZ12_status,
      NF1 = input$filter_NF1_status,
      PRC2 = input$filter_PRC2_status,
      TP53 = input$filter_TP53_status
    )
    
    for (f in names(filters)) {
      if (!is.null(filters[[f]]) && length(filters[[f]]) > 0) {
        filtering_data$samples <- filtering_data$samples[filtering_data$samples[[f]] %in% filters[[f]], ]
      }
    }
    
    if (!is.null(input$sample_select) && length(input$sample_select) > 0) {
      filtering_data$samples <- filtering_data$samples[rownames(filtering_data$samples) %in% input$sample_select, ]
    }
    
    filtering_data$counts <- filtering_data$counts[, rownames(filtering_data$samples)]
    filtering_data$norm_counts <- filtering_data$norm_counts[, rownames(filtering_data$samples)]
    filtering_data$species <- loaded_data_rv()$species
    
    filtered_data_rv(filtering_data)
    
    if (!is.null(dds_rv())) {
      filtered_dds <- dds_rv()[, rownames(filtering_data$samples)]
      filtered_dds_rv(filtered_dds)
    }
  })
  
  observeEvent(input$deselect_all, {
    req(loaded_data_rv(), filtered_data_rv(), dds_rv(), filtered_dds_rv())
    samples <- loaded_data_rv()$samples
    filtered_data_rv(loaded_data_rv())
    filtered_dds_rv(dds_rv())
    updateSelectInput(session, "filter_sample_type", choices = unique(samples$sample_type))
    updateSelectInput(session, "filter_source", choices = unique(samples$source))
    updateSelectInput(session, "filter_EED_status", choices = unique(samples$EED_status))
    updateSelectInput(session, "filter_data_type", choices = unique(samples$data_type))
    updateSelectInput(session, "filter_NF1_status", choices = unique(samples$NF1_status))
    updateSelectInput(session, "filter_SUZ12_status", choices = unique(samples$SUZ12_status))
    updateSelectInput(session, "filter_PRC2_status", choices = unique(samples$PRC2_status))
    updateSelectInput(session, "filter_TP53_status", choices = unique(samples$TP53_status))
    
    updateSelectInput(session, "sample_select", choices = rownames(samples), selected = NULL)
  })
  
  observeEvent(input$select_all, {
    req(loaded_data_rv(), dds_rv())
    data <- loaded_data_rv()
    filtered_data_rv(list(
      counts = data$counts,
      samples = data$samples,
      norm_counts = data$norm_counts,
      species = data$species
    ))
    filtered_dds_rv(dds_rv())
  })
  
  output$filteredDataTable <- renderDT({
    req(filtered_data_rv())
    tryCatch({
      datatable(filtered_data_rv()$samples)
    }, error = function(e) {
      showNotification(paste("Error rendering table:", e$message), type = "error")
      NULL
    })
  })
  
  output$download_sample_table <- downloadHandler(
    filename = function() "selected_sample_table.csv",
    content = function(file) {
      write.csv(filtered_data_rv()$samples, file, row.names = FALSE)
    }
  )
}