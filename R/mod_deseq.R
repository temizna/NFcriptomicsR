# === Module: mod_differential_expression ===

#' Differential Expression Module
#'
#' @description Handles DESeq2-based differential expression, renders results table, and provides a clustered heatmap of top genes
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param filtered_data_rv a reactive list with counts, samples, norm_counts, species
#' @param filtered_dds_rv a reactive valaue holding DESeqDataSet
#' @param res_reactive ReactiveVal to store DESeq2 result
#' @return None (outputs rendered)
#' @import DESeq2
#' @importFrom dplyr mutate
#' @importFrom DT renderDT datatable
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification downloadHandler renderPrint
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom graphics abline
#' @export
mod_differential_expression <- function(input, output, session, filtered_data_rv, filtered_dds_rv, res_reactive) {

  observe({
    req(filtered_data_rv())
    filtered_data<-filtered_data_rv()
    cols <- setdiff(colnames(filtered_data$samples),"batch")
    updateSelectInput(session, "metadata_column", choices = cols)
  })

  observeEvent(input$metadata_column, {
    req(filtered_data_rv())
    filtered_data<-filtered_data_rv()
    col <- input$metadata_column
    levels <- unique(as.character(filtered_data$samples[[col]]))
    updateSelectInput(session, "reference_condition", choices = levels)
    updateSelectInput(session, "test_condition", choices = levels)
  })

  observeEvent(input$run_de, {
    req(input$metadata_column, input$reference_condition, input$test_condition)
    req(filtered_dds_rv())
    # 🚨 Check for identical reference and test
    if (input$reference_condition == input$test_condition) {
      showNotification("Reference and Test conditions must be different.", type = "error")
      return()
    }
    dds <- filtered_dds_rv()
    dds[[input$metadata_column]] <- relevel(factor(dds[[input$metadata_column]]), ref = input$reference_condition)
    design(dds) <- as.formula(paste("~", input$metadata_column))

    dds <- DESeq(dds)
    res <- suppressMessages({results(dds, contrast = c(input$metadata_column, input$test_condition, input$reference_condition))})

    res$symbol <- rownames(res)

    if (is_ensembl_id(rownames(res))) {
      conv <- convert_ensembl_to_symbol(rownames(res), filtered_data$species)
      res$symbol <- conv[rownames(res)]
    }

    res_df <- as.data.frame(res)
    res_df$symbol[is.na(res_df$symbol)] <- rownames(res_df)[is.na(res_df$symbol)]
    res_df <- res_df[order(res_df$padj), ]

    res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
    res_df <- res_df[abs(res_df$log2FoldChange) >= input$lfc_threshold & res_df$padj <= input$padj_threshold, ]

    res_reactive(res_df)
    
    # Update the filtered_dds_rv with the final dds object (if you need to use it elsewhere)
    filtered_dds_rv(dds)  # This will update the filtered_dds_rv with the new DESeq2 results
    #print(resultsNames(dds)) 
    
  })

  output$deTable <- renderDT({
    req(res_reactive())
    datatable(res_reactive(), options = list(scrollX = TRUE))
  })

  generate_heatmap_plot <- function(filtered_data, top_n, res_data, metadata_column, cluster_columns) {
    # Ensure 'res_data' and 'filtered_data' are available
    top_genes <- head(res_data[order(res_data$padj), ], top_n)
    # Extract and log-transform expression data for the top genes
    expr <- filtered_data$norm_counts[rownames(top_genes), ]
    expr <- log2(expr + 1)
    # Get the metadata for grouping
    group_values <- filtered_data$samples[[metadata_column]]
    group_levels <- unique(group_values)
    
    # Define color palette for the groups
    palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(group_levels))
    
    # Create Heatmap Annotation
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = data.frame(Group = group_values),
      col = list(Group = setNames(palette, group_levels))
    )
    
    # Cluster columns based on user input or default to TRUE
    cluster_cols <- if (!is.null(cluster_columns)) cluster_columns else TRUE
    
    # Create and return the heatmap
    heatmap_plot <- ComplexHeatmap::Heatmap(expr,
                                            name = "log2(norm counts)",
                                            top_annotation = ha,
                                            cluster_rows = TRUE,
                                            cluster_columns = cluster_cols,
                                            show_column_names = FALSE,
                                            show_row_names = TRUE,
                                            row_names_gp = gpar(fontsize=4, fontface="bold"),
                                           # width = unit(8, "cm"), height = unit(8, "cm"),
                                            column_title = paste("Top", top_n, "Diff Genes"),
                                            column_title_gp = grid::gpar(fontface = "bold",fontsize = 6))
    
    return(heatmap_plot)
  }

  output$heatmapPlot <- renderPlot({
    # include input$cluster_columns so it triggers reactivity even when FALSE
    cluster <- isTRUE(input$cluster_columns)
    req(filtered_data_rv(), res_reactive())
    
    filtered_data <- filtered_data_rv()
    top_n <- input$num_genes
    metadata_column <- input$metadata_column
    
    top_genes <- head(res_reactive()[order(res_reactive()$padj), ], top_n)
    expr <- log2(filtered_data$norm_counts[rownames(top_genes), , drop = FALSE] + 1)
    
    group_values <- filtered_data$samples[[metadata_column]]
    group_levels <- unique(group_values)
    palette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(group_levels))
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = data.frame(Group = group_values),
      col = list(Group = setNames(palette, group_levels))
    )
    
    ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
      expr,
      name = "log2(norm counts)",
      top_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = cluster,
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 4, fontface = "bold"),
      column_title = paste("Top", top_n, "Diff Genes"),
      column_title_gp = grid::gpar(fontface = "bold")
    ))
  })
  
  
  
  # Download handler for heatmap plot
  output$download_heatmap <- downloadHandler(
    filename = function() { paste0("heatmap_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file)
      req(filtered_data_rv(), res_reactive())
      filtered_data <- filtered_data_rv()
      top_n <- input$num_genes
      metadata_column <- input$metadata_column
      cluster_columns <- isTRUE(input$cluster_columns)
      
      # Call the function to generate the heatmap plot
      heatmap_plot <- generate_heatmap_plot(
        filtered_data = filtered_data,
        top_n = top_n,
        res_data = res_reactive(),
        metadata_column = metadata_column,
        cluster_columns = cluster_columns
      )
      
      # Save the plot as a PDF
      dev.off()
    }
  )
  
  
  
  output$download_de_table <- downloadHandler(
    filename = function() "differential_expression_results.csv",
    content = function(file) {
      write.csv(res_reactive(), file, row.names = FALSE)
    }
  )
}
 
# === Register in server ===
# mod_differential_expression(input, output, session, filtered_data, filtered_dds_rv, res_reactive)

