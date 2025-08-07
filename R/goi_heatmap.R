#' Genes of Interest Heatmap Module
#'
#' Upload gene list and visualize as heatmap using user-selected metadata columns.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param filtered_data_rv reactiveValues with norm_counts, samples, species
#' @import ComplexHeatmap
#' @importFrom shiny req fileInput renderPlot downloadHandler selectInput checkboxInput uiOutput
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices pdf dev.off
#' @importFrom utils read.table
#' @export
mod_goi_heatmap <- function(input, output, session, filtered_data_rv) {
  goi_heatmap_data <- reactive({
    req(input$goi_file, filtered_data_rv()$norm_counts, input$goi_metadata_columns)
    
    gene_file <- input$goi_file$datapath
    gene_list <- tryCatch({
      ext <- tools::file_ext(input$goi_file$name)
      if (tolower(ext) == "csv") {
        read.csv(gene_file, header = TRUE)[[1]]
      } else if (tolower(ext) == "txt") {
        read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
      } else stop("Unsupported file format.")
    }, error = function(e) {
      showNotification(paste("Failed to read gene list:", e$message), type = "error")
      return(NULL)
    })
    
    selected_genes <- trimws(gene_list)
    available_genes <- rownames(filtered_data_rv()$norm_counts)
    species <- filtered_data_rv()$species
    
    if (is_ensembl_id(available_genes)) {
      symbol_to_ens <- convert_symbol_to_ensembl(selected_genes, species)
      found_genes <- unname(symbol_to_ens)
      found_genes <- found_genes[!is.na(found_genes) & found_genes %in% available_genes]
    } else {
      found_genes <- selected_genes[selected_genes %in% available_genes]
    }
    
    if (length(found_genes) == 0) {
      showNotification("None of the uploaded genes matched the dataset.", type = "error")
      return(NULL)
    }
    
    expr <- log2(filtered_data_rv()$norm_counts[found_genes, , drop = FALSE] + 1)
    
    sample_meta <- filtered_data_rv()$samples
    selected_cols <- input$goi_metadata_columns
    ann_colors <- list()
    ann_df <- data.frame(row.names = rownames(sample_meta))
    
    for (col in selected_cols) {
      if (all(is.na(sample_meta[[col]])) || length(unique(sample_meta[[col]])) == 0) next
      col_data <- as.factor(sample_meta[[col]])
      levels_col <- levels(col_data)
      if (length(levels_col) <= 8) {
        palette <- RColorBrewer::brewer.pal(8, "Dark2")[seq_along(levels_col)]
      } else {
        palette <- colorspace::rainbow_hcl(length(levels_col))
      }
      ann_df[[col]] <- col_data
      ann_colors[[col]] <- setNames(palette, levels_col)
    }
    
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = ann_df,
      col = ann_colors,
      which = "column"
    )
    
    list(expr = expr, ha = ha)
  })
  
  output$goi_heatmap <- renderPlot({
    heatmap_data <- goi_heatmap_data()
    req(heatmap_data)
    ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
      heatmap_data$expr,
      name = "log2(norm counts)",
      top_annotation = heatmap_data$ha,
      cluster_rows = TRUE,
      cluster_columns = isTRUE(input$cluster_columns),
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
      column_title = "Genes of Interest",
      column_title_gp = grid::gpar(fontface = "bold", fontsize = 6)
    ))
  })
  
  output$download_goi_heatmap <- downloadHandler(
    filename = function() {
      ext <- tolower(tools::file_ext(input$goi_heatmap_filename))
      if (ext %in% c("pdf", "svg")) input$goi_heatmap_filename else paste0("GOI_heatmap_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      heatmap_data <- goi_heatmap_data()
      req(heatmap_data)
      ext <- tools::file_ext(file)
      if (ext == "svg") svg(file) else if (ext == "pdf") pdf(file) else stop("Unsupported file type")
      ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
        heatmap_data$expr,
        name = "log2(norm counts)",
        top_annotation = heatmap_data$ha,
        cluster_rows = TRUE,
        cluster_columns = isTRUE(input$cluster_columns),
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 4, fontface = "bold"),
        column_title = "Genes of Interest",
        column_title_gp = grid::gpar(fontface = "bold", fontsize = 6)
      ))
      grDevices::dev.off()
    }
  )
}
