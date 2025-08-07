# === Module: mod_cancer_gene_census ===

#' Module to visualize overlapping genes from DE analysis with Cancer Gene Census
#'
#' This module compares differentially expressed genes (DEGs) from the DE results with a
#' Cancer Gene Census dataset, visualizes the overlap using a Venn diagram, and lists
#' the overlapping genes in a table.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param res_reactive Reactive expression containing DE results (e.g., logFC, padj, etc.)
#' @importFrom shiny req renderPlot renderUI actionButton showNotification downloadHandler 
#' @importFrom VennDiagram venn.diagram
#' @importFrom utils read.csv
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer
#' @importFrom shiny validate need
#' @importFrom utils install.packages
#' @export
mod_cancer_gene_census <- function(input, output, session, res_reactive) {
  
  # Load the Cancer Gene Census data from the app's stored extdata folder
  census_data <- reactive({
    # Read the CSV file from the app's extdata folder
    census_file_path <- system.file("extdata", "Census_allWed_Apr_23_2025.csv", package = "NFcriptomicsR")
    req(file.exists(census_file_path))
    census <- read.csv(census_file_path, stringsAsFactors = FALSE)
    return(census)
  })
  
  # Check if genes are in Ensembl format or symbols and convert to symbols if needed
  observeEvent(input$run_cancer_gene_census, {
    req(res_reactive(), census_data())  # Ensure DE results and Census data are available
    
    # Extract DE results and separate into upregulated and downregulated genes
    de_results <- res_reactive()
    de_up <- rownames(de_results[de_results$log2FoldChange > 1 & de_results$padj < 0.05, ])
    de_down <- rownames(de_results[de_results$log2FoldChange < -1 & de_results$padj < 0.05, ])
    
    # If the DE genes are in Ensembl format, convert them to symbols
    species <- "Homo sapiens"  # Assuming the species is Homo sapiens, adjust if needed
    if (!is_symbol(de_up)) {
      de_up <- convert_ensembl_to_symbol(de_up, species)
    }
    if (!is_symbol(de_down)) {
      de_down <- convert_ensembl_to_symbol(de_down, species)
    }
    
    # Extract Cancer Gene Census gene symbols
    census_genes <- census_data()$Gene_Symbol
    
    # Find overlaps between DE genes and Cancer Gene Census genes
    overlap_up <- intersect(de_up, census_genes)
    overlap_down <- intersect(de_down, census_genes)
    
    # Generate the Venn diagram for overlapping DEGs with Cancer Gene Census
    output$vennPlot <- renderPlot({
      venn_diagram <- venn.diagram(
        x = list("Upregulated DE Genes" = de_up, 
                 "Downregulated DE Genes" = de_down, 
                 "Cancer Gene Census" = census_genes),
        filename = NULL,
        fill = c("#A2C2E1", "#F1A7C7", "#A9E5A9"),
        alpha = 0.5,
        label.col = "black",
        main = "Venn Diagram: DE Genes vs Cancer Gene Census",
        cat.col = c("#A2C2E1", "#F1A7C7", "#A9E5A9"),
        cat.cex = 1.5,
        cat.pos = c(0, 0, 0),  # Adjust positions of category labels if needed
        main.pos = c(0.5, 1.05),  # Move title a bit towards the center
        sub.cex = 1.5,  # Title size of the subtitle
        sub.pos = c(0.5, 1.05),  # 
        disable.logging = TRUE,
        main.cex = 2
      )
      grid::grid.draw(venn_diagram)  # Ensure the plot is drawn on the grid
    })
    
    # Prepare overlapping genes table
    overlap_genes <- union(overlap_up, overlap_down)
    overlap_table <- census_data() %>%
      filter(Gene_Symbol %in% overlap_genes)
    
    # Render the table of overlapping genes
    output$overlappingGenesTable <- renderDT({
      DT::datatable(overlap_table, options = list(pageLength = 10))
    })
  })
  
  # Download the overlapping genes table
  output$download_cancer_gene_table <- downloadHandler(
    filename = function() {
      paste("overlapping_genes_table.csv")
    },
    content = function(file) {
      req(census_data(), res_reactive())
      # Prepare the overlapping genes table
      de_results <- res_reactive()
      de_up <- rownames(de_results[de_results$log2FoldChange > 1 & de_results$padj < 0.05, ])
      de_down <- rownames(de_results[de_results$log2FoldChange < -1 & de_results$padj < 0.05, ])
      census_genes <- census_data()$Gene_Symbol
      if (!is_symbol(de_up)) {
        de_up <- convert_ensembl_to_symbol(de_up, "Homo sapiens")
      }
      if (!is_symbol(de_down)) {
        de_down <- convert_ensembl_to_symbol(de_down, "Homo sapiens")
      }
      overlap_up <- intersect(de_up, census_genes)
      overlap_down <- intersect(de_down, census_genes)
      overlap_genes <- union(overlap_up, overlap_down)
      overlap_table <- census_data() %>%
        filter(Gene_Symbol %in% overlap_genes)
      write.csv(overlap_table, file)
    }
  )
  
  # Download the Venn diagram as PDF
  output$download_cgc_venn_plot <- downloadHandler(
    filename = function() {
      paste("_",input$test_condition,"_vs_",input$reference_condition,"_venn_diagram.pdf")
    },
    content = function(file) {
      pdf(file)
      # Ensure the Venn diagram is re-rendered on the PDF
      venn_diagram <- venn.diagram(
        x = list("Upregulated DE Genes" = de_up, 
                 "Downregulated DE Genes" = de_down, 
                 "Cancer Gene Census" = census_genes),
        filename = NULL,
        fill = c("#A2C2E1", "#F1A7C7", "#A9E5A9"),
        alpha = 0.5,
        label.col = "black",
        main = "Venn Diagram: DE Genes vs Cancer Gene Census",
        cat.col = c("#A2C2E1", "#F1A7C7", "#A9E5A9"),
        cat.cex = 1.5,
        cat.pos = c(0, 0, 0),  # Adjust positions of category labels if needed
        main.pos = c(0.5, 1.05),  # Move title a bit towards the center
        sub.cex = 1.5,  # Title size of the subtitle
        sub.pos = c(0.5, 1.05),  # 
        disable.logging = TRUE,
        main.cex = 2
      )
      grid::grid.draw(venn_diagram)
      dev.off()
    }
  )
}
