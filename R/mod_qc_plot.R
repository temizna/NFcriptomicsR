# === Module: mod_qc_plots ===

#' Quality Control Plot Module
#'
#' Generates various QC plots such as PCA, sample distance heatmap, mean-variance plots,
#' and gene expression variance histograms from normalized RNA-seq count data.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param filtered_data_rv a reactive list containing counts, samples, norm_counts, and species
#' @return None. Outputs are rendered to UI.
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs theme_minimal scale_color_brewer scale_y_continuous
#' @importFrom stringr str_split
#' @importFrom DT datatable
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var cor na.omit
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification reactiveVal downloadHandler
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @export
mod_qc_plot <- function(input, output, session, filtered_data_rv) {

  # Update QC group selection dropdown based on design metadata
  observe({
    req(filtered_data_rv()$samples)
    cols <- setdiff(colnames(filtered_data_rv()$samples), "batch")
    updateSelectInput(session, "group_select_qc", choices = cols)
  })
  generate_qc_plot <- function() {
    req(input$qc_plot_type, filtered_data_rv())
    
    filtered_data <- filtered_data_rv()
    group_col <- input$group_select_qc
    
    if (!group_col %in% colnames(filtered_data$samples)) {
      showNotification("Please select a valid grouping column for QC.", type = "error")
      return(NULL)
    }
    
    group_factor <- factor(filtered_data$samples[[group_col]])
    
    # Generate Plot Based on QC Plot Type
    if (input$qc_plot_type == "PCA") {
      expr <- log2(filtered_data$norm_counts + 1)
      expr <- expr[apply(expr, 1, function(x) var(x, na.rm = TRUE) > 0), , drop = FALSE]
      pca <- prcomp(t(expr), scale. = TRUE) 
      percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
      df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = group_factor, Sample = rownames(filtered_data$samples))
      
      p <- ggplot(df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
        geom_point(size = 3) +
        xlab(paste0("PC1 (", percentVar[1], "%)")) +
        ylab(paste0("PC2 (", percentVar[2], "%)")) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold")) +
        scale_color_brewer(palette = "Set2")
      return(p)
      
    } else if (input$qc_plot_type == "Sample Distance") {
      dist_matrix <- dist(t(log2(filtered_data$norm_counts + 1)))
      mat <- as.matrix(dist_matrix)
      rownames(mat) <- colnames(filtered_data$norm_counts)
      colnames(mat) <- colnames(filtered_data$norm_counts)
      
      p <- ComplexHeatmap::Heatmap(mat, name = "Distance", show_row_names=F,show_column_names = T, column_names_gp=gpar(size=3))
      return(p)
      
    } else if (input$qc_plot_type == "Mean-Variance") {
      means <- rowMeans(filtered_data$norm_counts)
      vars <- apply(filtered_data$norm_counts, 1, var)
      df <- data.frame(Mean = means, Variance = vars)
      
      p <- ggplot(df, aes(x = Mean, y = Variance)) +
        geom_point(alpha = 0.5) +
        scale_x_log10() +
        scale_y_log10() +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold")) +
        labs(title = "Mean-Variance Plot", x = "Mean Expression", y = "Variance")
      return(p)
      
    } else if (input$qc_plot_type == "Variance Histogram") {
      gene_vars <- apply(filtered_data$norm_counts, 1, var)
      df <- data.frame(Variance = gene_vars, Group = rep(group_factor, each = length(gene_vars)))
      
      p <- ggplot(df, aes(x = Variance, fill = Group)) +
        geom_histogram(bins = 100, position = "identity", alpha = 0.6, color = "black") +
        xlim(0, quantile(gene_vars, 0.99)) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold")) +
        labs(title = "Variance Histogram by Group", x = "Variance", y = "Gene Count") +
        scale_fill_brewer(palette = "Set1")
      return(p)
    }
  }
  
  output$qcPlot <- renderPlot({
    suppressMessages({generate_qc_plot()})
  })
  output$qcPlotDescription <- renderText({
    switch(input$qc_plot_type,
           "PCA" = "Principal Component Analysis (PCA) highlights the variation in gene expression across samples.",
           "Sample Distance" = "This heatmap shows pairwise distances between samples based on log2-normalized expression. 
           It is a measure of similarity between samples.",
           "Mean-Variance" = "This plot depicts the relationship between mean expression and variance across genes.",
           "Variance Histogram" = "This histogram displays the distribution of gene expression variances by group.",
           "Select a QC plot to view its description."
    )
  })

  # Download handler for QC plot
  output$download_qc_plot <- downloadHandler(
    filename = function() paste0(input$qc_plot_filename),
    content = function(file) {
      pdf(file)  # Open a PDF device
      print(generate_qc_plot())  # Print the plot to the PDF file
      dev.off() 
    }
  )
}

# === Register in server ===
# mod_qc_plots(input, output, session, data)

