# === Module: mod_volcano_ma_plot ===

#' Volcano and MA Plot Module
#'
#' Visualizes differential expression results as a volcano plot and MA plot,
#' includes labeling of top genes and download options
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param res_reactive ReactiveVal containing DE results
#' @param filtered_data_rv (used for species info)
#' @return None. Renders plots and download handlers
#' @import ggplot2 ggrepel
#' @importFrom dplyr arrange
#' @importFrom stringr str_split
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs theme_minimal scale_color_brewer scale_y_continuous
#' @importFrom stringr str_split
#' @importFrom DT datatable
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification downloadHandler renderPrint
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @export
mod_volcano_ma_plot <- function(input, output, session, res_reactive, filtered_data_rv) {
  
  # output$volcanoPlot <- renderPlot({
  #   req(res_reactive())
  #   res <- res_reactive()
  #   top_genes <- input$volcano_gene_label
  #   label_genes <- unlist(stringr::str_split(input$volcano_select, "[\\s,]+"))
  # 
  #   if (is_ensembl_id(res$symbol)) {
  #     label_genes <- convert_ensembl_to_symbol(label_genes, filtered_data_rv()$species)
  #   }
  # 
  #   res$color <- "NS"
  #   res$color[res$log2FoldChange >= input$volcano_lfc & res$padj <= input$volcano_padj] <- "Up"
  #   res$color[res$log2FoldChange <= -input$volcano_lfc & res$padj <= input$volcano_padj] <- "Down"
  # 
  #   res$label <- NA
  #   top_by_padj <- head(arrange(res, padj), top_genes)
  #   res$label[rownames(res) %in% rownames(top_by_padj)] <- res$symbol[rownames(res) %in% rownames(top_by_padj)]
  #   res$label[res$symbol %in% label_genes] <- res$symbol[res$symbol %in% label_genes]
  # 
  #   ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  #     geom_point(alpha = 0.6) +
  #     scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey")) +
  #     geom_text_repel(aes(label = label), max.overlaps = 50) +
  #     theme_minimal() +
  #     theme(axis.title = element_text(face = "bold")) +
  #     labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(Adjusted P-value)")
  # })
  # 
  # output$maPlot <- renderPlot({
  #   req(res_reactive())
  #   res <- res_reactive()
  # 
  #   res$color <- "NS"
  #   res$color[res$log2FoldChange >= input$volcano_lfc & res$padj <= input$volcano_padj] <- "Up"
  #   res$color[res$log2FoldChange <= -input$volcano_lfc & res$padj <= input$volcano_padj] <- "Down"
  # 
  #   ggplot(res, aes(x = baseMean, y = log2FoldChange, color = color)) +
  #     geom_point(alpha = 0.5) +
  #     scale_x_log10() +
  #     scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey")) +
  #     theme_minimal() +
  #     theme(axis.title = element_text(face = "bold")) +
  #     labs(title = "MA Plot", x = "Mean Expression (log10)", y = "Log2 Fold Change")
  # })

  # Volcano Plot Function
  generate_volcano_plot <- function(res, label_genes, top_genes, volcano_lfc, volcano_padj) {
    res$color <- "NS"
    res$color[res$log2FoldChange >= volcano_lfc & res$padj <= volcano_padj] <- "Up"
    res$color[res$log2FoldChange <= -volcano_lfc & res$padj <= volcano_padj] <- "Down"
    
    res$label <- NA
    top_by_padj <- head(arrange(res, padj), top_genes)
    res$label[rownames(res) %in% rownames(top_by_padj)] <- res$symbol[rownames(res) %in% rownames(top_by_padj)]
    res$label[res$symbol %in% label_genes] <- res$symbol[res$symbol %in% label_genes]
    
    p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey")) +
      geom_text_repel(aes(label = label, color = "black"), max.overlaps = 50, fontface = "bold") +
      theme_minimal() +
      theme(axis.title = element_text(face = "bold")) +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(Adjusted P-value)")
    
    return(p)
  }
  
  # MA Plot Function
  generate_ma_plot <- function(res, volcano_lfc, volcano_padj) {
    res$color <- "NS"
    res$color[res$log2FoldChange >= volcano_lfc & res$padj <= volcano_padj] <- "Up"
    res$color[res$log2FoldChange <= -volcano_lfc & res$padj <= volcano_padj] <- "Down"
    p <- ggplot(res, aes(x = baseMean, y = log2FoldChange, color = color)) +
      geom_point(alpha = 0.5) +
      scale_x_log10() +
      scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey")) +
      theme_minimal() +
      theme(axis.title = element_text(face = "bold")) +
      labs(title = "MA Plot", x = "Mean Expression (log10)", y = "Log2 Fold Change")
    
    return(p)
  }
  
  # Render Volcano Plot
  output$volcanoPlot <- renderPlot({
    req(res_reactive())
    res <- res_reactive()
    top_genes <- input$volcano_gene_label
    label_genes <- unlist(stringr::str_split(input$volcano_select, "[\\s,]+"))
    
    if (is_ensembl_id(res$symbol)) {
      label_genes <- convert_ensembl_to_symbol(label_genes, filtered_data_rv()$species)
    }
    
    p <- suppressMessages({generate_volcano_plot(res, label_genes, top_genes, input$volcano_lfc, input$volcano_padj)})
    print(suppressMessages({p}))
  })
  
  # Render MA Plot
  output$maPlot <- renderPlot({
    req(res_reactive())
    res <- res_reactive()
    p <- suppressMessages({generate_ma_plot(res, input$volcano_lfc, input$volcano_padj)})
    print(suppressMessages({p}))
  })
  
  # Volcano Plot Download Handler
  output$download_volcano_plot <- downloadHandler(
    filename = function() { paste0("volcano_plot_", Sys.Date(), ".pdf") },
    content = function(file) {
      tryCatch({
        message("Starting volcano plot download...")
        pdf(file)
        req(res_reactive(), input$volcano_gene_label, input$volcano_lfc, input$volcano_padj, input$volcano_select)
        res <- res_reactive()
        top_genes <- input$volcano_gene_label
        label_genes <- unlist(stringr::str_split(input$volcano_select, "[\\s,]+"))
        if (is_ensembl_id(res$symbol)) {
          label_genes <- convert_ensembl_to_symbol(label_genes, filtered_data_rv()$species)
        }
        
        # Convert to numeric explicitly
        res$padj <- as.numeric(res$padj)
        res$log2FoldChange <- as.numeric(res$log2FoldChange)
        p <- generate_volcano_plot(res, label_genes, top_genes, input$volcano_lfc, input$volcano_padj)
        print(p)
        dev.off()
        message("Download complete.")
      }, error = function(e) {
        message("Download error: ", e$message)
        showNotification(paste("Download error:", e$message), type = "error")
      })
    }
  )
  
  
  # MA Plot Download Handler
  output$download_ma_plot <- downloadHandler(
    filename = function() { paste0("ma_plot_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file)
      req(res_reactive())
      res <- res_reactive()
      p <- generate_ma_plot(res, input$volcano_lfc, input$volcano_padj)
      print(p)
      dev.off()
    }
  )
  
}

# === Register in server ===
# mod_volcano_ma_plot(input, output, session, res_reactive, data)

