# === Module: mod_cross_plot ===

#' Cross Plot and Venn Diagram Module
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param filtered_dds_rv A reactive value containing the DESeq2 object.
#' @param filtered_data_rv reactive data containing counts, samples, norm_counts, and species
#' @return None. Updates Cross Plot and Venn Diagram and download handlers.
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_vline geom_hline theme_minimal theme
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr inner_join select
#' @importFrom VennDiagram draw.pairwise.venn
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var cor
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification downloadHandler renderPrint
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom gridExtra grid.arrange
#' @export
mod_cross_plot <- function(input, output, session, filtered_data_rv, filtered_dds_rv) {
  crossplot_data <-reactiveVal(NULL)
  # Ensure that the dropdowns are populated once the data is available
  observe({
    req(filtered_data_rv()$samples)
    samples<-filtered_data_rv()$samples
    cols <- setdiff(colnames(samples),"batch")
    #print(cols)
    updateSelectInput(session, "metadata_column_x", choices = cols)
    updateSelectInput(session, "metadata_column_y", choices = cols)
  })

  # Update the test and reference conditions based on the selected metadata columns
  observeEvent(input$metadata_column_x, {
    req(filtered_data_rv()$samples)
    samples<-filtered_data_rv()$samples
    col <- input$metadata_column_x
    levels <- unique(as.character(samples[[col]]))
    updateSelectInput(session, "reference_condition_x", choices = levels)
    updateSelectInput(session, "test_condition_x", choices = levels)
  })

  observeEvent(input$metadata_column_y, {
    req(filtered_data_rv()$samples)
    samples<-filtered_data_rv()$samples
    col <- input$metadata_column_y
    levels <- unique(as.character(samples[[col]]))
    updateSelectInput(session, "reference_condition_y", choices = levels)
    updateSelectInput(session, "test_condition_y", choices = levels)
  })
  observe({
    req(input$metadata_column_y, input$metadata_column_x)
    # Only show the notification once when entering the "Cross Plot" tab
    if (is.null(input$run_crossplot)) return(NULL)
    
    # showNotification(
    #   "Please make sure you have entered your desired genes (space separated) to be displayed on the plots before you hit run crossplot button.",
    #   type = "warning",  # Type of notification, can be "warning", "error", "message", "success"
    #   duration = NULL    # Duration NULL means it will stay until dismissed
    # )
  })
  # Observing when the cross plot button is clicked
  observeEvent(input$run_crossplot, {
    req(filtered_dds_rv(),filtered_data_rv(), input$metadata_column_x, input$test_condition_x, input$reference_condition_x,
        input$metadata_column_y, input$test_condition_y, input$reference_condition_y)
    showNotification("Starting Cross Plot Analysis", type = "message")
    filtered_data<-filtered_data_rv()
    # 🚨 Check for identical reference and test
    if (input$reference_condition_x == input$test_condition_x) {
      showNotification("Reference and Test conditions must be different.", type = "error")
      return()
    }
    if (input$reference_condition_y == input$test_condition_y) {
      showNotification("Reference and Test conditions must be different.", type = "error")
      return()
    }
    # Create contrasts for X and Y axes based on user inputs
    x_comp <- c(input$metadata_column_x, input$test_condition_x, input$reference_condition_x)
    y_comp <- c(input$metadata_column_y, input$test_condition_y, input$reference_condition_y)
    
    # Check if results for the X and Y comparisons already exist
    contrast_exists_x <- paste(input$metadata_column_x, input$test_condition_x, input$reference_condition_x, sep = "_") %in% resultsNames(filtered_dds_rv())
    contrast_exists_y <- paste(input$metadata_column_y, input$test_condition_y, input$reference_condition_y, sep = "_") %in% resultsNames(filtered_dds_rv())
    
    # If either contrast doesn't exist, run DESeq2 analysis for the respective contrast
    if (!contrast_exists_x) {
      # Run DESeq2 for the X contrast
      design_formula_x <- as.formula(paste("~", input$metadata_column_x))
      dds_x <- suppressMessages({DESeq2::DESeqDataSetFromMatrix(countData = filtered_data$counts, colData = filtered_data$samples, design = design_formula_x)})
      dds_x <- suppressMessages({DESeq2::DESeq(dds_x) }) # Run DESeq2 analysis for X contrast
    } else { 
      dds_x<-filtered_dds_rv()
      }
    if (!contrast_exists_y) {
      # Run DESeq2 for the Y contrast
      design_formula_y <- as.formula(paste("~", input$metadata_column_y))
      dds_y <- suppressMessages({DESeq2::DESeqDataSetFromMatrix(countData = filtered_data$counts, colData = filtered_data$samples, design = design_formula_y)})
      dds_y <- suppressMessages({DESeq2::DESeq(dds_y) }) # Run DESeq2 analysis for Y contrast
    } else {
      dds_y<-filtered_dds_rv()
    }
    
    # Extract DESeq2 results for X and Y comparisons
    res_x_data <- results(dds_x, contrast = x_comp)
    res_y_data <- results(dds_y, contrast = y_comp)
    res_x_data <- res_x_data[!is.na(res_x_data$padj) & !is.na(res_x_data$log2FoldChange), ]
    res_y_data <- res_y_data[!is.na(res_y_data$padj) & !is.na(res_y_data$log2FoldChange), ]
    # Store results in reactive values
    crossplot_data({
      df_x <- data.frame(
        gene = rownames(res_x_data),
        log2FoldChange_x = res_x_data$log2FoldChange,
        padj_x = res_x_data$padj
      )
      df_y <- data.frame(
        gene = rownames(res_y_data),
        log2FoldChange_y = res_y_data$log2FoldChange,
        padj_y = res_y_data$padj
      )
      # Merge results based on gene names, keeping all genes from both datasets
      merged_df <- merge(df_x, df_y, by = "gene", all = TRUE)
      
      # Return the combined dataframe
      merged_df
    })
    #print(head(crossplot_data()))
    df <- crossplot_data()
    # Top N significant genes by combined adjusted p-values
    df$combined_padj <- pmax(df$padj_x, df$padj_y, na.rm = TRUE)
    top_genes <- head(df[order(df$combined_padj), "gene"], input$crossplot_topgenes)
    #print(top_genes)
    # Process highlight genes from input
    highlight_genes <- input$crossplot_gene_label
    highlight_genes <- if (nchar(trimws(highlight_genes)) > 0) {
      unlist(stringr::str_split(highlight_genes, "[\\s,]+"))
    } else {
      character(0)
    }
    if (length(highlight_genes) == 1 && highlight_genes == "") highlight_genes <- character(0)
    # Combine top genes and highlighted genes
    label_genes <- union(top_genes, highlight_genes)
    #print(label_genes)
    # Label genes based on this combined list
    df$label <- ifelse(df$gene %in% label_genes, df$gene, NA)
    df<-df[order(df$padj_x), ]
    # Update crossplot data
    crossplot_data(df)  # Update the reactive value with the new data
    #print("Updated Data with Labels:")
    #print(head(df))
    
  })

  # Crossplot rendering
  generate_cross_plot <- function(df, crossplot_gene_count, test_condition_x, reference_condition_x, test_condition_y, reference_condition_y) {
    # Top N significant genes by combined adjusted p-values
    df$combined_padj <- pmax(df$padj_x, df$padj_y, na.rm = TRUE)
    
    # Define categories based on log2FoldChange for X and Y comparisons
    df$category <- "Other"
    df$category[df$log2FoldChange_x >= 1 & df$log2FoldChange_y >= 1] <- "Up-Up"
    df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y <= -1] <- "Down-Down"
    df$category[abs(df$log2FoldChange_x) >= 1 & df$log2FoldChange_y < 1 & df$log2FoldChange_y >-1] <- "Comp1-only"
    df$category[df$log2FoldChange_x < 1 & df$log2FoldChange_x > -1 & abs(df$log2FoldChange_y) >= 1] <- "Comp2-only"
    df$category[df$log2FoldChange_x >= 1 & df$log2FoldChange_y <= -1] <- "Up-Down"
    df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y >= 1] <- "Down-Up"

    
    # Limit the dataset to the selected number of genes for plotting
    df <- head(df, crossplot_gene_count)
    
    # Calculate correlation values (Pearson and Spearman)
    pearson_r <- round(cor(df$log2FoldChange_x, df$log2FoldChange_y, method = 'pearson', use = 'complete.obs'), 3)
    spearman_rho <- round(cor(df$log2FoldChange_x, df$log2FoldChange_y, method = 'spearman', use = 'complete.obs'), 3)
    
    # Create the cross plot
    p <- ggplot(df, aes(x = log2FoldChange_x, y = log2FoldChange_y, color = category)) +
      annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.1,
               label = paste0("Pearson r = ", pearson_r, "   Spearman r = ", spearman_rho),
               size = 4) +
      geom_point(alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
      geom_vline(xintercept = c(-1, 1), linetype = "solid", color = "blue") +
      geom_hline(yintercept = c(-1, 1), linetype = "solid", color = "blue") +
      ggrepel::geom_text_repel(data = subset(df, !is.na(label)), aes(label = label,color = "black"), 
                               max.overlaps = Inf, size = 5, fontface="bold", box.padding = 0.1) +
      scale_color_manual(values = c(
        "Up-Up" = "firebrick",
        "Down-Down" = "royalblue",
        "Up-Down" = "goldenrod",
        "Down-up" = "purple",
        "Comp1-only" = "darkorange",
        "Comp2-only" = "darkgreen",
        
        "Other" = "gray"
      )) +
      theme_minimal() +
      labs(title = "Cross Plot: log2FC vs log2FC",
           x = paste("log2FC:", test_condition_x, "vs", reference_condition_x),
           y = paste("log2FC:", test_condition_y, "vs", reference_condition_y),
           color = "Regulation")
    
    return(p)
  }
  
  # Render the cross plot
  output$crossPlot <- renderPlot({
    req(crossplot_data())
    df <- crossplot_data()
    #print(head(df$label))
    p <- suppressMessages({generate_cross_plot(df, input$crossplot_gene_count, 
                             input$test_condition_x, input$reference_condition_x, 
                             input$test_condition_y, input$reference_condition_y)})
    print(suppressMessages({p}))
  })
  
  # Download handler for the cross plot
  output$download_cross_plot <- downloadHandler(
    filename = function() { paste0("cross_plot_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file)
      req(crossplot_data())
      df <- crossplot_data()
      p <- generate_cross_plot(df, input$crossplot_gene_count, 
                               input$test_condition_x, input$reference_condition_x, 
                               input$test_condition_y, input$reference_condition_y)
      print(p)
      dev.off()
    }
  )
  
  
  # Download handlers for both plots
  generate_cross_pathway_plot <- function(crossplot_data, filtered_data_rv) {
    df <- crossplot_data
    species <- filtered_data_rv()$species
    orgdb <- get_orgdb(species)
    
    # Top N significant genes by combined adjusted p-values
    df$combined_padj <- pmax(df$padj_x, df$padj_y, na.rm = TRUE)
    df$category <- "Other"
    df$category[df$log2FoldChange_x >= 1 & df$log2FoldChange_y >= 1] <- "Up-Up"
    df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y <= -1] <- "Down-Down"
    df$category[df$log2FoldChange_x > 1 & df$log2FoldChange_y < 1 &  df$log2FoldChange_y > -1] <- "Comp1-only"
    df$category[df$log2FoldChange_x < 1 & df$log2FoldChange_x > -1 & df$log2FoldChange_y >= 1] <- "Comp2-only"
    df$category[df$log2FoldChange_x >= 1 & df$log2FoldChange_y <= -1] <- "Up-Down"
    df$category[df$log2FoldChange_x <= -1 & df$log2FoldChange_y >= 1] <- "Down-Up"
    
    df=df[!(df$category=="Other"),]
    # Convert gene symbols to ENTREZ IDs
    if (is_symbol(df$gene)) {
      df_ids <- suppressMessages({bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)})
    } else {
      df_ids <- suppressMessages({bitr(df$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgdb)})
    }
    
    # Merge the data with the ENTREZ IDs
    df_merged <- merge(df, df_ids, by.x = "gene", by.y = 1)
    
    # Remove duplicate rows based on ENTREZID
    df_merged <- dplyr::distinct(df_merged, ENTREZID, .keep_all = TRUE)
    
    # Filter rows where log2FoldChange_x is greater than or equal to 1
    df_merged <- df_merged[abs(df_merged$log2FoldChange_x) >= 1, ]
    
    # Select relevant columns
    df_merged2 <- df_merged[, c("ENTREZID", "log2FoldChange_x", "category")]
    
    # Remove rows with NA values in log2FoldChange_x or ENTREZID
    df_merged2 <- df_merged2 %>%
      filter(!is.na(log2FoldChange_x) & !is.na(ENTREZID))
    
    # Perform pathway analysis using compareCluster
    # Assume orgdb has already been defined with get_orgdb(species)
    method <- input$cp_pathway_db
    if (method == "enrichGO") {
      formula_res <- clusterProfiler::compareCluster(
        ENTREZID ~ category,
        data = df_merged2,
        fun = "enrichGO",
        OrgDb = orgdb,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    } else if (method == "groupGO") {
      formula_res <- clusterProfiler::compareCluster(
        ENTREZID ~ category,
        data = df_merged2,
        fun = "groupGO",
        OrgDb = orgdb,
        ont = "BP",
        keyType = "ENTREZID",
        readable = TRUE
      )
    } else if (method == "enrichKEGG") {
      kegg_sp <- if (species == "Homo sapiens") "hsa" else "mmu"
      formula_res <- clusterProfiler::compareCluster(
        ENTREZID ~ category,
        data = df_merged2,
        fun = "enrichKEGG",
        organism = kegg_sp,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
    } else if (method == "enrichPathway") {
      reactome_code <- get_reactome_code(species)
      formula_res <- clusterProfiler::compareCluster(
        ENTREZID ~ category,
        data = df_merged2,
        fun = ReactomePA::enrichPathway,
        organism = reactome_code,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    } else {
      showNotification("Unsupported enrichment method selected.", type = "error")
      return(NULL)
    }
    
    
    # Check if formula_res is NULL or empty
    if (is.null(formula_res) || nrow(formula_res) == 0) {
      showNotification("No pathway enrichment results available.", type = "error")
      return(NULL)
    } else {
      return(enrichplot::dotplot(formula_res, x = "category"))
    }
  }
  
  # Render the Cross Pathway Plot
  output$crosspathplot <- renderPlot({
    req(crossplot_data(), filtered_data_rv())
    df <- crossplot_data()
    p <- generate_cross_pathway_plot(df, filtered_data_rv())
    if (!is.null(p)) {
      print(p)
    }
  })
  
  # Download Handler for Cross Pathway Plot
  output$download_cross_pathway_plot <- downloadHandler(
    filename = function() { paste0("cross_pathway_plot_", Sys.Date(), "_",input$test_condition_x,"_",input$test_condition_y,".pdf") },
    content = function(file) {
      pdf(file)
      df <- crossplot_data()
      p <- generate_cross_pathway_plot(df, filtered_data_rv())
      if (!is.null(p)) {
        print(p)
      }
      dev.off()
    }
  )
  # Generate Cross Venn Plot Function
  generate_cross_venn_plot <- function(crossplot_data) {
    df <- crossplot_data
    
    # Set the cutoffs for log2FoldChange and adjusted p-value
    lfc_cutoff <- 1
    padj_cutoff <- 0.05
    
    # Identify upregulated and downregulated genes for X and Y comparisons
    up_x <- rownames(df[df$log2FoldChange_x > lfc_cutoff & df$padj_x < padj_cutoff, ])
    up_y <- rownames(df[df$log2FoldChange_y > lfc_cutoff & df$padj_y < padj_cutoff, ])
    down_x <- rownames(df[df$log2FoldChange_x < -lfc_cutoff & df$padj_x < padj_cutoff, ])
    down_y <- rownames(df[df$log2FoldChange_y < -lfc_cutoff & df$padj_y < padj_cutoff, ])
    
    # Create the upregulated genes Venn plot with logging turned off
    up_venn <- VennDiagram::venn.diagram(
      x = list(X_Up = up_x, Y_Up = up_y),
      filename = NULL,
      fill = c("darkorange", "darkgreen"),
      alpha = 0.5,
      main = "Upregulated Genes Venn",
      disable.logging = TRUE
    )
    
    # Create the downregulated genes Venn plot with logging turned off
    down_venn <- VennDiagram::venn.diagram(
      x = list(X_Down = down_x, Y_Down = down_y),
      filename = NULL,
      fill = c("royalblue", "purple"),
      alpha = 0.5,
      main = "Downregulated Genes Venn",
      disable.logging = TRUE  # Ensure logging is turned off
    )
    # Return both plots as a list for further manipulation
    return(list(up_venn = up_venn, down_venn = down_venn))
  }
  overlapping_genes <- reactive({
    req(crossplot_data())
    df <- crossplot_data()
    lfc_cutoff <- 1
    padj_cutoff <- 0.05
    
    up_x <- df$gene[df$log2FoldChange_x > lfc_cutoff & df$padj_x < padj_cutoff]
    up_y <- df$gene[df$log2FoldChange_y > lfc_cutoff & df$padj_y < padj_cutoff]
    down_x <- df$gene[df$log2FoldChange_x < -lfc_cutoff & df$padj_x < padj_cutoff]
    down_y <- df$gene[df$log2FoldChange_y < -lfc_cutoff & df$padj_y < padj_cutoff]
    
    list(
      up_overlap = intersect(up_x, up_y),
      down_overlap = intersect(down_x, down_y)
    )
  })
  # Render the Cross Venn Plot
  output$crossVennPlot <- renderPlot({
    req(crossplot_data())
    df <- crossplot_data()
    venn_plots <- generate_cross_venn_plot(df)
    gridExtra::grid.arrange(
      grid::grid.grabExpr(grid::grid.draw(venn_plots$up_venn)),
      grid::grid.grabExpr(grid::grid.draw(venn_plots$down_venn)),
      ncol = 2
    )
  })
  
  # Download Handler for Cross Venn Plot
  output$download_cross_venn_plot <- downloadHandler(
    filename = function() { paste0("cross_venn_plot_", Sys.Date(),"_",input$test_condition_x,"_",input$test_condition_y,".pdf") },
    content = function(file) {
      pdf(file)
      df <- crossplot_data()
      venn_plots <- generate_cross_venn_plot(df)
      gridExtra::grid.arrange(
        grid::grid.grabExpr(grid::grid.draw(venn_plots$up_venn)),
        grid::grid.grabExpr(grid::grid.draw(venn_plots$down_venn)),
        ncol = 2
      )
      dev.off()
    }
  )
  

  output$download_overlap_genes <- downloadHandler(
    filename = function() {
      paste0("cross_plot_overlap_genes_", Sys.Date(), "_", input$test_condition_x, "_", input$test_condition_y, ".csv")
    },
    content = local(function(file) {
      overlaps <- overlapping_genes()
      up <- data.frame(Gene = overlaps$up_overlap, Direction = "Upregulated")
      down <- data.frame(Gene = overlaps$down_overlap, Direction = "Downregulated")
      write.csv(rbind(up, down), file, row.names = FALSE)
    })
  )
  
}
#mod_cross_plot(input, output, session, filtered_dds_rv)

