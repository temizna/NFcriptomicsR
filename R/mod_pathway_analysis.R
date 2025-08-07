#' Pathway Analysis Module
#'
#' This module handles pathway enrichment analysis using clusterProfiler and visualizations using enrichplot and pathfindR.
#'
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @param filtered_data_rv a list containing counts, norm_counts, samples, species
#' @param geneList_rv ReactiveVal for log2FC vector
#' @param kegg_pathway_results ReactiveVal for KEGG pathway analysis result
#' @param d1_merged_rv Reactive Data frame containing gene symbols, logFC, and padj
#' @param res_reactive ReactiveVal holding DESeq2 results
#' @param pathway_input_rv Reactive pathfindR input data frame
#' @param pathway_result_rv Reactive pathway result data frame
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var cor na.omit
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom clusterProfiler bitr
#' @importFrom clusterProfiler enrichKEGG enrichGO
#' @importFrom ReactomePA enrichPathway
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @importFrom DOSE enrichDO
#' @export
mod_pathway_analysis <- function(input, output, session, filtered_data_rv, res_reactive, geneList_rv, kegg_pathway_results, d1_merged_rv, pathway_input_rv, pathway_result_rv) {

  observeEvent(input$run_pathway, {
    req(filtered_data_rv())
    filtered_data<-filtered_data_rv()
    showNotification("Starting Pathway Analysis", type = "message")

    species <- filtered_data$species
    orgdb <- get_orgdb(species)
    res <- isolate(res_reactive())
    direction <- input$pathway_direction

    res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
    d1 <- res[, c("log2FoldChange", "padj")]
    d1$gene <- rownames(res)

    if (is_symbol(d1$gene)) {
      d1_ids <- bitr(d1$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
    } else {
      d1_ids <- bitr(d1$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgdb)
    }

    d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = 1)
    d1_merged <- d1_merged[!duplicated(d1_merged$ENTREZID), ]
    d1_merged_rv(d1_merged)

    gene_vector <- switch(direction,
      "Up" = d1_merged[d1_merged$log2FoldChange >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ],
      "Down" = d1_merged[d1_merged$log2FoldChange <= -input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ],
      d1_merged[abs(d1_merged$log2FoldChange) >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ]
    )

    gene_vector <- gene_vector[!is.na(gene_vector$log2FoldChange) & !is.na(gene_vector$ENTREZID), ]
    geneList <- gene_vector$log2FoldChange
    names(geneList) <- gene_vector$ENTREZID
    geneList <- geneList[!duplicated(names(geneList))]
    geneList <- sort(geneList, decreasing = TRUE)
    max_genes <- if (!is.null(input$max_genes)) input$max_genes else 1000
    if (length(geneList) > max_genes) geneList <- head(geneList, max_genes)

    selected_genes <- names(geneList)
    geneList_rv(geneList)

    if (length(selected_genes) < 10) {
      showNotification("Too few mapped genes for pathway analysis.", type = "error")
      return()
    }
    print("Starting Pathway Analysis")
    if (input$pathway_db == "GO") {
      pathway_result <- clusterProfiler::enrichGO(gene = selected_genes, OrgDb = orgdb, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = input$padj_threshold,
                                 qvalueCutoff  = input$pathway.qval,readable  = TRUE)
    } else if (input$pathway_db == "DOSE") {
      pathway_result <- DOSE::enrichDO(
        gene = selected_genes,
        ont = "DO",
        pvalueCutoff = input$padj_threshold,
        qvalueCutoff = input$pathway.qval,
        readable = TRUE)
    } else if (input$pathway_db == "KEGG") {
      kegg_sp <- if( filtered_data$species == "Homo sapiens") "hsa" else "mmu"
      x <- clusterProfiler::enrichKEGG(gene = selected_genes, organism = kegg_sp, pvalueCutoff = input$padj_threshold,qvalueCutoff =  input$pathway.qval)
      pathway_result <- setReadable(x, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      # pathfindR_input <- data.frame(
      # Gene.symbol = d1_merged$gene[match(selected_genes, d1_merged$ENTREZID)],
      # logFC = geneList[selected_genes],
      # adj.P.Val = d1_merged$padj[match(selected_genes, d1_merged$ENTREZID)]
   # )
    # pathfindR_input <- pathfindR_input[order(pathfindR_input$adj.P.Val), ]
    # if (nrow(pathfindR_input) > 1000) pathfindR_input <- head(pathfindR_input, 1000)
    # pathway_input_rv(pathfindR_input)
    # kegg_pathway_results(pathfindR::run_pathfindR(
    #   input = pathfindR_input,
    #   gene_sets = "KEGG",
    #   output_dir = file.path(getwd(), "kegg_pathview_outputs"),
    #   plot_enrichment_chart = FALSE
    # )) 
    } else {
      pathway_result <- tryCatch({
        x <- enrichPathway(gene = selected_genes, organism = get_reactome_code(filtered_data$species), pvalueCutoff = input$padj_threshold,
                      qvalueCutoff = input$pathway.qval, readable = TRUE)
        # enrichPathway ends here
        setReadable(x, OrgDb = orgdb, keyType = "ENTREZID")
      }, error = function(e) {
        showNotification(paste("Reactome pathway analysis failed:", e$message), type = "error")
        return(NULL)
      })
    } 
    if (!is.null(pathway_result) && nrow(pathway_result@result) > 0) {
      pathway_result <- setReadable(pathway_result, OrgDb = orgdb, keyType = "ENTREZID")
    } else {
      showNotification("No enriched KEGG pathways found.", type = "warning")
    }
    result_df <- as.data.frame(pathway_result@result)
    print("Columns in result_df:")
    print(colnames(result_df))
    print("Preview of result_df:")
    print(head(result_df))
    if (is.null(result_df) ||
        !all(c("ID", "geneID") %in% colnames(result_df)) ||
        nrow(result_df) < 2 ||
        anyNA(result_df$ID) ||
        anyNA(result_df$geneID)) {
      showNotification("Too few enriched terms to calculate term similarity for plots.", type = "warning")
      return()
    } else pathway_result <- pairwise_termsim(pathway_result) 
    pathway_result_rv(pathway_result)
      # KEGG Term-Gene Heatmap
  # output$keggHeatmapPlot <- renderPlot({
  #    req(kegg_pathway_results(), geneList_rv())
  #    selected_genes <- names(geneList_rv())
  #    d1_merged <- d1_merged_rv()
  #    gene_syms <- d1_merged$gene[match(selected_genes, d1_merged$ENTREZID)]
  #    padjs <- d1_merged$padj[match(selected_genes, d1_merged$ENTREZID)]
  #    valid_idx <- !is.na(gene_syms) & !is.na(padjs) & padjs >= 0 & padjs <= 1
  # 
  #    pathfindR_input <- data.frame(
  #      Gene.symbol = gene_syms[valid_idx],
  #      logFC = geneList_rv()[selected_genes][valid_idx],
  #      adj.P.Val = padjs[valid_idx]
  #    )
  # 
  #    pathfindR::term_gene_heatmap(
  #      result_df = kegg_pathway_results(),
  #      genes_df = pathfindR_input
  #    )
  # })

  # KEGG Pathway Image Rendering
  # output$keggPathwayImage <- renderImage({
  #   req(kegg_pathway_results(), geneList_rv(),d1_merged_rv())
  #   
  #   # Get the selected genes from the reactive gene list
  #   selected_genes <- names(geneList_rv())
  #   d1_merged <- d1_merged_rv()
  #   gene_syms <- d1_merged$gene[match(selected_genes, d1_merged$ENTREZID)]
  #   padjs <- d1_merged$padj[match(selected_genes, d1_merged$ENTREZID)]
  #   
  #   # Check if any gene symbols are actually Ensembl IDs and need conversion
  #   ensembl_ids <- gene_syms[grepl("^ENS", gene_syms)]  # Check for Ensembl IDs
  #   
  #   if (length(ensembl_ids) > 0) {
  #     # Convert Ensembl IDs to gene symbols using the provided utility function
  #     gene_syms_converted <- convert_ensembl_to_symbol(ensembl_ids, species = "Homo sapiens")
  #     # Update the gene symbols with converted values
  #     gene_syms[gene_syms %in% ensembl_ids] <- gene_syms_converted
  #   }
  #   
  #   # Filter the valid genes based on padj values (ensure both gene symbol and padj are valid)
  #   valid_idx <- !is.na(gene_syms) & !is.na(padjs) & padjs >= 0 & padjs <= 1
  #   
  #   # Prepare the data frame for pathfindR input
  #   pathfindR_input <- data.frame(
  #     Gene.symbol = gene_syms[valid_idx],
  #     logFC = geneList_rv()[selected_genes][valid_idx],
  #     adj.P.Val = padjs[valid_idx]
  #   )
  #   input_processed<-pathfindR::input_processing(
  #     pathfindR_input,
  #     p_val_threshold = 0.05,
  #     pin_name_path = "Biogrid",
  #     convert2alias = TRUE
  #   )
  #   # Visualize KEGG pathway terms using pathfindR's visualize_terms function
  #   numterm=min(nrow(kegg_pathway_results()), 10)
  #   pathway_visualization <- pathfindR::visualize_terms(kegg_pathway_results()[1:numterm,], input_processed)
  #   
  #   # Save the pathway visualization as a PNG image to a folder
  #   image_output_dir <- file.path(getwd(), "kegg_pathview_outputs", "term_visualizations")
  #   if (!dir.exists(image_output_dir)) dir.create(image_output_dir, recursive = TRUE)
  #   
  #   # Save the plot to the output directory (adjust filename if necessary)
  #   output_image_path <- file.path(image_output_dir, "kegg_pathway_visualization.png")
  #   ggsave(output_image_path, plot = pathway_visualization, width = 10, height = 8)
  #   
  #   # Retrieve the image file and render it in the Shiny app
  #   image_path <- list.files(
  #     path = image_output_dir,
  #     pattern = "_pathfindR\\.png$",
  #     full.names = TRUE
  #   )
  #   
  #   req(length(image_path) > 0 && file.exists(image_path[1]))
  #   
  #   # Return the image output for rendering in Shiny
  #   list(
  #     src = normalizePath(image_path[1]),
  #     contentType = "image/png",
  #     width = "100%",
  #     alt = "KEGG Pathway Visualization",
  #     deleteFile = FALSE
  #   )
  # }, deleteFile = FALSE)
  


 output$dotPlot <- renderPlot({
    #req(pathway_result)
    enrichplot::dotplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
  })

  output$download_dot_plot <- downloadHandler(
    filename = function() { paste0("Pathway_","_",input$test_condition,"_vs_",input$reference_condition,input$pathway_db,"_",input$pathway_direction,"_dot_plot.pdf", sep="")  },
    content = function(file) {
      pdf(file)
      print(enrichplot::dotplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
      dev.off()
    }
  )
  
  output$pathheatmapPlot <- renderPlot({
    #req(pathway_result)
    enrichplot::heatplot(pathway_result, foldChange=geneList, showCategory = 5) # + theme(axis.text.y = element_text(size = 6, face = "bold"))
  })
  
  output$download_heatmap_plot <- downloadHandler(
    filename = function() { paste0("Pathway_","_",input$test_condition,"_vs_",input$reference_condition,input$pathway_db,"_",input$pathway_direction,"_heatmap_plot.pdf", sep="") },
    content = function(file) {
      pdf(file)
      print(enrichplot::heatplot(pathway_result, foldChange=geneList , showCategory = 5)) #+ theme(axis.text.y = element_text(size = 6, face = "bold")))
      dev.off()
    }
  )
  
  output$treePlot <- renderPlot({
    #req(pathway_result)
    enrichplot::treeplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
  })
  
  output$download_tree_plot <- downloadHandler(
    filename = function() { paste0("Pathway_",input$test_condition,"_vs_",input$reference_condition,"_",input$pathway_db,"_",input$pathway_direction,"_tree_plot.pdf", sep="") },
    content = function(file) {
      pdf(file)
      print(enrichplot::treeplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
      dev.off()
    }
  )
  
  output$upsetPlot <- renderPlot({
    #req(pathway_result)
    enrichplot::upsetplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
  })
  
  output$download_upset_plot <- downloadHandler(
    filename = function() { paste0("Pathway_",input$test_condition,"_vs_",input$reference_condition,"_",input$pathway_db,"_",input$pathway_direction,"_upset_plot.pdf", sep="") },
    content = function(file) {
      pdf(file)
      print(enrichplot::upsetplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
      dev.off()
    }
  )

  output$emapPlot <- renderPlot({
    #req(pathway_result)
    #pathway_result_filtered <- subset(pathway_result, padj < 0.05)
    enrichplot::emapplot(pathway_result,showCategory = 10)
  })

  output$download_emap_plot <- downloadHandler(
    filename = function() { paste0("Pathway_",input$test_condition,"_vs_",input$reference_condition,"_",input$pathway_db,"_",input$pathway_direction,"_emap_plot.pdf", sep="") },
    content = function(file) {
      req(pathway_result)
      #pathway_result_filtered <- subset(pathway_result, padj < 0.05)
      pdf(file)
      print(enrichplot::emapplot(pathway_result,showCategory = 10))
      dev.off()
    }
  )

  output$cnetPlot <- renderPlot({
    #req(pathway_result)
    #pathway_result_filtered <- subset(pathway_result, padj < 0.05)
    enrichplot::cnetplot(pathway_result,showCategory = 10)
  })

  output$download_cnet_plot <- downloadHandler(
    filename = function() { paste0("Pathway_",input$test_condition,"_vs_",input$reference_condition,"_",input$pathway_db,"_",input$pathway_direction,"_cnet_plot.pdf", sep="") },
    content = function(file) {
      pdf(file)
      print(enrichplot::cnetplot(pathway_result))
      dev.off()
    }
  )

  output$circularPlot <- renderPlot({
    req(pathway_result,geneList_rv())
    validate(need(nrow(pathway_result@result) > 0, "No enriched terms to show circular plot."))
    enrichplot::cnetplot(pathway_result, layout = input$circular_layout, foldChange=geneList_rv(),
                         showCategory = 5,circular = TRUE,colorEdge = TRUE)
  })

  output$download_circular_plot <- downloadHandler(
    filename = function() { paste0("Pathway_",input$test_condition,"_vs_",input$reference_condition,"-",input$pathway_db,"_",input$pathway_direction,"_",input$circular_layout,"_circular_plot.pdf", sep="") },
    content = function(file) {
      req(pathway_result, geneList_rv())
      pdf(file)
      print(enrichplot::cnetplot(pathway_result, layout = input$circular_layout,  foldChange=geneList_rv(),
                                 showCategory = 5,circular = TRUE, colorEdge = TRUE))
      dev.off()
    }
  )

  output$pathwayTable <- renderDT({
    req(pathway_result)
    as.data.frame(pathway_result)
  })

  output$download_pathway_table <- downloadHandler(
    filename = function() { paste0("Pathway_",input$test_condition,"_vs_",input$reference_condition,"_",input$pathway_db,"_",input$pathway_direction,"_results.csv", sep="")  },
    content = function(file) {
      write.csv(as.data.frame(pathway_result), file, row.names = FALSE)
    }
  )
  })

}
# === Register in server ===
# mod_pathway_analysis(input, output, session, filtered_data, res_reactive, kegg_pathway_results,d1_merged_rv, pathway_input_rv)

