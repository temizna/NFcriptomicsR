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
mod_pathway_analysis <- function(input, output, session,
                                 filtered_data_rv, res_reactive,
                                 geneList_rv, kegg_pathway_results,
                                 d1_merged_rv, pathway_input_rv,
                                 pathway_result_rv) {
  
  observeEvent(input$run_pathway, {
    req(filtered_data_rv())
    filtered_data <- filtered_data_rv()
    
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
    
    pathway_result <- NULL
    
    tryCatch({
      if (input$pathway_db == "GO") {
        go_data <- GOSemSim::godata('org.Hs.eg.db', ont = "BP", computeIC = FALSE)
        pathway_result <- clusterProfiler::enrichGO(
          gene = selected_genes,
          OrgDb = orgdb,
          keyType = "ENTREZID",
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval,
          readable = TRUE
        )
      } else if (input$pathway_db == "DOSE") {
        pathway_result <- DOSE::enrichDO(
          gene = selected_genes,
          ont = "DO",
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval,
          readable = TRUE
        )
      } else if (input$pathway_db == "KEGG") {
        kegg_sp <- if (species == "Homo sapiens") "hsa" else "mmu"
        pathway_result <- clusterProfiler::enrichKEGG(
          gene = selected_genes,
          organism = kegg_sp,
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval
        )
        pathway_result <- clusterProfiler::setReadable(pathway_result, OrgDb = orgdb, keyType = "ENTREZID")
      } else {
        pathway_result <- ReactomePA::enrichPathway(
          gene = selected_genes,
          organism = get_reactome_code(species),
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval,
          readable = TRUE
        )
        pathway_result <- clusterProfiler::setReadable(pathway_result, OrgDb = orgdb, keyType = "ENTREZID")
      }
    }, error = function(e) {
      showNotification(paste("❌ Pathway enrichment failed:", e$message), type = "error")
      return()
    })
    
    if (is.null(pathway_result) || nrow(pathway_result@result) == 0) {
      showNotification("⚠️ No enriched terms found.", type = "warning")
      return()
    }
    
    result_df <- as.data.frame(pathway_result@result)
    message("Preview of pathway result:")
    print(head(result_df, 6))
    print(str(result_df))
    
    if (!all(c("ID", "geneID") %in% colnames(result_df)) ||
        nrow(result_df) < 2 ||
        anyNA(result_df$ID) ||
        anyNA(result_df$geneID)) {
      showNotification("⚠️ Too few enriched terms to compute term similarity for plots.", type = "warning")
      pathway_result_rv(pathway_result)
      return()
    }
    
    # === Compute Term Similarity (if applicable) ===
    if (inherits(pathway_result, "enrichResult")) {
      pathway_result <- tryCatch({
        if (input$pathway_db == "GO") {
          semData <- GOSemSim::godata(annoDb = orgdb, ont = "BP")
          enrichplot::pairwise_termsim(pathway_result, semData = semData)
        } else {
          enrichplot::pairwise_termsim(pathway_result)
        }
      }, error = function(e) {
        showNotification(paste("⚠️ Failed to compute term similarity:", e$message), type = "warning")
        return(pathway_result)
      })
    }
    
    if (!is.null(pathway_result@termsim)) {
      print("Term similarity matrix dimensions:")
      print(dim(pathway_result@termsim))
    }
    
    # Store the result
    pathway_result_rv(pathway_result)
  })
  
  # === Render Plots ===
  
  output$dotPlot <- renderPlot({
    req(pathway_result_rv())
    enrichplot::dotplot(pathway_result_rv()) + theme(axis.text.y = element_text(size = 6, face = "bold"))
  })
  
  output$pathheatmapPlot <- renderPlot({
    req(pathway_result_rv(), geneList_rv())
    enrichplot::heatplot(pathway_result_rv(), foldChange = geneList_rv(), showCategory = 5)
  })
  
  output$treePlot <- renderPlot({
    req(pathway_result_rv())
    enrichplot::treeplot(pathway_result_rv()) + theme(axis.text.y = element_text(size = 6, face = "bold"))
  })
  
  output$upsetPlot <- renderPlot({
    req(pathway_result_rv())
    enrichplot::upsetplot(pathway_result_rv()) + theme(axis.text.y = element_text(size = 6, face = "bold"))
  })
  
  output$emapPlot <- renderPlot({
    req(pathway_result_rv())
    enrichplot::emapplot(pathway_result_rv(), showCategory = 10)
  })
  
  output$cnetPlot <- renderPlot({
    req(pathway_result_rv())
    enrichplot::cnetplot(pathway_result_rv(), showCategory = 10)
  })
  
  output$circularPlot <- renderPlot({
    req(pathway_result_rv(), geneList_rv())
    enrichplot::cnetplot(pathway_result_rv(),
                         layout = input$circular_layout,
                         foldChange = geneList_rv(),
                         showCategory = 5,
                         circular = TRUE,
                         colorEdge = TRUE)
  })
  
  # === Render Table ===
  
  output$pathwayTable <- renderDT({
    req(pathway_result_rv())
    df <- as.data.frame(pathway_result_rv())
    validate(need(nrow(df) > 0, "No pathway results available."))
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # === Downloads ===
  
  output$download_pathway_table <- downloadHandler(
    filename = function() {
      paste0("Pathway_", input$test_condition, "_vs_", input$reference_condition, "_", input$pathway_db, "_", input$pathway_direction, "_results.csv")
    },
    content = function(file) {
      write.csv(as.data.frame(pathway_result_rv()), file, row.names = FALSE)
    }
  )
  
  output$download_dot_plot <- downloadHandler(
    filename = function() { paste0("dotplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file); print(enrichplot::dotplot(pathway_result_rv())); dev.off()
    }
  )
  
  output$download_heatmap_plot <- downloadHandler(
    filename = function() { paste0("heatplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file); print(enrichplot::heatplot(pathway_result_rv(), foldChange = geneList_rv())); dev.off()
    }
  )
  
  output$download_tree_plot <- downloadHandler(
    filename = function() { paste0("treeplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file); print(enrichplot::treeplot(pathway_result_rv())); dev.off()
    }
  )
  
  output$download_upset_plot <- downloadHandler(
    filename = function() { paste0("upsetplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file); print(enrichplot::upsetplot(pathway_result_rv())); dev.off()
    }
  )
  
  output$download_emap_plot <- downloadHandler(
    filename = function() { paste0("emapplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file); print(enrichplot::emapplot(pathway_result_rv(), showCategory = 10)); dev.off()
    }
  )
  
  output$download_cnet_plot <- downloadHandler(
    filename = function() { paste0("cnetplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file); print(enrichplot::cnetplot(pathway_result_rv())); dev.off()
    }
  )
  
  output$download_circular_plot <- downloadHandler(
    filename = function() { paste0("circularplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      pdf(file); print(enrichplot::cnetplot(pathway_result_rv(), layout = input$circular_layout, foldChange = geneList_rv(), circular = TRUE, colorEdge = TRUE)); dev.off()
    }
  )
}
