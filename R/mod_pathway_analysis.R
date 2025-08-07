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
    # === PATHWAY ANALYSIS BASED ON DB SELECTION ===
    pathway_result <- NULL
    
    tryCatch({
      if (input$pathway_db == "GO") {
        # GO requires semantic data available before computing similarity
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
        kegg_sp <- if (filtered_data$species == "Homo sapiens") "hsa" else "mmu"
        pathway_result <- clusterProfiler::enrichKEGG(
          gene = selected_genes,
          organism = kegg_sp,
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval
        )
        pathway_result <- setReadable(pathway_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        
      } else {  # Reactome
        pathway_result <- enrichPathway(
          gene = selected_genes,
          organism = get_reactome_code(filtered_data$species),
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval,
          readable = TRUE
        )
        pathway_result <- setReadable(pathway_result, OrgDb = orgdb, keyType = "ENTREZID")
      }
    }, error = function(e) {
      showNotification(paste("❌ Pathway enrichment failed:", e$message), type = "error")
      return()
    })
    
    # === CHECK AND PROCESS RESULTS ===
    if (is.null(pathway_result) || nrow(pathway_result@result) == 0) {
      showNotification("⚠️ No enriched terms found.", type = "warning")
      return()
    }
    
    result_df <- as.data.frame(pathway_result@result)
    
    # Diagnostics for dev/debug
    message("Preview of pathway result:")
    print(utils::head(result_df, 6))
    print(str(result_df))
    
    if (!all(c("ID", "geneID") %in% colnames(result_df)) ||
        nrow(result_df) < 2 ||
        anyNA(result_df$ID) ||
        anyNA(result_df$geneID)) {
      showNotification("⚠️ Too few enriched terms to compute term similarity for plots.", type = "warning")
      pathway_result_rv(pathway_result)
      return()
    }
    
    # === PAIRWISE TERM SIMILARITY ===
    # Only call if there's enough enriched terms
    if (inherits(pathway_result, "enrichResult") &&
        nrow(pathway_result@result) >= 2 &&
        all(c("ID", "geneID") %in% colnames(pathway_result@result))) {
      
      pathway_result <- tryCatch({
        # Use semantic data only for GO
        if (input$pathway_db == "GO") {
          semData <- GOSemSim::godata(OrgDb = orgdb, ont = "BP")
          enrichplot::pairwise_termsim(pathway_result, semData = semData)
        } else {
          enrichplot::pairwise_termsim(pathway_result)
        }
      }, error = function(e) {
        showNotification(paste("Failed to compute term similarity:", e$message), type = "warning")
        return(pathway_result)
      })
    }
    
    if (!is.null(pathway_result@termsim)) {
      print("Term similarity matrix dimensions:")
      print(dim(pathway_result@termsim))
    } else {
      print("No termsim slot available.")
    }
    
    
    
    # Store result
    pathway_result_rv(pathway_result)
  
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

