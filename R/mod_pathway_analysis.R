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
    species <- filtered_data$species
    orgdb <- get_orgdb(species)
    showNotification("Starting Pathway Analysis", type = "message")
    
    res <- isolate(res_reactive())
    direction <- input$pathway_direction
    
    res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
    d1 <- res[, c("log2FoldChange", "padj")]
    d1$gene <- rownames(res)
    
    if (is_symbol(d1$gene)) {
      d1_ids <- clusterProfiler::bitr(d1$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
    } else {
      d1_ids <- clusterProfiler::bitr(d1$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgdb)
    }
    
    d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = 1)
    d1_merged <- d1_merged[!duplicated(d1_merged$ENTREZID), ]
    d1_merged_rv(d1_merged)
    
    gene_vector <- switch(direction,
                          "Up"   = d1_merged[d1_merged$log2FoldChange >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ],
                          "Down" = d1_merged[d1_merged$log2FoldChange <= -input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ],
                          d1_merged[abs(d1_merged$log2FoldChange) >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ])
    
    gene_vector <- gene_vector[!is.na(gene_vector$log2FoldChange) & !is.na(gene_vector$ENTREZID), ]
    geneList <- gene_vector$log2FoldChange
    names(geneList) <- gene_vector$ENTREZID
    geneList <- sort(geneList[!duplicated(names(geneList))], decreasing = TRUE)
    
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
        pathway_result <- clusterProfiler::enrichGO(
          gene          = selected_genes,
          OrgDb         = orgdb,
          keyType       = "ENTREZID",
          ont           = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff  = input$padj_threshold,
          qvalueCutoff  = input$pathway.qval,
          readable      = TRUE
        )
      } else if (input$pathway_db == "DOSE") {
        pathway_result <- DOSE::enrichDO(
          gene         = selected_genes,
          ont          = "DO",
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval,
          readable     = TRUE
        )
      } else if (input$pathway_db == "KEGG") {
        org_code <- if (species == "Homo sapiens") "hsa" else "mmu"
        pathway_result <- clusterProfiler::enrichKEGG(
          gene         = selected_genes,
          organism     = org_code,
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval
        )
        pathway_result <- clusterProfiler::setReadable(pathway_result, OrgDb = orgdb, keyType = "ENTREZID")
      } else {
        pathway_result <- ReactomePA::enrichPathway(
          gene         = selected_genes,
          organism     = get_reactome_code(species),
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval,
          readable     = TRUE
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
    
    # Term similarity
    if (inherits(pathway_result, "enrichResult") &&
        all(c("ID", "geneID") %in% colnames(pathway_result@result)) &&
        nrow(pathway_result@result) >= 2) {
      pathway_result <- tryCatch({
        if (input$pathway_db == "GO") {
          semData <- GOSemSim::godata(annoDb = orgdb, ont = "BP", computeIC = FALSE)
          enrichplot::pairwise_termsim(pathway_result, semData = semData)
        } else {
          enrichplot::pairwise_termsim(pathway_result)
        }
      }, error = function(e) {
        showNotification(paste("⚠️ Term similarity error:", e$message), type = "warning")
        return(pathway_result)
      })
    }
    
    pathway_result_rv(pathway_result)
  })
  
  # === Plots ===
  output$dotPlot <- renderPlot({
    req(pathway_result_rv())
    tryCatch({
      enrichplot::dotplot(pathway_result_rv()) +
        theme(axis.text.y = element_text(size = 6, face = "bold"))
    }, error = function(e) {
      showNotification("❌ Failed to render dot plot", type = "error")
      plot.new(); text(0.5, 0.5, "Dot plot failed.")
    })
  })
  
  output$pathheatmapPlot <- renderPlot({
    req(pathway_result_rv(), geneList_rv())
    tryCatch({
      enrichplot::heatplot(pathway_result_rv(), foldChange = geneList_rv(), showCategory = 5)
    }, error = function(e) {
      showNotification("❌ Failed to render heatmap", type = "error")
      plot.new(); text(0.5, 0.5, "Heatmap failed.")
    })
  })
  
  output$treePlot <- renderPlot({
    req(pathway_result_rv())
    tryCatch({
      enrichplot::treeplot(pathway_result_rv()) +
        theme(axis.text.y = element_text(size = 6, face = "bold"))
    }, error = function(e) {
      showNotification("❌ Failed to render tree plot", type = "error")
      plot.new(); text(0.5, 0.5, "Tree plot failed.")
    })
  })
  
  output$upsetPlot <- renderPlot({
    req(pathway_result_rv())
    tryCatch({
      enrichplot::upsetplot(pathway_result_rv()) +
        theme(axis.text.y = element_text(size = 6, face = "bold"))
    }, error = function(e) {
      showNotification("❌ Failed to render upset plot", type = "error")
      plot.new(); text(0.5, 0.5, "Upset plot failed.")
    })
  })
  
  output$emapPlot <- renderPlot({
    req(pathway_result_rv())
    tryCatch({
      enrichplot::emapplot(pathway_result_rv(), showCategory = 10)
    }, error = function(e) {
      showNotification("❌ Failed to render emapplot", type = "error")
      plot.new(); text(0.5, 0.5, "Emapplot failed.")
    })
  })
  
  output$cnetPlot <- renderPlot({
    req(pathway_result_rv())
    tryCatch({
      enrichplot::cnetplot(pathway_result_rv(), showCategory = 10)
    }, error = function(e) {
      showNotification("❌ Failed to render cnetplot", type = "error")
      plot.new(); text(0.5, 0.5, "Cnetplot failed.")
    })
  })
  
  output$circularPlot <- renderPlot({
    req(pathway_result_rv(), geneList_rv())
    tryCatch({
      enrichplot::cnetplot(pathway_result_rv(),
                           layout = input$circular_layout,
                           foldChange = geneList_rv(),
                           showCategory = 5,
                           circular = TRUE,
                           colorEdge = TRUE)
    }, error = function(e) {
      showNotification("❌ Failed to render circular plot", type = "error")
      plot.new(); text(0.5, 0.5, "Circular plot failed.")
    })
  })
  
  
  # === Table ===
  output$pathwayTable <- renderDT({
    req(pathway_result_rv())
    DT::datatable(as.data.frame(pathway_result_rv()), options = list(scrollX = TRUE))
  })
  
  # === Download Handlers ===
  output$download_dot_plot <- downloadHandler(
    filename = function() { paste0("dotplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      tryCatch({
        pdf(file)
        print(enrichplot::dotplot(pathway_result_rv()))
        dev.off()
      }, error = function(e) {
        showNotification("❌ Failed to download dot plot", type = "error")
        pdf(file)
        plot.new(); text(0.5, 0.5, "Download failed.")
        dev.off()
      })
    }
  )
  
  
  output$download_heatmap_plot <- downloadHandler(
    filename = function() { paste0("heatplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      tryCatch({
        pdf(file)
        print(enrichplot::heatplot(pathway_result_rv(), foldChange = geneList_rv()))
        dev.off()
      }, error = function(e) {
        showNotification("❌ Failed to download heatmap plot", type = "error")
        pdf(file); plot.new(); text(0.5, 0.5, "Download failed."); dev.off()
      })
    }
  )
  
  output$download_tree_plot <- downloadHandler(
    filename = function() { paste0("treeplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      tryCatch({
        pdf(file)
        print(enrichplot::treeplot(pathway_result_rv()))
        dev.off()
      }, error = function(e) {
        showNotification("❌ Failed to download tree plot", type = "error")
        pdf(file); plot.new(); text(0.5, 0.5, "Download failed."); dev.off()
      })
    }
  )
  
  output$download_upset_plot <- downloadHandler(
    filename = function() { paste0("upsetplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      tryCatch({
        pdf(file)
        print(enrichplot::upsetplot(pathway_result_rv()))
        dev.off()
      }, error = function(e) {
        showNotification("❌ Failed to download upset plot", type = "error")
        pdf(file); plot.new(); text(0.5, 0.5, "Download failed."); dev.off()
      })
    }
  )
  
  output$download_emap_plot <- downloadHandler(
    filename = function() { paste0("emapplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      tryCatch({
        pdf(file)
        print(enrichplot::emapplot(pathway_result_rv(), showCategory = 10))
        dev.off()
      }, error = function(e) {
        showNotification("❌ Failed to download emapplot", type = "error")
        pdf(file); plot.new(); text(0.5, 0.5, "Download failed."); dev.off()
      })
    }
  )
  
  output$download_cnet_plot <- downloadHandler(
    filename = function() { paste0("cnetplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      tryCatch({
        pdf(file)
        print(enrichplot::cnetplot(pathway_result_rv(), showCategory = 10))
        dev.off()
      }, error = function(e) {
        showNotification("❌ Failed to download cnetplot", type = "error")
        pdf(file); plot.new(); text(0.5, 0.5, "Download failed."); dev.off()
      })
    }
  )
  
  output$download_circular_plot <- downloadHandler(
    filename = function() { paste0("circularplot_", Sys.Date(), ".pdf") },
    content = function(file) {
      tryCatch({
        pdf(file)
        print(enrichplot::cnetplot(pathway_result_rv(),
                                   layout = input$circular_layout,
                                   foldChange = geneList_rv(),
                                   showCategory = 5,
                                   circular = TRUE,
                                   colorEdge = TRUE))
        dev.off()
      }, error = function(e) {
        showNotification("❌ Failed to download circular plot", type = "error")
        pdf(file); plot.new(); text(0.5, 0.5, "Download failed."); dev.off()
      })
    }
  )
  
}
