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
                                 filtered_data_rv,
                                 res_reactive,
                                 geneList_rv,
                                 kegg_pathway_results,
                                 d1_merged_rv,
                                 pathway_result_rv) {
  
  observeEvent(input$run_pathway, {
    req(filtered_data_rv())
    filtered_data <- filtered_data_rv()
    showNotification("Starting Pathway Analysis...", type = "message")
    
    species <- filtered_data$species
    orgdb <- get_orgdb(species)
    res <- isolate(res_reactive())
    direction <- input$pathway_direction
    
    # Clean and prepare DESeq2 results
    res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
    d1 <- res[, c("log2FoldChange", "padj")]
    d1$gene <- rownames(res)
    
    # Map gene IDs
    if (is_symbol(d1$gene)) {
      d1_ids <- bitr(d1$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
    } else {
      d1_ids <- bitr(d1$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgdb)
    }
    
    d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = 1)
    d1_merged <- d1_merged[!duplicated(d1_merged$ENTREZID), ]
    d1_merged_rv(d1_merged)
    
    # Filter genes by direction
    gene_vector <- switch(
      direction,
      "Up"   = d1_merged[d1_merged$log2FoldChange >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ],
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
    
    # Run enrichment
    showNotification("Running enrichment analysis...", type = "message")
    pathway_result <- NULL
    if (input$pathway_db == "GO") {
      pathway_result <- clusterProfiler::enrichGO(
        gene = selected_genes, OrgDb = orgdb, keyType = "ENTREZID",
        ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = input$padj_threshold, qvalueCutoff = input$pathway.qval,
        readable = TRUE
      )
    } else if (input$pathway_db == "DOSE") {
      pathway_result <- DOSE::enrichDO(
        gene = selected_genes,
        pvalueCutoff = input$padj_threshold,
        qvalueCutoff = input$pathway.qval,
        readable = TRUE
      )
    } else if (input$pathway_db == "KEGG") {
      kegg_sp <- if (filtered_data$species == "Homo sapiens") "hsa" else "mmu"
      x <- clusterProfiler::enrichKEGG(
        gene = selected_genes, organism = kegg_sp,
        pvalueCutoff = input$padj_threshold, qvalueCutoff = input$pathway.qval
      )
      pathway_result <- setReadable(x, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    } else {
      pathway_result <- tryCatch({
        x <- ReactomePA::enrichPathway(
          gene = selected_genes,
          organism = get_reactome_code(filtered_data$species),
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval,
          readable = TRUE
        )
        setReadable(x, OrgDb = orgdb, keyType = "ENTREZID")
      }, error = function(e) {
        showNotification(paste("Reactome analysis failed:", e$message), type = "error")
        return(NULL)
      })
    }
    
    # Validate enrichment result
    if (is.null(pathway_result) || nrow(pathway_result@result) == 0) {
      showNotification("No enriched pathways found.", type = "warning")
      return()
    }
    
    result_df <- as.data.frame(pathway_result@result)
    print("Columns in result_df:")
    print(colnames(result_df))
    print(head(result_df))
    
    # Safe computation of pairwise similarity
    sig_terms <- result_df[result_df$p.adjust <= input$pathway.qval, ]
    if (nrow(sig_terms) >= 2) {
      pathway_result <- tryCatch({
        pairwise_termsim(pathway_result)
      }, error = function(e) {
        showNotification(paste("Term similarity computation skipped:", e$message), type = "warning")
        return(pathway_result)
      })
    } else {
      showNotification("Too few significant terms to compute term similarity.", type = "warning")
    }
    
    pathway_result_rv(pathway_result)
    
    # --- Visualization Outputs ---
    
    output$dotPlot <- renderPlot({
      req(pathway_result)
      if (is.null(pathway_result@result) || nrow(pathway_result@result) == 0) {
        validate("No enriched terms available for dot plot.")
      }
      tryCatch({
        enrichplot::dotplot(pathway_result) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Dot plot unavailable (insufficient terms).")
      })
    })
    
    output$pathheatmapPlot <- renderPlot({
      req(pathway_result)
      if (is.null(pathway_result@result) || nrow(pathway_result@result) == 0) {
        validate("No enriched terms for heatmap.")
      }
      tryCatch({
        enrichplot::heatplot(pathway_result, foldChange = geneList, showCategory = 5)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Heatmap unavailable (insufficient gene-term mapping).")
      })
    })
    
    output$treePlot <- renderPlot({
      req(pathway_result)
      if (is.null(pathway_result@result) || nrow(pathway_result@result) < 2) {
        validate("Too few enriched terms for tree plot.")
      }
      tryCatch({
        enrichplot::treeplot(pathway_result) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6, face = "bold"))
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Tree plot unavailable.")
      })
    })
    
    output$emapPlot <- renderPlot({
      req(pathway_result)
      if (is.null(pathway_result@result) || nrow(pathway_result@result) < 2) {
        validate("Too few terms to draw enrichment map.")
      }
      tryCatch({
        enrichplot::emapplot(pathway_result, showCategory = 10)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Enrichment map unavailable.")
      })
    })
    
    output$cnetPlot <- renderPlot({
      req(pathway_result)
      if (is.null(pathway_result@result) || nrow(pathway_result@result) < 2) {
        validate("Too few enriched terms for cnet plot.")
      }
      tryCatch({
        enrichplot::cnetplot(pathway_result, showCategory = 10)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Cnet plot unavailable.")
      })
    })
    
    output$circularPlot <- renderPlot({
      req(pathway_result, geneList_rv())
      if (is.null(pathway_result@result) || nrow(pathway_result@result) < 2) {
        validate("No enriched terms for circular plot.")
      }
      tryCatch({
        enrichplot::cnetplot(
          pathway_result,
          layout = input$circular_layout,
          foldChange = geneList_rv(),
          showCategory = 5,
          circular = TRUE,
          colorEdge = TRUE
        )
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Circular plot unavailable.")
      })
    })
    
    
    output$pathwayTable <- DT::renderDT({
      as.data.frame(pathway_result)
    })
    
    showNotification("Pathway analysis completed successfully!", type = "message")
  })
}
