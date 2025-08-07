# === Module: mod_gsea_analysis ===

#' GSEA Module
#'
#' Runs Gene Set Enrichment Analysis (GSEA) using clusterProfiler and msigdbr,
#' supports multiple databases (GO, KEGG, Reactome, Hallmark) and visualizes
#' the top enriched pathways with dot plots and enrichment scores.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param filtered_data_rv a list containing species and normalized expression
#' @param res_reactive ReactiveVal containing DE results
#' @import clusterProfiler msigdbr enrichplot
#' @importFrom utils head read.csv write.csv str
#' @importFrom stats as.formula dist model.matrix prcomp quantile relevel var
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom grid gpar
#' @importFrom shiny isolate req renderPlot renderImage showNotification
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom shinythemes shinytheme
#' @export
mod_gsea_analysis <- function(input, output, session, filtered_data_rv, res_reactive) {

  observeEvent(input$run_gsea, {
    req(res_reactive(), input$gsea_db, input$lfc_threshold, input$padj_threshold, input$gsea_pvalue)
    showNotification("Starting GSEA Analysis", type = "message")
    res <- res_reactive()
    res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
    res$gene <- rownames(res)

    if (!is_symbol(res$gene)) {
      conv <- convert_ensembl_to_symbol(res$gene, filtered_data_rv()$species)
      res$gene <- unname(conv[rownames(res)])
    }

    d1 <- res[, c("log2FoldChange", "padj")]
    d1$gene <- rownames(res)
    gene_vector <- d1[abs(d1$log2FoldChange) >= input$lfc_threshold & d1$padj <= input$padj_threshold, ]
    geneList <- gene_vector$log2FoldChange
    names(geneList) <- gene_vector$gene
    geneList <- sort(geneList, decreasing = TRUE)

    species_name <- if (filtered_data_rv()$species == "Mus musculus") "Mus musculus" else "Homo sapiens"
    db_collection <- switch(input$gsea_db,
      "GO" = "C5",
      "KEGG" = "C2",
      "Reactome" = "C2",
      "Hallmark" = "H",
      "Cancer Cell Atlas" = "C4",
      "Cancer Gene Neighbourhoods" ="C4",
      "Cancer Modules" = "C4",
      "Txn Factor Targets" = "C3"
    )
    db_subcollection <- switch(input$gsea_db,
      "GO" = "GO:BP",
      "KEGG" = "CP:KEGG_MEDICUS",
      "Reactome" = "CP:REACTOME",
      "Hallmark" = NULL,
      "Cancer Cell Atlas" = "3CA",
      "Cancer Gene Neighbourhoods" ="CGN",
      "Cancer Modules" = "CM",
      "Txn Factor Targets" = "TFT:GTRD"
    )

    m_df <- msigdbr(species = species_name, collection = db_collection, subcollection = db_subcollection)
    term2gene <- m_df[, c("gs_name", "gene_symbol")]
    colnames(term2gene) <- c("ID", "gene")

    gsea_result <- clusterProfiler::GSEA(
      geneList = geneList,
      TERM2GENE = term2gene,
      pvalueCutoff = input$gsea_pvalue,
      verbose = FALSE
    )

    if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) < 1) {
      showNotification("No enriched terms found in GSEA.", type = "warning")
      return()
    }
    updateSelectInput(session, "gsea_selected_pathway", choices = gsea_result@result$ID)
    output$gseaDotPlot <- renderPlot({
      req(gsea_result)
      df <- gsea_result@result
      df$.sign <- ifelse(df$NES > 0, "Activated", "Inhibited")
      color_by <- input$gsea_color_scale
      top_n <- input$gsea_top_n
      
      df <- head(df[order(df$p.adjust), ], top_n)
      
      if (color_by == "NES") {
        df$coreSize <- sapply(strsplit(df$core_enrichment, "/"), length)
        df$geneRatio <- df$coreSize / df$setSize  # GeneRatio equivalent
        
        p <- ggplot(df, aes(x = geneRatio, y = reorder(Description, geneRatio), size = coreSize, color = NES)) +
          geom_point() +
          scale_color_gradient2(low = "blue", high = "red", midpoint = 0) +
          labs(
            x = "Gene Ratio",
            y = NULL,
            color = "NES",
            size = "Gene Count"
          ) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 6, face = "bold"))
        if (input$gsea_split_dotplot && any(df$NES < 0, na.rm = TRUE)) {
          p <- p + facet_wrap(~.sign)
        }
        return(p)
      } else {
        if (input$gsea_split_dotplot && any(df$NES < 0, na.rm = TRUE)) {
          return(
            enrichplot::dotplot(gsea_result, showCategory = top_n, split = ".sign", color = color_by) +
              facet_grid(. ~ .sign) +
              theme(axis.text.y = element_text(size = 6, face = "bold"))
          )
        } else {
          return(
            enrichplot::dotplot(gsea_result, showCategory = top_n, color = color_by) +
              theme(axis.text.y = element_text(size = 6, face = "bold"))
          )
        }
      }
    })

    output$download_gsea_dot_plot <- downloadHandler(
      filename = function() { paste0("GSEA_","_",input$test_condition,"_vs_",input$reference_condition,"_",input$gsea_db,"_dot_plot.pdf", sep="")  },
      content = function(file) {
        pdf(file)
        gsea_result@result$.sign <- ifelse(gsea_result@result$NES > 0, "Activated", "Inhibited")
        color_by <- input$gsea_color_scale

        p <- if (input$gsea_split_dotplot && min(gsea_result@result$NES, na.rm = TRUE) < 0) {
          enrichplot::dotplot(gsea_result, showCategory = input$gsea_top_n, split = ".sign", color = color_by) +
            facet_grid(. ~ .sign) +
            theme(axis.text.y = element_text(size = 6, face = "bold"))
        } else {
          enrichplot::dotplot(gsea_result, showCategory = input$gsea_top_n, color = color_by) +
            theme(axis.text.y = element_text(size = 6, face = "bold"))
        }
        print(p)
        dev.off()
      }
    )

    output$gseaEnrichmentPlot <- renderPlot({
      req(input$gsea_selected_pathway, gsea_result)
      if (!(input$gsea_selected_pathway %in% gsea_result@result$ID)) {
        showNotification("Selected pathway not found in GSEA results.", type = "error")
        return(NULL)
      }
      enrichplot::gseaplot2(gsea_result, geneSetID = input$gsea_selected_pathway)
    })
    
    output$download_gsea_enrichment_plot <- downloadHandler(
      filename = function() {
        paste0("gsea_enrichment_plot_","_",input$test_condition,"_vs_",input$reference_condition,"_", gsub("\\W+", "_", input$gsea_selected_pathway), ".pdf")
      },
      content = function(file) {
        pdf(file, width = 8, height = 6)
        print(enrichplot::gseaplot2(gsea_result, geneSetID = input$gsea_selected_pathway))
        dev.off()
      }
    )
    
    output$GSEAupsetPlot <- renderPlot({
      #req(pathway_result)
      enrichplot::upsetplot(gsea_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
    })
    
    output$download_gsea_upset_plot <- downloadHandler(
      filename = function() { paste0("GSEA_",input$gsea_db,"_",input$test_condition,"_vs_",input$reference_condition,"_upset_plot.pdf", sep="") },
      content = function(file) {
        pdf(file)
        print(enrichplot::upsetplot(gsea_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
        dev.off()
      }
    )
    output$gseaTable <- renderDT({
      datatable(as.data.frame(gsea_result), options = list(scrollX = TRUE))
    })
    output$download_gsea_table <- downloadHandler(
      filename = function() { paste0("GSEA_",input$gsea_db,"_",input$test_condition,"_vs_",input$reference_condition,"_result.csv", sep="")  },
      content = function(file) {
        write.csv(as.data.frame(gsea_result), file, row.names = FALSE)
      }
    )
  })
}

# === Register in server ===
# mod_gsea_analysis(input, output, session, data, res_reactive)

