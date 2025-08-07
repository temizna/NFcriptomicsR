#' Non-overlapping Pathway Analysis Module
#'
#' This module handles the pathway enrichment analysis for non-overlapping genes. It filters out the genes that are already included in the initial pathway results and reruns the pathway analysis using clusterProfiler.
#'
#' @param input Shiny input object containing the user interface input values.
#' @param output Shiny output object to render the pathway analysis results and visualizations.
#' @param session Shiny session object for managing the session and UI updates.
#' @param filtered_data_rv Reactive object containing the filtered data, which includes normalized counts, samples, and species information.
#' @param res_reactive Reactive object containing the DESeq2 results (log2 fold change, p-value, and adjusted p-value).
#' @param geneList_rv Reactive object containing the initial list of genes, typically with log2 fold change values and corresponding Entrez IDs, used for the first pathway analysis.
#' @param geneListU_rv Reactive object containing the unique list of genes, with log2 fold change values and corresponding Entrez IDs, used for the non-overlapping pathway analysis.
#' @param pathway_result_rv Reactive object containing the initial pathway analysis results from clusterProfiler, including gene IDs and enrichment terms. 
#' @importFrom shiny req renderPlot renderUI downloadHandler showNotification
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom enrichplot dotplot heatplot treeplot cnetplot upsetplot emapplot
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom grDevices dev.off
#' @export
mod_pathway_analysis_non_overlap <- function(input, output, session, geneList_rv, filtered_data_rv, pathway_result_rv, res_reactive, geneListU_rv) {
  
  # Filter out the common genes between pathway results and original gene list
  observeEvent(input$run_non_overlap_pathway, {
    req(geneList_rv(), pathway_result_rv())  # Ensure we have initial results and gene list
    
    # Get the initial results and gene list
    pathway_results <- pathway_result_rv()
    initial_gene_list <- geneList_rv()
    
    # Extract the genes that are common between pathway_results and initial gene list
    pathway_genes <- unique(unlist(strsplit(pathway_results$geneID, "/")))
    common_genes <- intersect(pathway_genes, names(initial_gene_list))
    
    # Update the gene list by removing common genes
    new_gene_list <- initial_gene_list[!names(initial_gene_list) %in% common_genes]
    geneListU_rv(new_gene_list)  # Store the updated gene list in a reactive value
    
    # Now rerun the pathway analysis with the updated gene list
    species <- filtered_data_rv()$species
    orgdb <- get_orgdb(species)
    
    res <- res_reactive()
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
    
    selected_genes <- names(new_gene_list)
    
    # Run the pathway analysis with the updated gene list
    if (input$pathway_db == "GO") {
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
    } else if (input$pathway_db == "KEGG") {
      kegg_sp <- if( filtered_data_rv()$species == "Homo sapiens") "hsa" else "mmu"
      x <- clusterProfiler::enrichKEGG(
        gene = selected_genes, 
        organism = kegg_sp, 
        pvalueCutoff = input$padj_threshold,
        qvalueCutoff =  input$pathway.qval
      )
      pathway_result <- setReadable(x, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    } else {
      pathway_result <- tryCatch({
        x <- enrichPathway(
          gene = selected_genes, 
          organism = get_reactome_code(filtered_data_rv()$species), 
          pvalueCutoff = input$padj_threshold,
          qvalueCutoff = input$pathway.qval, 
          readable = TRUE
        )
        setReadable(x, OrgDb = orgdb, keyType = "ENTREZID")
      }, error = function(e) {
        showNotification(paste("Reactome pathway analysis failed:", e$message), type = "error")
        return(NULL)
      })
    }
    
    if (!is.null(pathway_result) && nrow(pathway_result@result) > 0) {
      pathway_result <- setReadable(pathway_result, OrgDb = orgdb, keyType = "ENTREZID")
    } else {
      showNotification("No enriched pathways found.", type = "warning")
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
    } else {
      pathway_result <- pairwise_termsim(pathway_result) 
    }
    
    # Proceed with rendering visualizations
    output$dotPlot_nonOL <- renderPlot({
      enrichplot::dotplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
    })
    
    output$download_nonOL_dot_plot <- downloadHandler(
      filename = function() { paste0("Pathway_",input$pathway_db,"_",input$pathway_direction,"_",input$test_condition,"_vs_",input$reference_condition,"_nonOL_dot_plot.pdf", sep="")  },
      content = function(file) {
        pdf(file)
        print(enrichplot::dotplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
        dev.off()
      }
    )
    
    output$heatmapPlot_nonOL <- renderPlot({
      enrichplot::heatplot(pathway_result, foldChange=geneListU_rv(), showCategory = 5)
    })
    
    output$download_nonOL_heatmap_plot <- downloadHandler(
      filename = function() { paste0("Pathway_",input$pathway_db,"_",input$pathway_direction,"_",input$test_condition,"_vs_",input$reference_condition,"_nonOL_heatmap_plot.pdf", sep="") },
      content = function(file) {
        pdf(file)
        print(enrichplot::heatplot(pathway_result, foldChange=geneListU_rv() , showCategory = 5))
        dev.off()
      }
    )
    
    output$treePlot_nonOL <- renderPlot({
      enrichplot::treeplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
    })
    
    output$download_nonOL_tree_plot <- downloadHandler(
      filename = function() { paste0("Pathway_",input$pathway_db,"_",input$pathway_direction,"_",input$test_condition,"_vs_",input$reference_condition,"_nonOL_tree_plot.pdf", sep="") },
      content = function(file) {
        pdf(file)
        print(enrichplot::treeplot(pathway_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
        dev.off()
      }
    )
    
    output$pathwayTable_nonOL <- renderDT({
      req(pathway_result)
      as.data.frame(pathway_result)
    })
    
    output$download_pathway_nonOL_table <- downloadHandler(
      filename = function() { paste0("Pathway_",input$pathway_db,"_",input$pathway_direction,"_",input$test_condition,"_vs_",input$reference_condition,"_nonOL_results.csv", sep="")  },
      content = function(file) {
        write.csv(as.data.frame(pathway_result), file, row.names = FALSE)
      }
    )
  })
}

