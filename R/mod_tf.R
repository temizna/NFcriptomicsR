# === Module: mod_tf_enrichment_analysis ===

#' Transcription Factor Enrichment Analysis Module
#'
#' This module allows users to perform transcription factor enrichment analysis
#' based on selected datasets stored locally (gzipped) in inst/extdata.
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param res_reactive Reactive value containing the gene list (log2FC vector) for enrichment
#' @param filtered_data_rv reactive list containing data
#' @importFrom gson read.gmt
#' @importFrom dplyr left_join
#' @importFrom vroom vroom
#' @export
mod_tf_enrichment_analysis <- function(input, output, session,res_reactive,filtered_data_rv) {
  
  load_tf_data <- function(tf_data_source) {
    tf_data_file_map <- list(
      "TRANSFAC and JASPAR PWMs" = "TRANSFAC_and_JASPAR_PWMs.h.gmt.gz",
      "ENCODE and ChEA Consensus TFs from ChIP-X" = "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt.gz",
      "TRRUST_Transcription_Factors_2019" = "TRRUST_Transcription_Factors_2019h.gmt.gz",
      "TF_Perturbations_Followed_by_Expression" = "TF_Perturbations_Followed_by_Expression_Human.gmt.gz",
      "hTFtarget" = "TF-Target-information.txt.gz",
      "TFLink" = "TFLink_Homo_sapiens_interactions_All_GMT_proteinName_v1.0.gmt.gz"
    )
    
    tf_data_file <- tf_data_file_map[[tf_data_source]]
    if (is.null(tf_data_file)) stop("Invalid TF data source selected.")
    
    tf_file_path <- system.file("extdata", tf_data_file, package = "TranscriptoPathR")
    if (tf_file_path == "") stop("TF data file not found.")
    
    if (grepl(".gmt", tf_data_file)) {
      tf_data <- gson::read.gmt(tf_file_path)
    } else {
      tf_data <- vroom::vroom(tf_file_path, delim = "\t", show_col_types = FALSE)
    }
    
    return(tf_data)
  }
  
  # Reactive value to store enrichment results
  processed_tf_data <- reactiveVal()
  
  # Observe button click
  observeEvent(input$run_tf_enrichment, {
    req(res_reactive(),filtered_data_rv()$species)
    print("Running TF Enrichment Analysis...")
    tf_data_source <- input$tf_data_source
    print(paste("Selected TF data source:", tf_data_source))
    
    tryCatch({
      # Load TF data using helper function
      tf_data <- load_tf_data(tf_data_source)
      
      #print("TF data loaded:")
      #print(head(tf_data))
      
      # Process TF data: always map target genes to ENTREZ IDs
      tf_data <- tf_data %>%
        rename(TF_name = 1, target_gene = 2)
      
      tf_data_entrez <- suppressMessages({bitr(
        tf_data$target_gene,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      )})
      
      tf_data_gmt <- tf_data %>%
        dplyr::left_join(tf_data_entrez, by = c("target_gene" = "SYMBOL")) %>%
        dplyr::select(TF_name, ENTREZID) %>%
        dplyr::filter(!is.na(ENTREZID))
      
      species <- filtered_data_rv()$species
      orgdb <- get_orgdb(species)
      res <- isolate(res_reactive())
      direction <- input$gene_direction
      
      res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
      d1 <- res[, c("log2FoldChange", "padj")]
      d1$gene <- toupper(rownames(res)) # converting mouse to human assuming similar gene names for ortologs
      
      if (is_symbol(d1$gene)) {
        d1_ids <- suppressMessages({bitr(d1$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)})
      } else {
        d1_ids <- suppressMessages({bitr(d1$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgdb)})
      }
      
      d1_merged <- merge(d1, d1_ids, by.x = "gene", by.y = 1)
      d1_merged <- d1_merged[!duplicated(d1_merged$ENTREZID), ]
      gene_vector <- switch(direction,
                            "Up" = d1_merged[d1_merged$log2FoldChange >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ],
                            "Down" = d1_merged[d1_merged$log2FoldChange <= -input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ],
                            "Both"= d1_merged[abs(d1_merged$log2FoldChange) >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ]
      )
      #gene_vector <- d1_merged[abs(d1_merged$log2FoldChange) >= input$lfc_threshold & d1_merged$padj <= input$padj_threshold, ]
      gene_vector <- gene_vector[!is.na(gene_vector$log2FoldChange) & !is.na(gene_vector$ENTREZID), ]
      geneList <- gene_vector$log2FoldChange
      names(geneList) <- gene_vector$ENTREZID
      geneList <- geneList[!duplicated(names(geneList))]
      geneList <- sort(geneList, decreasing = TRUE)
      max_genes <- 1000
      if (length(geneList) > max_genes) geneList <- head(geneList, max_genes)
      
      selected_genes <- names(geneList)
      if(input$enrichment_method=="Over_representation"){
        tf_result <- clusterProfiler::enricher(
          gene = selected_genes,
          TERM2GENE = tf_data_gmt,
          pvalueCutoff = 0.05,
          qvalueCutoff = input$tf.qval
        )
      } else {
        tf_result <- tryCatch({
          clusterProfiler::GSEA(
            geneList = geneList,
            TERM2GENE = tf_data_gmt,
            pvalueCutoff = input$tf.qval,
            verbose = FALSE
          )
        }, error = function(e) {
          showNotification(paste("TF:GSEA failed:", e$message), type = "error")
          return(NULL)
        })
      }


      if (is.null(tf_result) || nrow(as.data.frame(tf_result)) < 1) {
        showNotification("No enriched terms found in TF enrichment.", type = "warning")
        return()
      }
      
      output$tf_dotplot <- renderPlot({
        enrichplot::dotplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
      })
      
      output$download_tf_dotplot <- downloadHandler(
        filename = function() paste0("TF_Enrichment_", input$tf_data_source, "_dotplot.pdf"),
        content = function(file) {
          p<-enrichplot::dotplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
          ggsave(file, p,  width = 10, height = 8, units = "in")
        }
      )
      if(input$enrichment_method=="GSEA") {
        output$tf_ridgeplot <- renderPlot({
          enrichplot::ridgeplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
        })
        
        output$download_tf_ridgeplot <- downloadHandler(
          filename = function() paste0("TF_Enrichment_", input$tf_data_source, "_ridgeplot.pdf"),
          content = function(file) {
            p<-enrichplot::ridgeplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
            ggsave(file, p,  width = 10, height = 8, units = "in")
          }
        )
      }  else {
        showNotification("No ridge plot in over representation analysis.", type = "warning")
        return()
      }

      output$tf_results_table <- renderDT({
        as.data.frame(tf_result@result)
      })
      
      output$download_tf_results_table <- downloadHandler(
        filename = function() paste0("TF_Enrichment_", input$tf_data_source, "_results.csv"),
        content = function(file) {
          write.csv(as.data.frame(tf_result@result), file, row.names = FALSE)
        }
      )
      
    }, error = function(e) {
      showNotification(paste("Error loading or processing TF data:", e$message), type = "error")
    })
  })
}

