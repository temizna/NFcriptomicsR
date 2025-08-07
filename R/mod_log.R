#' Logger Module Server
#'
#' @param id Shiny module id
#' @param input_all The full Shiny input object (passed from main server)
#' @return A reactive expression with log entries
#' @export
# Logger Module Server (fixed)
mod_logger_server <- function(id, input_all) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    user_log <- reactiveVal(character())
    
    # === Helper to log params ===
    log_inputs <- function(params) {
      paste(
        sapply(params, function(pn) {
          val <- input_all[[pn]]
          val_str <- if (length(val) > 1) paste(val, collapse = ", ") else as.character(val)
          paste0(pn, ": ", val_str)
        }),
        collapse = "; "
      )
    }
    
    # === Explicit observeEvent for each action button ===
    
    
    # Sample Select
    observeEvent(input_all$select_all, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      entry <- paste0("[", timestamp, "] User clicked sample_select-select_all")
      user_log(c(user_log(), entry))
    })
    
    observeEvent(input_all$run_filter, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      params<-c("filter_dim","filter_treatment","filter_batch","filter_cellline","filter_PR","filter_ER","sample_select")
      entry <- paste0("[", timestamp, "] User clicked sample_select-run_filter\nParameters: ", log_inputs(params))
      user_log(c(user_log(), entry))
    })
    
    observeEvent(input_all$deselect_all, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      entry <- paste0("[", timestamp, "] User clicked sample_select-deselect_all")
      user_log(c(user_log(), entry))
    })
    
    # Differential Expression
    observeEvent(input_all$run_de, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      params <- c("metadata_column", "reference_condition", "test_condition",
                  "lfc_threshold", "padj_threshold", "num_genes", "cluster_columns")
      entry <- paste0("[", timestamp, "] User clicked run_de\nParameters: ", log_inputs(params))
      user_log(c(user_log(), entry))
    })
    
    # Cross Plot
    observeEvent(input_all$run_crossplot, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      params <- c("crossplot_gene_label", "metadata_column_x", "reference_condition_x",
                  "test_condition_x", "metadata_column_y", "reference_condition_y",
                  "test_condition_y", "crossplot_gene_count", "crossplot_topgenes")
      entry <- paste0("[", timestamp, "] User clicked run_crossplot\nParameters: ", log_inputs(params))
      user_log(c(user_log(), entry))
    })
    
    # GSEA
    observeEvent(input_all$run_gsea, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      params <- c("gsea_split_dotplot", "gsea_color_scale", "gsea_db",
                  "gsea_top_n", "lfc_threshold", "padj_threshold", "gsea_pvalue")
      entry <- paste0("[", timestamp, "] User clicked run_gsea\nParameters: ", log_inputs(params))
      user_log(c(user_log(), entry))
    })
    
    # Pathway Analysis
    observeEvent(input_all$run_pathway, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      params <- c("pathway_db", "pathway_direction", "circular_layout",
                  "lfc_threshold", "padj_threshold", "pathway.qval", "max_genes")
      entry <- paste0("[", timestamp, "] User clicked run_pathway\nParameters: ", log_inputs(params))
      user_log(c(user_log(), entry))
    })
    
    # Non-overlap Pathway
    observeEvent(input_all$run_non_overlap_pathway, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      entry <- paste0("[", timestamp, "] User clicked run_non_overlap_pathway")
      user_log(c(user_log(), entry))
    })
    
    # TF Enrichment
    observeEvent(input_all$run_tf_enrichment, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      params <- c("tf_data_source", "tf.qval")
      entry <- paste0("[", timestamp, "] User clicked run_tf_enrichment\nParameters: ", log_inputs(params))
      user_log(c(user_log(), entry))
    })
    
    # Cancer Gene Census
    observeEvent(input_all$run_cancer_gene_census, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      entry <- paste0("[", timestamp, "] User clicked run_cancer_gene_census")
      user_log(c(user_log(), entry))
    })
    
    # Consensus Cluster
    observeEvent(input_all$run_pca, {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      params <- c("max_pc", "variance_threshold","similarity_threshold","max_genes", "min_genes", "pca_enrich_method")
      entry <- paste0("[", timestamp, "] User clicked run_PCA\nParameters: ", log_inputs(params))
      user_log(c(user_log(), entry))
    })
    
    # === Outputs ===
    output$log_text <- renderText({ paste(user_log(), collapse = "\n\n") })
    
    output$download_log <- downloadHandler(
      filename = function() paste0("user_session_log_", Sys.Date(), ".txt"),
      content = function(file) { writeLines(user_log(), con = file) }
    )
    
    session$onSessionEnded(function() {
      log_dir <- file.path(getwd(), "logs")
      if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
      file_path <- file.path(log_dir, paste0("user_session_log_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".txt"))
      
      # freeze log before accessing outside reactive context
      log_snapshot <- isolate(user_log())
      writeLines(log_snapshot, file_path)
    })
    return(user_log)
  })
}


#' Logger Module UI
#'
#' @param id Shiny module id
#' @return A Shiny UI element to display the log and a download button
#' @importFrom shiny  h3 div verbatimTextOutput downloadButton moduleServer reactive reactiveValues renderText
#' @export
mod_logger_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("User Session Log"),
    div(
      style = "overflow-y: auto; height: 500px; border: 1px solid #ccc; padding: 10px; white-space: pre-wrap; background-color: #f9f9f9; font-family: monospace; font-size: 12px;",
      verbatimTextOutput(ns("log_text"))
    ),
    downloadButton(ns("download_log"), "Download Log")
  )
}

