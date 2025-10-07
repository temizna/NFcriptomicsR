# R/logger.R

#' Append a message to the app log (safe no-op fallback)
#'
#' Prints a timestamped message to the console.
#'
#' @param msg   Message to log (character)
#' @param level Log level (INFO, WARN, ERROR)
#' @export
append_log <- function(msg, level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %-5s %s\n", ts, toupper(level), msg))
}

#' Logger Module UI
#'
#' @param id Shiny module id
#' @return A Shiny UI element to display the log and a download button
#' @importFrom shiny NS tagList h3 div verbatimTextOutput downloadButton
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

#' Logger Module Server (fixed & robust)
#'
#' Call from your main server: `mod_logger_server("logger", input_all = input)`
#'
#' @param id Shiny module id
#' @param input_all The full Shiny input object (passed from main server)
#' @return A reactive expression with log entries
#' @importFrom shiny moduleServer observeEvent reactiveVal renderText downloadHandler
#' @export
mod_logger_server <- function(id, input_all) {
  moduleServer(id, function(input, output, session) {
    user_log <- reactiveVal(character())
    
    # format list of parameter values pulled from input_all
    log_inputs <- function(params) {
      if (is.null(params) || !length(params)) return("")
      paste(
        sapply(params, function(pn) {
          val <- input_all[[pn]]
          if (length(val) > 1) val <- paste(val, collapse = ", ")
          sprintf("%s: %s", pn, if (length(val)) as.character(val) else "NULL")
        }),
        collapse = "; "
      )
    }
    
    # single helper: append to UI buffer + console
    add_entry <- function(title, params = NULL, level = "INFO") {
      ts    <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      param <- log_inputs(params)
      line  <- if (nzchar(param))
        sprintf("[%s] %s\nParameters: %s", ts, title, param)
      else
        sprintf("[%s] %s", ts, title)
      
      # UI buffer
      old <- user_log()
      user_log(c(old, line))
      
      # console
      append_log(line, level)
      invisible(line)
    }
    
    # --------------------------
    # Sample Select
    # --------------------------
    observeEvent(input_all$select_all, {
      add_entry("User clicked sample_select-select_all")
    }, ignoreInit = TRUE)
    
    observeEvent(input_all$run_filter, {
      add_entry(
        "User clicked sample_select-run_filter",
        params = c("filter_dim","filter_treatment","filter_batch","filter_cellline","filter_PR","filter_ER","sample_select")
      )
    }, ignoreInit = TRUE)
    
    observeEvent(input_all$deselect_all, {
      add_entry("User clicked sample_select-deselect_all")
    }, ignoreInit = TRUE)
    
    # --------------------------
    # Differential Expression
    # --------------------------
    observeEvent(input_all$run_de, {
      add_entry(
        "User clicked run_de",
        params = c("metadata_column", "reference_condition", "test_condition",
                   "lfc_threshold", "padj_threshold", "num_genes", "cluster_columns")
      )
    }, ignoreInit = TRUE)
    
    # --------------------------
    # Cross Plot
    # --------------------------
    observeEvent(input_all$run_crossplot, {
      add_entry(
        "User clicked run_crossplot",
        params = c("crossplot_gene_label", "metadata_column_x", "reference_condition_x",
                   "test_condition_x", "metadata_column_y", "reference_condition_y",
                   "test_condition_y", "crossplot_gene_count", "crossplot_topgenes")
      )
    }, ignoreInit = TRUE)
    
    # --------------------------
    # GSEA
    # --------------------------
    observeEvent(input_all$run_gsea, {
      add_entry(
        "User clicked run_gsea",
        params = c("gsea_split_dotplot", "gsea_color_scale", "gsea_db",
                   "gsea_top_n", "lfc_threshold", "padj_threshold", "gsea_pvalue")
      )
    }, ignoreInit = TRUE)
    
    # --------------------------
    # GSVA / ssGSEA (module id: gsva)
    # These inputs exist because your GSVA module is instantiated with id = "gsva"
    # --------------------------
    observeEvent(input_all[["gsva-run"]], {
      add_entry(
        "Clicked gsva-run",
        params = c("gsva-gset_source","gsva-gsva_db","gsva-gsva_method","gsva-mx_opt",
                   "gsva-min_gs_size","gsva-max_gs_size","gsva-use_de_defaults",
                   "gsva-gsva_metadata_column","gsva-gsva_reference_condition","gsva-gsva_test_condition")
      )
    }, ignoreInit = TRUE)
    
    observeEvent(input_all[["gsva-dl_table"]],    { add_entry("Clicked gsva-dl_table") },    ignoreInit = TRUE)
    observeEvent(input_all[["gsva-dl_scores"]],   { add_entry("Clicked gsva-dl_scores") },   ignoreInit = TRUE)
    observeEvent(input_all[["gsva-dl_volcano"]],  { add_entry("Clicked gsva-dl_volcano") },  ignoreInit = TRUE)
    observeEvent(input_all[["gsva-dl_heatmap"]],  { add_entry("Clicked gsva-dl_heatmap") },  ignoreInit = TRUE)
    observeEvent(input_all[["gsva-dl_boxplots"]], { add_entry("Clicked gsva-dl_boxplots") }, ignoreInit = TRUE)
    
    # --------------------------
    # Pathway Analysis
    # --------------------------
    observeEvent(input_all$run_pathway, {
      add_entry(
        "User clicked run_pathway",
        params = c("pathway_db", "pathway_direction", "circular_layout",
                   "lfc_threshold", "padj_threshold", "pathway.qval", "max_genes")
      )
    }, ignoreInit = TRUE)
    
    # Non-overlap Pathway
    observeEvent(input_all$run_non_overlap_pathway, {
      add_entry("User clicked run_non_overlap_pathway")
    }, ignoreInit = TRUE)
    
    # TF Enrichment
    observeEvent(input_all$run_tf_enrichment, {
      add_entry("User clicked run_tf_enrichment", params = c("tf_data_source", "tf.qval"))
    }, ignoreInit = TRUE)
    
    # Cancer Gene Census
    observeEvent(input_all$run_cancer_gene_census, {
      add_entry("User clicked run_cancer_gene_census")
    }, ignoreInit = TRUE)
    
    # PCA
    observeEvent(input_all$run_pca, {
      add_entry(
        "User clicked run_PCA",
        params = c("max_pc", "variance_threshold","similarity_threshold","max_genes", "min_genes", "pca_enrich_method")
      )
    }, ignoreInit = TRUE)
    
    # -------- Outputs --------
    output$log_text <- renderText({
      paste(user_log(), collapse = "\n\n")
    })
    
    output$download_log <- downloadHandler(
      filename = function() paste0("user_session_log_", Sys.Date(), ".txt"),
      content  = function(file) writeLines(user_log(), con = file)
    )
    
    session$onSessionEnded(function() {
      log_dir <- file.path(getwd(), "logs")
      if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
      file_path <- file.path(
        log_dir,
        paste0("user_session_log_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".txt")
      )
      writeLines(isolate(user_log()), file_path)
    })
    
    return(user_log)
  })
}
