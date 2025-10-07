# R/mod_gsva.R

#' GSVA / ssGSEA module UI (uses DE tab selections by default)
#' @param id Shiny module ID
#' @export
mod_gsva_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # Gene set source
        shiny::radioButtons(
          ns("gset_source"), "Gene Set Source",
          choices = c("Built-in MSigDB" = "msigdb", "Upload GMT/RDS" = "upload"),
          selected = "msigdb"
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'msigdb'", ns("gset_source")),
          shiny::selectInput(
            ns("gsva_db"), "Select Database:",
            choices = c("GO","KEGG","Reactome","Hallmark",
                        "Cancer Cell Atlas","Cancer Gene Neighbourhoods",
                        "Cancer Modules","Txn Factor Targets"),
            selected = "Hallmark"
          ),
          shiny::helpText("Species uses filtered_data_rv()$species. Gene IDs are matched to SYMBOL.")
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'upload'", ns("gset_source")),
          shiny::fileInput(ns("gene_sets"), "Upload Gene Sets (.gmt or .rds/.RDS)",
                           accept = c(".gmt", ".rds", ".RDS"))
        ),
        
        # Method + params (default to ssGSEA)
        shiny::selectInput(ns("gsva_method"), "Method",
                           choices = c("ssgsea","gsva"), selected = "ssgsea"),
        shiny::checkboxInput(ns("mx_opt"),
                             "GSVA: mx.diff | ssGSEA: normalize scores", value = TRUE),
        shiny::numericInput(ns("min_gs_size"), "Min gene set size", value = 10, min = 1, step = 1),
        shiny::numericInput(ns("max_gs_size"), "Max gene set size", value = 500, min = 10, step = 10),
        
        shiny::hr(),
        # Use DE tab selections by default; show overrides only if unchecked
        shiny::checkboxInput(ns("use_de_defaults"),
                             "Use selections from Differential Expression tab", value = TRUE),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == false", ns("use_de_defaults")),
          shiny::selectInput(ns("gsva_metadata_column"), "Variable to test:", choices = NULL),
          shiny::selectInput(ns("gsva_reference_condition"), "Reference Condition:", choices = NULL),
          shiny::selectInput(ns("gsva_test_condition"), "Test Condition:", choices = NULL)
        ),
        
        shiny::hr(),
        # Dynamic multi-select for extra column annotations
        shiny::uiOutput(ns("ann_cols_ui")),
        
        shiny::actionButton(ns("run"), "Run GSVA + Differential"),
        
        shiny::hr(),
        shiny::downloadButton(ns("dl_table"),   "Download Results (CSV)"),
        shiny::downloadButton(ns("dl_scores"),  "Download GSVA Scores (CSV)"),
        shiny::downloadButton(ns("dl_volcano"), "Download Volcano (PDF)"),
        shiny::downloadButton(ns("dl_heatmap"), "Download Heatmap (PDF)"),
        shiny::downloadButton(ns("dl_boxplots"),"Download Boxplots (PDF)")
      ),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Results Table", DT::dataTableOutput(ns("tbl"))),
          shiny::tabPanel("Volcano", shiny::plotOutput(ns("volcano"), height = "520px")),
          shiny::tabPanel("Heatmap", shiny::plotOutput(ns("heatmap"), height = "820px")),
          shiny::tabPanel("Top Pathway Boxplots", shiny::plotOutput(ns("boxplots"), height = "600px"))
        )
      )
    )
  )
}

#' GSVA / ssGSEA module Server
#'
#' @param id Shiny module ID
#' @param dds_rv reactiveVal returning a DESeqDataSet
#' @param filtered_data_rv() ReactiveValues with at least $species
#' @param group_var_react optional reactive() for DE tab "Variable to test"
#' @param ref_level_react optional reactive() for DE tab "Reference Condition"
#' @param test_level_react optional reactive() for DE tab "Test Condition"
#' @importFrom SummarizedExperiment assay colData
#' @importFrom DESeq2 vst
#' @importFrom GSVA gsva gsvaParam ssgseaParam
#' @importFrom GSEABase getGmt geneIds
#' @importFrom msigdbr msigdbr msigdbr_species
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom DT renderDT datatable
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal geom_boxplot facet_wrap
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.newpage grid.text
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr arrange
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom scales hue_pal
#' @importFrom stats model.matrix relevel aggregate
#' @importFrom utils write.csv
#' @importFrom grDevices pdf dev.off
#' @importFrom BiocParallel SerialParam
#' @importFrom shiny showNotification
#' @export
mod_gsva_server <- function(
    id,
    dds_rv,
    filtered_data_rv,
    group_var_react = NULL,
    ref_level_react = NULL,
    test_level_react = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # --- Populate local override controls (only used if use_de_defaults == FALSE)
    shiny::observe({
      shiny::req(dds_rv())
      meta <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
      fact_cols <- names(meta)[sapply(meta, function(x) is.factor(x) || (is.character(x) && length(unique(x)) <= 50))]
      shiny::updateSelectInput(session, "gsva_metadata_column",
                               choices = fact_cols,
                               selected = if (length(fact_cols)) fact_cols[1] else character(0))
    })
    shiny::observeEvent(input$gsva_metadata_column, {
      shiny::req(dds_rv(), input$gsva_metadata_column)
      meta <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
      lv <- unique(as.character(meta[[input$gsva_metadata_column]]))
      shiny::updateSelectInput(session, "gsva_reference_condition",
                               choices = lv, selected = if (length(lv)) lv[1] else character(0))
      shiny::updateSelectInput(session, "gsva_test_condition",
                               choices = lv, selected = if (length(lv) >= 2) lv[2] else if (length(lv)) lv[1] else character(0))
    }, ignoreInit = TRUE)
    
    # creates the multi-select dynamically (extra column annotations)
    output$ann_cols_ui <- shiny::renderUI({
      shiny::req(dds_rv())
      meta <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
      shiny::selectInput(
        ns("ann_cols"),
        "Additional column annotations:",
        choices = colnames(meta),
        multiple = TRUE
      )
    })
    
    # Helper to decide which selections to use (DE tab by default)
    gv_selected <- shiny::reactive({
      if (isTRUE(input$use_de_defaults) && !is.null(group_var_react)) group_var_react() else input$gsva_metadata_column
    })
    ref_selected <- shiny::reactive({
      if (isTRUE(input$use_de_defaults) && !is.null(ref_level_react)) ref_level_react() else input$gsva_reference_condition
    })
    test_selected <- shiny::reactive({
      if (isTRUE(input$use_de_defaults) && !is.null(test_level_react)) test_level_react() else input$gsva_test_condition
    })
    
    # ---- analysis
    results_react <- shiny::eventReactive(input$run, {
      shiny::req(dds_rv(), filtered_data_rv()$species)
      
      # Pull group/levels (DE tab by default)
      gv    <- gv_selected()
      refL0 <- ref_selected()
      testL0<- test_selected()
      shiny::req(gv, refL0, testL0)
      
      # 1) VST transform (genes x samples)
      vsd  <- DESeq2::vst(dds_rv(), blind = TRUE)
      expr <- SummarizedExperiment::assay(vsd)
      
      # 1b) Match gene IDs to SYMBOL (like your GSEA module)
      if (!is_symbol(rownames(expr))) {
        conv <- convert_ensembl_to_symbol(rownames(expr), filtered_data_rv()$species)
        keep <- !is.na(conv)
        expr <- expr[keep, , drop = FALSE]
        rownames(expr) <- conv[keep]
        # collapse duplicate symbols by mean
        df <- as.data.frame(expr)
        df$gene <- rownames(df)
        df <- stats::aggregate(. ~ gene, data = df, FUN = mean)
        rownames(df) <- df$gene; df$gene <- NULL
        expr <- as.matrix(df)
      }
      rownames(expr) <- make.unique(rownames(expr))
      
      # 2) Species mapping (exactly like your GSEA)
      species_map <- list(
        "Homo sapiens" = "Homo sapiens",
        "Mus musculus" = "Mus musculus",
        "Rattus norvegicus" = "Rattus norvegicus",
        "Canis familiaris" = "dog",
        "Saccharomyces cerevisiae" = "Saccharomyces cerevisiae"
      )
      species_name <- species_map[[filtered_data_rv()$species]]
      shiny::validate(shiny::need(!is.null(species_name), "Selected species not supported by msigdbr."))
      
      # 3) Gene sets (MSigDB or upload) - SYMBOL-based
      if (identical(input$gset_source, "msigdb")) {
        supported <- msigdbr::msigdbr_species()$species_name
        shiny::validate(shiny::need(species_name %in% supported,
                                    paste0("Species '", species_name, "' not in msigdbr_species().")))
        coll <- .map_db_to_collection(input$gsva_db)
        subc <- .map_db_to_subcollection(input$gsva_db)
        m_df <- msigdbr::msigdbr(species = species_name, collection = coll, subcollection = subc)
        shiny::validate(shiny::need(nrow(m_df) > 0, "No gene sets returned from msigdbr for this selection."))
        gsets <- split(m_df$gene_symbol, m_df$gs_name)
      } else {
        shiny::req(input$gene_sets)
        gsets <- .read_gene_sets(input$gene_sets$datapath)
        if (inherits(gsets, "GeneSetCollection")) gsets <- GSEABase::geneIds(gsets)
        shiny::validate(shiny::need(is.list(gsets), "Uploaded RDS must be a list or GeneSetCollection."))
      }
      
      # Optional sanity: ensure overlap
      overlap_n <- length(intersect(rownames(expr), unique(unlist(gsets))))
      if (overlap_n < 25) shiny::showNotification(
        paste0("Warning: low symbol overlap between expression and gene sets (n=", overlap_n, ")."),
        type = "warning"
      )
      
      # 4) GSVA / ssGSEA - support legacy GSVA and newer Param API
      gsva_ns   <- asNamespace("GSVA")
      gsva_fun  <- get("gsva", gsva_ns)
      exports   <- getNamespaceExports("GSVA")
      has_param_api <- all(c("gsvaParam","ssgseaParam") %in% exports)
      
      minSize <- input$min_gs_size; maxSize <- input$max_gs_size
      if (has_param_api) {
        if (identical(input$gsva_method, "gsva")) {
          param <- GSVA::gsvaParam(
            exprData = expr,
            geneSets = gsets,
            minSize  = minSize,
            maxSize  = maxSize,
            kcdf     = "Gaussian",
            maxDiff  = isTRUE(input$mx_opt)
          )
        } else {
          param <- GSVA::ssgseaParam(
            exprData  = expr,
            geneSets  = gsets,
            minSize   = minSize,
            maxSize   = maxSize,
            alpha     = 0.25,
            normalize = isTRUE(input$mx_opt)
          )
        }
        gsva_mat <- GSVA::gsva(param, BPPARAM = BiocParallel::SerialParam())
      } else {
        formals_gsva <- names(formals(gsva_fun))
        gsva_args <- list(
          expr          = expr,
          gset.idx.list = gsets,
          method        = input$gsva_method,
          min.sz        = minSize,
          max.sz        = maxSize,
          kcdf          = "Gaussian",
          verbose       = FALSE
        )
        if ("BPPARAM" %in% formals_gsva) gsva_args$BPPARAM <- BiocParallel::SerialParam()
        if ("parallel.sz" %in% formals_gsva) gsva_args$parallel.sz <- 1L
        if (identical(input$gsva_method, "gsva"))   gsva_args$mx.diff     <- isTRUE(input$mx_opt)
        if (identical(input$gsva_method, "ssgsea")) gsva_args$ssgsea.norm <- isTRUE(input$mx_opt)
        gsva_mat <- do.call(gsva_fun, gsva_args)
      }
      
      # 5) limma differential - robust to non-syntactic level names
      meta <- as.data.frame(SummarizedExperiment::colData(dds_rv()))
      shiny::validate(shiny::need(gv %in% colnames(meta), sprintf("Grouping column '%s' not in colData(dds).", gv)))
      group_orig <- factor(meta[[gv]])
      shiny::validate(shiny::need(all(c(refL0, testL0) %in% levels(group_orig)),
                                  "Selected levels not found in the grouping column."))
      
      # Keep display labels, but use syntactic labels in design/contrasts
      lvl_disp <- levels(group_orig)
      lvl_safe <- make.names(lvl_disp)
      group_safe <- factor(group_orig, levels = lvl_disp, labels = lvl_safe)
      
      design <- stats::model.matrix(~ 0 + group_safe)
      colnames(design) <- lvl_safe
      
      refL  <- lvl_safe[match(refL0,  lvl_disp)]
      testL <- lvl_safe[match(testL0, lvl_disp)]
      
      fit <- limma::lmFit(gsva_mat, design)
      contr_str <- paste(testL, "-", refL, sep = "")
      fit2 <- limma::eBayes(limma::contrasts.fit(fit, limma::makeContrasts(contrasts = contr_str, levels = design)))
      tt <- tibble::rownames_to_column(limma::topTable(fit2, number = Inf, sort.by = "P"), "Pathway")
      
      list(
        gsva_mat = gsva_mat,
        tt = tt,
        meta = meta,
        group_disp = factor(group_orig),   # display labels
        group_safe = group_safe,           # safe labels used in modeling
        contrast_label = paste0(testL0, "_vs_", refL0),
        group_var = gv
      )
    })
    
    # ---- outputs
    output$tbl <- DT::renderDT({
      shiny::req(results_react())
      dplyr::arrange(results_react()$tt, .data$adj.P.Val) |>
        DT::datatable(options = list(pageLength = 25, scrollX = TRUE))
    })
    
    output$volcano <- shiny::renderPlot({
      shiny::req(results_react())
      df <- results_react()$tt
      df$neglog10FDR <- -log10(pmax(df$adj.P.Val, 1e-300))
      sig <- df$adj.P.Val < 0.05 & abs(df$logFC) > 0.2
      ggplot2::ggplot(df, ggplot2::aes(x = .data$logFC, y = .data$neglog10FDR, color = sig)) +
        ggplot2::geom_point(alpha = 0.65, size = 1.8) +
        ggplot2::scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D81B60")) +
        ggplot2::labs(title = paste0("Differential Pathway Activity: ", results_react()$contrast_label),
                      x = "logFC (GSVA score)", y = "-log10(FDR)") +
        ggplot2::theme_minimal(base_size = 12)
    })
    
    # ---- Heatmap helpers (annotation + z-score)
    .build_col_annotation <- function(meta, group_disp, ann_cols = NULL) {
      df <- data.frame(Group = group_disp, check.names = FALSE)
      col_list <- list()
      
      # Group colors
      g_lvls <- levels(group_disp)
      col_list$Group <- setNames(scales::hue_pal()(length(g_lvls)), g_lvls)
      
      # Add extra columns
      if (!is.null(ann_cols) && length(ann_cols)) {
        for (cn in ann_cols) {
          v <- meta[[cn]]
          if (is.null(v)) next
          if (is.numeric(v)) {
            rng <- stats::quantile(v, probs = c(0.02, 0.98), na.rm = TRUE)
            if (any(is.na(rng)) || diff(rng) == 0) rng <- range(v, na.rm = TRUE)
            col_list[[cn]] <- circlize::colorRamp2(
              c(rng[1], mean(rng, na.rm = TRUE), rng[2]),
              c("#08306b", "#f7f7f7", "#7f0000")
            )
            df[[cn]] <- v
          } else {
            f <- factor(as.character(v))
            lv <- levels(f)
            n  <- length(lv)
            pal <- if (n <= 8) RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)]
            else scales::hue_pal()(n)
            col_list[[cn]] <- setNames(pal, lv)
            df[[cn]] <- f
          }
        }
      }
      
      ComplexHeatmap::HeatmapAnnotation(
        df = df,
        col = col_list,
        annotation_name_side = "left",
        annotation_legend_param = list(nrow = 2)
      )
    }
    
    .row_zscore <- function(mat) {
      mat <- as.matrix(mat)
      z <- t(scale(t(mat))); z[is.na(z)] <- 0; z
    }
    
    # ---- Heatmap (screen)
    output$heatmap <- shiny::renderPlot({
      shiny::req(results_react())
      gsva_mat <- results_react()$gsva_mat
      tt <- results_react()$tt
      topN <- head(tt$Pathway, 50)
      sub  <- gsva_mat[intersect(rownames(gsva_mat), topN), , drop = FALSE]
      shiny::validate(shiny::need(nrow(sub) > 1, "No pathways available for heatmap."))
      
      zmat <- .row_zscore(sub)
      ann  <- .build_col_annotation(results_react()$meta,
                                    results_react()$group_disp,
                                    input$ann_cols)
      
      hm_col <- circlize::colorRamp2(c(-2, 0, 2), c("#2b8cbe", "white", "#e34a33"))
      
      grid::grid.newpage()
      ComplexHeatmap::draw(
        ComplexHeatmap::Heatmap(
          zmat,
          name = "Z-score",
          col  = hm_col,
          top_annotation = ann,
          show_row_names = TRUE,
          show_column_names = FALSE,
          cluster_rows = TRUE,
          cluster_columns = TRUE
        ),
        heatmap_legend_side = "right",
        annotation_legend_side = "bottom"
      )
    })
    
    # ---- Boxplots (screen)
    output$boxplots <- shiny::renderPlot({
      shiny::req(results_react())
      gsva_mat <- results_react()$gsva_mat
      tt <- results_react()$tt
      meta <- results_react()$meta
      gv <- results_react()$group_var
      top6 <- head(tt$Pathway, 6)
      if (!length(top6)) return(NULL)
      
      df <- t(gsva_mat[top6, , drop = FALSE]) |> as.data.frame()
      df[[gv]] <- meta[[gv]]
      df_long <- df |>
        tibble::rownames_to_column("Sample") |>
        tidyr::pivot_longer(cols = tidyselect::all_of(top6), names_to = "Pathway", values_to = "Score")
      ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[gv]], y = .data$Score, fill = .data[[gv]])) +
        ggplot2::geom_boxplot(outlier.alpha = 0.25) +
        ggplot2::facet_wrap(~ Pathway, scales = "free_y", ncol = 3) +
        ggplot2::labs(x = gv, y = "GSVA score", title = "Top Differential Pathways") +
        ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "none")+
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
        )
    })
    
    # ---- downloads (PDF/CSV)
    output$dl_table <- shiny::downloadHandler(
      filename = function() paste0("gsva_differential_", results_react()$contrast_label, ".csv"),
      content  = function(file) utils::write.csv(results_react()$tt, file, row.names = FALSE)
    )
    output$dl_scores <- shiny::downloadHandler(
      filename = function() "gsva_scores_matrix.csv",
      content  = function(file) utils::write.csv(results_react()$gsva_mat, file)
    )
    output$dl_volcano <- shiny::downloadHandler(
      filename = function() paste0("gsva_volcano_", results_react()$contrast_label, ".pdf"),
      content  = function(file) {
        df <- results_react()$tt
        df$neglog10FDR <- -log10(pmax(df$adj.P.Val, 1e-300))
        sig <- df$adj.P.Val < 0.05 & abs(df$logFC) > 0.2
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$logFC, y = .data$neglog10FDR, color = sig)) +
          ggplot2::geom_point(alpha = 0.65, size = 1.8) +
          ggplot2::scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D81B60")) +
          ggplot2::labs(title = paste0("Differential Pathway Activity: ", results_react()$contrast_label),
                        x = "logFC (GSVA score)", y = "-log10(FDR)") +
          ggplot2::theme_minimal(base_size = 12)
        grDevices::pdf(file, width = 7.5, height = 5); print(p); grDevices::dev.off()
      }
    )
    output$dl_heatmap <- shiny::downloadHandler(
      filename = function() paste0("gsva_heatmap_top50_", results_react()$contrast_label, ".pdf"),
      content  = function(file) {
        gsva_mat <- results_react()$gsva_mat
        tt <- results_react()$tt
        topN <- head(tt$Pathway, 50)
        sub  <- gsva_mat[intersect(rownames(gsva_mat), topN), , drop = FALSE]
        
        grDevices::pdf(file, width = 8.5, height = 11)
        if (nrow(sub) < 2) {
          grid::grid.newpage(); grid::grid.text("No pathways available for heatmap.")
          grDevices::dev.off(); return()
        }
        zmat <- .row_zscore(sub)
        ann  <- .build_col_annotation(results_react()$meta,
                                      results_react()$group_disp,
                                      input$ann_cols)
        hm_col <- circlize::colorRamp2(c(-2, 0, 2), c("#2b8cbe", "white", "#e34a33"))
        ComplexHeatmap::draw(
          ComplexHeatmap::Heatmap(
            zmat,
            name = "Z-score",
            col  = hm_col,
            top_annotation = ann,
            show_row_names = TRUE,
            show_column_names = FALSE,
            cluster_rows = TRUE,
            cluster_columns = TRUE
          ),
          heatmap_legend_side = "bottom",
          annotation_legend_side = "bottom"
        )
        grDevices::dev.off()
      }
    )
    output$dl_boxplots <- shiny::downloadHandler(
      filename = function() paste0("gsva_top6_boxplots_", results_react()$contrast_label, ".pdf"),
      content  = function(file) {
        gsva_mat <- results_react()$gsva_mat
        tt <- results_react()$tt
        meta <- results_react()$meta
        gv <- results_react()$group_var
        top6 <- head(tt$Pathway, 6)
        
        grDevices::pdf(file, width = 10, height = 7)
        if (!length(top6)) {
          grid::grid.newpage(); grid::grid.text("No pathways available for boxplots.")
          grDevices::dev.off(); return()
        }
        df <- t(gsva_mat[top6, , drop = FALSE]) |> as.data.frame()
        df[[gv]] <- meta[[gv]]
        df_long <- df |>
          tibble::rownames_to_column("Sample") |>
          tidyr::pivot_longer(cols = tidyselect::all_of(top6), names_to = "Pathway", values_to = "Score")
        p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[gv]], y = .data$Score, fill = .data[[gv]])) +
          ggplot2::geom_boxplot(outlier.alpha = 0.25) +
          ggplot2::facet_wrap(~ Pathway, scales = "free_y", ncol = 3) +
          ggplot2::labs(x = gv, y = "GSVA score", title = "Top Differential Pathways") +
          ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "none")+
          ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
          )
        print(p)
        grDevices::dev.off()
      }
    )
    
    invisible(list(results_react = results_react))
  })
}

# ---- helpers --------------------------------------------------------------

.map_db_to_collection <- function(db) {
  switch(db,
         "GO"="C5","KEGG"="C2","Reactome"="C2","Hallmark"="H",
         "Cancer Cell Atlas"="C4","Cancer Gene Neighbourhoods"="C4",
         "Cancer Modules"="C4","Txn Factor Targets"="C3","C2")
}

.map_db_to_subcollection <- function(db) {
  switch(db,
         "GO"="GO:BP","KEGG"="CP:KEGG_MEDICUS","Reactome"="CP:REACTOME",
         "Hallmark"=NULL,"Cancer Cell Atlas"="3CA","Cancer Gene Neighbourhoods"="CGN",
         "Cancer Modules"="CM","Txn Factor Targets"="TFT:GTRD",NULL)
}

.read_gene_sets <- function(path) {
  if (grepl("\\.gmt$", basename(path), ignore.case = TRUE)) {
    GSEABase::geneIds(GSEABase::getGmt(path))
  } else {
    obj <- readRDS(path)
    if (inherits(obj, "GeneSetCollection")) GSEABase::geneIds(obj) else obj
  }
}
