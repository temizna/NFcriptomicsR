Package: NFcriptomicsR 
Title: A Comprehensive Shiny App for RNA-Seq Analysis of  Neurofibroma/MPNST/NF1/SUZ12/EED cell lines and primary samples
Version: 0.5.0
Authors@R: 
    person(Nuri Alpay, Temiz, email = temizna@umn.edu, role = c("aut", "cre"))
Description: 
    NFcriptomicsR provides an interactive Shiny interface for complete RNA-Seq analysis of 219 breast cancer cell lines,
    including data normalization, quality control, differential expression,
    gene expression visualization, pathway analysis, and GSEA using CRAN and Bioconductor tools.
License: MIT
Encoding: UTF-8
LazyData: true
RoxygenNote: 7.2.3

Depends:
    R (>= 4.1.0)

Imports:
    shiny,
    shinythemes,
    ggplot2,
    dplyr,
    tidyr,
    stringr,
    readr,
    tibble,
    reshape2,
    DT,
    RColorBrewer,
    rDESeq2,
    Biobase,
    clusterProfiler,
    enrichplot,
    ComplexHeatmap,
    msigdbr,
    ReactomePA,
    pathview,
    AnnotationDbi,
    org.Hs.eg.db,
    org.Mm.eg.db

Suggests:
    testthat,
    knitr,
    rmarkdown

VignetteBuilder: knitr

---

# NFcriptomicsR 

**NFcriptomicsR** is a comprehensive Shiny-based application for end-to-end RNA-Seq data analysis of Neurofibroma. The data comes preloaded. It provides an interactive GUI for both novice and advanced users to perform quality control, normalization, differential expression, pathway analysis, and gene set enrichment analysis (GSEA) with no coding experience.
## ✨ Features

- **Data Input**  
  - Pre-loaded RNA-Seq data.  

- **Sample Selection**  
  - Interactive filtering of samples for downstream analysis

- **Gene Expression Visualization**  
  - Boxplots of selected genes grouped by metadata categories
  - Table of one-way ANOVA comparison between groups for selected genes followed by Tukey's HSD test.

- **Quality Control (QC)**  
  - PCA, sample distance heatmaps, mean-variance plots, variance histograms

- **Genes of interest Heatmap Visualization**
  - Gene expression heatmap of uploaded genes of interest grouped by metadata categories
- **Differential Expression Analysis**  
  - DESeq2-based analysis with customizable thresholds and conditions

- **Heatmap**  
  - Top DE genes with hierarchical clustering and group annotations

- **Volcano and MA Plots**  
  - Visualization of differential expression results with gene labeling

- **Cross Plot**  
  - Compare DE results between two conditions or experiments
  - Cross plot of DE genes, venn diagram and heatmap of overlaps and table of overlaps
- **Pathway Analysis**  
  - GO, KEGG, Reactome enrichment using clusterProfiler  
  - Visualizations: dot, cnet, circular, emap, heatmap, tree, upset plots  
  - KEGG pathway download (no-rendering)  

- **GSEA (Gene Set Enrichment Analysis)**  
  - Supports MSigDB collections: Hallmark, GO, KEGG, Reactome, Cancer Cell Atlas, Cancer Modules, Txn factor Targets  
  - Dot plots, enrichment plot, upset plot  and enrichment tables
- **Non-overlap Pathway Analysis**
  - Overlapping genes from the initial selected pathway analysis is removed and same pathway analysis is re-run. 
  - Visualizations: Dot, tree, heatmap plots and pathway tables 
- **Transcription Factor Enrichment**
  - Supports TRANSFAC and JASPAR, ENCODE, TRRUST, TF Perturbations followed by Expression, hFTarget, TFLink
  - supports GSEA or over representation analysis
  - Enrichment dot plots and enrich term tables , ridgeplot with GSEA option
- **PCA Analysis (Simple Consensus Clustering)**
  - Principle components covering user selected variance is used to create sample clusters, select gene features and comparison pathway analysis
  - Visualizations: Contributing genes heatmap, sample correlation heatmap, dot plot of enrichment comparison between components, enrichment and contributing genes tables
- **Session Log**
  - Live session log that logs all variable selections in order of clicks and their related action buttons (i.e. run differential expression)
## 🧬 Supported Species

- Homo sapiens


## 📦 Installation

Install required CRAN and Bioconductor dependencies before installing NFcriptomicsR:

```r
# Install Bioconductor manager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c(
  "DESeq2", "GEOquery", "Biobase", "clusterProfiler", "enrichplot",
  "ComplexHeatmap", "msigdbr", "ReactomePA", "pathview",
  "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi"
))

#install devtools if not already installed
install.packages("devtools")
library("devtools")
# Clone or download this repo, then from root directory:
devtools::install_github("temizna/NFcriptomicsR", force=TRUE)
```

## 🚀 Running the App

After installation, launch the app from R or RStudio:

```r
library(NFcriptomicsR)
NFcriptomicsR::run_app()
```

This will open the app in your default web browser.

## 📂 File Requirements
Pre-loaded data. The data is normalized as below:
countData <- read.csv(paste(datadir,"subread_counts_gene_symbol2.csv",sep=""), row.names=1,check.names = F)
coldata<-read.csv(paste(datadir,"design_stemshiny.csv",sep=""),header=T, row.names = 1,check.names = F)
countData2<-countData
temp<-rownames(coldata)[rownames(coldata) %in% colnames(countData)]
coldata2=coldata[temp,]
countData=countData2[,rownames(coldata2)]
batch=coldata[colnames(countData),3]
rms<-which(is.na(batch))
#countData=countData[,-rms]
temp <- countData[rowSums(countData) >= 10, ]
countData <- temp[apply(temp, 1, var) > 0.1, ]

batch=coldata[colnames(countData),3]
group=coldata[colnames(countData),2]
adjusted <- ComBat_seq(as.matrix(countData), batch=batch, group=group, full_mod = T)
dds <- DESeqDataSetFromMatrix(countData = adjusted, colData = coldata2, design = ~ treatment+dimension+cellline+PR+ER)
dds <- DESeq(dds)
norm_data_vst <- DESeq2::vst(dds, blind = FALSE)
norm_data<-assay(norm_data_vst)
save(adjusted, metadata, norm_data, file = "inst/extdata/processed_cellline_data.rda")
saveRDS(dds, file = "inst/extdata/processed_dds.rds", compress = "xz")
## 🛠 Development

All modules are structured as reusable Shiny modules. Additional features can be added easily. Utility functions are separated into a `utils.R` file.

To ensure all required Bioconductor packages are installed:

```r
# R/zzz.R
.onLoad <- function(libname, pkgname) {
  required_bioc <- c(
    "DESeq2", "GEOquery", "Biobase", "clusterProfiler", "enrichplot",
    "ComplexHeatmap", "msigdbr", "ReactomePA", "pathview",
    "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi"
  )

  missing_pkgs <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    message("Installing missing Bioconductor packages: ", paste(missing_pkgs, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
  }
}
```

Also ensure the following are loaded in your `run_app()` function:

```r
run_app <- function() {
  options(shiny.maxRequestSize = 250 * 1024^2)

  library(shiny)
  library(shinythemes)
  library(DESeq2)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(GEOquery)
  library(pathview)
  library(ComplexHeatmap)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(msigdbr)

  shiny::shinyApp(
    ui = NFcriptomicsR::app_ui(),
    server = NFcriptomicsR::app_server
  )
}
```


## 📝 License

MIT © Nuri Alpay Temiz

