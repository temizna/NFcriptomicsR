.onLoad <- function(libname, pkgname) {
  bioc_required <- c(
    "DESeq2", "org.Hs.eg.db","clusterProfiler", 
    "enrichplot",  
    "ComplexHeatmap", "msigdbr", "ReactomePA"
  )

  cran_required <- c("ggplot2", "reshape2", "dplyr", "tibble", "stringr", 
                     "readr", "DT", "RColorBrewer", "shiny", "shinythemes")

  install_if_missing <- function(pkgs, installer) {
    for (pkg in pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        installer(pkg)
      }
    }
  }

  # Install CRAN packages if needed
  install_if_missing(cran_required, utils::install.packages)

  # Install Bioconductor packages if needed
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  install_if_missing(bioc_required, BiocManager::install)
}

