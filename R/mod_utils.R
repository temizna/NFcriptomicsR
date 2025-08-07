# === Module: mod_utils ===

#' Convert Ensembl IDs to Gene Symbols
#'
#' @description Converts Ensembl IDs to gene symbols using the appropriate OrgDb.
#' 
#' @param gene_ids Character vector of Ensembl IDs
#' @param species Character string: either "Homo sapiens" or "Mus musculus"
#'
#' @return A named character vector with Ensembl IDs as names and gene symbols as values
#' @importFrom stats setNames
#' @importFrom shinythemes shinytheme
#' @export
convert_ensembl_to_symbol <- function(gene_ids, species = "Homo sapiens") {
  org_db <- get_orgdb(species)
  id_type <- "ENSEMBL"
  res <- tryCatch({
    AnnotationDbi::select(org_db, keys = gene_ids, columns = c("SYMBOL"), keytype = id_type)
  }, error = function(e) {
    warning("Ensembl to Symbol conversion failed:", e$message)
    return(NULL)
  })
  out <- setNames(res$SYMBOL, res$ENSEMBL)
  return(out[!duplicated(names(out))])
}
#' Convert Gene Symbols to Ensembl IDs
#'
#' @param gene_symbols Character vector of gene symbols
#' @param species Character string
#' @return Named character vector
#' @export
convert_symbol_to_ensembl <- function(gene_symbols, species = "Homo sapiens") {
  org_db <- get_orgdb(species)
  id_type <- "SYMBOL"
  res <- tryCatch({
    AnnotationDbi::select(org_db, keys = gene_symbols, columns = "ENSEMBL", keytype = id_type)
  }, error = function(e) {
    warning("Symbol to Ensembl conversion failed: ", e$message)
    return(NULL)
  })
  out <- setNames(res$ENSEMBL, res$SYMBOL)
  return(out[!duplicated(names(out))])
}
# Suppress global variable notes for NSE variables
utils::globalVariables(c(
  "log2FoldChange", "padj", "label", "baseMean", "color",
  "log2FoldChange_x", "log2FoldChange_y",
  "PC1", "PC2", "Sample", "Group", "Mean", "Variance",
  "Gene_Symbol","overlapping_genes","Gene","NES","NS", "TF_name",
  "geneRatio","Description","coreSize","pc","tagList"
))

#' Check if IDs are Ensembl
#'
#' @description Checks whether a vector of gene IDs are mostly Ensembl IDs.
#'
#' @param ids Character vector of gene identifiers
#' @return Logical indicating whether majority of IDs resemble Ensembl format
#' @export
is_ensembl_id <- function(ids) {
  mean(grepl("^ENSG\\d+|^ENSMUSG\\d+", ids)) > 0.5
}

#' Get OrgDb object from species name
#'
#' @description Returns the appropriate OrgDb object for the given species.
#'
#' @param species Character string: either "Homo sapiens" or "Mus musculus"
#' @return Corresponding OrgDb object
#' @export
get_orgdb <- function(species) {
  if (species == "Homo sapiens") {
    requireNamespace("org.Hs.eg.db")
    return(org.Hs.eg.db)
  } else {
    stop("Unsupported species")
  }
}

#' Get KEGG Code from species name
#'
#' @description Returns the KEGG organism code for a given species.
#'
#' @param species Character string: either "Homo sapiens" or "Mus musculus"
#' @return Character KEGG organism code
#' @export
get_kegg_code <- function(species) {
  if (species == "Homo sapiens") {
    return("hsa")
  } else if (species == "Mus musculus") {
    return("mmu")
  } else {
    stop("Unsupported species")
  }
}
#' Get Reactome Code from species name
#'
#' @description Returns the Reactome organism code for a given species.
#'
#' @param species Character string: either "Homo sapiens" or "Mus musculus"
#' @return Character Reactome organism code
#' @export
get_reactome_code <- function(species) {
  if (species == "Homo sapiens") {
    return("human")  # Use "human" for Reactome
  } else if (species == "Mus musculus") {
    return("mouse")  # Use "mouse" for Reactome
  } else {
    return("human")  # Default to human if species is not recognized
  }
}
#' Check if IDs are gene symbols
#'
#' @description Heuristically determines whether gene IDs appear to be symbols rather than Ensembl.
#'
#' @param ids Character vector of gene identifiers
#' @return Logical indicating whether IDs are mostly gene symbols
#' @export
is_symbol <- function(ids) {
  mean(!grepl("^ENSG|^ENSMUSG", ids)) > 0.5
}

