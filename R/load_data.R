# mod_loaded_data_loader.R
#' Load Preloaded loaded_data Module
#'
#' This module loads the preprocessed RNASeq loaded_data (counts, samples, normalized counts, and DESeq2 object)
#' from `.rda` files located in the `inst/extloaded_data/` directory. It assigns the loaded loaded_data to the reactive values.
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param dds_rv A reactiveVal container for storing the DESeq2 loaded_dataset.
#' @param loaded_data_rv A list  where the loaded_data will be stored. It should contain the following elements:
#'   - `counts`: Raw count loaded_data.
#'   - `samples`: Sample metaloaded_data.
#'   - `norm_counts`: Normalized counts (rlog-transformed).
#'   - `dds`: The DESeq2 loaded_dataset object.
#' @return A list with elements: counts, samples, norm_counts, dds
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#' @importFrom shiny observeEvent observe reactive
#' @export
load_preloaded_data <- function(input, output, session, loaded_data_rv, dds_rv) {
  # Load preprocessed loaded_data from the inst/extloaded_data folder
  rda_path <- system.file("extdata", "processed_cellline_data.rda", package = "NFcriptomicsR")
#  dds_path <- system.file("extdata", "processed_dds.rds", package = "NFcriptomicsR")
  
  # Check if the paths exist
  if (!file.exists(rda_path)) {
    stop("Processed loaded_data not found!")
  }
 # if (!file.exists(dds_path)) {
#    stop("Processed DESeq2 object not found!")
 # }
  
  # Load the loaded_data and DESeq2 object outside the observer (only once)
  load(rda_path, envir = .GlobalEnv)
 # print(ls())  # Debugging check: check if the loaded_data was loaded successfully
  

 # print(paste("Loaded DESeq2 object:", ls()))  # Debugging check for dds
  
  # Extract loaded_data
  raw_counts <- get("adjusted", envir = .GlobalEnv)
  metadata <- get("metadata", envir = .GlobalEnv)
  norm_data <- get("norm_data", envir = .GlobalEnv)
  #dds <- get("dds", envir = .GlobalEnv)  # The DESeq2 object
  # Create DDS
  cat("Row names of metadata:\n")
  print(head(rownames(metadata)))
  cat("Column names of raw_counts:\n")
  print(head(colnames(raw_counts)))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design =  ~ sample_type)
  dds <- DESeq2::estimateSizeFactors(dds)
  # Ensure 'dds' object is loaded correctly
  if (is.null(dds)) {
    stop("The 'dds' object was not loaded properly.")
  }
  
  # Store the DESeq2 object in the reactive value container
  dds_rv(dds)
  #print(head(metadata))
  # Assign the loaded loaded_data to the loaded_data object
  loaded_data <- list(
    counts = raw_counts,
    samples = metadata,
    norm_counts = norm_data,
    species = "Homo sapiens"
  )
  
  # Check the structure of loaded_data
  #print(str(loaded_data))
  loaded_data_rv(loaded_data)
}
