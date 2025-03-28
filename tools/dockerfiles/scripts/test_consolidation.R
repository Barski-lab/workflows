#!/usr/bin/env Rscript
#
# Test script to validate that consolidation of common functions is working correctly
#

# Set options
options(warn = -1)
options("width" = 400)

# Load required libraries
suppressMessages({
  library(tidyverse)
  library(glue)
})

cat("============================================\n")
cat("Testing consolidation of common functions...\n")
cat("============================================\n\n")

# Helper function to source files with fallback paths
source_with_fallback <- function(filepath, absolute_path = NULL) {
  # Try absolute path first if provided
  if (!is.null(absolute_path) && file.exists(absolute_path)) {
    message(paste("Sourcing from absolute path:", absolute_path))
    return(source(absolute_path))
  }
  
  # Try relative path
  if (file.exists(filepath)) {
    message(paste("Sourcing from relative path:", filepath))
    return(source(filepath))
  }
  
  # If we get here, try a standard Docker path
  docker_path <- file.path("/usr/local/bin", filepath)
  if (file.exists(docker_path)) {
    message(paste("Sourcing from Docker path:", docker_path))
    return(source(docker_path))
  }
  
  # If all fails, error
  stop(paste("Could not find file to source:", filepath))
}

# Load common functions
cat("1. Loading common utility functions...\n")
tryCatch({
  source_with_fallback("functions/common/utilities.R", "/usr/local/bin/functions/common/utilities.R")
  cat("✓ Utilities loaded successfully\n")
}, error = function(e) {
  cat("❌ Error loading utilities.R:", e$message, "\n")
})

cat("\n2. Loading common visualization functions...\n")
tryCatch({
  source_with_fallback("functions/common/visualization.R", "/usr/local/bin/functions/common/visualization.R")
  cat("✓ Visualization functions loaded successfully\n")
}, error = function(e) {
  cat("❌ Error loading visualization.R:", e$message, "\n")
})

cat("\n3. Loading common clustering functions...\n")
tryCatch({
  source_with_fallback("functions/common/clustering.R", "/usr/local/bin/functions/common/clustering.R")
  cat("✓ Clustering functions loaded successfully\n")
}, error = function(e) {
  cat("❌ Error loading clustering.R:", e$message, "\n")
})

cat("\n4. Loading common export functions...\n")
tryCatch({
  source_with_fallback("functions/common/export_functions.R", "/usr/local/bin/functions/common/export_functions.R")
  cat("✓ Export functions loaded successfully\n")
}, error = function(e) {
  cat("❌ Error loading export_functions.R:", e$message, "\n")
})

# Test function existence
cat("\n5. Checking for existence of key functions...\n")
common_functions <- c(
  "validate_metadata", 
  "apply_batch_correction", 
  "get_clustered_data", 
  "create_mds_plot", 
  "export_gct_data"
)

for (func in common_functions) {
  if (exists(func)) {
    cat(glue("✓ Function {func} exists\n"))
  } else {
    cat(glue("❌ Function {func} NOT found\n"))
  }
}

cat("\n============================================\n")
cat("Consolidation test completed!\n")
cat("============================================\n") 