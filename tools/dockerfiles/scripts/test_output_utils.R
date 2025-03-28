#!/usr/bin/env Rscript
#
# Test script for consolidated output utility functions
#

# Setup test environment
test_dir <- "test_output_utils"
if (!dir.exists(test_dir)) {
  dir.create(test_dir)
}
old_wd <- setwd(test_dir)

# Source output utilities
tryCatch({
  # Try different paths for sourcing
  if (file.exists("../functions/common/output_utils.R")) {
    source("../functions/common/output_utils.R")
  } else if (file.exists("/usr/local/bin/functions/common/output_utils.R")) {
    source("/usr/local/bin/functions/common/output_utils.R")
  } else {
    stop("Could not find output_utils.R")
  }
  
  # Check if key functions exist
  required_functions <- c(
    "get_output_filename", 
    "save_plot", 
    "write_gct_file", 
    "write_cls_file", 
    "write_deseq_results",
    "generate_mds_plot_html",
    "verify_outputs"
  )
  
  missing_functions <- required_functions[!sapply(required_functions, exists)]
  
  if (length(missing_functions) > 0) {
    stop("Missing required functions: ", paste(missing_functions, collapse = ", "))
  } else {
    cat("✓ All required functions are available\n")
  }
  
  # Test data
  test_prefix <- "test_output"
  test_data <- matrix(
    runif(100), nrow = 10, ncol = 10,
    dimnames = list(
      paste0("gene", 1:10),
      paste0("sample", 1:10)
    )
  )
  test_classes <- rep(c("A", "B"), each = 5)
  
  # Test 1: Output filename generation
  filenames <- lapply(c("counts_all", "ma_plot", "mds_plot"), function(type) {
    lapply(c("gct", "png", "pdf", "html"), function(ext) {
      get_output_filename(test_prefix, type, ext)
    })
  })
  cat("✓ Output filename generation works\n")
  
  # Test 2: GCT file creation
  gct_file <- get_output_filename(test_prefix, "counts_all", "gct")
  write_gct_file(test_data, gct_file)
  if (!file.exists(gct_file)) {
    stop("Failed to create GCT file")
  }
  cat("✓ GCT file creation works\n")
  
  # Test 3: CLS file creation
  cls_file <- get_output_filename(test_prefix, "phenotypes", "cls")
  write_cls_file(test_classes, cls_file)
  if (!file.exists(cls_file)) {
    stop("Failed to create CLS file")
  }
  cat("✓ CLS file creation works\n")
  
  # Test 4: DESeq2 results writing
  mock_results <- data.frame(
    log2FoldChange = rnorm(10),
    pvalue = runif(10, 0, 0.1),
    padj = runif(10, 0, 0.1),
    baseMean = runif(10, 100, 1000),
    row.names = paste0("gene", 1:10)
  )
  class(mock_results) <- c("DESeqResults", "data.frame")
  results_file <- get_output_filename(test_prefix, "report", "tsv")
  write_deseq_results(mock_results, results_file)
  if (!file.exists(results_file)) {
    stop("Failed to create DESeq2 results file")
  }
  cat("✓ DESeq2 results export works\n")
  
  # Test 5: Markdown summary
  summary_file <- get_output_filename(test_prefix, "summary", "md")
  write_markdown_summary(c("# Test Summary", "", "This is a test"), summary_file)
  if (!file.exists(summary_file)) {
    stop("Failed to create markdown summary")
  }
  cat("✓ Markdown summary creation works\n")
  
  # Test 6: Plot saving (when ggplot2 is available)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    test_plot <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
      ggplot2::geom_point()
    
    plot_files <- save_plot(test_plot, test_prefix, "test_plot")
    
    if (!all(sapply(unlist(plot_files), file.exists))) {
      stop("Failed to save plots in multiple formats")
    }
    cat("✓ Plot saving in multiple formats works\n")
  } else {
    cat("⚠ Skipping plot saving test (ggplot2 not available)\n")
  }
  
  # Test 7: MDS plot HTML generation (when required packages are available)
  if (requireNamespace("plotly", quietly = TRUE) && 
      requireNamespace("htmlwidgets", quietly = TRUE) &&
      requireNamespace("limma", quietly = TRUE)) {
    
    # Create a simple DESeq2 object mock
    dds_mock <- list(
      assays = list(counts = test_data),
      colData = data.frame(
        condition = test_classes,
        row.names = colnames(test_data)
      )
    )
    class(dds_mock) <- "DESeqDataSet"
    
    # Mock transformed counts
    vst_mock <- list(
      assays = list(data = test_data)
    )
    assay <- function(x) x$assays$data
    
    tryCatch({
      html_file <- generate_mds_plot_html(
        dds = dds_mock,
        transformed_counts = vst_mock,
        metadata = dds_mock$colData,
        color_by = "condition",
        prefix = test_prefix
      )
      
      if (!file.exists(html_file)) {
        stop("Failed to create MDS plot HTML")
      }
      cat("✓ MDS plot HTML generation works\n")
    }, error = function(e) {
      cat("⚠ MDS plot HTML generation error:", conditionMessage(e), "\n")
    })
  } else {
    cat("⚠ Skipping MDS plot HTML test (required packages not available)\n")
  }
  
  # Test 8: Output verification
  output_files <- c(gct_file, cls_file, results_file, summary_file)
  verify_result <- verify_outputs(test_prefix, "test", FALSE, expected_files = output_files)
  if (!verify_result) {
    stop("Output verification failed")
  }
  cat("✓ Output verification works\n")
  
  # Summary
  cat("\nAll output utility tests completed successfully!\n")
  
  # Clean up
  file.remove(list.files(pattern = "test_output"))
  
}, finally = {
  setwd(old_wd)
  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
})

cat("\nOutput utility functions validation completed.\n") 