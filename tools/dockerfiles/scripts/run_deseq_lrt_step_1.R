#!/usr/bin/env Rscript
options(warn = -1)
options("width" = 400)
# Set a global error handler that aborts the script on any error
options(error = function() {
  message("An unexpected error occurred. Aborting script.")
  quit(save = "no", status = 1, runLast = FALSE)
})


suppressMessages(library(argparse)) # For argument parsing
suppressMessages(library(BiocParallel)) # For parallel processing
suppressMessages(library(pheatmap)) # For heatmap
suppressMessages(library(DESeq2)) # For DESeq2
suppressMessages(library(tidyverse)) # For data manipulation
suppressMessages(library(Glimma)) # For glimmaMDS
suppressMessages(library(cmapR)) # For write_gct
suppressMessages(library(sva)) # For ComBat_seq
suppressMessages(library(limma)) # For removeBatchEffect
suppressMessages(library(hopach)) # For clustering
suppressMessages(library(rlang)) # For traceback handling
suppressMessages(library(stringr)) # For string manipulation
suppressMessages(library(RColorBrewer)) # For color palettes
suppressMessages(library(gridExtra)) # For combining plots
suppressMessages(library(glue)) # For string interpolation

options(rlang_backtrace_on_error = "full")

mutate <- dplyr::mutate
filter <- dplyr::filter
group_by <- dplyr::group_by
slice <- dplyr::slice
rename <- dplyr::rename
select <- dplyr::select
`%>%` <- magrittr::`%>%`
`%in%` <- base::`%in%`

###### v0.0.6 ######
#
# Changes:
# - Added function to save error_message, work to catch the following commonplace errors right now:
#   - Mismatch between sample names and metadata sample names [testing]
#   - Non-full-rank model provided [testing]
#   - Missing "batch" column [testing]
#   - Wrong design-formula syntax [testing]
#   - 0 genes after RPKM filtering [testing]
#   - Incorrect metadata file format [testing]
#   - Wrong cluster parametrs
#   - Insufficient Biological Replicates
#   - Mismatch between design-formula covars and metadata columns [testing]
# 
###### v0.0.5 ######
#
# Changes:
# - 'combatseq' batch correction is applied in Step 1 if specified.
# - 'model' is noted and passed to Step 2 for application.
# - Batch correction options are handled appropriately.
# - Necessary information is saved for Step 2.
#
###### v0.0.4 ######
#
# Changes:
# - Adjusted batch correction to be optional.
# - Batch correction with limma is only applied if specified.
#
###### v0.0.3 ######
#
# Changes:
# - Reads preprocessed data from RDS files.
# - Accepts `contrast_indices`, `fdr`, `lfcThreshold`, and `regulation` as inputs.
# - Applies batch correction using `limma::removeBatchEffect` if specified.
# - Uses `results` function with appropriate parameters.
#
###### v0.0.2 ######
#
# Changes:
# - Save DESeq2 dataset objects (dds) and LRT results (dsq_lrt_res) into RDS files.
# - Save expression data (expression_data_df) into an RDS file.
# - Save metadata (metadata_df) into an RDS file.
# - Adjust outputs to include these RDS files.
#
###### v0.0.1 ######
#
# Note: at least two biological replicates are required for every compared category.
#
# All input CSV/TSV files should have the following header (case-sensitive)
# <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
# <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV
#
# Format of the input files is identified based on file's extension
# *.csv - CSV
# *.tsv - TSV
# Otherwise used CSV by default
#
# The output file's rows order corresponds to the rows order of the first CSV/TSV file.
# Output file is always saved in TSV format
#
# Output file includes only intersected rows from all input files. Intersected by
# RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand
#
# Additionally we calculate -LOG10(pval) and -LOG10(padj)
#
# Example of CSV metadata file set with --meta
#
# ,time,condition
# DH1,day5,WT
# DH2,day5,KO
# DH3,day7,WT
# DH4,day7,KO
# DH5,day7,KO
#
# where time, condition, day5, day7, WT, KO should be a single words (without spaces)
# and DH1, DH2, DH3, DH4, DH5 correspond to the --names (spaces are allowed)
#
# --contrast should be set based on your metadata file in a form of Factor Numerator Denominator
# where Factor      - columns name from metadata file
#       Numerator   - category from metadata file to be used as numerator in fold change calculation
#       Denominator - category from metadata file to be used as denominator in fold change calculation
# for example condition WT KO
# if --contrast is set as a single string "condition WT KO" then is will be splitted by space
#

# --- Added error reporting function ---
report_error <- function(msg, details = NULL, recommendations = NULL) {
  # Format the timestamp in a friendly, readable way
  timestamp <- format(Sys.time(), "%A, %B %d, %Y %I:%M:%S %p")
  
  # Construct the Markdown report
  md_msg <- paste0(
    "# ❌ Error Report\n\n",
    "**Timestamp:** ", timestamp, "\n\n",
    "**Error Message:** ", msg, "\n\n"
  )
  
  # Add details section if provided
  if (!is.null(details)) {
    md_msg <- paste0(
      md_msg,
      "## Details\n\n",
      "```\n",
      paste(details, collapse = "\n"),
      "\n```\n\n"
    )
  }
  
  # Add recommendations section if provided
  if (!is.null(recommendations)) {
    recs <- paste0("- ", paste(recommendations, collapse = "\n- "))
    md_msg <- paste0(
      md_msg,
      "## Recommendations\n\n",
      recs, "\n"
    )
  }
  
  # Write the Markdown message to a file
  writeLines(md_msg, con = "error_report.txt")
  
  # Stop execution by raising an error with the given message
  stop(msg)
}


READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")

get_file_type <- function(filename) {
  ext <- tools::file_ext(filename)
  separator <- ","
  if (ext == "tsv") {
    separator <- "\t"
  }
  return(separator)
}

load_expression_data <- function(filenames,
                                 prefixes,
                                 read_colname,
                                 rpkm_colname,
                                 intersect_by) {
  collected_isoforms <- NULL
  for (i in 1:length(filenames)) {
    isoforms <- read.table(
      filenames[i],
      sep = get_file_type(filenames[i]),
      header = TRUE,
      stringsAsFactors = FALSE
    )
    # FOR TEST ONLY TO REDUCE NUMBER OF ROWS AND SPEED UP TESTING
    if (!is.null(args$test_mode) && args$test_mode) {
      print("Test mode is ON, each sample will be limited to 500 rows")
      isoforms <- isoforms %>%
        dplyr::slice_head(n = 500)
    } else {
      print("Test mode is OFF, processing all rows")
    }

    print(paste0("Load ", nrow(isoforms), " rows from ", filenames[i]))
    colnames(isoforms)[colnames(isoforms) == read_colname] <- paste(prefixes[i], read_colname, sep = " ")
    colnames(isoforms)[colnames(isoforms) == rpkm_colname] <- paste(prefixes[i], rpkm_colname, sep = " ")
    if (is.null(collected_isoforms)) {
      collected_isoforms <- isoforms
    } else {
      collected_isoforms <- dplyr::full_join(collected_isoforms, isoforms, by = intersect_by)
    }
  }
  print(paste0("Number of rows common for all loaded files ", nrow(collected_isoforms)))

  return(collected_isoforms)
}

assert_args <- function(args) {
  if (length(args$input) != length(args$name)) {
    print("Exiting: --input and --name have different number of values")
    quit(
      save = "no",
      status = 1,
      runLast = FALSE
    )
  }
  tryCatch(
    expr = {
      # Try to load design formula
      design_formula <- as.formula(tolower(args$design))
    },
    error = function(e) {
      print(paste0("Exiting: failed to load --design ", tolower(args$design), " as formula"))
      quit(
        save = "no",
        status = 1,
        runLast = FALSE
      )
    }
  )
  tryCatch(
    expr = {
      # Try to load reduced formula
      reduced_formula <- as.formula(tolower(args$reduced))
    },
    error = function(e) {
      print(paste0("Exiting: failed to load --reduced ", tolower(args$reduced), " as formula"))
      quit(
        save = "no",
        status = 1,
        runLast = FALSE
      )
    }
  )
  return(args)
}


create_releveled_deseq_object <- function(dds, relevel_list) {
  # Rebuild DESeqDataSet to ensure a proper copy
  dds_copy <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts(dds),
    colData   = as.data.frame(colData(dds)),
    design    = design(dds)
  )
  # Apply releveling for each factor as specified in relevel_list
  for (fac in names(relevel_list)) {
    new_ref <- relevel_list[[fac]]
    dds_copy[[fac]] <- relevel(dds_copy[[fac]], ref = new_ref)
  }
  return(dds_copy)
}

# Function to check if the design-formula syntax correct
check_design_formula <- function(formula_string, type = "Design") {
  log_message(glue::glue("Validating {type} formula syntax: '{formula_string}'"), "STEP")
  
  # Try to parse the formula
  tryCatch(
    {
      f <- as.formula(formula_string)
      
      # Additional validation for formula terms
      terms <- attr(terms(f), "term.labels")
      debug_log(glue::glue("Formula terms: {paste(terms, collapse=', ')}"))
      
      # Check for interactions
      has_interactions <- any(grepl(":", terms, fixed = TRUE))
      log_message(glue::glue("Formula contains interaction terms: {has_interactions}"), "INFO")
      
      return(f)
    },
    error = function(e) {
      log_message(glue::glue("Formula parsing error: {e$message}"), "ERROR")
      report_error(
        msg = paste(type, "Formula Parsing Error"),
        details = c(
          paste0("Failed to parse the ", type, " formula: '", formula_string, "'."),
          "Ensure the formula starts with '~' and uses '+' (as independent), ':' (for interaction) or '*' (as + and : at once) to add covariates. E.g. '~factor1 + factor2 + factor1:factor2' (or just '~factor1 * factor2'). "
        ),
        recommendations = c("Please check the syntax of your formula and try again.")
      )
    }
  )
}

# Function to check if the design-formula covariates are the same as colnames in metdata
check_design_covariates <- function(design_formula, metadata_df) {
  # Extract covariate names from the design formula
  covariates <- all.vars(design_formula)
  
  # Identify missing covariates
  missing_covariates <- setdiff(covariates, colnames(metadata_df))
  
  if (length(missing_covariates) > 0) {
    report_error(
      msg = "Missing Covariates in Metadata",
      details = paste("The following covariates are missing in metadata:", paste(missing_covariates, collapse = ", ")),
      recommendations = c("Ensure metadata contains all covariates specified in the design formula.")
    )
  } else {
    message("All covariates in the design formula are present in the metadata.")
  }
}

# Function to apply batch correction using ComBat_seq
apply_combat_seq <- function(count_data, design, column_data) {

  start_time <- proc.time()

  print("Applying ComBat_seq batch correction")

  print("Checking full-rank of the model...")
  validate_full_rank_design(design, column_data)

  corrected_counts <- sva::ComBat_seq(
    counts = as.matrix(count_data),
    batch = column_data$batch,
    covar_mod = model.matrix(design, data = column_data)
  )

  end_time <- proc.time() - start_time

  print("Total time of apply_combat_seq execution: ")
  print(end_time) # This shows user, system, and elapsed times

  return(corrected_counts)
}

# Enhanced batch correction handling with detailed logging
handle_batch_correction <- function(args, countData, metadata_df, design_formula) {
  log_message("Handling batch correction setup...", "STEP")
  
  # Store original values
  original_counts <- NULL
  corrected_counts <- NULL
  batch_warning <- NULL
  original_design_formula <- design_formula
  
  # Check if batch correction is requested
  if (args$batchcorrection == "none") {
    log_message("No batch correction requested", "INFO")
    return(list(
      countData = countData,
      design_formula = design_formula,
      original_counts = NULL,
      corrected_counts = NULL,
      batch_warning = NULL
    ))
  }
  
  log_message(glue::glue("Batch correction method requested: {args$batchcorrection}"), "INFO")
  
  # Check if 'batch' column exists
  if (!("batch" %in% colnames(metadata_df))) {
    batch_warning <- "⚠️ You specified batch correction but there's no 'batch' column in your metadata. Proceeding without batch correction."
    log_message(batch_warning, "WARNING")
    return(list(
      countData = countData,
      design_formula = design_formula,
      original_counts = NULL,
      corrected_counts = NULL,
      batch_warning = batch_warning
    ))
  }
  
  # Ensure 'batch' is a factor
  metadata_df$batch <- as.factor(metadata_df$batch)
  batch_levels <- levels(metadata_df$batch)
  log_message(glue::glue("Batch factor has {length(batch_levels)} levels: {paste(batch_levels, collapse=', ')}"), "INFO")
  
  # Check if we have sufficient samples in each batch
  batch_counts <- table(metadata_df$batch)
  debug_log("Samples per batch:", batch_counts)
  
  # Apply appropriate batch correction
  if (args$batchcorrection == "combatseq") {
    log_message("Applying ComBat-seq batch correction to count data...", "STEP")
    
    # Check if any batch has only one sample
    if (any(batch_counts < 2)) {
      log_message("Warning: Some batches have fewer than 2 samples, which might affect batch correction", "WARNING")
    }
    
    original_counts <- countData
    
    # Timing the ComBat-seq operation
    start_time <- proc.time()
    corrected_counts <- apply_combat_seq(countData, design_formula, metadata_df)
    end_time <- proc.time() - start_time
    
    log_message(glue::glue("ComBat-seq batch correction completed in {round(end_time['elapsed'], 2)} seconds"), "SUCCESS")
    
    # Compare before and after statistics
    debug_log("Original counts summary statistics:")
    print(summary(as.vector(as.matrix(original_counts))))
    
    debug_log("Corrected counts summary statistics:")
    print(summary(as.vector(as.matrix(corrected_counts))))
    
    countData <- corrected_counts
    log_message("ComBat-seq batch correction applied. Original design formula maintained.", "SUCCESS")
    
  } else if (args$batchcorrection == "model") {
    log_message("Including 'batch' in the design formula for limma batch correction...", "STEP")
    
    # Show the original design formula
    debug_log(glue::glue("Original design formula: {deparse(design_formula)}"))
    
    # Remove '~' from the original design formula string and add batch
    original_design <- substring(deparse(design_formula), 2)
    new_design_formula <- as.formula(paste("~ batch +", original_design))
    design_formula <- new_design_formula
    
    log_message(glue::glue("Modified design formula to include batch: {deparse(design_formula)}"), "SUCCESS")
  }
  
  return(list(
    countData = countData,
    design_formula = design_formula,
    original_counts = original_counts,
    corrected_counts = corrected_counts,
    batch_warning = batch_warning
  ))
}

get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq2 for multi-factor analysis using LRT (likelihood ratio or chi-squared test)")
  parser$add_argument(
    "-i",
    "--input",
    help = "Grouped by gene / TSS/ isoform expression files, formatted as CSV/TSV",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-n",
    "--name",
    help = "Unique names for input files, no special characters, spaces are allowed. Number and order corresponds to --input",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-m",
    "--meta",
    help = "Metadata file to describe relation between samples, where first column corresponds to --name, formatted as CSV/TSV",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-d",
    "--design",
    help = "Design formula. Should start with ~, like ~condition+celltype+condition:celltype. See DESeq2 manual for details",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-r",
    "--reduced",
    help = "Reduced formula to compare against with the term(s) of interest removed. Should start with ~. See DESeq2 manual for details",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "--batchcorrection",
    help = paste(
      "Specifies the batch correction method to be applied.",
      "- 'combatseq' applies ComBat_seq at the beginning of the analysis, removing batch effects from the counts before differential expression analysis.",
      "- 'model' applies removeBatchEffect from the limma package after differential expression analysis.",
      "- Default: none"
    ),
    type = "character",
    choices = c("none", "combatseq", "model"),
    default = "none"
  )
  parser$add_argument(
    "--fdr",
    help = paste(
      "In the exploratory visualization part of the analysis output only features",
      "with adjusted p-value (FDR) not bigger than this value. Also the significance",
      "cutoff used for optimizing the independent filtering. Default: 0.1."
    ),
    type = "double",
    default = 0.1
  )
  parser$add_argument(
    "--lfcthreshold",
    help = paste(
      "Log2 fold change threshold for determining significant differential expression.",
      "Genes with absolute log2 fold change greater than this threshold will be considered.",
      "Default: 0.59 (about 1.5 fold change)"
    ),
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = paste(
      "Flag to indicate whether to use lfcthreshold as the null hypothesis value in the results function call.",
      "If TRUE, lfcthreshold is used in the hypothesis test (i.e., genes are tested against this threshold).",
      "If FALSE, the null hypothesis is set to 0, and lfcthreshold is used only as a downstream filter.",
      "Default: FALSE"
    ),
    action = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--rpkm_cutoff",
    help = paste(
      "RPKM cutoff for filtering genes. Genes with RPKM values below this threshold will be excluded from the analysis.",
      "Default: NULL (no filtering)"
    ),
    type    = "integer",
    default = NULL
  )
  parser$add_argument(
    "--scaling_type",
    help = paste(
      "Specifies the type of scaling to be applied to the expression data.",
      "- 'minmax' applies Min-Max scaling, normalizing values to a range of [-2, 2].",
      "- 'zscore' applies Z-score standardization, centering data to mean = 0 and standard deviation = 1.",
      "- Default: none (no scaling applied)."
    ),
    type = "character",
    choices = c("minmax", "zscore"),
    default = "zscore"
  )
  parser$add_argument(
    "--rowdist",
    help = paste(
      "Distance metric for HOPACH row clustering. Ignored if --cluster is not",
      "provided. Default: cosangle"
    ),
    type = "character",
    default = "cosangle",
    choices = c(
      "cosangle",
      "euclid",
      "abseuclid",
      "cor",
      "abscor"
    )
  )
  parser$add_argument(
    "--columndist",
    help = paste(
      "Distance metric for HOPACH column clustering. Ignored if --cluster is not",
      "provided. Default: euclid"
    ),
    type = "character",
    default = "euclid",
    choices = c(
      "cosangle",
      "euclid",
      "abseuclid",
      "cor",
      "abscor"
    )
  )
  parser$add_argument(
    "--cluster",
    help = paste(
      "Hopach clustering method to be run on normalized read counts for the",
      "exploratory visualization part of the analysis. Default: do not run",
      "clustering"
    ),
    type = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
  )
  parser$add_argument(
    "--k",
    help    = "Number of levels (depth) for Hopach clustering: min - 1, max - 15. Default: 3.",
    type    = "integer",
    default = 3
  )
  parser$add_argument(
    "--kmax",
    help    = "Maximum number of clusters at each level for Hopach clustering: min - 2, max - 9. Default: 5.",
    type    = "integer",
    default = 5
  )
  parser$add_argument(
    "-o",
    "--output",
    help = "Output prefix for generated files",
    type = "character",
    default = "./deseq"
  )
  parser$add_argument(
    "-p",
    "--threads",
    help = "Threads number",
    type = "integer",
    default = 1
  )
  parser$add_argument(
    "--lrt_only_mode",
    help    = "Run LRT only, no contrasts",
    action  = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--test_mode",
    help = "Run for test, only first 500 rows, clustering minimised",
    action = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--clean_logs",
    help = "Clean existing log files before starting new analysis",
    action = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--debug",
    help = "Enable debug mode for more verbose logging",
    action = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--save_rdata",
    help = "Save intermediate RData files for debugging",
    action = "store_true",
    default = FALSE
  )
  args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
  return(args)
}

generate_lrt_md <- function(deseq_results, full_formula, reduced_formula, output_file, alpha = args$fdr,
                            batch_warning = NULL, rpkm_filtered_count = NULL) {
  # Initialize the markdown content
  md_content <- ""

  # Add batch warning if present
  if (!is.null(batch_warning)) {
    md_content <- paste0(md_content, "# **Warning**\n\n", batch_warning, "\n\n---\n\n")
  }

  # Add DESeq LRT results summary based on padj value only
  if (!is.null(deseq_results)) {
    # Start summarizing the LRT results
    md_content <- paste0(md_content, "# Likelihood Ratio Test (LRT) Results\n\n---\n\n")

    # Describe the full and reduced formulas
    md_content <- paste0(
      md_content,
      "Based on your **full formula**: `", full_formula, "` and **reduced formula**: `", reduced_formula, "`, ",
      "this LRT analysis tests whether removing the interaction term (or terms) significantly affects gene expression. ",
      "The test uses only the **FDR adjusted p-value** (padj) to determine significance, as Log Fold Change (LFC) is irrelevant in the context of LRT.\n\n"
    )

    # Calculate the number of significant genes based on padj value
    significant_genes <- sum(deseq_results$padj < alpha, na.rm = TRUE)
    total_genes <- sum(!is.na(deseq_results$padj))

    # Extract the outliers and low counts from the DESeq results summary
    summary_output <- capture.output(summary(deseq_results))

    outliers <- gsub(".*: ", "", summary_output[6]) # Outliers
    low_counts <- gsub(".*: ", "", summary_output[7]) # Low counts (independent filtering)
    mean_count <- gsub("[^0-9]", "", summary_output[8])

    # Add a summary of significant genes
    lrt_summary <- paste0(
      "### Results Summary\n\n",
      "From this LRT analysis, **", significant_genes, " genes** (out of ", total_genes, " tested) are identified as significant with a padj value < ", alpha, ".\n\n"
    )

    # Add information about outliers and low counts
    lrt_summary <- paste0(
      lrt_summary,
      "**Outliers**<sup>1</sup>: ", outliers, " of genes were detected as outliers and excluded from analysis.\n\n",
      "**Low counts**<sup>2</sup>: ", low_counts, " of genes were removed due to low counts (mean <", mean_count, ") and independent filtering.\n\n",
      "Arguments of ?DESeq2::results():   \n<sup>1</sup> - see 'cooksCutoff',\n<sup>2</sup> - see 'independentFiltering'\n\n"
    )

    # Add this summary to the markdown content
    md_content <- paste0(md_content, lrt_summary)

    # Add explanation for the next steps
    next_steps <- paste0(
      "---\n\n",
      "### Next Steps\n\n",
      "If the number of significant genes is substantial, consider including the interaction term in your design formula ",
      "for a more detailed analysis of differential gene expression.\n\n",
      "For further insights and to explore detailed contrasts using the Wald test for the complex design formula, ",
      "please visit the **Complex Interaction Analysis** tab for more information.\n\n"
    )

    md_content <- paste0(md_content, next_steps)
  }

  if (!is.null(rpkm_filtered_count)) {
    md_content <- glue::glue("{md_content}### RPKM Filtering Summary\n\n")
    md_content <- glue::glue("{md_content}**Genes filtered out by RPKM threshold**: {rpkm_filtered_count}\n\n")
  }
  # Write the content to the output file
  writeLines(md_content, con = output_file)
}

# Enhanced logging function with file output and debug levels
log_message <- function(message, type = "INFO", log_file = "deseq_analysis.log") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  prefix <- switch(type,
                  "INFO" = "ℹ️",
                  "WARNING" = "⚠️",
                  "ERROR" = "❌",
                  "SUCCESS" = "✅",
                  "DEBUG" = "🔍",
                  "STEP" = "▶️",
                  "▶️")
  
  formatted_message <- glue::glue("[{timestamp}] {prefix} {message}")
  message(formatted_message)
  
  # Always write to log file
  write(formatted_message, file = log_file, append = TRUE)
}

# Setup debug mode handling
enable_debug_mode <- function(args) {
  # Create a global debug flag
  assign("DEBUG_MODE", args$debug, envir = .GlobalEnv)
  
  if (args$debug) {
    log_message("Debug mode enabled - verbose logging will be shown", "DEBUG")
    # Set R options for more verbose error reporting
    options(error = function() {
      traceback(3)
      log_message("Error encountered, see traceback above", "ERROR")
      if (!interactive()) quit(status = 1)
    })
  }
  
  # Clean up any existing log file at the start
  if (args$clean_logs) {
    file.remove("deseq_analysis.log")
    log_message("Started new log file", "INFO")
  }
}

# Debug-specific logging that only appears when debug mode is on
debug_log <- function(message, data = NULL) {
  if (exists("DEBUG_MODE") && DEBUG_MODE) {
    log_message(message, "DEBUG")
    
    # If data object is provided, print its summary
    if (!is.null(data)) {
      if (is.data.frame(data)) {
        log_message(glue::glue("Data dimensions: {nrow(data)} rows x {ncol(data)} columns"), "DEBUG")
        if (nrow(data) > 0) {
          class_info <- sapply(data, class)
          log_message("Column types: ", "DEBUG")
          print(class_info)
          
          # Print a few rows for inspection
          log_message("Preview of data:", "DEBUG")
          print(head(data, 5))
        }
      } else if (is.matrix(data)) {
        log_message(glue::glue("Matrix dimensions: {nrow(data)} rows x {ncol(data)} columns"), "DEBUG")
        log_message("Matrix preview:", "DEBUG")
        print(head(data, 5))
      } else if (is.list(data)) {
        log_message(glue::glue("List with {length(data)} elements"), "DEBUG")
        log_message(glue::glue("List names: {paste(names(data), collapse=', ')}"), "DEBUG")
      } else {
        log_message("Object summary:", "DEBUG")
        print(summary(data))
      }
    }
  }
}

# Main execution function to organize workflow
run_deseq_analysis <- function(args) {
  # Start time tracking for the entire analysis
  total_start_time <- proc.time()
  
  # Setup and validation
  message("=== DESeq2 Analysis Pipeline ===")
  message(glue::glue("Analysis started at {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"))
  
  # 1. Load and validate metadata
  metadata_df <- load_and_validate_metadata(args)
  
  # 2. Load and validate expression data
  expression_data <- load_and_validate_expression_data(args, metadata_df)
  
  # 3. Process count data and apply batch correction
  count_data_results <- process_count_data(args, expression_data$counts, metadata_df, 
                                         expression_data$design_formula)
  
  # 4. Run DESeq2 analysis
  deseq_results <- run_deseq2(count_data_results$countData, 
                            metadata_df, 
                            count_data_results$design_formula,
                            args$reduced,
                            args)
  
  # 5. Process and export results
  export_results(deseq_results, expression_data$expression_df, 
                metadata_df, args, count_data_results$batch_warning,
                expression_data$rpkm_filtered_count)
  
  # Report total elapsed time
  total_elapsed <- proc.time() - total_start_time
  message(glue::glue("=== Analysis completed in {round(total_elapsed['elapsed']/60, 1)} minutes ==="))
  message(glue::glue("Finished at {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"))
}

# Function to load and validate metadata
load_and_validate_metadata <- function(args) {
  message("Loading metadata...")
  
  # Check if metadata file is formatted correctly
  check_file_delimiter(args$meta)
  
  # Load metadata
  metadata_df <- read.table(
    args$meta,
    sep = get_file_type(args$meta),
    header = TRUE,
    stringsAsFactors = FALSE,
    row.names = 1
  )
  
  # Clean metadata column and row names
  colnames(metadata_df) <- clean_sample_names(colnames(metadata_df))
  rownames(metadata_df) <- clean_sample_names(rownames(metadata_df))
  
  message(glue::glue("Loaded metadata for {nrow(metadata_df)} samples with {ncol(metadata_df)} covariates"))
  
  # Validate metadata
  check_batch_column(metadata_df, args$batchcorrection)
  
  # Check design formulas
  design_formula <- check_design_formula(args$design, "Full Design-Formula")
  reduced_formula <- check_design_formula(args$reduced, "Reduced Design-Formula")
  
  # Verify covariates exist in metadata
  check_design_covariates(design_formula, metadata_df)
  
  # Make sure there are sufficient replicates for all factors
  check_replicates(metadata_df, all.vars(design_formula))
  
  # Add formulas to metadata for convenience
  attr(metadata_df, "design_formula") <- design_formula
  attr(metadata_df, "reduced_formula") <- reduced_formula
  
  return(metadata_df)
}

# Function to load and validate expression data
load_and_validate_expression_data <- function(args, metadata_df) {
  message("Loading expression data...")
  
  # Clean sample names for consistency
  clean_names <- clean_sample_names(args$name)
  
  # Load expression data
  expression_data_df <- load_expression_data(args$input, clean_names, READ_COL, RPKM_COL, INTERSECT_BY)
  message(glue::glue("Loaded expression data for {nrow(expression_data_df)} genes from {length(args$input)} files"))
  
  # Apply RPKM filtering if specified
  rpkm_filtered_count <- NULL
  if (!is.null(args$rpkm_cutoff)) {
    message(glue::glue("Applying RPKM cutoff of {args$rpkm_cutoff}..."))
    initial_gene_count <- nrow(expression_data_df)
    
    expression_data_df <- filter_rpkm(expression_data_df, args$rpkm_cutoff)
    
    # Ensure at least some genes remain
    if (nrow(expression_data_df) == 0) {
      report_error(
        "No genes remaining after RPKM filtering!",
        details = "All genes were removed due to the RPKM cutoff being too high.",
        recommendations = c("Reduce the RPKM threshold to retain more genes.")
      )
    }
    
    # Calculate how many genes were removed
    rpkm_filtered_count <- initial_gene_count - nrow(expression_data_df)
    
    message(glue::glue("{rpkm_filtered_count} genes removed, {nrow(expression_data_df)} genes retained"))
  }
  
  # Process count data
  read_counts_columns <- grep(
    paste(READ_COL, sep = ""),
    colnames(expression_data_df),
    value = TRUE,
    ignore.case = TRUE
  )
  
  read_counts_data_df <- expression_data_df %>%
    dplyr::mutate_at("GeneId", toupper) %>%
    dplyr::distinct(GeneId, .keep_all = TRUE) %>%
    remove_rownames() %>%
    column_to_rownames("GeneId")
  
  read_counts_data_df <- read_counts_data_df[read_counts_columns]
  
  # Clean count data column names
  colnames(read_counts_data_df) <- lapply(colnames(read_counts_data_df), function(s) {
    paste(head(unlist(strsplit(s, " ", fixed = TRUE)), -1), collapse = " ")
  })
  colnames(read_counts_data_df) <- clean_sample_names(colnames(read_counts_data_df))
  
  # Verify sample name consistency between metadata and counts
  validate_sample_consistency(metadata_df, read_counts_data_df)
  
  # Reorder count data columns to match metadata
  read_counts_data_df <- read_counts_data_df[, rownames(metadata_df)]
  
  # Return all results
  return(list(
    expression_df = expression_data_df,
    counts = read_counts_data_df,
    design_formula = attr(metadata_df, "design_formula"),
    reduced_formula = attr(metadata_df, "reduced_formula"),
    rpkm_filtered_count = rpkm_filtered_count
  ))
}

# Function to validate sample names consistency
validate_sample_consistency <- function(metadata_df, count_data_df) {
  log_message("Validating sample consistency between metadata and expression data...", "STEP")
  
  debug_log("Metadata sample names:", rownames(metadata_df))
  debug_log("Count data column names:", colnames(count_data_df))
  
  # Check for exact matches
  samples_in_metadata_not_in_data <- setdiff(rownames(metadata_df), colnames(count_data_df))
  samples_in_data_not_in_metadata <- setdiff(colnames(count_data_df), rownames(metadata_df))
  
  # More advanced debugging - check for case-insensitive matches or partial matches
  if (length(samples_in_metadata_not_in_data) > 0) {
    debug_log("Checking for case-insensitive matches for missing samples...")
    
    # Try case-insensitive matching
    lower_count_names <- tolower(colnames(count_data_df))
    for (sample in samples_in_metadata_not_in_data) {
      matches <- grep(tolower(sample), lower_count_names, value = TRUE)
      if (length(matches) > 0) {
        debug_log(glue::glue("Sample '{sample}' has potential case-insensitive matches: {paste(matches, collapse=', ')}"))
      }
    }
  }
  
  # Report results
  if (length(samples_in_metadata_not_in_data) > 0) {
    report_error(
      "Missing Samples in Expression Data",
      details = c(
        "### Issue:",
        "The following samples are present in the metadata but missing from the expression data:",
        "",
        paste(paste0("- ", samples_in_metadata_not_in_data), collapse = "\n")
      ),
      recommendations = c(
        "Ensure that sample names in the metadata file match exactly those in the expression data.",
        "Check for typos, extra spaces, or differences in letter case.",
        "Verify that the correct expression data file was loaded."
      )
    )
  }
  
  if (length(samples_in_data_not_in_metadata) > 0) {
    report_error(
      "Missing Samples in Metadata",
      details = c(
        "### Issue:",
        "The following samples are present in the expression data but missing from the metadata:",
        "",
        paste(paste0("- ", samples_in_data_not_in_metadata), collapse = "\n")
      ),
      recommendations = c(
        "Ensure that all samples in the expression data are listed in the metadata file.",
        "Check for typos, extra spaces, or differences in letter case.",
        "Confirm that the metadata file corresponds to the current expression dataset."
      )
    )
  }
  
  log_message("✓ Sample consistency validation passed", "SUCCESS")
}

# Enhanced DESeq2 analysis with detailed progress reporting
run_deseq2 <- function(countData, metadata_df, design_formula, reduced_formula, args) {
  log_message("Starting DESeq2 analysis...", "STEP")
  
  # Verify input data dimensions
  log_message(glue::glue("Count data dimensions: {nrow(countData)} genes x {ncol(countData)} samples"), "INFO")
  log_message(glue::glue("Metadata dimensions: {nrow(metadata_df)} samples x {ncol(metadata_df)} variables"), "INFO")
  
  # Create DESeq2 dataset
  log_message("Creating DESeqDataSet object...", "INFO")
  dse <- tryCatch({
    DESeqDataSetFromMatrix(
      countData = countData,
      colData = metadata_df,
      design = design_formula
    )
  }, error = function(e) {
    log_message(glue::glue("Error creating DESeqDataSet: {e$message}"), "ERROR")
    stop(e$message)
  })
  
  log_message(glue::glue("Created DESeqDataSet with {nrow(dse)} genes and {ncol(dse)} samples"), "SUCCESS")
  
  # Save the raw DESeqDataSet if debugging is enabled
  if (args$save_rdata) {
    saveRDS(dse, file = "debug_dse_raw.rds")
    log_message("Saved raw DESeqDataSet to debug_dse_raw.rds", "DEBUG")
  }
  
  # Run DESeq2 with Wald test for contrasts
  log_message("Running DESeq2 with Wald test (this may take a while)...", "STEP")
  start_time <- proc.time()
  
  dsq_wald <- tryCatch({
    DESeq(
      dse,
      test = "Wald",
      quiet = FALSE,
      parallel = TRUE
    )
  }, error = function(e) {
    log_message(glue::glue("Error in DESeq Wald test: {e$message}"), "ERROR")
    stop(e$message)
  })
  
  end_time <- proc.time() - start_time
  log_message(glue::glue("DESeq2 Wald test completed in {round(end_time['elapsed']/60, 2)} minutes"), "SUCCESS")
  
  # Save Wald results if debugging
  if (args$save_rdata) {
    saveRDS(dsq_wald, file = "debug_dsq_wald.rds")
    log_message("Saved Wald test results to debug_dsq_wald.rds", "DEBUG")
  }
  
  # Normalize counts
  log_message("Normalizing counts with rlog transformation...", "STEP")
  start_time <- proc.time()
  
  if (args$batchcorrection == "model" && "batch" %in% colnames(metadata_df)) {
    log_message("Applying rlog transformation and batch effect removal...", "INFO")
    
    rlog_transformed <- tryCatch({
      rlog(dsq_wald, blind = FALSE)
    }, error = function(e) {
      log_message(glue::glue("Error in rlog transformation: {e$message}"), "ERROR")
      stop(e$message)
    })
    
    rlog_counts <- assay(rlog_transformed)
    
    # Prepare design matrix without 'batch' for removeBatchEffect
    design_formula_no_batch <- as.formula(paste0("~", str_remove(as.character(dsq_wald@design)[2], " \\+ batch")))
    design_matrix <- model.matrix(design_formula_no_batch, data = metadata_df)
    
    debug_log("Design matrix for batch correction:", head(design_matrix))
    
    # Apply removeBatchEffect
    log_message("Removing batch effects from normalized counts...", "INFO")
    corrected_counts <- tryCatch({
      limma::removeBatchEffect(rlog_counts, batch = metadata_df$batch, design = design_matrix)
    }, error = function(e) {
      log_message(glue::glue("Error in batch effect removal: {e$message}"), "ERROR")
      stop(e$message)
    })
    
    normCounts <- corrected_counts
    
    # Compare before/after batch correction
    debug_log("Before batch correction (summary):")
    print(summary(as.vector(rlog_counts)))
    debug_log("After batch correction (summary):")
    print(summary(as.vector(normCounts)))
    
  } else {
    log_message("Applying standard rlog transformation without batch correction...", "INFO")
    rlog_transformed <- tryCatch({
      rlog(dsq_wald, blind = FALSE)
    }, error = function(e) {
      log_message(glue::glue("Error in rlog transformation: {e$message}"), "ERROR")
      stop(e$message)
    })
    
    normCounts <- assay(rlog_transformed)
  }
  
  end_time <- proc.time() - start_time
  log_message(glue::glue("Normalization completed in {round(end_time['elapsed'], 2)} seconds"), "SUCCESS")
  
  # Run DESeq2 with LRT
  log_message("Running DESeq2 with Likelihood Ratio Test...", "STEP")
  start_time <- proc.time()
  
  dsq_lrt <- tryCatch({
    DESeq(
      dse,
      test = "LRT",
      reduced = reduced_formula,
      quiet = FALSE,
      parallel = TRUE
    )
  }, error = function(e) {
    log_message(glue::glue("Error in DESeq LRT test: {e$message}"), "ERROR")
    stop(e$message)
  })
  
  end_time <- proc.time() - start_time
  log_message(glue::glue("DESeq2 LRT test completed in {round(end_time['elapsed']/60, 2)} minutes"), "SUCCESS")
  
  # Save LRT results if debugging
  if (args$save_rdata) {
    saveRDS(dsq_lrt, file = "debug_dsq_lrt.rds")
    log_message("Saved LRT test results to debug_dsq_lrt.rds", "DEBUG")
  }
  
  # Get LRT results
  log_message("Extracting LRT test results...", "INFO")
  dsq_lrt_res <- tryCatch({
    results(
      dsq_lrt,
      alpha = args$fdr,
      independentFiltering = TRUE
    )
  }, error = function(e) {
    log_message(glue::glue("Error extracting LRT results: {e$message}"), "ERROR")
    stop(e$message)
  })
  
  # Log summary of LRT results
  significant_count <- sum(dsq_lrt_res$padj < args$fdr, na.rm = TRUE)
  total_tested <- sum(!is.na(dsq_lrt_res$padj))
  log_message(glue::glue("LRT results: {significant_count} significant genes of {total_tested} tested (FDR < {args$fdr})"), "INFO")
  
  # Generate contrasts if not in LRT-only mode
  contrast_df <- data.frame()
  if (!args$lrt_only_mode) {
    log_message("Generating contrasts...", "STEP")
    contrast_df <- generate_contrasts(dsq_wald, args)
    log_message(glue::glue("Generated {nrow(contrast_df)} contrast entries"), "INFO")
  } else {
    log_message("Skipping contrast generation (LRT-only mode)", "INFO")
  }
  
  # Return all results
  log_message("DESeq2 analysis completed successfully", "SUCCESS")
  return(list(
    dsq_wald = dsq_wald,
    dsq_lrt = dsq_lrt,
    dsq_lrt_res = dsq_lrt_res,
    normCounts = normCounts,
    contrast_df = contrast_df
  ))
}

# Function to export all results
export_results <- function(deseq_results, expression_data_df, metadata_df, args, 
                          batch_warning, rpkm_filtered_count) {
  message("Exporting analysis results...")
  output_prefix <- args$output
  
  # Export LRT contrasts with significant gene counts
  contrasts_df <- generate_effect_contrasts(metadata_df, as.formula(args$design), deseq_results$dsq_lrt, args$fdr)
  
  # Write contrasts to TSV
  contrasts_filename <- paste(output_prefix, "_contrasts_table.tsv", sep = "")
  write.table(contrasts_df, file = contrasts_filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  message(glue::glue("✓ Exported contrasts to {contrasts_filename}"))
  
  # Annotate expression data with LRT results
  annotated_expression_df <- expression_data_df %>%
    bind_cols(as.data.frame(deseq_results$dsq_lrt_res) %>% select(pvalue, padj)) %>%
    mutate(
      `-LOG10(pval)` = -log10(as.numeric(pvalue)),
      `-LOG10(padj)` = -log10(as.numeric(padj))
    )
  
  # Export LRT report
  lrt_report_filename <- paste0(output_prefix, "_lrt_result.md")
  generate_lrt_md(
    deseq_results$dsq_lrt_res, 
    args$design, 
    args$reduced, 
    lrt_report_filename, 
    alpha = args$fdr, 
    batch_warning = batch_warning, 
    rpkm_filtered_count = rpkm_filtered_count
  )
  message(glue::glue("✓ Exported LRT markdown report to {lrt_report_filename}"))
  
  # Export full expression table
  results_filename <- paste0(output_prefix, "_gene_exp_table.tsv")
  write.table(
    annotated_expression_df,
    file = results_filename,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  message(glue::glue("✓ Exported results table to {results_filename}"))
  
  # Clean up graphics
  graphics.off()
  
  # Process sample columns for the annotations
  sample_order <- colnames(deseq_results$dsq_lrt)
  
  colnames(annotated_expression_df) <- gsub(" Rpkm", "", colnames(annotated_expression_df))
  annotated_expression_df <- annotated_expression_df %>%
    select(-one_of(sample_order)) # remove all sample columns
  
  colnames(annotated_expression_df) <- gsub(" TotalReads", "", colnames(annotated_expression_df))
  
  annotated_expression_df <- annotated_expression_df %>%
    select(-one_of(sample_order)) %>% # remove all sample columns
    distinct(GeneId, .keep_all = TRUE) %>%
    remove_rownames() %>%
    column_to_rownames("GeneId")
  
  # Export charts and GCT files
  export_charts(
    deseq_results$dsq_lrt_res,
    annotated_expression_df, 
    metadata_df, 
    deseq_results$normCounts, 
    output_prefix, 
    args
  )
  message("✓ Exported charts and GCT files")
}

# Export GCT with subset of significant genes only
subset_significant_genes <- function(lrt_results, alpha = 0.1) {
  log_message("=== Debugging before subsetting ===")
  # Get genes with padj < alpha
  sig_genes <- rownames(lrt_results)[which(lrt_results$padj < alpha)]
  log_message(paste("Checking if all row names in row_metadata_filtered are present in normCounts:"))
  log_message(paste("All row names present:", !is.null(sig_genes)))
  log_message(paste("Number of rows after subsetting:", length(sig_genes)))
  return(sig_genes)
}

# Create a function to standardize plot settings
configure_plot_theme <- function() {
  # Set default theme for all plots
  ggplot2::theme_set(
    ggplot2::theme_classic() + 
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "right"
    )
  )
}

# Call this at the start of your script
configure_plot_theme()

# Main script execution with enhanced error handling
main <- function() {
  # Set up error handling
  tryCatch({
    # Parse arguments
    args <- get_args()
    
    # Setup debugging and logging
    if ("debug" %in% names(args) && args$debug) {
      enable_debug_mode(args)
    }
    
    # Configure BiocParallel
    register(MulticoreParam(args$threads))
    log_message(glue::glue("Using {args$threads} CPU threads for parallel processing"), "INFO")
    
    # Run the analysis
    run_deseq_analysis(args)
    
    # Report successful completion
    log_message("DESeq2 analysis completed successfully! ✓", "SUCCESS")
    
  }, error = function(e) {
    # Capture full error information
    error_msg <- paste("Analysis failed:", e$message)
    error_trace <- paste(capture.output(traceback()), collapse="\n")
    
    log_message(error_msg, "ERROR")
    log_message("Error traceback:", "ERROR")
    message(error_trace)
    
    # Check if error_report.txt exists
    if (file.exists("error_report.txt")) {
      log_message("See error_report.txt for detailed error information", "ERROR")
    } else {
      # Create basic error report if none exists
      writeLines(
        paste0(
          "# ❌ Error Report\n\n",
          "**Timestamp:** ", format(Sys.time(), "%A, %B %d, %Y %I:%M:%S %p"), "\n\n",
          "**Error Message:** ", e$message, "\n\n",
          "## Error Traceback\n\n```\n", error_trace, "\n```\n"
        ),
        con = "error_report.txt"
      )
      log_message("Basic error report written to error_report.txt", "ERROR")
    }
    
    # Exit with error code
    quit(status = 1)
  })
}

# Execute main function
main()

# Generate all possible contrasts from DESeq2 results object
generate_contrasts <- function(dds, args) {
  log_message("Generating contrasts for differential expression analysis...", "STEP")
  
  # Get results names from DESeq2 object
  result_names <- resultsNames(dds)
  debug_log(glue::glue("Available result names from DESeq2: {paste(result_names, collapse=', ')}"))
  
  # Skip intercept
  result_names <- result_names[result_names != "Intercept"]
  
  # Initialize contrast dataframe  
  contrast_df <- data.frame(
    factor = character(),
    contrast_name = character(),
    numerator = character(),
    denominator = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # If no valid results, return empty dataframe
  if (length(result_names) == 0) {
    log_message("No valid contrasts found in DESeq2 results.", "WARNING")
    return(contrast_df)
  }
  
  # Process main effects
  main_effect_contrasts <- generate_main_effect_contrasts(dds, result_names)
  if (nrow(main_effect_contrasts) > 0) {
    contrast_df <- rbind(contrast_df, main_effect_contrasts)
    log_message(glue::glue("Generated {nrow(main_effect_contrasts)} main effect contrasts"), "SUCCESS")
  } else {
    log_message("No main effect contrasts generated", "INFO")
  }
  
  # Process interaction effects
  interaction_effect_contrasts <- generate_interaction_effect_contrasts(dds, result_names)
  if (nrow(interaction_effect_contrasts) > 0) {
    contrast_df <- rbind(contrast_df, interaction_effect_contrasts)
    log_message(glue::glue("Generated {nrow(interaction_effect_contrasts)} interaction effect contrasts"), "SUCCESS")
  } else {
    log_message("No interaction effect contrasts generated", "INFO")
  }
  
  # Final validation
  if (nrow(contrast_df) == 0) {
    log_message("Warning: No contrasts could be generated. Check your design formula and metadata.", "WARNING")
  } else {
    log_message(glue::glue("Successfully generated {nrow(contrast_df)} total contrasts"), "SUCCESS")
  }
  
  return(contrast_df)
}

# Generate contrasts for main effects
generate_main_effect_contrasts <- function(dds, result_names) {
  log_message("Processing main effect contrasts...", "INFO")
  
  # Initialize dataframe
  main_effect_contrasts <- data.frame(
    factor = character(),
    contrast_name = character(),
    numerator = character(),
    denominator = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # Get all factor columns from colData
  col_data <- as.data.frame(colData(dds))
  factor_cols <- sapply(col_data, is.factor)
  factor_names <- names(factor_cols)[factor_cols]
  
  # Skip 'batch' column if it exists and is a factor
  if ("batch" %in% factor_names) {
    factor_names <- factor_names[factor_names != "batch"]
    debug_log("Excluding 'batch' from contrast generation")
  }
  
  debug_log(glue::glue("Processing factor columns: {paste(factor_names, collapse=', ')}"))
  
  # Process each factor
  for (factor_name in factor_names) {
    # Get levels for this factor
    factor_levels <- levels(col_data[[factor_name]])
    
    if (length(factor_levels) < 2) {
      debug_log(glue::glue("Skipping factor '{factor_name}' - fewer than 2 levels"))
      next
    }
    
    debug_log(glue::glue("Processing factor '{factor_name}' with levels: {paste(factor_levels, collapse=', ')}"))
    
    # Process factor according to number of levels
    if (length(factor_levels) == 2) {
      # Binary factor - simple case
      log_message(glue::glue("Processing binary factor: {factor_name}"), "INFO")
      
      # Check if the contrast exists in result_names
      contrast_name <- paste0(factor_name, factor_levels[2])
      
      if (contrast_name %in% result_names) {
        main_effect_contrasts <- rbind(
          main_effect_contrasts,
          data.frame(
            factor = factor_name,
            contrast_name = contrast_name,
            numerator = factor_levels[2],
            denominator = factor_levels[1],
            type = "main_effect",
            stringsAsFactors = FALSE
          )
        )
        debug_log(glue::glue("Added binary contrast: {contrast_name}"))
      } else {
        # Try alternate naming pattern
        alt_contrast_name <- find_contrast_name(result_names, factor_name, factor_levels[2])
        
        if (!is.null(alt_contrast_name)) {
          main_effect_contrasts <- rbind(
            main_effect_contrasts,
            data.frame(
              factor = factor_name,
              contrast_name = alt_contrast_name,
              numerator = factor_levels[2],
              denominator = factor_levels[1],
              type = "main_effect",
              stringsAsFactors = FALSE
            )
          )
          debug_log(glue::glue("Added binary contrast with alternate name: {alt_contrast_name}"))
        } else {
          log_message(glue::glue("Warning: Could not find contrast for binary factor '{factor_name}'"), "WARNING")
        }
      }
    } else {
      # Multi-level factor
      log_message(glue::glue("Processing multi-level factor: {factor_name} with {length(factor_levels)} levels"), "INFO")
      
      # Reference level is the first one
      ref_level <- factor_levels[1]
      
      # For each non-reference level
      for (i in 2:length(factor_levels)) {
        level <- factor_levels[i]
        contrast_name <- paste0(factor_name, level)
        
        if (contrast_name %in% result_names) {
          main_effect_contrasts <- rbind(
            main_effect_contrasts,
            data.frame(
              factor = factor_name,
              contrast_name = contrast_name,
              numerator = level,
              denominator = ref_level,
              type = "main_effect",
              stringsAsFactors = FALSE
            )
          )
          debug_log(glue::glue("Added multi-level contrast: {contrast_name}"))
        } else {
          # Try alternate naming pattern
          alt_contrast_name <- find_contrast_name(result_names, factor_name, level)
          
          if (!is.null(alt_contrast_name)) {
            main_effect_contrasts <- rbind(
              main_effect_contrasts,
              data.frame(
                factor = factor_name,
                contrast_name = alt_contrast_name,
                numerator = level,
                denominator = ref_level,
                type = "main_effect",
                stringsAsFactors = FALSE
              )
            )
            debug_log(glue::glue("Added multi-level contrast with alternate name: {alt_contrast_name}"))
          } else {
            log_message(glue::glue("Warning: Could not find contrast for level '{level}' of factor '{factor_name}'"), "WARNING")
          }
        }
      }
    }
  }
  
  log_message(glue::glue("Completed main effect contrast generation with {nrow(main_effect_contrasts)} contrasts"), "INFO")
  return(main_effect_contrasts)
}

# Helper function to find contrast name using various patterns
find_contrast_name <- function(result_names, factor_name, level) {
  # Common patterns for contrast names in DESeq2
  patterns <- c(
    paste0(factor_name, "_", level),
    paste0(factor_name, level),
    paste0(factor_name, ".", level),
    paste0(tolower(factor_name), "_", tolower(level)),
    paste0(tolower(factor_name), tolower(level)),
    paste0(tolower(factor_name), ".", tolower(level))
  )
  
  debug_log(glue::glue("Searching for contrast patterns: {paste(patterns, collapse=', ')}"))
  
  # Try exact matches first
  for (pattern in patterns) {
    if (pattern %in% result_names) {
      return(pattern)
    }
  }
  
  # Try partial/fuzzy matches
  for (pattern in patterns) {
    matches <- grep(pattern, result_names, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      debug_log(glue::glue("Found fuzzy match: {matches[1]} for pattern {pattern}"))
      return(matches[1])
    }
  }
  
  # No match found
  return(NULL)
}

# Generate contrasts for interaction effects
generate_interaction_effect_contrasts <- function(dds, result_names) {
  log_message("Processing interaction effect contrasts...", "INFO")
  
  # Initialize dataframe
  interaction_effect_contrasts <- data.frame(
    factor = character(),
    contrast_name = character(),
    numerator = character(),
    denominator = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # Identify interaction terms in result_names
  interaction_terms <- result_names[grepl(":", result_names, fixed = TRUE)]
  
  if (length(interaction_terms) == 0) {
    debug_log("No interaction terms found in results names")
    return(interaction_effect_contrasts)
  }
  
  debug_log(glue::glue("Found interaction terms: {paste(interaction_terms, collapse=', ')}"))
  
  # Process each interaction term
  for (term in interaction_terms) {
    # Split by colon to get factors and levels
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    
    if (length(parts) < 2) {
      debug_log(glue::glue("Skipping malformed interaction term: {term}"))
      next
    }
    
    # Extract factor names and levels (this is complex and may need customization)
    # This is a simplified approach and may need to be adapted
    factors <- character(0)
    levels <- character(0)
    
    for (part in parts) {
      # Try to separate factor name from level
      # This assumes a pattern like "factor1level1"
      # Extract alphabetic prefix (factor) and the rest (level)
      match <- regexpr("[A-Za-z]+", part)
      if (match > 0) {
        factor_name <- substr(part, match, attr(match, "match.length"))
        level <- substr(part, attr(match, "match.length") + 1, nchar(part))
        
        factors <- c(factors, factor_name)
        levels <- c(levels, level)
      }
    }
    
    if (length(factors) != length(parts) || length(levels) != length(parts)) {
      debug_log(glue::glue("Could not parse interaction term components for: {term}"))
      next
    }
    
    # Create a descriptive factor name for the interaction
    factor_name <- paste(factors, collapse = ":")
    
    # Create descriptive numerator/denominator
    numerator <- paste(levels, collapse = ":")
    denominator <- "reference levels"
    
    # Add to results
    interaction_effect_contrasts <- rbind(
      interaction_effect_contrasts,
      data.frame(
        factor = factor_name,
        contrast_name = term,
        numerator = numerator,
        denominator = denominator,
        type = "interaction_effect",
        stringsAsFactors = FALSE
      )
    )
    
    debug_log(glue::glue("Added interaction contrast: {term}"))
  }
  
  log_message(glue::glue("Completed interaction effect contrast generation with {nrow(interaction_effect_contrasts)} contrasts"), "INFO")
  return(interaction_effect_contrasts)
}

# Function to validate the design matrix has full rank
validate_full_rank_design <- function(design_formula, metadata) {
  log_message("Validating design matrix rank...", "STEP")
  
  tryCatch({
    # Create design matrix
    design_matrix <- model.matrix(design_formula, data = metadata)
    
    # Check matrix rank
    matrix_rank <- qr(design_matrix)$rank
    matrix_cols <- ncol(design_matrix)
    
    debug_log(glue::glue("Design matrix: {matrix_cols} columns, rank {matrix_rank}"))
    
    if (matrix_rank < matrix_cols) {
      # Matrix is not full rank
      report_error(
        "Non-Full-Rank Design Matrix",
        details = c(
          glue::glue("Your design matrix has {matrix_cols} columns but rank {matrix_rank}."),
          "This indicates linear dependencies in your design, which will cause model fitting problems."
        ),
        recommendations = c(
          "Check for factors with only one sample per group",
          "Ensure batch variables don't perfectly correlate with experimental factors",
          "Consider simplifying your design formula",
          "Remove confounding factors that perfectly correlate with other factors"
        )
      )
    }
    
    log_message("Design matrix has full rank ✓", "SUCCESS")
    return(TRUE)
    
  }, error = function(e) {
    report_error(
      "Design Matrix Error",
      details = c(
        "Failed to create or validate the design matrix.",
        glue::glue("Error message: {e$message}")
      ),
      recommendations = c(
        "Check that all variables in your design formula exist in the metadata",
        "Ensure your metadata contains the correct data types (factors vs. numeric)",
        "Verify there are no missing values in the metadata"
      )
    )
  })
}

# Function to validate reduced model is properly nested in full model
validate_model_nesting <- function(full_formula, reduced_formula) {
  log_message("Validating model nesting relationship...", "STEP")
  
  # Extract terms from both formulas
  full_terms <- attr(terms(full_formula), "term.labels")
  reduced_terms <- attr(terms(reduced_formula), "term.labels")
  
  # Check if all reduced terms are in full model
  missing_terms <- setdiff(reduced_terms, full_terms)
  
  if (length(missing_terms) > 0) {
    report_error(
      "Invalid Reduced Model",
      details = c(
        "The reduced model contains terms not present in the full model:",
        paste(missing_terms, collapse=", ")
      ),
      recommendations = c(
        "The reduced model must be nested within the full model",
        "Remove the extra terms from the reduced model",
        "Ensure both models use the same variable names"
      )
    )
  }
  
  # Check if full model has at least one more term than reduced
  if (length(full_terms) <= length(reduced_terms)) {
    report_error(
      "Invalid Model Comparison",
      details = c(
        "The full model does not have more terms than the reduced model.",
        glue::glue("Full model terms: {paste(full_terms, collapse=', ')}"),
        glue::glue("Reduced model terms: {paste(reduced_terms, collapse=', ')}")
      ),
      recommendations = c(
        "The full model must contain additional terms compared to the reduced model",
        "Add terms of interest to the full model, or",
        "Remove terms from the reduced model"
      )
    )
  }
  
  # Check if there are interaction terms in the full model but not in reduced
  full_interactions <- full_terms[grepl(":", full_terms, fixed=TRUE)]
  reduced_interactions <- reduced_terms[grepl(":", reduced_terms, fixed=TRUE)]
  
  if (length(full_interactions) > 0 && length(setdiff(full_interactions, reduced_interactions)) == 0) {
    log_message("Warning: Full model has interactions but they are all present in reduced model. LRT may not test interactions.", "WARNING")
  }
  
  log_message("Model nesting validation passed ✓", "SUCCESS")
  return(TRUE)
}

# Function to check for sufficient replicates for all factors
check_replicates <- function(metadata_df, factor_names) {
  log_message("Checking for sufficient biological replicates...", "STEP")
  
  # Minimum number of replicates required
  min_replicates <- 2
  
  # Track problematic factors
  insufficient_replicates <- list()
  
  # Process only factor columns (skip batch)
  for (factor_name in factor_names) {
    if (factor_name == "batch" || !factor_name %in% colnames(metadata_df)) {
      next
    }
    
    # Convert to factor if not already
    if (!is.factor(metadata_df[[factor_name]])) {
      debug_log(glue::glue("Converting {factor_name} to factor"))
      metadata_df[[factor_name]] <- as.factor(metadata_df[[factor_name]])
    }
    
    # Count samples per level
    level_counts <- table(metadata_df[[factor_name]])
    debug_log(glue::glue("Factor {factor_name} level counts:"))
    print(level_counts)
    
    # Check if any level has fewer than minimum replicates
    problem_levels <- names(level_counts[level_counts < min_replicates])
    
    if (length(problem_levels) > 0) {
      insufficient_replicates[[factor_name]] <- problem_levels
      log_message(
        glue::glue("Warning: Factor '{factor_name}' has levels with fewer than {min_replicates} replicates: {paste(problem_levels, collapse=', ')}"),
        "WARNING"
      )
    }
  }
  
  # If problems found, report error
  if (length(insufficient_replicates) > 0) {
    # Construct details message
    details <- c("The following factors have levels with insufficient replicates:")
    for (factor_name in names(insufficient_replicates)) {
      levels <- insufficient_replicates[[factor_name]]
      details <- c(details, glue::glue("  - {factor_name}: {paste(levels, collapse=', ')}"))
    }
    
    report_error(
      "Insufficient Biological Replicates",
      details = details,
      recommendations = c(
        "DESeq2 requires at least 2 biological replicates per condition for reliable results.",
        "Consider adding more samples or combining factor levels with few samples.",
        "If this is unavoidable, consider using an alternative approach like edgeR."
      )
    )
  }
  
  log_message("All factors have sufficient replicates ✓", "SUCCESS")
  return(TRUE)
}

# Function to check file delimiter
check_file_delimiter <- function(file_path) {
  log_message(glue::glue("Checking format of file: {file_path}"), "INFO")
  
  # Try to read the first few lines
  tryCatch({
    lines <- readLines(file_path, n = 10)
    
    # Determine expected delimiter based on file extension
    expected_delimiter <- get_file_type(file_path)
    expected_char <- ifelse(expected_delimiter == ",", "comma", "tab")
    
    # Count delimiters in header line
    header <- lines[1]
    comma_count <- str_count(header, ",")
    tab_count <- str_count(header, "\t")
    
    # Determine actual delimiter
    if (comma_count > tab_count) {
      actual_char <- "comma"
    } else if (tab_count > comma_count) {
      actual_char <- "tab"
    } else {
      # If no delimiters found, report error
      report_error(
        "Invalid File Format",
        details = c(
          glue::glue("Could not detect a valid delimiter in file: {file_path}"),
          "The header line doesn't contain commas or tabs."
        ),
        recommendations = c(
          "Check that your file is properly formatted as CSV or TSV",
          "Ensure the file has a header line with column names",
          "Verify that the file is not empty or corrupted"
        )
      )
    }
    
    # If delimiter doesn't match file extension, warn user
    if (actual_char != expected_char) {
      log_message(
        glue::glue("Warning: File extension suggests {expected_char}-separated values, but content appears to be {actual_char}-separated"),
        "WARNING"
      )
    }
    
    log_message(glue::glue("File format check passed: {actual_char}-separated values"), "SUCCESS")
    
  }, error = function(e) {
    report_error(
      "File Reading Error",
      details = c(
        glue::glue("Failed to read file: {file_path}"),
        glue::glue("Error message: {e$message}")
      ),
      recommendations = c(
        "Check that the file exists and is accessible",
        "Verify the file is not corrupted or empty",
        "Ensure you have the necessary permissions to read the file"
      )
    )
  })
}

# Function to filter by RPKM
filter_rpkm <- function(expression_data_df, rpkm_cutoff) {
  log_message(glue::glue("Filtering genes by RPKM cutoff of {rpkm_cutoff}..."), "STEP")
  
  # Find RPKM columns
  rpkm_columns <- grep(" Rpkm$", colnames(expression_data_df), value = TRUE)
  
  if (length(rpkm_columns) == 0) {
    report_error(
      "RPKM Filtering Error",
      details = "Could not find any RPKM columns in the expression data",
      recommendations = c(
        "Check that the input files contain RPKM values",
        "Verify the column naming convention (' Rpkm' suffix)",
        "If your data doesn't include RPKM values, disable RPKM filtering"
      )
    )
  }
  
  debug_log(glue::glue("Found {length(rpkm_columns)} RPKM columns: {paste(head(rpkm_columns, 3), collapse=', ')}..."))
  
  # Count genes before filtering
  initial_count <- nrow(expression_data_df)
  log_message(glue::glue("Starting with {initial_count} genes"), "INFO")
  
  # Filter genes where at least one sample has RPKM >= cutoff
  filtered_df <- expression_data_df %>%
    filter(across(all_of(rpkm_columns), ~ any(. >= rpkm_cutoff)))
  
  # Count genes after filtering
  final_count <- nrow(filtered_df)
  filtered_count <- initial_count - final_count
  
  log_message(glue::glue("Removed {filtered_count} genes with RPKM < {rpkm_cutoff} in all samples"), "INFO")
  log_message(glue::glue("Retained {final_count} genes ({round(final_count/initial_count*100, 1)}% of original)"), "SUCCESS")
  
  if (final_count == 0) {
    report_error(
      "No Genes Remaining After RPKM Filtering",
      details = glue::glue("All {initial_count} genes were filtered out with RPKM cutoff of {rpkm_cutoff}"),
      recommendations = c(
        "Lower the RPKM cutoff value",
        "Check your expression data for potential issues",
        "Consider using raw counts without RPKM filtering"
      )
    )
  }
  
  return(filtered_df)
}

# Function to process count data and apply batch correction
process_count_data <- function(args, countData, metadata_df, design_formula) {
  log_message("Processing count data...", "STEP")
  
  # Quick sanity check on count data
  debug_log("Count data statistics:", summary(as.matrix(countData)))
  
  # Check for negative values in count data
  if (any(countData < 0, na.rm = TRUE)) {
    report_error(
      "Negative Values in Count Data",
      details = "Found negative values in count data. RNA-seq counts must be non-negative integers.",
      recommendations = c(
        "Check your input files for data processing errors",
        "Ensure your data represents raw counts, not normalized values",
        "If using pre-processed data, verify the processing steps"
      )
    )
  }
  
  # Check for non-integer values if not in test mode
  if (!args$test_mode && !all(countData == floor(countData), na.rm = TRUE)) {
    log_message("Warning: Count data contains non-integer values. DESeq2 expects integer counts for RNA-seq data.", "WARNING")
  }
  
  # Apply batch correction if specified
  batch_correction_results <- handle_batch_correction(args, countData, metadata_df, design_formula)
  
  # Validate the final design matrix has full rank
  validate_full_rank_design(batch_correction_results$design_formula, metadata_df)
  
  # Validate model nesting if using LRT
  if (!args$lrt_only_mode) {
    reduced_formula <- as.formula(args$reduced)
    validate_model_nesting(batch_correction_results$design_formula, reduced_formula)
  }
  
  log_message("Count data processing completed successfully", "SUCCESS")
  return(batch_correction_results)
}

# Improved diagnosis of interaction terms
diagnose_interaction_terms <- function(design_formula, results_names) {
  log_message("Diagnosing interaction terms in design formula...", "STEP")
  
  # Extract the formula string
  formula_str <- as.character(design_formula)[2]
  debug_log(glue::glue("Formula: {formula_str}"), "DEBUG")
  
  # Check for explicit interaction terms
  has_colon <- grepl(":", formula_str, fixed = TRUE)
  has_asterisk <- grepl("*", formula_str, fixed = TRUE, useBytes = TRUE)
  
  # Extract all terms from the design formula
  model_terms <- attr(terms(design_formula), "term.labels")
  interaction_terms <- grep(":", model_terms, value = TRUE, fixed = TRUE)
  
  if (length(interaction_terms) == 0) {
    if (has_colon || has_asterisk) {
      log_message("⚠️ Formula contains ':' or '*' operators but no interaction terms were identified in the model", "WARNING")
    } else {
      log_message("No interaction terms found in design formula", "INFO")
    }
    return(list(has_interactions = FALSE, terms = character(0)))
  }
  
  log_message(glue::glue("Found {length(interaction_terms)} interaction terms: {paste(interaction_terms, collapse=', ')}"), "INFO")
  
  # Check if interaction terms appear in results
  interaction_results <- grep(":", results_names, value = TRUE, fixed = TRUE)
  
  if (length(interaction_results) == 0) {
    log_message("⚠️ Interaction terms in formula but none found in DESeq2 results - check factor levels", "WARNING")
    return(list(has_interactions = TRUE, terms = interaction_terms, in_results = FALSE))
  }
  
  # Match formula interactions to results
  matched_interactions <- intersect(interaction_terms, results_names)
  unmatched_interactions <- setdiff(interaction_terms, results_names)
  
  if (length(matched_interactions) > 0) {
    log_message(glue::glue("Found {length(matched_interactions)} interaction terms in DESeq2 results"), "SUCCESS")
  }
  
  if (length(unmatched_interactions) > 0) {
    log_message(glue::glue("⚠️ {length(unmatched_interactions)} interaction terms from formula not found in results: {paste(unmatched_interactions, collapse=', ')}"), "WARNING")
  }
  
  return(list(
    has_interactions = TRUE, 
    terms = interaction_terms, 
    in_results = (length(interaction_results) > 0),
    matched = matched_interactions,
    unmatched = unmatched_interactions
  ))
}

# Enhanced sample name cleaning function
clean_sample_names <- function(names, preserve_case = FALSE) {
  log_message("Cleaning sample names...", "INFO")
  
  if (length(names) == 0) {
    return(names)
  }
  
  debug_log(glue::glue("Original names: {paste(head(names, 3), collapse=', ')}..."))
  
  # Function to clean individual name
  clean_name <- function(name) {
    # Replace problematic characters with underscore
    clean <- gsub("[^a-zA-Z0-9_]", "_", name)
    
    # Ensure name starts with a letter or underscore (R requirement)
    if (!grepl("^[a-zA-Z_]", clean)) {
      clean <- paste0("X", clean)
    }
    
    # Convert to lowercase if not preserving case
    if (!preserve_case) {
      clean <- tolower(clean)
    }
    
    return(clean)
  }
  
  # Apply cleaning to all names
  cleaned_names <- sapply(names, clean_name)
  
  # Check for duplicates
  if (any(duplicated(cleaned_names))) {
    log_message("Warning: Duplicate names found after cleaning. Adding unique suffixes.", "WARNING")
    
    # Find duplicates
    dups <- cleaned_names[duplicated(cleaned_names)]
    
    # Add suffixes to duplicates
    for (dup in unique(dups)) {
      idx <- which(cleaned_names == dup)
      if (length(idx) > 1) {
        for (i in 1:length(idx)) {
          cleaned_names[idx[i]] <- paste0(cleaned_names[idx[i]], "_", i)
        }
      }
    }
  }
  
  debug_log(glue::glue("Cleaned names: {paste(head(cleaned_names, 3), collapse=', ')}..."))
  log_message(glue::glue("Cleaned {length(names)} sample names"), "SUCCESS")
  
  return(cleaned_names)
}

# Function to check batch column if batch correction is requested
check_batch_column <- function(metadata_df, batch_correction) {
  log_message("Checking batch correction settings...", "INFO")
  
  if (batch_correction == "none") {
    log_message("No batch correction requested", "INFO")
    return(TRUE)
  }
  
  log_message(glue::glue("Batch correction method requested: {batch_correction}"), "INFO")
  
  # Check if 'batch' column exists
  if (!("batch" %in% colnames(metadata_df))) {
    report_error(
      "Missing Batch Column",
      details = c(
        "You requested batch correction, but there is no 'batch' column in your metadata.",
        glue::glue("Requested batch correction method: {batch_correction}"),
        glue::glue("Available metadata columns: {paste(colnames(metadata_df), collapse=', ')}")
      ),
      recommendations = c(
        "Add a 'batch' column to your metadata file",
        "Turn off batch correction by setting --batchcorrection=none",
        "Ensure the batch column is named exactly 'batch' (case-sensitive)"
      )
    )
  }
  
  # Check if batch column has at least 2 levels
  batch_values <- unique(metadata_df$batch)
  if (length(batch_values) < 2) {
    report_error(
      "Insufficient Batch Levels",
      details = c(
        glue::glue("The 'batch' column contains only {length(batch_values)} unique value(s): {paste(batch_values, collapse=', ')}"),
        "Batch correction requires at least 2 different batch values."
      ),
      recommendations = c(
        "Check your metadata file to ensure batches are correctly labeled",
        "Turn off batch correction if you don't have multiple batches"
      )
    )
  }
  
  # If using ComBat-seq, check that no batch has only one sample
  if (batch_correction == "combatseq") {
    batch_counts <- table(metadata_df$batch)
    single_sample_batches <- names(batch_counts[batch_counts == 1])
    
    if (length(single_sample_batches) > 0) {
      report_error(
        "Single-Sample Batches",
        details = c(
          "ComBat-seq requires multiple samples per batch for reliable correction.",
          glue::glue("The following batch(es) have only one sample: {paste(single_sample_batches, collapse=', ')}")
        ),
        recommendations = c(
          "Ensure each batch has at least 2 samples",
          "Consider using 'model' batch correction instead",
          "Or merge small batches if biologically appropriate"
        )
      )
    }
  }
  
  log_message("Batch column validation passed ✓", "SUCCESS")
  return(TRUE)
}

# Function to export visualization charts and GCT files
export_charts <- function(lrt_results, annotated_expression_df, metadata_df, normCounts, output_prefix, args) {
  log_message("Generating visualizations and GCT files...", "STEP")
  
  # Create output directories if they don't exist
  plots_dir <- file.path(dirname(output_prefix), "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
    log_message(glue::glue("Created plots directory: {plots_dir}"), "INFO")
  }
  
  # Export normalized counts as GCT file
  export_gct_data(normCounts, annotated_expression_df, metadata_df, output_prefix, args$fdr)
  
  # Generate PCA plot
  generate_pca_plot(normCounts, metadata_df, file.path(plots_dir, paste0(basename(output_prefix), "_pca.pdf")))
  
  # Generate sample distance heatmap
  generate_sample_distance_heatmap(normCounts, metadata_df, file.path(plots_dir, paste0(basename(output_prefix), "_sample_distance.pdf")))
  
  # Generate MDS plot using Glimma
  generate_mds_plot(normCounts, metadata_df, file.path(plots_dir, paste0(basename(output_prefix), "_mds")))
  
  # Generate MA plot for LRT results
  generate_ma_plot(lrt_results, file.path(plots_dir, paste0(basename(output_prefix), "_ma_plot.pdf")), args$fdr)
  
  log_message("Visualizations and GCT files generated successfully", "SUCCESS")
}

# Export data as GCT file for downstream analysis
export_gct_data <- function(normCounts, annotated_expression_df, metadata_df, output_prefix, fdr_threshold = 0.1) {
  log_message("Exporting normalized counts as GCT file...", "STEP")
  
  # Create GCT file path
  gct_file <- paste0(output_prefix, "_normalized_counts.gct")
  
  # Filter for significant genes if we have p-values
  if ("padj" %in% colnames(annotated_expression_df)) {
    log_message(glue::glue("Filtering for genes with FDR < {fdr_threshold}"), "INFO")
    sig_genes <- rownames(annotated_expression_df)[annotated_expression_df$padj < fdr_threshold]
    
    if (length(sig_genes) > 0) {
      # Keep only significant genes in the normalized counts
      filtered_counts <- normCounts[rownames(normCounts) %in% sig_genes, ]
      log_message(glue::glue("Exporting {nrow(filtered_counts)} significant genes to GCT"), "INFO")
      
      # Write filtered data to GCT format
      tryCatch({
        cmapR::write_gct(
          ds = cmapR::mat2ds(filtered_counts),
          ofile = gct_file
        )
        log_message(glue::glue("Filtered GCT file written to: {gct_file}"), "SUCCESS")
      }, error = function(e) {
        log_message(glue::glue("Error writing filtered GCT file: {e$message}"), "WARNING")
        # Try writing without filtering
        cmapR::write_gct(
          ds = cmapR::mat2ds(normCounts),
          ofile = gct_file
        )
        log_message("Wrote unfiltered GCT file as fallback", "INFO")
      })
    } else {
      log_message("No genes passed FDR threshold, writing all normalized counts", "WARNING")
      cmapR::write_gct(
        ds = cmapR::mat2ds(normCounts),
        ofile = gct_file
      )
    }
  } else {
    # No filtering, write all normalized counts
    log_message("No FDR information available, writing all normalized counts", "INFO")
    cmapR::write_gct(
      ds = cmapR::mat2ds(normCounts),
      ofile = gct_file
    )
  }
  
  # Also write a non-filtered version
  full_gct_file <- paste0(output_prefix, "_all_normalized_counts.gct")
  cmapR::write_gct(
    ds = cmapR::mat2ds(normCounts),
    ofile = full_gct_file
  )
  
  log_message(glue::glue("Full normalized counts written to: {full_gct_file}"), "INFO")
}

# Generate PCA plot for samples
generate_pca_plot <- function(normCounts, metadata_df, output_file) {
  log_message("Generating PCA plot...", "STEP")
  
  # Transpose for PCA
  pca_data <- t(normCounts)
  
  # Calculate PCA
  pca_result <- prcomp(pca_data, scale. = TRUE)
  
  # Extract variance explained
  variance_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  
  # Create data frame for plotting
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$sample <- rownames(pca_df)
  
  # Match with metadata
  pca_df <- merge(pca_df, metadata_df, by.x = "sample", by.y = "row.names")
  
  # Identify factor columns for coloring
  factor_cols <- sapply(metadata_df, is.factor)
  factor_names <- names(factor_cols)[factor_cols]
  
  # Pick first factor for coloring (excluding batch)
  color_by <- factor_names[!factor_names %in% "batch"]
  if (length(color_by) > 0) {
    color_by <- color_by[1]
  } else {
    color_by <- NULL
  }
  
  # Create PCA plot
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes_string(color = color_by), size = 3) +
    xlab(paste0("PC1: ", variance_explained[1], "% variance")) +
    ylab(paste0("PC2: ", variance_explained[2], "% variance")) +
    theme_minimal() +
    theme(legend.position = "right") +
    ggtitle("PCA Plot of Samples")
  
  # If we have batch information, create a plot with batch highlighted
  if ("batch" %in% colnames(metadata_df)) {
    batch_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = batch), size = 3) +
      xlab(paste0("PC1: ", variance_explained[1], "% variance")) +
      ylab(paste0("PC2: ", variance_explained[2], "% variance")) +
      theme_minimal() +
      theme(legend.position = "right") +
      ggtitle("PCA Plot by Batch")
    
    # Combine plots
    combined_plot <- gridExtra::grid.arrange(pca_plot, batch_plot, ncol = 1)
    
    # Save combined plot
    ggsave(output_file, combined_plot, width = 10, height = 14, units = "in")
  } else {
    # Save single plot
    ggsave(output_file, pca_plot, width = 10, height = 7, units = "in")
  }
  
  log_message(glue::glue("PCA plot saved to: {output_file}"), "SUCCESS")
}

# Generate sample distance heatmap
generate_sample_distance_heatmap <- function(normCounts, metadata_df, output_file) {
  log_message("Generating sample distance heatmap...", "STEP")
  
  # Calculate sample distances
  sample_dists <- dist(t(normCounts))
  sample_dist_matrix <- as.matrix(sample_dists)
  
  # Get annotation data for samples
  ann_data <- metadata_df
  rownames(ann_data) <- rownames(metadata_df)
  
  # Remove columns that have too many unique values (likely not useful for annotation)
  ann_cols <- sapply(ann_data, function(x) length(unique(x)) <= 10)
  ann_data <- ann_data[, ann_cols, drop = FALSE]
  
  # Create color palettes for annotations
  ann_colors <- list()
  for (col_name in colnames(ann_data)) {
    if (is.factor(ann_data[[col_name]]) || is.character(ann_data[[col_name]])) {
      levels <- unique(ann_data[[col_name]])
      colors <- colorRampPalette(RColorBrewer::brewer.pal(min(9, length(levels)), "Set1"))(length(levels))
      ann_colors[[col_name]] <- setNames(colors, levels)
    }
  }
  
  # Generate heatmap
  pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    annotation_col = ann_data,
    annotation_colors = ann_colors,
    main = "Sample Distances",
    filename = output_file,
    width = 10,
    height = 8
  )
  
  log_message(glue::glue("Sample distance heatmap saved to: {output_file}"), "SUCCESS")
}

# Generate MDS plot using Glimma
generate_mds_plot <- function(normCounts, metadata_df, output_prefix) {
  log_message("Generating interactive MDS plot with Glimma...", "STEP")
  
  # Create interactive MDS plot
  tryCatch({
    # Create a copy of metadata with row names as a column
    meta_for_mds <- metadata_df
    meta_for_mds$sample <- rownames(meta_for_mds)
    
    # Run Glimma MDS plot
    glimmaMDS(
      normCounts,
      groups = meta_for_mds,
      path = dirname(output_prefix),
      html = paste0(basename(output_prefix), ".html"),
      launch = FALSE
    )
    
    log_message(glue::glue("Interactive MDS plot saved to: {output_prefix}.html"), "SUCCESS")
  }, error = function(e) {
    log_message(glue::glue("Error generating MDS plot: {e$message}"), "WARNING")
  })
}

# Generate MA plot for LRT results
generate_ma_plot <- function(lrt_results, output_file, fdr_threshold = 0.1) {
  log_message("Generating MA plot for LRT results...", "STEP")
  
  # Create the MA plot
  tryCatch({
    pdf(output_file, width = 10, height = 8)
    DESeq2::plotMA(lrt_results, alpha = fdr_threshold, main = "MA Plot - LRT Results")
    dev.off()
    
    log_message(glue::glue("MA plot saved to: {output_file}"), "SUCCESS")
  }, error = function(e) {
    log_message(glue::glue("Error generating MA plot: {e$message}"), "WARNING")
  })
}

# Function to generate contrast dataframe for LRT
generate_contrast_dataframe <- function(metadata_df, reduced_formula, full_formula) {
  # Ensure formulas are formula objects
  if (is.character(reduced_formula)) {
    reduced_formula <- as.formula(reduced_formula)
  }
  if (is.character(full_formula)) {
    full_formula <- as.formula(full_formula)
  }
  
  # Get terms in full and reduced formulas
  full_terms <- attr(terms(full_formula), "term.labels")
  reduced_terms <- attr(terms(reduced_formula), "term.labels")
  
  # Find terms that are in full but not in reduced model
  contrast_terms <- setdiff(full_terms, reduced_terms)
  
  # If no difference, can't perform LRT test
  if (length(contrast_terms) == 0) {
    stop("No difference between full and reduced models. Cannot perform LRT.")
  }
  
  # Create design matrix for contrast extraction
  design <- model.matrix(full_formula, metadata_df)
  
  # Identify the coefficients corresponding to contrast terms
  coef_names <- colnames(design)
  
  # For interaction terms, we need to match partial names
  contrast_coefs <- c()
  contrast_labels <- c()
  
  for (term in contrast_terms) {
    # Replace : with .* for regex matching of interaction terms
    regex_term <- gsub(":", ".*", term)
    # Find matching coefficients
    matches <- grep(regex_term, coef_names, value = TRUE)
    contrast_coefs <- c(contrast_coefs, matches)
    # For each match, add the corresponding contrast term
    contrast_labels <- c(contrast_labels, rep(term, length(matches)))
  }
  
  # Create contrast dataframe
  contrast_df <- data.frame(
    coefficient = contrast_coefs,
    contrast = paste("LRT", contrast_labels, sep = "_"),
    stringsAsFactors = FALSE
  )
  
  return(contrast_df)
}

# Function to process the LRT results
process_lrt_results <- function(dds, metadata_df, full_formula, reduced_formula, alpha = 0.1) {
  log_message("Generate contrasts")
  
  # Run LRT test
  dds_lrt <- DESeq(dds, test = "LRT", reduced = reduced_formula, parallel = TRUE)
  
  # Get LRT results
  res_lrt <- results(dds_lrt, alpha = alpha)
  
  # Description of results
  log_message("LRT Results description")
  print(mcols(res_lrt, use.names = TRUE))
  
  # Print LRT summary
  print(summary(res_lrt))
  
  # Generate contrast dataframe for reference
  contrast_df <- generate_contrast_dataframe(metadata_df, reduced_formula, full_formula)
  
  # Return results
  list(
    dds_lrt = dds_lrt,
    res_lrt = res_lrt,
    contrast_df = contrast_df
  )
}

# Main function to export contrast results
export_contrast_results <- function(lrt_results, output_prefix) {
  # Extract results from the LRT analysis
  res_lrt <- lrt_results$res_lrt
  contrast_df <- lrt_results$contrast_df
  
  # Write contrast dataframe to file
  contrast_output_file <- paste0(output_prefix, "_contrasts.csv")
  write.csv(contrast_df, file = contrast_output_file, row.names = FALSE, quote = FALSE)
  log_message(paste("Exported contrast table to", contrast_output_file))
  
  # Write LRT results to file
  results_output_file <- paste0(output_prefix, "_lrt_results.csv")
  results_df <- as.data.frame(res_lrt)
  results_df$gene <- rownames(results_df)
  write.csv(results_df, file = results_output_file, row.names = FALSE, quote = FALSE)
  log_message(paste("Exported LRT results to", results_output_file))
  
  # Return file paths
  list(
    contrast_file = contrast_output_file,
    results_file = results_output_file
  )
}

# Separate generation of main effect and interaction contrasts
generate_effect_contrasts <- function(metadata_df, full_formula, dds_lrt = NULL, alpha = 0.1) {
  # Ensure formula is a formula object
  if (is.character(full_formula)) {
    full_formula <- as.formula(full_formula)
  }
  
  # Start timing
  start_time <- Sys.time()
  log_message("🔍 Debug: Starting contrast generation...", type = "DEBUG")
  
  # Get all terms in the formula
  all_terms <- attr(terms(full_formula), "term.labels")
  
  # Get unique factors from the design formula
  factors <- all.vars(full_formula)
  factors <- factors[factors %in% colnames(metadata_df)]
  
  log_message(paste("Factors identified for contrast analysis:", paste(factors, collapse = ", ")), type = "DEBUG")
  
  # Track all contrasts
  all_contrasts <- list()
  main_contrast_count <- 0
  interaction_contrast_count <- 0
  
  # For each factor, get the levels
  factor_levels <- lapply(factors, function(f) {
    levels(metadata_df[[f]])
  })
  names(factor_levels) <- factors
  
  # Log factor levels for debugging
  level_string <- paste(sapply(names(factor_levels), function(f) {
    paste0(f, ":", paste(factor_levels[[f]], collapse = ","))
  }), collapse = " | ")
  log_message(paste("Factor levels:", level_string), type = "DEBUG")
  
  # Generate main effect contrasts - for each factor within each level of other factors
  log_message("🔹 Generating main effect contrasts...", type = "DEBUG")
  main_start_time <- Sys.time()
  
  main_contrasts <- list()
  
  # For each factor, generate main effects contrasts
  for (factor_name in factors) {
    factor_values <- factor_levels[[factor_name]]
    if (length(factor_values) < 2) {
      next  # Need at least 2 levels for contrasts
    }
    
    # Generate pairwise comparisons for this factor
    if (length(factors) == 1) {
      # Simple case - only one factor
      contrasts <- expand.grid(
        numerator = factor_values, 
        denominator = factor_values, 
        stringsAsFactors = FALSE
      )
      contrasts <- contrasts[contrasts$numerator != contrasts$denominator, ]
      
      # Add metadata
      contrasts$effect <- "main"
      contrasts$specificity_group <- factor_name
      contrasts$contrast <- paste0(factor_name, "_", contrasts$numerator, "_vs_", contrasts$denominator)
      
      main_contrasts[[length(main_contrasts) + 1]] <- contrasts
    } else {
      # Complex case - multiple factors, create contrasts within each level of other factors
      
      # Get all other factors
      other_factors <- factors[factors != factor_name]
      
      # For each combination of other factor levels
      other_factor_combinations <- expand.grid(lapply(other_factors, function(f) factor_levels[[f]]))
      colnames(other_factor_combinations) <- other_factors
      
      for (i in 1:nrow(other_factor_combinations)) {
        combo <- other_factor_combinations[i, , drop = FALSE]
        
        # Create specificity group based on the combination
        spec_group_parts <- c()
        for (of in other_factors) {
          spec_group_parts <- c(spec_group_parts, paste0(of, "_", combo[[of]]))
        }
        specificity_group <- paste(spec_group_parts, collapse = "_")
        
        # Create contrasts for this factor within this combination
        contrasts <- expand.grid(
          numerator = factor_values, 
          denominator = factor_values, 
          stringsAsFactors = FALSE
        )
        contrasts <- contrasts[contrasts$numerator != contrasts$denominator, ]
        
        # Add metadata
        contrasts$effect <- "main"
        contrasts$specificity_group <- specificity_group
        contrasts$contrast <- paste0(factor_name, "_", contrasts$numerator, "_vs_", contrasts$denominator)
        
        main_contrasts[[length(main_contrasts) + 1]] <- contrasts
      }
    }
  }
  
  # Combine all main contrasts
  main_contrasts_df <- do.call(rbind, main_contrasts)
  if (!is.null(main_contrasts_df) && nrow(main_contrasts_df) > 0) {
    main_contrast_count <- nrow(main_contrasts_df)
    all_contrasts$main <- main_contrasts_df
  }
  
  log_message(paste("✅ Main effect contrasts generated in", round(difftime(Sys.time(), main_start_time, units = "secs"), 3), "sec"), type = "DEBUG")
  
  # Generate interaction contrasts if there's more than one factor
  if (length(factors) > 1) {
    log_message("🔹 Generating interaction effect contrasts...", type = "DEBUG")
    interaction_start_time <- Sys.time()
    
    # Look for interaction terms in the formula
    interaction_terms <- all_terms[grepl(":", all_terms)]
    
    if (length(interaction_terms) > 0) {
      all_interaction_contrasts <- list()
      
      # For each interaction term
      for (term in interaction_terms) {
        # Split into constituent factors
        int_factors <- strsplit(term, ":")[[1]]
        
        # Check if all factors exist in metadata
        if (!all(int_factors %in% colnames(metadata_df))) {
          next
        }
        
        # For each pair of interacting factors
        for (i in 1:length(int_factors)) {
          for (j in 1:length(int_factors)) {
            if (i == j) next
            
            factor1 <- int_factors[i]
            factor2 <- int_factors[j]
            
            # Get levels for each factor
            levels1 <- factor_levels[[factor1]]
            levels2 <- factor_levels[[factor2]]
            
            # Generate contrasts for factor1 at each level of factor2
            for (level2 in levels2) {
              # Contrasts between levels of factor1 at this level of factor2
              for (level1a in levels1) {
                for (level1b in levels1) {
                  if (level1a == level1b) next
                  
                  # Create the interaction contrast
                  interaction_contrasts <- data.frame(
                    effect = "interaction",
                    specificity_group = paste0(factor1, "_", level1a, "_vs_", level1b),
                    numerator = paste0(factor2, level2),
                    denominator = paste0(factor2, levels2[1]), # Use first level as reference
                    contrast = paste0(factor1, level1a, ".", factor2, level2),
                    stringsAsFactors = FALSE
                  )
                  
                  all_interaction_contrasts[[length(all_interaction_contrasts) + 1]] <- interaction_contrasts
                }
              }
            }
            
            # Generate contrasts for factor2 at each level of factor1
            for (level1 in levels1) {
              # Contrasts between levels of factor2 at this level of factor1
              for (level2a in levels2) {
                for (level2b in levels2) {
                  if (level2a == level2b) next
                  
                  # Create the interaction contrast
                  interaction_contrasts <- data.frame(
                    effect = "interaction",
                    specificity_group = paste0(factor2, "_", level2a, "_vs_", level2b),
                    numerator = paste0(factor1, level1),
                    denominator = paste0(factor1, levels1[1]), # Use first level as reference
                    contrast = paste0(factor1, level1, ".", factor2, level2a),
                    stringsAsFactors = FALSE
                  )
                  
                  all_interaction_contrasts[[length(all_interaction_contrasts) + 1]] <- interaction_contrasts
                }
              }
            }
          }
        }
      }
      
      # Combine all interaction contrasts
      interaction_contrasts_df <- do.call(rbind, all_interaction_contrasts)
      if (!is.null(interaction_contrasts_df) && nrow(interaction_contrasts_df) > 0) {
        interaction_contrast_count <- nrow(interaction_contrasts_df)
        all_contrasts$interaction <- interaction_contrasts_df
      }
    }
    
    log_message(paste("✅ Interaction contrasts generated in", round(difftime(Sys.time(), interaction_start_time, units = "secs"), 3), "sec"), type = "DEBUG")
  }
  
  # Log summary
  log_message(paste("Total main contrasts generated:", main_contrast_count), type = "DEBUG")
  log_message(paste("Total interaction contrasts generated:", interaction_contrast_count), type = "DEBUG")
  log_message(paste("✅ Debug: Contrast generation completed in", round(difftime(Sys.time(), start_time, units = "secs"), 2), "sec"), type = "DEBUG")
  
  # Return contrast dataframe
  all_contrasts_df <- NULL
  if (!is.null(all_contrasts$main) && nrow(all_contrasts$main) > 0) {
    all_contrasts_df <- all_contrasts$main
  }
  
  if (!is.null(all_contrasts$interaction) && nrow(all_contrasts$interaction) > 0) {
    if (is.null(all_contrasts_df)) {
      all_contrasts_df <- all_contrasts$interaction
    } else {
      all_contrasts_df <- rbind(all_contrasts_df, all_contrasts$interaction)
    }
  }
  
  # Add contrast number
  if (!is.null(all_contrasts_df)) {
    all_contrasts_df$contrast_number <- seq_len(nrow(all_contrasts_df))
    
    # Count significant genes for each contrast if dds_lrt is provided
    if (!is.null(dds_lrt)) {
      all_contrasts_df$significant_genes <- 0
      
      # Get LRT results
      res_lrt <- results(dds_lrt, alpha = alpha)
      
      # Count significant genes
      sig_count <- sum(!is.na(res_lrt$padj) & res_lrt$padj < alpha, na.rm = TRUE)
      
      # Assign same count to all contrasts for LRT results
      all_contrasts_df$significant_genes <- sig_count
    } else {
      all_contrasts_df$significant_genes <- NA
    }
    
    # Reorder columns to match expected format
    all_contrasts_df <- all_contrasts_df[, c("contrast_number", "effect", "contrast", 
                                           "specificity_group", "numerator", "denominator", 
                                           "significant_genes")]
  }
  
  return(all_contrasts_df)
}

# Function to export contrast results with differential expression data
export_contrasts_with_de_data <- function(lrt_results, dds, output_prefix, alpha = 0.1) {
  # Extract components
  res_lrt <- lrt_results$res_lrt
  dds_lrt <- lrt_results$dds_lrt
  
  # Get the contrast information
  contrasts_df <- generate_effect_contrasts(colData(dds), design(dds))
  
  if (is.null(contrasts_df) || nrow(contrasts_df) == 0) {
    log_message("No contrasts generated. Skipping export.", type = "WARNING")
    return(NULL)
  }
  
  # Save basic contrast information
  contrast_output_file <- paste0(output_prefix, "_contrasts.csv")
  write.csv(contrasts_df, file = contrast_output_file, row.names = FALSE, quote = FALSE)
  log_message(paste("Exported contrast table to", contrast_output_file))
  
  # Get normalized counts for export
  norm_counts <- counts(dds_lrt, normalized = TRUE)
  
  # Convert LRT results to data frame
  lrt_results_df <- as.data.frame(res_lrt)
  lrt_results_df$gene <- rownames(lrt_results_df)
  
  # Write LRT results
  lrt_output_file <- paste0(output_prefix, "_lrt_results.csv")
  write.csv(lrt_results_df, file = lrt_output_file, row.names = FALSE, quote = FALSE)
  log_message(paste("Exported LRT results to", lrt_output_file))
  
  # Return file paths
  return(list(
    contrast_file = contrast_output_file,
    lrt_results_file = lrt_output_file
  ))
}

