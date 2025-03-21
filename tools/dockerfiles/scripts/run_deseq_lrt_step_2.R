#!/usr/bin/env Rscript
options(warn = -1)
options("width" = 300)

suppressMessages(library(argparse))
suppressMessages(library(BiocParallel))
suppressMessages(library(pheatmap))
suppressMessages(library(Glimma))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(hopach))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(DESeq2))
suppressMessages(library(limma)) # For removeBatchEffect
suppressMessages(library(cmapR))
suppressMessages(library(data.table))
suppressMessages(library(rlang)) # For traceback handling

options(rlang_backtrace_on_error = "full")

mutate <- dplyr::mutate
filter <- dplyr::filter
group_by <- dplyr::group_by
slice <- dplyr::slice
rename <- dplyr::rename
select <- dplyr::select
arrange <- dplyr::arrange
distinct <- dplyr::distinct
`%>%` <- magrittr::`%>%`
`%in%` <- base::`%in%`
`%/%` <- base::`%/%`
`%%` <- base::`%%`

##########################################################################################
#
# v0.0.5 LRT Step 2
#
# Changes:
# - Applied 'model' batch correction in Step 2 if specified.
# - Used batch information from Step 1.
# - Ensured necessary data is passed between steps.
#
################################################################################################
# v0.0.3 - Modified to save all possible information into RDS files and handle batch correction
#
# Changes:
# - Added batch correction options (`CombatSeq` and `model`).
# - Modified design formula based on batch correction method.
# - Saved all necessary data (dds, contrasts, expression data, metadata) into RDS files.
#
##########################################################################################
#
# v0.0.2 LRT Step 2
#
# Changes:
# - Modified to read RDS files generated from the first script.
# - Inputs adjusted to accept RDS files and contrast indices.
# - Removed redundant code; now uses precomputed data.
#
##########################################################################################
#
# v0.0.1 LRT Step 2
#
# - Implemented **contrast generation** for complex experimental designs:
#   - Supports **main effects** (e.g., factor1 levels) and **interaction effects** (factor1 vs factor2).
#   - Automatically generates contrasts based on user-defined factors and levels.
#   - Designed for comparing multiple conditions and extracting meaningful differential gene expression patterns.
#
# - Added **log2 fold-change threshold** (default: 0.59, ~1.5 fold change) for filtering significant genes.
#
# - Integrated **batch correction** with two options:
#   - `combatseq`: Corrects for batch effects prior to differential analysis.
#   - `model`: Removes batch effects post-analysis using limma.
#
# - Visual outputs include:
#   - **Heatmaps** displaying top 30 most variable genes.
#   - **MDS plots** for visualizing sample distances.
#   - **MA plots** for visualizing differential expression.
#
# - **Clustering**: Supports HOPACH clustering with min-max scaling (row/column/both).
#
# - Exported **GCT files** compatible with GSEA, for both filtered and unfiltered data.
#
# - Multi-threading support through BiocParallel for faster processing.
#
# - Handles input in CSV/TSV formats and outputs standardized files.
#
##########################################################################################
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
# The output file's rows order corresponds to the rows order of the first CSV/TSV file in
# the untreated group. Output is always saved in TSV format
#
# Output file includes only intersected rows from all input files. Intersected by
# RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand
#
# DESeq/DESeq2 always compares untreated_vs_treated groups (condition-1-vs-condition-2)
#
# Additionally we calculate -LOG10(pval) and -LOG10(padj)
#
# Use -un and -tn to set custom names for treated and untreated conditions
#
# Use -ua and -ta to set aliases for input expression files. Should be unique
# Exports GCT files to be used by GSEA. GCT files is always with uppercase GeneId
#
##########################################################################################


READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")

get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq2 analysis using contrasts from previous LRT step")
  parser$add_argument(
    "-dsq",
    "--dsq_obj_data",
    help = "RDS file containing contrassts and expression data from step 1",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-cd",
    "--contrast_df",
    help     = "TSV file containing contrasts data",
    type     = "character",
    required = TRUE
  )
  parser$add_argument(
    "-bcm",
    "--batchcorrection",
    help = "RDS file containing the batch correction method used in step 1",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-ci",
    "--contrast_indices",
    help = "Comma-separated list of integers representing contrast indices (e.g., 1,2,3)",
    type = "character",
    required = TRUE
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
    action = "store_true"
  )
  parser$add_argument(
    "--regulation",
    help = paste(
      "Direction of differential expression comparison. β is the log2 fold change.",
      "'both' for both up and downregulated genes (|β| > lfcThreshold); ",
      "'up' for upregulated genes (β > lfcThreshold); ",
      "'down' for downregulated genes (β < -lfcThreshold). ",
      "Default: both"
    ),
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
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
    help = "Output prefix. Default: deseq",
    type = "character",
    default = "./deseq"
  )
  parser$add_argument(
    "-p",
    "--threads",
    help = "Threads",
    type = "integer",
    default = 1
  )
  parser$add_argument(
    "--test_mode",
    help = "Run for test, only first 500 rows, clustering minimised",
    action = "store_true",
    default = FALSE
  )
  args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
  return(args)
}

# Load arguments
args <- get_args()

# Set threads
register(MulticoreParam(args$threads))

# Function to log messages with timestamps
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

rebuild_dds <- function(dds, factor_ref_list) {
  # 1) Convert colData to a plain data.frame so we can manipulate it
  colData_df <- as.data.frame(colData(dds))

  # 2) Re-level factors as specified in factor_ref_list, e.g.
  #    factor_ref_list might be list(grouped_ref = c("grouped" = "0H_GFP_N"))
  for (f_name in names(factor_ref_list)) {
    # factor_ref_list[[f_name]] is the new reference level for factor f_name
    ref_level            <- factor_ref_list[[f_name]]
    colData_df[[f_name]] <- relevel(as.factor(colData_df[[f_name]]), ref = ref_level)
  }

  # 3) Build a new DESeqDataSet
  dds_new <- DESeqDataSetFromMatrix(
    countData = counts(dds),
    colData   = colData_df,
    design    = design(dds)
  )
  return(dds_new)
}


# Load contrasts ta
contrast_df <- data.table::fread(args$contrast_df)


# Load RDS files with detailed logging
log_message(paste("Loading contrasts from", args$dsq_obj_data))
all_contrasts <- readRDS(args$dsq_obj_data)
log_message("Contrasts Loaded:")
print(str(all_contrasts, max.level = 0))
print(all_contrasts)
log_message("Structure of Contrasts:")
glimpse(all_contrasts)


contrasts_values <- sapply(all_contrasts, function(x) x$contrast)
print("Contrasts:")
print(contrasts_values)
spec_group_values <- sapply(all_contrasts, function(x) x$specificity_group)
print("Specificity Groups:")
print(spec_group_values)


log_message(paste("Loading expression data from", args$dsq_obj_data))
expression_data_df <- all_contrasts$expression_data_df
all_contrasts <- purrr::list_modify(all_contrasts, expression_data_df = NULL)

log_message("Expression Data Loaded:")
print(head(expression_data_df))
log_message("Structure of Expression Data:")
glimpse(expression_data_df)

if (!is.null(args$rpkm_cutoff)) {
  print("Using RPKM cutoff for additional filtering of expression data:")
  print(args$rpkm_cutoff)
  expression_data_df <- filter_rpkm(expression_data_df, args$rpkm_cutoff)
  print("Expression data after RPKM filtering: ")
  print(head(expression_data_df))
  print(dim(expression_data_df))
}

log_message(paste("Loading DESeq2 object from Contrasts"))
dds <- all_contrasts$deseq_obj
all_contrasts <- purrr::list_modify(all_contrasts, deseq_obj = NULL)

log_message("DESeq2 Object Loaded:")
print(dds)
log_message("Sample Names in DESeq2 Object:")
print(colnames(dds))

log_message(paste("Loading metadata from DESeq2 object"))
metadata_df <- colData(dds)
log_message("Metadata Loaded:")
print(head(metadata_df))
log_message("Structure of Metadata:")
glimpse(metadata_df)

batch_correction_method <- args$batchcorrection
log_message(paste("Batch correction method used:", batch_correction_method))

# Convert the comma-separated string into a vector of integers
contrast_vector <- as.integer(unlist(strsplit(args$contrast_indices, ",")))
log_message(paste("Contrast Indices Parsed:", paste(contrast_vector, collapse = ", ")))

# Ensure that the contrasts indices are valid
if (any(contrast_vector > length(all_contrasts))) {
  stop("One or more contrast indices are out of bounds.")
}

# Harmonize sample names
log_message("Starting harmonization of sample names...")

# Define READ_COL correctly based on your actual column naming
READ_COL <- "TotalReads" # Ensure this matches your read counts column pattern

# Extract read counts columns
read_counts_columns <- grep(
  paste(READ_COL, sep = ""),
  colnames(expression_data_df),
  value = TRUE,
  ignore.case = TRUE
)

log_message(paste("Columns matching READ_COL pattern ('", READ_COL, "'):", paste(read_counts_columns, collapse = ", ")))

read_counts_data_df <- expression_data_df[read_counts_columns]
log_message("Read Counts Data Extracted:")
print(head(read_counts_data_df))
log_message("Structure of Read Counts Data:")
glimpse(read_counts_data_df)

# Clean and standardize column names
log_message("Cleaning and standardizing column names...")
original_colnames <- colnames(read_counts_data_df)
log_message("Original Column Names:")
print(original_colnames)

cleaned_colnames <- sapply(original_colnames, function(s) {
  # Remove the last word if separated by space
  words <- unlist(strsplit(s, " ", fixed = TRUE))
  if (length(words) > 1) {
    paste(head(words, -1), collapse = "_")
  } else {
    s
  }
})

# Convert to lowercase and trim whitespace
cleaned_colnames <- tolower(cleaned_colnames)
cleaned_colnames <- trimws(cleaned_colnames)
cleaned_colnames <- gsub("[^[:alnum:]_]", "", cleaned_colnames)

log_message("Cleaned Column Names:")
print(cleaned_colnames)

# Assign cleaned column names back to the dataframe
colnames(read_counts_data_df) <- cleaned_colnames
log_message("Read Counts Data with Cleaned Column Names:")
print(head(read_counts_data_df))
log_message("Sample Names in Read Counts Data:")
print(colnames(read_counts_data_df))

log_message("Sample Names in DESeq2 Object After Cleaning:")
print(colnames(dds))
print("Checking for exact match between DESeq2 object and read counts data...")
print(dds)

# Check for exact match (regardless of order)
if (!setequal(colnames(dds), colnames(read_counts_data_df))) {
  # Identify mismatched samples
  missing_in_read_counts <- setdiff(colnames(dds), colnames(read_counts_data_df))
  missing_in_deseq <- setdiff(colnames(read_counts_data_df), colnames(dds))

  if (length(missing_in_read_counts) > 0) {
    log_message("Samples in DESeq2 object but not in read counts data:")
    print(missing_in_read_counts)
  }

  if (length(missing_in_deseq) > 0) {
    log_message("Samples in read counts data but not in DESeq2 object:")
    print(missing_in_deseq)
  }

  stop("Sample names in DESeq2 object and read counts data do not match.")
} else {
  # Reorder read counts columns to match DESeq2 object
  read_counts_data_df <- read_counts_data_df[, colnames(dds)]
  log_message("Reordered read counts data to match DESeq2 object sample names.")

  # Verify the reordering
  if (!all(colnames(dds) == colnames(read_counts_data_df))) {
    stop("Reordering failed: Sample names still do not match.")
  } else {
    log_message("Read counts data successfully reordered to match DESeq2 object.")
  }
}

# Apply limma batch correction if specified
if (batch_correction_method == "model" && "batch" %in% colnames(metadata_df)) {
  # Case 3: After DESeq2, apply rlog transformation and remove batch effects using limma
  print("Applying rlog transformation and limma batch effect removal")
  rlog_transformed <- rlog(dds, blind = FALSE)
  rlog_counts <- assay(rlog_transformed)

  # Prepare design matrix without 'batch' for removeBatchEffect
  design_formula <- as.formula(paste0("~", str_remove(as.character(dds@design)[2], " \\+ batch")))

  design_matrix <- model.matrix(design_formula, data = metadata_df)

  # Apply removeBatchEffect on rlog normalized counts
  normCounts <- limma::removeBatchEffect(rlog_counts, batch = metadata_df$batch, design = design_matrix)
} else {
  # If no batch correction is needed, just perform rlog
  normCounts <- assay(rlog(dds, blind = FALSE))
}

# Updated get_contrast_res function that handles main effects (single & multi)
# as well as interaction contrasts.
get_contrast_res <- function(dds_base, contrast_row) {
  # Global args are assumed to be defined: args$fdr, args$use_lfc_thresh, args$lfcthreshold, args$regulation
  altHypothesis <- switch(args$regulation,
                          "up"   = "greater",
                          "down" = "less",
                          "both" = "greaterAbs",
                          "greaterAbs")  # default

  lfcThreshold <- if (args$use_lfc_thresh) args$lfcthreshold else 0

  if (contrast_row$effect_type == "main") {
    # --- MAIN EFFECT CONTRAST ---
    # Check if this is a single-factor or multi-factor scenario.
    col_names <- colnames(colData(dds_base))
    if (contrast_row$specificity_group %in% col_names) {
      # SINGLE-FACTOR: specificity_group is the factor to re-level.
      factor_ref_list <- setNames(list(contrast_row$denominator), contrast_row$specificity_group)
    } else {
      # MULTI-FACTOR: specificity_group is composite, e.g. "otherFactor_otherLevel".
      parts <- strsplit(contrast_row$specificity_group, "_")[[1]]
      if (length(parts) < 2) {
        stop("Unexpected format for specificity_group in multi-factor contrast")
      }
      # The first part is the name of the additional factor.
      other_factor <- parts[1]
      other_level  <- paste(parts[-1], collapse = "_")

      # Extract main factor from the contrast name (assumed format: "mainFactor_lvl_vs_ref")
      main_parts <- strsplit(contrast_row$contrast, "_")[[1]]
      main_factor <- main_parts[1]

      factor_ref_list <- c(
        setNames(list(contrast_row$denominator), main_factor),
        setNames(list(other_level), other_factor)
      )
    }

    # Rebuild the DESeqDataSet with proper reference levels and run DESeq2.
    dds_subset <- rebuild_dds(dds_base, factor_ref_list)
    dds_subset <- DESeq(dds_subset, test = "Wald")

    message("DESeq object subset for main effect contrast: ", contrast_row$contrast)
    print(dds_subset)
    print(resultsNames(dds_subset))

    res <- results(dds_subset,
                   name                 = contrast_row$contrast,
                   alpha                = args$fdr,
                   lfcThreshold         = lfcThreshold,
                   independentFiltering = TRUE,
                   altHypothesis        = altHypothesis)

    message("DESeq2 subset results for main effect contrast obtained.")
    return(res)

  } else if (contrast_row$effect_type == "interaction") {
    # --- INTERACTION CONTRAST ---
    # Use the helper function to extract factor names and their levels.
    factors_levels <- extract_factors_and_levels(contrast_row$contrast)
    factor1 <- factors_levels$factor1
    level1  <- factors_levels$level1
    factor2 <- factors_levels$factor2
    level2  <- factors_levels$level2

    # Determine which factor was re-leveled during generation by examining specificity_group.
    # In generate_interaction_effect_contrasts:
    #   - In the first branch, specificity_group is built as: paste(factor2, level2, "vs", ref_level2)
    #     and re-leveling is applied to factor1.
    #   - In the second branch, specificity_group is: paste(factor1, level1, "vs", ref_level1)
    #     and re-leveling is applied to factor2.
    sg_parts <- strsplit(contrast_row$specificity_group, "_")[[1]]
    if (sg_parts[1] == factor2) {
      # Then the intended re-leveling is for factor1.
      relevel_factor <- factor1
      # The stored denominator is of the form paste0(factor1, ref_level1);
      # remove the factor name prefix to get the actual ref level.
      ref_level <- sub(paste0("^", factor1), "", contrast_row$denominator)
    } else if (sg_parts[1] == factor1) {
      # Then re-level factor2.
      relevel_factor <- factor2
      ref_level <- sub(paste0("^", factor2), "", contrast_row$denominator)
    } else {
      stop("Could not determine releveling factor for interaction contrast based on specificity_group")
    }

    factor_ref_list <- setNames(list(ref_level), relevel_factor)

    # Optionally, you might consider subsetting the DESeqDataSet based on the other factor's level
    # (as commented in your generation code). For now we simply rebuild the DESeqDataSet.
    dds_subset <- rebuild_dds(dds_base, factor_ref_list)
    dds_subset <- DESeq(dds_subset, test = "Wald")

    message("DESeq object subset for interaction contrast: ", contrast_row$contrast)
    print(dds_subset)
    print(resultsNames(dds_subset))

    res <- results(dds_subset,
                   name                 = contrast_row$contrast,
                   alpha                = args$fdr,
                   lfcThreshold         = lfcThreshold,
                   independentFiltering = TRUE,
                   altHypothesis        = altHypothesis)

    message("DESeq2 subset results for interaction contrast obtained.")
    return(res)

  } else {
    stop("Unknown effect_type in contrast_row")
  }
}

# Function to export MDS plot
export_mds_html_plot <- function(norm_counts_data, location) {
  tryCatch(
    expr = {
      htmlwidgets::saveWidget(
        glimmaMDS(
          x = norm_counts_data,
          groups = as.data.frame(metadata_df),
          labels = rownames(metadata_df)
        ),
        file = location
      )
    },
    error = function(e) {
      print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
    }
  )
}

# Function to export CLS files
export_cls <- function(categories, location) {
  tryCatch(
    expr = {
      output_stream <- file(location, "w")
      on.exit(close(output_stream), add = TRUE)
      cat(
        paste(length(categories),
              length(levels(categories)),
              "1",
              sep = "\t"
        ),
        paste("#", paste(
          unique(as.character(categories)),
          collapse = "\t"
        ), sep = "\t"),
        paste(paste(as.character(categories),
                    collapse = "\t"
        ), sep = "\t"),
        file = output_stream,
        sep = "\n"
      )
      print(paste("Exporting CLS data to", location,
                  sep = " "
      ))
    },
    error = function(e) {
      print(paste("Failed to export CLS data to ", location, "with error - ", e,
                  sep = ""
      ))
    }
  )
}

# Function to generate clusters
get_clustered_data <- function(expression_data, by = "row", k = 3, kmax = 5, dist = "cosangle", scaling_type = "zscore") {

  start_time <- proc.time()

  if (!(by %in% c("row", "col"))) {
    stop("Invalid value for 'by'. Choose either 'row' or 'col'.")
  }

  # If clustering by columns, transpose so that columns become rows for scaling.
  if (by == "col") {
    print("Transposing expression data to scale columns")
    expression_data <- t(expression_data)
  }

  # Apply scaling per row (which is the desired unit, either original rows or columns)
  expression_data <- switch(scaling_type,
    "minmax" = t(apply(expression_data, 1, scale_min_max)),
    "zscore" = {
      scaled_data <- t(scale(t(expression_data), center = TRUE, scale = TRUE))
      scaled_data[is.na(scaled_data)] <- 0  # Handle zero-variance rows
      scaled_data
    },
    stop("Invalid scaling type. Choose 'minmax' or 'zscore'.")
  )

  # If data was transposed for column scaling, transpose back to original orientation.
  if (by == "col") {
    print("Transposing expression data back to original orientation")
    expression_data <- t(expression_data)
  }

  print(paste0("Running HOPACH for ", nrow(expression_data), "  features"))
  hopach_results <- hopach::hopach(expression_data,
                                   verbose = TRUE,
                                   K       = k,
                                   kmax    = kmax,
                                   khigh   = kmax,
                                   d = dist
  )

  print("Parsing cluster labels")
  # hopach_results$clustering$labels gives final cluster labels as integers
  # hopach returns them as numeric without the "c" prefix, so we add "c" ourselves or just rely on numeric.
  # Actually hopach by default returns numeric labels (1,11,12...). We'll add "c" prefix ourselves for consistency.

  # Final labels (no prefix 'c' in hopach by default)
  final_labels      <- hopach_results$clustering$labels[hopach_results$clustering$order]
  # Convert to character
  final_labels_char <- as.character(final_labels)

  # Add "c" prefix if desired (optional)
  # final_labels_char <- paste0("c", final_labels_char)

  # Each digit in the label corresponds to a level.
  # Determine max number of levels
  max_levels <- max(nchar(final_labels_char))

  # Create a data frame for clusters
  clusters           <- data.frame(Label = final_labels_char, stringsAsFactors = FALSE)
  rownames(clusters) <- rownames(expression_data)[hopach_results$clustering$order]

  # Split labels into levels
  # For each label, we split into characters and assign to new columns
  level_data <- do.call(rbind, lapply(clusters$Label, function(lbl) {
    # Split into individual characters
    chars <- unlist(strsplit(lbl, split = ""))
    # If shorter than max_levels, pad with NA
    if (length(chars) < max_levels) {
      chars <- c(chars, rep(NA, max_levels - length(chars)))
    }
    return(chars)
  }))

  # Name the columns
  colnames(level_data) <- paste0("Cluster_Level_", seq_len(max_levels))

  # Combine into clusters
  clusters <- cbind(clusters, level_data)

  # Optionally remove the original 'Label' column if not needed
  # Or rename it to something else
  clusters$HCL   <- paste0("c", clusters$Label)
  clusters$Label <- NULL

  end_time <- proc.time() - start_time

  print("Total time of get_clustered_data function execution: ")
  print(end_time)

  return(list(
    order      = hopach_results$clustering$order,
    expression = expression_data,
    clusters   = clusters
  ))
}

# Function for min-max scaling
scale_min_max <- function(x,
                          min_range = -2,
                          max_range = 2) {
  min_val <- min(x)
  max_val <- max(x)
  scaled_x <-
    (x - min_val) / (max_val - min_val) * (max_range - min_range) + min_range
  return(scaled_x)
}


filter_rpkm <- function(expression_df, n) {
  expression_df %>%
    filter(if_any(contains("Rpkm"), ~. > n))
}

# Function to cluster data and re-order based on clustering results
cluster_and_reorder <- function(normCounts, col_metadata, row_metadata, args) {

  start_time <- proc.time()

  if (args$cluster != "none") {
    if (args$test_mode) {
        k    <- 2
        kmax <- 2
      } else {
        k    <- args$k
        kmax <- args$kmax
      }
    # Column clustering if requested
    if (args$cluster == "column" || args$cluster == "both") {
      clustered_data_cols <- get_clustered_data(normCounts, by = "col", k = k, kmax = kmax, scaling_type = args$scaling_type, dist = args$columndist)
      normCounts          <- normCounts[, clustered_data_cols$order, drop = FALSE]
      col_metadata        <- col_metadata[clustered_data_cols$order, , drop = FALSE]
      # After reordering, cbind cluster info
      col_metadata        <- cbind(col_metadata, clustered_data_cols$clusters)
    }
    # Row clustering if requested
    if (args$cluster == "row" || args$cluster == "both") {
      clustered_data_rows <- get_clustered_data(normCounts, by = "row", k = k, kmax = kmax, scaling_type = args$scaling_type, dist = args$rowdist)
      normCounts          <- clustered_data_rows$expression[clustered_data_rows$order, , drop = FALSE]
      row_metadata        <- row_metadata[clustered_data_rows$order, , drop = FALSE]
      # After reordering rows, add cluster annotations
      row_metadata        <- cbind(row_metadata, clustered_data_rows$clusters)
    }
  } else {
    # No clustering
  }

  end_time <- proc.time() - start_time

  print("Total time of execution cluster_and_reorder function: ")
  print(end_time)

  print("Clustered data:")
  print(head(normCounts))
  print("Row metadata:")
  print(head(row_metadata))
  print("Column metadata:")
  print(head(col_metadata))

  return(list(normCounts = normCounts, col_metadata = col_metadata, row_metadata = row_metadata))
}

# TODO: the following funcion should be updated to handle the any contrast data with modifiee LFC and padj column
#  names into appropriate contrast like contrast_1_LFC and so on and the filtering should be then based on either of
#  them

# Function to add metadata columns and create the final results data frame
# Updated function: add contrast-specific metadata to results
add_metadata_to_results <- function(expression_data_df, deseq_result, read_count_cols, contrast_name) {
  res_df <- as.data.frame(deseq_result) %>%
    select(baseMean, log2FoldChange, pvalue, padj) %>%
    rename(
      !!paste0(contrast_name, "_LFC")    := log2FoldChange,
      !!paste0(contrast_name, "_pvalue") := pvalue,
      !!paste0(contrast_name, "_FDR")    := padj
    )

  print("DESeq2 results data frame:")
  print(head(res_df))

  # Merge with the main expression data
  expression_data_merged <- data.frame(
    cbind(expression_data_df[, !colnames(expression_data_df) %in% read_count_cols], res_df),
    check.names = FALSE,
    check.rows  = FALSE
  )
  return(expression_data_merged)
}

# Function to export DESeq2 report
export_deseq_report <- function(expression_data_df, output_prefix) {
  expression_data_df_filename <- paste0(output_prefix, "_gene_exp_table.tsv")
  write.table(
    expression_data_df,
    file = expression_data_df_filename,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  print(paste("Export DESeq report to", expression_data_df_filename))
}

# TODO: here filtering should be based on either of LFC or padj column for ANY of the contrast (like contrast_1_padj
#  <= fdr | contrast_2_padj <= fdr, etc.)
# Function to export normalized counts and filtered counts to GCT format
# Updated GCT export: filter if ANY FDR column meets the threshold
export_gct_data <- function(normCounts, row_metadata, col_metadata) {
  tryCatch({
    # Ensure col_metadata columns are vectors
    col_metadata <- col_metadata %>% mutate_all(as.vector)

    # Initialize col_order_vector (will be non-NULL if user specifies a column order)
    col_order_vector <- tolower(rownames(metadata_df))
    print("Order of columns for heatmap without clustering, based on metadata input:")
    print(col_order_vector)

    normCounts   <- normCounts[, col_order_vector, drop = FALSE]
    col_metadata <- col_metadata[col_order_vector, , drop = FALSE]

    gct_data     <- new("GCT", mat = normCounts, rdesc = row_metadata, cdesc = col_metadata)
    cmapR::write_gct(ds = gct_data, ofile = "counts_all.gct", appenddim = FALSE)
    print(paste("Exported GCT (all) to", "counts_all.gct"))

    # Filter rows by any FDR column
    fdr_cols <- grep("_FDR$", colnames(row_metadata), value = TRUE)
    print("FDR columns for filtering:")
    print(fdr_cols)
    lfc_cols <- grep("_LFC$", colnames(row_metadata), value = TRUE)
    print("LFC columns for filtering:")
    print(lfc_cols)
    if (length(fdr_cols) == 0) {
      warning("No FDR columns found. No filtering performed for counts_filtered.gct.")
      row_metadata_filtered <- row_metadata
    } else {
      row_metadata_filtered <- row_metadata %>%
        filter(if_any(all_of(fdr_cols), ~. <= args$fdr)) %>%
        filter(if_any(all_of(abs(lfc_cols)), ~. >= args$lfcthreshold))
    }

    if (nrow(row_metadata_filtered) == 0) {
      warning(paste("No genes passed the FDR/LFC thresholds of", args$fdr, "/", args$lfcthreshold, "and will be exported as empty GCT."))
    }

    filtered_normCounts <- normCounts[rownames(row_metadata_filtered), , drop = FALSE]

    # Cluster filtered data
    clustered_data <- cluster_and_reorder(filtered_normCounts, col_metadata, row_metadata_filtered, args)

    gct_data_filtered <- new("GCT",
                             mat   = clustered_data$normCounts,
                             rdesc = clustered_data$row_metadata,
                             cdesc = clustered_data$col_metadata)
    cmapR::write_gct(ds = gct_data_filtered, ofile = "counts_filtered.gct", appenddim = FALSE)
    print(paste("Exported GCT (filtered) to", "counts_filtered.gct"))
  }, error = function(e) {
    print(paste("Failed to export GCT data:", e$message))
  })
}

# Iterate through each index in contrast_vector and process the corresponding contrast
annotated_expression_combined_df <- NULL

print("Starting DESeq2 analysis for contrast_vector:")
print(contrast_vector)

for (contrast_index in contrast_vector) {

  contrast_df_subset <- filter(contrast_df, contrast_number == contrast_index)

  contrast_index <- which(contrasts_values == contrast_df_subset$contrast & spec_group_values ==
    contrast_df_subset$specificity_group)
  
  contrast_selected <- all_contrasts[[contrast_index]]

  print("Selected contrast structure:")
  print(contrast_selected)

  contrast_name <- paste0("contrast_", contrast_index, "_", contrast_selected$contrast, "_",
                          contrast_selected$specificity_group)
  print(paste("Processing contrast:", contrast_name))

  # Get DESeq2 results for the specific contrast
  deseq_result <- get_contrast_res(dds_base = dds, contrast_row = contrast_selected)

  print("DESeq2 results obtained.")
  print("Summary:")
  print(summary(deseq_result))
  print("Head:")
  print(head(deseq_result))

  expression_data_df <- expression_data_df %>%
    dplyr::mutate_at("GeneId", toupper) %>%
    dplyr::distinct(GeneId, .keep_all = TRUE)

  print("Expression data cleaned:")
  print(head(expression_data_df))

  # Add metadata to results
  annotated_expression_df <- add_metadata_to_results(
    expression_data_df = expression_data_df,
    deseq_result       = deseq_result,
    read_count_cols    = read_counts_columns,
    contrast_name = paste0(contrast_selected$contrast, "_", contrast_selected$specificity_group)
  )

  # Remove " Rpkm" from all column names
  colnames(annotated_expression_df) <- gsub(" Rpkm", "", colnames(annotated_expression_df))

  print("Metadata added to results:")
  print(head(annotated_expression_df))

  # Export DESeq2 report
  export_deseq_report(annotated_expression_df, contrast_name)
  print("Exported DESeq2 report.")

  # Merge into combined annotation
  # If it's the first iteration, just assign directly
  if (is.null(annotated_expression_combined_df) || nrow(annotated_expression_combined_df) == 0) {
    annotated_expression_combined_df <- annotated_expression_df
  } else {
    annotated_expression_combined_df <- annotated_expression_combined_df %>%
      full_join(annotated_expression_df)
  }
}

export_mds_html_plot(normCounts, "mds_plot.html")
print("Exported MDS plot.")

annotated_expression_combined_df <- annotated_expression_combined_df[, !duplicated(names
                                                                                   (annotated_expression_combined_df))]

# Specify the desired order of some key columns
annotation_cols <- c("GeneId", "RefseqId", "Chrom", "TxStart", "TxEnd", "Strand", "baseMean")

# Get the current sample order from the DESeq object
sample_order <- colnames(dds)

# Extract contrast columns: any column name that starts with 'contrast_'
contrast_cols <- grep("vs", colnames(annotated_expression_combined_df), value = TRUE)

print("Annotation contrast columns:")
print(contrast_cols)

# Now create the final order:
# 1. Annotation columns (in the specified order, if they exist)
# 2. The sample columns (in the exact order they appear in the DESeq object)
# 3. All contrast-related columns
final_cols <- c(
  intersect(annotation_cols, colnames(annotated_expression_combined_df)),
  intersect(sample_order, colnames(annotated_expression_combined_df)),
  intersect(contrast_cols, colnames(annotated_expression_combined_df))
)

print("Final column order:")
print(final_cols)

# Reorder the columns of the data frame
annotated_expression_combined_df <- annotated_expression_combined_df[, final_cols, drop = FALSE]

# After the loop, ensure that 'row_metadata' is correctly set
# Remove sample columns from row metadata
row_metadata <- annotated_expression_combined_df %>%
  select(-one_of(sample_order)) %>% # remove all sample columns
  distinct(GeneId, .keep_all = TRUE) %>%
  remove_rownames() %>%
  column_to_rownames("GeneId")

# TODO: here normCounts should be the entire data
# Export GCT data
export_gct_data(normCounts, row_metadata, as.data.frame(colData(dds)))

print("DESeq2 analysis complete.")
