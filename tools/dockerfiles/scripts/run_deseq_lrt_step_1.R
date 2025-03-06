#!/usr/bin/env Rscript
options(warn = -1)
options("width" = 400)

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

options(rlang_backtrace_on_error = "full")

mutate <- dplyr::mutate
filter <- dplyr::filter
group_by <- dplyr::group_by
slice <- dplyr::slice
rename <- dplyr::rename
select <- dplyr::select
`%>%` <- magrittr::`%>%`
`%in%` <- base::`%in%`

###### v0.0.5 ######
#
# Changes:
# - 'combatseq' batch correction is applied in Step 1 if specified.
# - 'limmaremovebatcheffect' is noted and passed to Step 2 for application.
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

###### run_deseq_lrt_step_1.R ######

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
      collected_isoforms <- merge(collected_isoforms,
        isoforms,
        by = intersect_by,
        sort = FALSE
      )
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

# Function to apply batch correction using ComBat_seq
apply_combat_seq <- function(count_data, design, column_data) {

  start_time <- proc.time()

  print("Applying ComBat_seq batch correction")
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

get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq2 for multi-factor analysis using LRT (likelihood ratio or chi-squared test)")
  parser$add_argument(
    "-i",
    "--input",
    help = "Grouped by gene / TSS/ isoform expression files, formatted as CSV/TSV",
    type = "character",
    required = "True",
    nargs = "+"
  )
  parser$add_argument(
    "-n",
    "--name",
    help = "Unique names for input files, no special characters, spaces are allowed. Number and order corresponds to --input",
    type = "character",
    required = "True",
    nargs = "+"
  )
  parser$add_argument(
    "-m",
    "--meta",
    help = "Metadata file to describe relation between samples, where first column corresponds to --name, formatted as CSV/TSV",
    type = "character",
    required = "True"
  )
  parser$add_argument(
    "-d",
    "--design",
    help = "Design formula. Should start with ~, like ~condition+celltype+condition:celltype. See DESeq2 manual for details",
    type = "character",
    required = "True"
  )
  parser$add_argument(
    "-r",
    "--reduced",
    help = "Reduced formula to compare against with the term(s) of interest removed. Should start with ~. See DESeq2 manual for details",
    type = "character",
    required = "True"
  )
  parser$add_argument(
    "--batchcorrection",
    help = paste(
      "Specifies the batch correction method to be applied.",
      "- 'combatseq' applies ComBat_seq at the beginning of the analysis, removing batch effects from the counts before differential expression analysis.",
      "- 'limmaremovebatcheffect' applies removeBatchEffect from the limma package after differential expression analysis.",
      "- Default: none"
    ),
    type = "character",
    choices = c("none", "combatseq", "limmaremovebatcheffect"),
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
      "abscosangle",
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
      "abscosangle",
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
  args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
  return(args)
}

generate_lrt_md <- function(deseq_results, full_formula, reduced_formula, output_file, alpha = args$fdr,
                            batch_warning = NULL) {
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

  # Write the content to the output file
  writeLines(md_content, con = output_file)
}

# Function to generate main effect contrasts with different reference levels
generate_main_effect_contrasts <- function(dds, factors, factor_levels) {

  start_time <- proc.time()
  contrasts  <- list()

  # Helper function to rebuild a DESeqDataSet from scratch
  rebuild_dds <- function(dds, factor_ref_list) {
    # 1) Convert colData to a plain data.frame so we can manipulate it
    colData_df <- as.data.frame(colData(dds))

    # 2) Re-level factors as specified in factor_ref_list, e.g.
    #    factor_ref_list might be list(grouped_ref = c("grouped" = "0H_GFP_N"))
    for (f_name in names(factor_ref_list)) {
      # factor_ref_list[[f_name]] is the new reference level for factor f_name
      ref_level <- factor_ref_list[[f_name]]
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

  # 1) SINGLE-FACTOR scenario
  #    (e.g. ~ grouped with multiple levels)
  if (length(factors) == 1) {
    factor <- factors[[1]]
    all_levels <- factor_levels[[factor]]

    # For each level, define it as reference, compare to others
    for (ref_level in all_levels) {

      # Build a new copy where 'ref_level' is the reference
      dds_temp <- rebuild_dds(
        dds,
        factor_ref_list = setNames(list(ref_level), factor)
      )

      # Run DESeq once for that new reference
      dds_temp <- DESeq(dds_temp, test = "Wald")

      # Compare each *other* level to ref_level
      for (lvl in all_levels) {
        if (lvl != ref_level) {
          contrast_name <- paste(factor, lvl, "vs", ref_level, sep = "_")

          # results() using name = "factor_lvl_vs_ref" style
          contrast_res <- results(
            dds_temp,
            name                = contrast_name,
            alpha               = args$fdr,
            lfcThreshold        = lfcthreshold,
            independentFiltering = TRUE
          )

          significant_genes <- nrow(subset(contrast_res, padj < args$fdr & abs(log2FoldChange) > args$lfcthreshold))

          contrasts <- append(contrasts, list(list(
            effect_type       = "main",
            specificity_group = factor,
            numerator         = lvl,
            denominator       = ref_level,
            contrast          = contrast_name,
            # subset            = dds_temp,
            contrast_res      = contrast_res,
            significant_genes = significant_genes
          )))
        }
      }
    }

    # 2) MULTI-FACTOR scenario
  } else {
    for (factor in factors) {
      other_factors <- setdiff(factors, factor)

      for (other_factor in other_factors) {
        other_levels <- factor_levels[[other_factor]]

        for (other_level in other_levels) {

          for (ref_level in factor_levels[[factor]]) {

            # Build a new copy where factor=ref_level and other_factor=other_level
            dds_temp <- rebuild_dds(
              dds,
              factor_ref_list = c(
                setNames(list(ref_level), factor),
                setNames(list(other_level), other_factor)
              )
            )

            # Run DESeq with these new references
            dds_temp <- DESeq(dds_temp, test = "Wald")

            # For each non-ref level of 'factor'
            for (lvl in factor_levels[[factor]]) {
              if (lvl != ref_level) {
                specificity_group <- paste(other_factor, other_level, sep = "_")
                contrast_name     <- paste(factor, lvl, "vs", ref_level, sep = "_")

                contrast_res <- results(
                  dds_temp,
                  name                = contrast_name,
                  alpha               = args$fdr,
                  lfcThreshold        = lfcthreshold,
                  independentFiltering = TRUE
                )

                significant_genes <- nrow(subset(contrast_res, padj < args$fdr & abs(log2FoldChange) > args$lfcthreshold))

                contrasts <- append(contrasts, list(list(
                  effect_type       = "main",
                  specificity_group = specificity_group,
                  numerator         = lvl,
                  denominator       = ref_level,
                  contrast          = contrast_name,
                  # subset            = dds_temp,
                  contrast_res      = contrast_res,
                  significant_genes = significant_genes
                )))
              }
            }
          }
        }
      }
    }
  }

  end_time <- proc.time() - start_time
  message("Total time of generate_main_effect_contrasts function execution: ")
  print(end_time)

  return(contrasts)
}

# Function to dynamically extract the factor and level names from interaction terms
extract_factors_and_levels <- function(term) {
  parts <- strsplit(term, "\\.")[[1]]
  factor1 <- sub("[0-9A-Z_]+$", "", parts[1])
  factor2 <- sub("[0-9A-Z_]+$", "", parts[2])
  level1 <- sub("^.*?([0-9A-Z_]+$)", "\\1", parts[1])
  level2 <- sub("^.*?([0-9A-Z_]+$)", "\\1", parts[2])
  list(factor1 = factor1, level1 = level1, factor2 = factor2, level2 = level2)
}

# Function to generate interaction effect contrasts
generate_interaction_effect_contrasts <- function(dds) {

  start_time <- proc.time()

  contrasts <- list()
  interaction_names <- resultsNames(dds)

  # Identify interaction terms in the results names
  interaction_terms <- grep("\\.", interaction_names, value = TRUE)

  for (interaction in interaction_terms) {
    factors_levels <- extract_factors_and_levels(interaction)
    factor1 <- factors_levels$factor1
    level1 <- factors_levels$level1
    factor2 <- factors_levels$factor2
    level2 <- factors_levels$level2

    levels1 <- levels(colData(dds)[[factor1]])
    levels2 <- levels(colData(dds)[[factor2]])

    # Generate contrasts for factor1 vs factor2
    for (ref_level1 in levels1) {
      for (ref_level2 in levels2) {
        if ((ref_level1 != level1 || ref_level2 != level2) && (ref_level2 != level2)) {
          specificity_group <- paste(factor2, level2, "vs", ref_level2, sep = "_")
          numerator <- paste0(factor1, level1)
          denominator <- paste0(factor1, ref_level1)

          if (numerator != denominator && specificity_group != paste(factor2, ref_level2, "vs", ref_level2, sep = "_")) {
            # dds_subset <- dds[colData(dds)[[factor2]] == ref_level2, ]
            dds_subset <- dds
            colData(dds_subset)[[factor1]] <- relevel(colData(dds_subset)[[factor1]], ref = ref_level1)

            dds_subset <- DESeq(dds_subset, test = "Wald")

            contrast_res <- results(dds_subset,
              name = interaction,
              alpha = args$fdr,
              lfcThreshold = lfcthreshold,
              independentFiltering = TRUE
            )

            significant_genes <- nrow(subset(contrast_res, padj < args$fdr & abs(log2FoldChange) > args$lfcthreshold))

            contrasts <- append(contrasts, list(list(
              effect_type = "interaction",
              specificity_group = specificity_group,
              numerator = numerator,
              denominator = denominator,
              contrast = interaction,
              # subset = dds_subset,
              contrast_res = contrast_res,
              significant_genes = significant_genes
            )))
          }
        }
      }
    }

    # Generate contrasts for factor2 vs factor1
    for (ref_level2 in levels2) {
      for (ref_level1 in levels1) {
        if ((ref_level2 != level2 || ref_level1 != level1) && (ref_level1 != level1)) {
          specificity_group <- paste(factor1, level1, "vs", ref_level1, sep = "_")
          numerator <- paste0(factor2, level2)
          denominator <- paste0(factor2, ref_level2)

          if (numerator != denominator && specificity_group != paste(factor1, ref_level1, "vs", ref_level1, sep = "_")) {
            # dds_subset <- dds[colData(dds)[[factor1]] == ref_level1, ]
            dds_subset <- dds
            colData(dds_subset)[[factor2]] <- relevel(colData(dds_subset)[[factor2]], ref = ref_level2)

            dds_subset <- DESeq(dds_subset, test = "Wald")

            contrast_res <- results(dds_subset,
              name = interaction,
              alpha = args$fdr,
              lfcThreshold = lfcthreshold,
              independentFiltering = TRUE
            )

            significant_genes <- nrow(subset(contrast_res, padj < args$fdr & abs(log2FoldChange) > args$lfcthreshold))

            contrasts <- append(contrasts, list(list(
              effect_type = "interaction",
              specificity_group = specificity_group,
              numerator = numerator,
              denominator = denominator,
              contrast = interaction,
              # subset = dds_subset,
              contrast_res = contrast_res,
              significant_genes = significant_genes
            )))
          }
        }
      }
    }
  }

  end_time <- proc.time() - start_time

  print("Total time of generate_interaction_effect_contrasts function execution: ")
  print(end_time)

  return(contrasts)
}

# Main function to generate the dataframe with all possible contrasts
generate_contrasts <- function(dds) {

  start_time <- proc.time()

  # Extract the design formula and model matrix
  design_formula <- design(dds)

  # Get the levels of each factor in the design
  factors <- all.vars(design_formula)
  factors <- factors[!str_detect(factors, "batch")]
  print("Analyzing the following factors:")
  print(factors)

  factor_levels <- lapply(factors, function(f) levels(colData(dds)[[f]]))
  names(factor_levels) <- factors
  print("Analyzing the following factor levels:")
  print(factor_levels)

  # Generate main and interaction effect contrasts
  main_contrasts <- generate_main_effect_contrasts(dds, factors, factor_levels)
  interaction_contrasts <- generate_interaction_effect_contrasts(dds)
  all_contrasts <- c(main_contrasts, interaction_contrasts)

  print("Head of all_contrasts:")
  print(head(all_contrasts))

  all_contrasts <- purrr::list_modify(all_contrasts, expression_data_df = expression_data_df)
  all_contrasts <- purrr::list_modify(all_contrasts, deseq_obj = dds)

  print(paste("Exporting contrasts list to", paste0(args$output, "_contrasts.rds"), sep = " "))
  saveRDS(all_contrasts, file = paste0(args$output, "_contrasts.rds"))

  # Create a dataframe to store the contrasts and number of significant genes
  contrast_df <- data.frame(
    effect = character(),
    specificity_group = character(),
    contrast = character(),
    numerator = character(),
    denominator = character(),
    significant_genes = integer(),
    stringsAsFactors = FALSE
  )

  # Run DESeq2 for each contrast and store the results
  count <- 1
  for (contrast in all_contrasts) {
    print(count)

    contrast_df <- rbind(contrast_df, data.frame(
      effect = contrast$effect_type,
      contrast = contrast$contrast,
      specificity_group = contrast$specificity_group,
      numerator = contrast$numerator,
      denominator = contrast$denominator,
      significant_genes = contrast$significant_genes,
      stringsAsFactors = FALSE
    ))
    count <- count + 1
  }

  # Remove duplicate contrasts
  contrast_df <- contrast_df %>%
    group_by(specificity_group) %>%
    distinct(contrast, .keep_all = TRUE) %>%
    arrange(specificity_group) %>%
    ungroup() %>%
    ungroup() %>%
    mutate(contrast_number = row_number()) %>%
    select(contrast_number, everything()) %>%
    arrange(desc(significant_genes))

  end_time <- proc.time() - start_time

  print("Total time of generate_contrasts function execution: ")
  print(end_time)

  # Sort the dataframe by the number of significant genes
  return(contrast_df)
}

# Function to clean sample names
clean_sample_names <- function(names) {
  names <- trimws(names) # Remove leading/trailing whitespace
  names <- tolower(names) # Convert to lowercase
  names <- gsub("\\s+", "_", names) # Replace spaces with underscores
  names <- gsub("[^[:alnum:]_]", "", names) # Remove non-alphanumeric characters except underscores
  names <- gsub("^_+|_+$", "", names) # Remove leading/trailing underscores
  return(names)
}

# Function to export MDS plot
export_mds_html_plot <- function(norm_counts_data, location) {
  tryCatch(
    expr = {
      htmlwidgets::saveWidget(
        Glimma::glimmaMDS(
          x      = norm_counts_data,
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

# Function to export normalized counts and filtered counts to GCT format
# TODO: this function should also return levels of HOPACH clustering in the row annotation
# Updated export_gct_data function:
# 1. Write full GCT file without clustering.
# 2. Filter rows by padj.
# 3. Cluster only the filtered subset.
# 4. Write filtered (and clustered) GCT file.

export_gct_data <- function(normCounts, row_metadata, col_metadata, output_prefix) {
  tryCatch({
    # Ensure col_metadata columns are vectors
    col_metadata <- col_metadata %>% mutate_all(as.vector)

    # Initialize col_order_vector (will be non-NULL if user specifies a column order)
    col_order_vector <- tolower(rownames(metadata_df))
    print("Order of columns for heatmap without clustering, based on metadata input:")
    print(col_order_vector)

    normCounts   <- normCounts[, col_order_vector, drop = FALSE]
    col_metadata <- col_metadata[col_order_vector, , drop = FALSE]

    # Create and export the initial GCT file
    gct_data <- new("GCT", mat = normCounts, rdesc = row_metadata, cdesc = col_metadata)
    cmapR::write_gct(ds = gct_data, ofile = paste0(output_prefix, "_counts_all.gct"), appenddim = FALSE)
    print(paste("Exporting GCT data to", paste0(output_prefix, "_counts_all.gct")))

    # Filter by padj if available
    if ('padj' %in% colnames(row_metadata)) {
      row_metadata_filtered <- row_metadata %>% filter(padj <= args$fdr)
      print("=== After Filtering by padj ===")
      print(paste("Number of rows in row_metadata_filtered:", nrow(row_metadata_filtered)))
      print("Sample row names in row_metadata_filtered:")
      print(head(rownames(row_metadata_filtered)))
    } else {
      row_metadata_filtered <- row_metadata
      print("=== No Filtering by padj ===")
      print(paste("Number of rows in row_metadata_filtered:", nrow(row_metadata_filtered)))
    }

    # Convert to uppercase and trim whitespace
    rownames(row_metadata_filtered) <- toupper(trimws(rownames(row_metadata_filtered)))
    rownames(normCounts)            <- toupper(trimws(rownames(normCounts)))

    # Check for duplicates in normCounts
    if (any(duplicated(rownames(normCounts)))) {
      print("Duplicate gene IDs found in normCounts. Removing duplicates.")
      normCounts <- normCounts[!duplicated(rownames(normCounts)),]
    }

    # Check for duplicates in row_metadata_filtered
    if (any(duplicated(rownames(row_metadata_filtered)))) {
      print("Duplicate gene IDs found in row_metadata_filtered. Removing duplicates.")
      row_metadata_filtered <- row_metadata_filtered[!duplicated(rownames(row_metadata_filtered)),]
    }

    # Check if any rows pass the filter
    if (nrow(row_metadata_filtered) == 0) {
      warning(paste("No genes passed the FDR threshold of", args$fdr, "or the log2 fold change threshold of", lfcthreshold))
      # Optionally, skip exporting filtered GCT
      return(NULL)
    }

    # Debugging before subsetting
    print("=== Debugging before subsetting ===")
    print("Checking if all row names in row_metadata_filtered are present in normCounts:")
    all_present <- all(rownames(row_metadata_filtered) %in% rownames(normCounts))
    print(paste("All row names present:", all_present))

    if (!all_present) {
      missing_rows <- setdiff(rownames(row_metadata_filtered), rownames(normCounts))
      print(paste("Number of missing row names:", length(missing_rows)))
      print("Examples of missing row names:")
      print(head(missing_rows))
    }

    # Subset only the matching row names to prevent "subscript out of bounds" errors
    common_rows           <- intersect(rownames(row_metadata_filtered), rownames(normCounts))
    filtered_normCounts   <- normCounts[common_rows, , drop = FALSE]
    row_metadata_filtered <- row_metadata_filtered[common_rows, , drop = FALSE]
    print(paste("Number of rows after subsetting:", nrow(filtered_normCounts)))

    # Existing debug statements
    print("=== Debugging inside export_gct_data ===")
    print("Dimensions of normCounts:")
    print(dim(normCounts))
    print("Dimensions of row_metadata:")
    print(dim(row_metadata))
    print("Columns of row_metadata:")
    print(colnames(row_metadata))
    print("A few row names of row_metadata:")
    print(head(rownames(row_metadata)))

    if ("padj" %in% colnames(row_metadata)) {
      print("padj column found. Summary of padj:")
      print(summary(row_metadata$padj))
    } else {
      print("padj column not found in row_metadata.")
    }

    # Now cluster the filtered data
    clustered_data <- cluster_and_reorder(filtered_normCounts, col_metadata, row_metadata_filtered, args)

    clustered_data$normCounts   <- clustered_data$normCounts[, col_order_vector, drop = FALSE]
    clustered_data$col_metadata <- clustered_data$col_metadata[col_order_vector, , drop = FALSE]

    # After clustering:
    print("Dimensions of row_metadata_filtered:")
    print(dim(row_metadata_filtered))
    print("A few row names of row_metadata_filtered:")
    print(head(rownames(row_metadata_filtered)))

    # Check intersection with normCounts rows:
    matching_rows <- intersect(rownames(normCounts), rownames(row_metadata_filtered))
    print("Number of matching row names between normCounts and row_metadata_filtered:")
    print(length(matching_rows))

    non_matching_rows <- setdiff(rownames(row_metadata_filtered), rownames(normCounts))
    print("Number of row names in row_metadata_filtered not found in normCounts:")
    print(length(non_matching_rows))
    if (length(non_matching_rows) > 0) {
      print("Examples of non-matching row names:")
      print(head(non_matching_rows))
    }

    # Create and export the filtered GCT file
    gct_data_filtered <- new("GCT",
                             mat   = clustered_data$normCounts,
                             rdesc = clustered_data$row_metadata,
                             cdesc = clustered_data$col_metadata)

    cmapR::write_gct(ds = gct_data_filtered, ofile = paste0(output_prefix, "_counts_filtered.gct"), appenddim = FALSE)
    print(paste("Exporting GCT data to", paste0(output_prefix, "_counts_filtered.gct")))

  }, error = function(e) {
    print(paste("Failed to export GCT data to", output_prefix, "with error -", e$message))
  })
}

# Removed clustering from here. Just export MDS and call export_gct_data().
export_charts <- function(res, annotated_expression_df, column_data, normCounts, output, args) {
  # Export MDS plot
  print("Exporting MDS plot")
  # In case we have original counts, we need to use them for MDS plot
  # It means that we have some sort of batch-correction function applied
  if (exists("original_counts")) {

    dse <- DESeqDataSetFromMatrix(
      countData = original_counts,
      colData   = metadata_df,
      design    = design_formula
    )

    rlog_original_counts <- assay(rlog(dse, blind = FALSE))

    print("Exporting MDS plot with original counts")

    export_mds_html_plot(rlog_original_counts, paste0(output, "_mds_plot.html"))

    print("Exporting MDS plot with batch correction")

    # normCounts is supposed to be batch-corrected here and rlog-normalised
    export_mds_html_plot(normCounts, paste0(output, "_mds_plot_corrected.html"))

  } else {
    export_mds_html_plot(normCounts, paste0(output, "_mds_plot.html"))
  }

  # Now just call export_gct_data directly without clustering here.
  export_gct_data(normCounts, annotated_expression_df, column_data, output)
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
      clustered_data_cols <- get_clustered_data(normCounts, transpose = TRUE, k = k, kmax = kmax, scaling_type = args$scaling_type, dist = args$columndist)
      normCounts          <- normCounts[, clustered_data_cols$order, drop = FALSE]
      col_metadata        <- col_metadata[clustered_data_cols$order, , drop = FALSE]
      # After reordering, cbind cluster info
      col_metadata        <- cbind(col_metadata, clustered_data_cols$clusters)
    }
    # Row clustering if requested
    if (args$cluster == "row" || args$cluster == "both") {
      clustered_data_rows <- get_clustered_data(normCounts, transpose = FALSE, k = k, kmax = kmax, scaling_type = args$scaling_type, dist = args$rowdist)
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


# Function to clean sample names
clean_sample_names <- function(names) {
  names <- trimws(names) # Remove leading/trailing whitespace
  names <- tolower(names) # Convert to lowercase
  names <- gsub("\\s+", "_", names) # Replace spaces with underscores
  names <- gsub("[^[:alnum:]_]", "", names) # Remove non-alphanumeric characters except underscores
  names <- gsub("^_+|_+$", "", names) # Remove leading/trailing underscores
  return(names)
}

# Function to export MDS plot
export_mds_html_plot <- function(norm_counts_data, location) {
  tryCatch(
    expr = {
      htmlwidgets::saveWidget(
        Glimma::glimmaMDS(
          x      = norm_counts_data,
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

# Function to export normalized counts and filtered counts to GCT format
# TODO: this function should also return levels of HOPACH clustering in the row annotation
# Updated export_gct_data function:
# 1. Write full GCT file without clustering.
# 2. Filter rows by padj.
# 3. Cluster only the filtered subset.
# 4. Write filtered (and clustered) GCT file.

export_gct_data <- function(normCounts, row_metadata, col_metadata, output_prefix) {
  tryCatch({
    # Ensure col_metadata columns are vectors
    col_metadata <- col_metadata %>% mutate_all(as.vector)

    # Create and export the initial GCT file
    gct_data <- new("GCT", mat = normCounts, rdesc = row_metadata, cdesc = col_metadata)
    cmapR::write_gct(ds = gct_data, ofile = paste0(output_prefix, "_counts_all.gct"), appenddim = FALSE)
    print(paste("Exporting GCT data to", paste0(output_prefix, "_counts_all.gct")))

    # Filter by padj if available
    if ('padj' %in% colnames(row_metadata)) {
      row_metadata_filtered <- row_metadata %>% filter(padj <= args$fdr)
      print("=== After Filtering by padj ===")
      print(paste("Number of rows in row_metadata_filtered:", nrow(row_metadata_filtered)))
      print("Sample row names in row_metadata_filtered:")
      print(head(rownames(row_metadata_filtered)))
    } else {
      row_metadata_filtered <- row_metadata
      print("=== No Filtering by padj ===")
      print(paste("Number of rows in row_metadata_filtered:", nrow(row_metadata_filtered)))
    }

    # Convert to uppercase and trim whitespace
    rownames(row_metadata_filtered) <- toupper(trimws(rownames(row_metadata_filtered)))
    rownames(normCounts)            <- toupper(trimws(rownames(normCounts)))

    # Check for duplicates in normCounts
    if (any(duplicated(rownames(normCounts)))) {
      print("Duplicate gene IDs found in normCounts. Removing duplicates.")
      normCounts <- normCounts[!duplicated(rownames(normCounts)),]
    }

    # Check for duplicates in row_metadata_filtered
    if (any(duplicated(rownames(row_metadata_filtered)))) {
      print("Duplicate gene IDs found in row_metadata_filtered. Removing duplicates.")
      row_metadata_filtered <- row_metadata_filtered[!duplicated(rownames(row_metadata_filtered)),]
    }

    # Check if any rows pass the filter
    if (nrow(row_metadata_filtered) == 0) {
      warning(paste("No genes passed the FDR threshold of", args$fdr))
      # Optionally, skip exporting filtered GCT
      return(NULL)
    }

    # **New Debug Statements Before Subsetting**
    print("=== Debugging before subsetting ===")
    print("Checking if all row names in row_metadata_filtered are present in normCounts:")
    all_present <- all(rownames(row_metadata_filtered) %in% rownames(normCounts))
    print(paste("All row names present:", all_present))

    if (!all_present) {
      missing_rows <- setdiff(rownames(row_metadata_filtered), rownames(normCounts))
      print(paste("Number of missing row names:", length(missing_rows)))
      print("Examples of missing row names:")
      print(head(missing_rows))
    }

    # Proceed to subset only the matching row names
    # This prevents the "subscript out of bounds" error
    common_rows           <- intersect(rownames(row_metadata_filtered), rownames(normCounts))
    filtered_normCounts <- normCounts[common_rows, , drop = FALSE]
    row_metadata_filtered <- row_metadata_filtered[common_rows, , drop = FALSE]  # This line ensures consistency
    print(paste("Number of rows after subsetting:", nrow(filtered_normCounts)))

    # **Existing Debug Statements Inside export_gct_data**
    print("=== Debugging inside export_gct_data ===")
    print("Dimensions of normCounts:")
    print(dim(normCounts))
    print("Dimensions of row_metadata:")
    print(dim(row_metadata))
    print("Columns of row_metadata:")
    print(colnames(row_metadata))
    print("A few row names of row_metadata:")
    print(head(rownames(row_metadata)))

    if ("padj" %in% colnames(row_metadata)) {
      print("padj column found. Summary of padj:")
      print(summary(row_metadata$padj))
    } else {
      print("padj column not found in row_metadata.")
    }

    # Now cluster the filtered data
    clustered_data <- cluster_and_reorder(filtered_normCounts, col_metadata, row_metadata_filtered, args)

    # After clustering:
    print("Dimensions of row_metadata_filtered:")
    print(dim(row_metadata_filtered))
    print("A few row names of row_metadata_filtered:")
    print(head(rownames(row_metadata_filtered)))

    # Check intersection with normCounts rows:
    matching_rows <- intersect(rownames(normCounts), rownames(row_metadata_filtered))
    print("Number of matching row names between normCounts and row_metadata_filtered:")
    print(length(matching_rows))

    non_matching_rows <- setdiff(rownames(row_metadata_filtered), rownames(normCounts))
    print("Number of row names in row_metadata_filtered not found in normCounts:")
    print(length(non_matching_rows))
    if (length(non_matching_rows) > 0) {
      print("Examples of non-matching row names:")
      print(head(non_matching_rows))
    }

    # Create and export the filtered GCT file
    gct_data_filtered <- new("GCT",
                             mat   = clustered_data$normCounts,
                             rdesc = clustered_data$row_metadata,
                             cdesc = clustered_data$col_metadata)

    cmapR::write_gct(ds = gct_data_filtered, ofile = paste0(output_prefix, "_counts_filtered.gct"), appenddim = FALSE)
    print(paste("Exporting GCT data to", paste0(output_prefix, "_counts_filtered.gct")))

  }, error = function(e) {
    print(paste("Failed to export GCT data to", output_prefix, "with error -", e$message))
  })
}

# Removed clustering from here. Just export MDS and call export_gct_data().
export_charts <- function(res, annotated_expression_df, column_data, normCounts, output, args) {
  # Export MDS plot
  print("Exporting MDS plot")
  # In case we have original counts, we need to use them for MDS plot
  # It means that we have some sort of batch-correction function applied
  if (exists("original_counts")) {

    dse <- DESeqDataSetFromMatrix(
      countData = original_counts,
      colData   = metadata_df,
      design    = design_formula
    )

    rlog_original_counts <- assay(rlog(dse, blind = FALSE))

    print("Exporting MDS plot with original counts")

    export_mds_html_plot(rlog_original_counts, paste0(output, "_mds_plot.html"))

    print("Exporting MDS plot with batch correction")

    # normCounts is supposed to be batch-corrected here and rlog-normalised
    export_mds_html_plot(normCounts, paste0(output, "_mds_plot_corrected.html"))

  } else {
    export_mds_html_plot(normCounts, paste0(output, "_mds_plot.html"))
  }

  # Now just call export_gct_data directly without clustering here.
  export_gct_data(normCounts, annotated_expression_df, column_data, output)
}

# Function to generate clusters
get_clustered_data <- function(expression_data, transpose = FALSE, k = 3, kmax = 5) {

  start_time <- proc.time()

  if (transpose) {
    print("Transposing expression data")
    expression_data <- t(expression_data)
  }

  # Apply scaling per row
  expression_data <- t(apply(expression_data, 1, scale_min_max))

  if (transpose) {
    print("Transposing expression data back")
    expression_data <- t(expression_data)
  }

  print(paste0("Running HOPACH for ", nrow(expression_data), "  features"))
  hopach_results <- hopach::hopach(expression_data,
                                   verbose = TRUE,
                                   K       = k,
                                   kmax    = kmax,
                                   khigh   = kmax
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
    order = hopach_results$clustering$order,
    expression = expression_data,
    clusters = clusters
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

cluster_and_reorder <- function(normCounts, col_metadata, row_metadata, args) {

  start_time <- proc.time()

  if (args$cluster != "none") {
    # Column clustering if requested
    if (args$cluster == "column" || args$cluster == "both") {
      clustered_data_cols <- get_clustered_data(normCounts, transpose = TRUE)
      # TODO: should we use drop T or F here? (it fails with F, but works with T), I'm not quite sure if is it
      #  correct way
      normCounts <- normCounts[, clustered_data_cols$order]
      col_metadata        <- col_metadata[clustered_data_cols$order, , drop = FALSE]
      # After reordering, cbind cluster info
      col_metadata        <- cbind(col_metadata, clustered_data_cols$clusters)
    }
    # Row clustering if requested
    if (args$cluster == "row" || args$cluster == "both") {
      if (args$test_mode) {
        k    <- 2
        kmax <- 2
      } else {
        k    <- args$k
        kmax <- args$kmax
      }
      clustered_data_rows <- get_clustered_data(normCounts, transpose = FALSE, k = k, kmax = kmax)
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


# Parse arguments
args <- get_args()

# Set threads
register(MulticoreParam(args$threads))

# Load metadata
metadata_df <- read.table(
  args$meta,
  sep = get_file_type(args$meta),
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1
)
print(paste("Load metadata from", args$meta, sep = " "))
print(metadata_df)

# Load design formula
design_formula <- as.formula(tolower(args$design))
print("Load design formula")
print(design_formula)

# Load reduced formula
reduced_formula <- as.formula(tolower(args$reduced))
print("Load reduced formula")
print(reduced_formula)

print("Using use_lfc_thresh argument as:")
print(args$use_lfc_thresh)
print("So, the lfc threshold is:")
lfcthreshold <- ifelse(args$use_lfc_thresh, args$lfcthreshold, 0)
print(lfcthreshold)

# Clean sample names
args$name <- clean_sample_names(args$name)
colnames(metadata_df) <- clean_sample_names(colnames(metadata_df))
rownames(metadata_df) <- clean_sample_names(rownames(metadata_df))

print("Processed args$name:")
print(args$name)

print("Processed rownames(metadata_df):")
print(rownames(metadata_df))

# Load expression data
expression_data_df <- load_expression_data(args$input, args$name, READ_COL, RPKM_COL, INTERSECT_BY)
print("Expression data to analyze: ")
print(head(expression_data_df))
print(dim(expression_data_df))

if (!is.null(args$rpkm_cutoff)) {
  print("Using RPKM cutoff for filtering:")
  print(args$rpkm_cutoff)
  expression_data_df <- filter_rpkm(expression_data_df, args$rpkm_cutoff)
  print("Expression data after RPKM filtering: ")
  print(head(expression_data_df))
  print(dim(expression_data_df))
}

# Select all columns with read counts data, reorder them based on the row names from metadata_df
read_counts_columns <- grep(
  paste(READ_COL, sep = ""),
  colnames(expression_data_df),
  value = TRUE,
  ignore.case = TRUE
)

# TODO: double-check if we should use GeneId or RefseqId here
read_counts_data_df <- expression_data_df %>%
  dplyr::mutate_at("GeneId", toupper) %>%
  dplyr::distinct(GeneId, .keep_all = TRUE) %>% # to prevent from failing when input files are not grouped by GeneId
  remove_rownames() %>%
  column_to_rownames("GeneId")

read_counts_data_df <- read_counts_data_df[read_counts_columns]

# Clean column names of read_counts_data_df
colnames(read_counts_data_df) <- lapply(colnames(read_counts_data_df), function(s) {
  paste(head(unlist(strsplit(s, " ", fixed = TRUE)), -1), collapse = " ")
})
colnames(read_counts_data_df) <- clean_sample_names(colnames(read_counts_data_df))

print("Processed colnames(read_counts_data_df):")
print(colnames(read_counts_data_df))

print("Read counts data")
print(head(read_counts_data_df))

# Check for mismatches in sample names
samples_in_metadata_not_in_data <- setdiff(rownames(metadata_df), colnames(read_counts_data_df))
samples_in_data_not_in_metadata <- setdiff(colnames(read_counts_data_df), rownames(metadata_df))

if (length(samples_in_metadata_not_in_data) > 0) {
  print("Error: The following samples are in the metadata but not in the expression data:")
  print(samples_in_metadata_not_in_data)
  stop("Sample names mismatch between metadata and expression data.")
}

if (length(samples_in_data_not_in_metadata) > 0) {
  print("Error: The following samples are in the expression data but not in the metadata:")
  print(samples_in_data_not_in_metadata)
  stop("Sample names mismatch between expression data and metadata.")
}

corrected_counts <- NULL

# Reorder read counts data based on metadata
read_counts_data_df <- read_counts_data_df[, rownames(metadata_df)]

# Initialize batch warning message
batch_warning <- NULL

countData <- read_counts_data_df

# Check for batch correction and handle accordingly
if (args$batchcorrection != "none") {
  if ("batch" %in% colnames(metadata_df)) {
    # Ensure 'batch' is a factor
    metadata_df$batch <- as.factor(metadata_df$batch)

    if (args$batchcorrection == "combatseq") {
      # Case 2: Apply ComBat_seq batch correction to raw counts
      print("Applying ComBat_seq batch correction")
      original_counts <- countData
      corrected_counts <- apply_combat_seq(countData, design_formula, metadata_df)
      countData  <- corrected_counts
      # Design formula remains as provided
    } else if (args$batchcorrection == "limmaremovebatcheffect") {
      # Case 3: Include 'batch' in the design formula
      print("Including 'batch' in the design formula for limma batch correction")
      # Remove '~' from the original design formula string
      original_design <- substring(tolower(args$design), 2)
      # Create new design formula with 'batch' included
      design_formula <- as.formula(paste("~ batch +", original_design))
      # Note: We do not modify counts at this stage
    }
  } else {
    batch_warning <- "You provided batch correction but there's no 'batch' column in your metadata."
    print(batch_warning)
    # Proceed without batch correction
    # countData remains as read_counts_data_df
  }
} else {
  # Case 1: No batch correction
  # countData remains as read_counts_data_df
  # Design formula remains as provided
  print("No batch correction provided; proceeding with default settings")
}

print("Starting DESeq object creation using:")
print("Coldata:")
print(metadata_df)
print("Countdata:")
print(head(countData))
print("Design formula:")
print(design_formula)

# Create DESeq2 dataset
dse <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = metadata_df,
  design = design_formula
)

# Extract the raw counts matrix
cts <- counts(dse) # a matrix of counts, rows = genes, columns = samples

# Sum counts per sample (sum columns to get total reads per sample)
sample_sums <- colSums(cts)

# Create a data frame for plotting
plot_df <- data.frame(
  Sample = names(sample_sums),
  Count  = sample_sums
)

# Now 'plot_df' has three columns: Sample, Count, and Statistic
# You can plot them using ggplot2, similar to the previous logic

stat_barchart <- ggplot(plot_df, aes(x = Sample, y = Count)) +
  geom_col(position = "dodge", fill = "royalblue") +
  labs(
    title = "Total Reads by Sample",
    x     = "Sample",
    y     = "Count"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("alignment_stats_barchart.png",
       plot   = stat_barchart,
       dpi    = 400,
       height = 3,
       width  = max(10, length(sample_sums) * 0.5),
       units  = "in"
)

print("Run DESeq2 using Wald")
dsq_wald <- DESeq(
  dse,
  test = "Wald",
  quiet = FALSE,
  parallel = TRUE
)

# Obtain normalized counts for downstream analyses
if (args$batchcorrection == "limmaremovebatcheffect" && "batch" %in% colnames(metadata_df)) {
  # Case 3: After DESeq2, apply rlog transformation and remove batch effects using limma
  print("Applying rlog transformation and limma batch effect removal")

  original_counts <- countData

  rlog_transformed <- rlog(dsq_wald, blind = FALSE)
  rlog_counts <- assay(rlog_transformed)

  # Prepare design matrix without 'batch' for removeBatchEffect
  design_formula <- as.formula(paste0("~", str_remove(as.character(dsq_wald@design)[2], " \\+ batch")))

  design_matrix <- model.matrix(design_formula, data = metadata_df)
  # Apply removeBatchEffect
  corrected_counts <- limma::removeBatchEffect(rlog_counts, batch = metadata_df$batch, design = design_matrix)
  normCounts <- corrected_counts
} else {
  # Case 1 and Case 2: Use rlog-transformed counts without additional batch correction
  print("Applying rlog transformation without additional batch effect removal")
  rlog_transformed <- rlog(dsq_wald, blind = FALSE)
  normCounts <- assay(rlog_transformed)
}

# Obtain normalized counts for downstream analyses
if (args$batchcorrection == "limmaremovebatcheffect" && "batch" %in% colnames(metadata_df)) {
  # Case 3: After DESeq2, apply rlog transformation and remove batch effects using limma
  print("Applying rlog transformation and limma batch effect removal")

  original_counts <- countData

  rlog_transformed <- rlog(dsq_wald, blind = FALSE)
  rlog_counts <- assay(rlog_transformed)

  # Prepare design matrix without 'batch' for removeBatchEffect
  design_formula <- as.formula(paste0("~", str_remove(as.character(dsq_wald@design)[2], " \\+ batch")))

  design_matrix <- model.matrix(design_formula, data = metadata_df)
  # Apply removeBatchEffect
  corrected_counts <- limma::removeBatchEffect(rlog_counts, batch = metadata_df$batch, design = design_matrix)
  normCounts <- corrected_counts
} else {
  # Case 1 and Case 2: Use rlog-transformed counts without additional batch correction
  print("Applying rlog transformation without additional batch effect removal")
  rlog_transformed <- rlog(dsq_wald, blind = FALSE)
  normCounts <- assay(rlog_transformed)
}

print("Run DESeq2 using LRT")
dsq_lrt <- DESeq(
  dse,
  test = "LRT",
  reduced = reduced_formula,
  quiet = FALSE,
  parallel = TRUE
)

if (!args$lrt_only_mode) {
  print("Generate contrasts")

  contrast_df <- generate_contrasts(dsq_wald)

  contrasts_filename <- paste(args$output, "_contrasts_table.tsv", sep = "")
  write.table(
    contrast_df,
    file      = contrasts_filename,
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )
}


dsq_lrt_res <- results(
  dsq_lrt,
  alpha = args$fdr,
  independentFiltering = TRUE
)

print("LRT Results description")
print(mcols(dsq_lrt_res))

# Export results to TSV file
annotated_expression_df <- expression_data_df %>%
  bind_cols(as.data.frame(dsq_lrt_res) %>% select(pvalue, padj)) %>%
  mutate(
    `-LOG10(pval)` = -log10(as.numeric(pvalue)),
    `-LOG10(padj)` = -log10(as.numeric(padj))
  )


lrt_report_filename <- paste0(args$output, "_lrt_result.md")
summary(dsq_lrt_res)
generate_lrt_md(dsq_lrt_res, args$design, args$reduced, lrt_report_filename, alpha = args$fdr, batch_warning = batch_warning)
print(paste("Export LRT markdown report to", lrt_report_filename, sep = " "))

results_filename <- paste0(args$output, "_gene_exp_table.tsv")
write.table(
  annotated_expression_df,
  file = results_filename,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
print(paste("Export results to", results_filename, sep = " "))
graphics.off()

sample_order <- colnames(dse)

colnames(annotated_expression_df) <- gsub(" Rpkm", "", colnames(annotated_expression_df))

annotated_expression_df <- annotated_expression_df %>%
  select(-one_of(sample_order)) # remove all sample columns

colnames(annotated_expression_df) <- gsub(" TotalReads", "", colnames(annotated_expression_df))

annotated_expression_df <- annotated_expression_df %>%
  select(-one_of(sample_order)) %>% # remove all sample columns
  distinct(GeneId, .keep_all = TRUE) %>%
  remove_rownames() %>%
  column_to_rownames("GeneId")

export_charts(dsq_lrt_res, annotated_expression_df, metadata_df, normCounts, args$output, args)
