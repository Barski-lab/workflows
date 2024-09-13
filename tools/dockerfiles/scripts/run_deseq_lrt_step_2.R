#!/usr/bin/env Rscript
options(warn = -1)
options("width" = 300)
options(
  error = function() {
    traceback(3)
    quit(
      save = "no",
      status = 1,
      runLast = FALSE
    )
  }
)


suppressMessages(library(argparse))
suppressMessages(library(BiocParallel))
suppressMessages(library(pheatmap))
suppressMessages(library(Glimma))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(hopach))
suppressMessages(library(cmapR))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(DESeq2))


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
#   - `limmaremovebatcheffect`: Removes batch effects post-analysis using limma.
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
RPKM_UNTREATED_ALIAS <- "RpkmCondition1"
RPKM_TREATED_ALIAS <- "RpkmCondition2"


get_file_type <- function(filename) {
  ext <- tools::file_ext(filename)
  separator <- ","
  if (ext == "tsv") {
    separator <- "\t"
  }
  return(separator)
}


export_mds_html_plot <- function(norm_counts_data, location) {
  tryCatch(
    expr = {
      htmlwidgets::saveWidget(
        glimmaMDS(
          x = assay(norm_counts_data),
          groups = as.data.frame(SummarizedExperiment::colData(norm_counts_data)),
          labels = rownames(SummarizedExperiment::colData(norm_counts_data))
        ),
        file = location
      )
    },
    error = function(e) {
      print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
    }
  )
}


export_gct <- function(counts_mat,
                       row_metadata,
                       col_metadata,
                       location) {
  tryCatch(
    expr = {
      row_metadata <- row_metadata %>%
        rownames_to_column("id") %>%
        mutate_at("id", as.vector)
      col_metadata <- col_metadata %>%
        rownames_to_column("id") %>%
        mutate_at("id", as.vector)
      gct_data <- new("GCT",
                      mat = counts_mat[row_metadata$id, col_metadata$id],
                      # to guarantee the order and number of row/columns
                      rdesc = row_metadata,
                      cdesc = col_metadata
      )
      write_gct(
        ds = gct_data,
        ofile = location,
        appenddim = FALSE
      )
      print(paste("Exporting GCT data to", location, sep = " "))
    },
    error = function(e) {
      print(paste("Failed to export GCT data to", location, sep = " "))
    }
  )
}


export_cls <- function(categories, location) {
  base::tryCatch(
    expr = {
      output_stream <- base::file(location, "w")
      on.exit(base::close(output_stream), add = TRUE) # can't put it in 'finally' as there is no access to output_stream variable
      base::cat(
        base::paste(length(categories), # number of datasets
                    length(base::levels(categories)), # number of different categories
                    "1", # should be always 1
                    sep = "\t"
        ),
        base::paste("#", base::paste(
          base::unique(as.character(categories)), # preserves the order, but removes duplicates
          collapse = "\t"
        ), sep = "\t"),
        base::paste(base::paste(as.character(categories),
                                collapse =
                                  "\t"
        ), sep = "\t"),
        file = output_stream,
        sep = "\n"
      )
      base::print(base::paste("Exporting CLS data to", location,
                              sep =
                                " "
      ))
    },
    error = function(e) {
      base::print(base::paste("Failed to export CLS data to ", location, "with error - ", e,
                              sep =
                                ""
      ))
    }
  )
}

get_clustered_data <- function(expression_data, transpose) {
  if (transpose) {
    print("Transposing expression data")
    expression_data <- t(expression_data)
  }

  expression_data <- apply(
    expression_data,
    1,
    FUN = function(x) {
      scale_min_max(x)
    }
  )

  print("Running HOPACH")
  hopach_results <- hopach(expression_data)

  if (transpose) {
    print("Transposing expression data")
    expression_data <- t(expression_data)
  }

  print("Parsing cluster labels")
  clusters <- as.data.frame(hopach_results$clustering$labels)
  colnames(clusters) <- "label"
  clusters <- cbind(clusters, "HCL" = outer(clusters$label, 10^c((nchar(
    trunc(clusters$label)
  )[1] - 1):0), function(a, b) {
    paste0("c", a %/% b %% 10)
  }))
  clusters <- clusters[, c(-1), drop = FALSE]

  return(list(
    order = as.vector(hopach_results$clustering$order),
    expression = expression_data,
    clusters = clusters
  ))
}

scale_min_max <- function(x,
                          min_range = -2,
                          max_range = 2) {
  min_val <- min(x)
  max_val <- max(x)
  scaled_x <-
    (x - min_val) / (max_val - min_val) * (max_range - min_range) + min_range
  return(scaled_x)
}


load_isoform_set <- function(filenames,
                             prefixes,
                             read_colname,
                             rpkm_colname,
                             rpkm_colname_alias,
                             conditions,
                             intersect_by,
                             digits,
                             batch_metadata,
                             collected_data = NULL) {
  for (i in 1:length(filenames)) {
    isoforms <- read.table(
      filenames[i],
      sep = get_file_type(filenames[i]),
      header = TRUE,
      stringsAsFactors = FALSE
    )
    new_read_colname <- paste(prefixes[i], conditions, sep = "_")
    colnames(isoforms)[colnames(isoforms) == read_colname] <- new_read_colname
    colnames(isoforms)[colnames(isoforms) == rpkm_colname] <- paste(conditions, i, rpkm_colname,
                                                                    sep =
                                                                      " "
    )
    if (!is.null(batch_metadata)) {
      batch <- batch_metadata[prefixes[i], "batch"]
      print(
        paste(
          "Load ",
          nrow(isoforms),
          " rows from '",
          filenames[i],
          "' as '",
          new_read_colname,
          "', batch '",
          batch,
          "'",
          sep = ""
        )
      )
      column_data_frame <- data.frame(conditions, batch,
                                      row.names =
                                        c(new_read_colname)
      )
    } else {
      print(
        paste(
          "Load ",
          nrow(isoforms),
          " rows from '",
          filenames[i],
          "' as '",
          new_read_colname,
          "'",
          sep = ""
        )
      )
      column_data_frame <- data.frame(conditions, row.names = c(new_read_colname))
    }
    if (is.null(collected_data)) {
      collected_data <- list(
        collected_isoforms = isoforms,
        read_colnames = c(new_read_colname),
        column_data = column_data_frame
      )
    } else {
      collected_data$collected_isoforms <- merge(collected_data$collected_isoforms,
                                                 isoforms,
                                                 by = intersect_by,
                                                 sort = FALSE
      )
      collected_data$read_colnames <- c(collected_data$read_colnames, new_read_colname)
      collected_data$column_data <- rbind(collected_data$column_data, column_data_frame)
    }
  }
  rpkm_columns <- grep(
    paste("^", conditions, " [0-9]+ ", rpkm_colname, sep = ""),
    colnames(collected_data$collected_isoforms),
    value = TRUE,
    ignore.case = TRUE
  )
  collected_data$collected_isoforms[rpkm_colname_alias] <- format(rowSums(collected_data$collected_isoforms[, rpkm_columns, drop = FALSE]) / length(filenames),
                                                                  digits = digits
  )
  collected_data$rpkm_colnames <- c(collected_data$rpkm_colnames, rpkm_colname_alias)
  collected_data$collected_isoforms <- collected_data$collected_isoforms[, !colnames(collected_data$collected_isoforms) %in% rpkm_columns]
  return(collected_data)
}


export_ma_plot <- function(data,
                           rootname,
                           width = 800,
                           height = 800,
                           resolution = 72) {
  tryCatch(
    expr = {
      png(
        filename = paste(rootname, ".png", sep = ""),
        width = width,
        height = height,
        res = resolution
      )
      plotMA(data)
      dev.off()

      pdf(
        file = paste(rootname, ".pdf", sep = ""),
        width = round(width / resolution),
        height = round(height / resolution)
      )
      plotMA(data)
      dev.off()

      cat(paste("\nExport MA-plot to ", rootname, ".(png/pdf)", "\n",
                sep =
                  ""
      ))
    },
    error = function(e) {
      dev.off()
      cat(paste(
        "\nFailed to export MA-plot to ",
        rootname,
        ".(png/pdf)",
        "\n",
        sep = ""
      ))
    }
  )
}

export_heatmap <- function(mat_data,
                           column_data,
                           rootname,
                           width = 800,
                           height = 800,
                           resolution = 72) {
  tryCatch(
    expr = {
      png(
        filename = paste(rootname, ".png", sep = ""),
        width = width,
        height = height,
        res = resolution
      )
      pheatmap(
        mat = mat_data,
        main = "Top 30 genes from VST normalized read counts",
        annotation_col = column_data,
        cluster_rows = FALSE,
        show_rownames = TRUE,
        cluster_cols = FALSE
      )
      dev.off()

      pdf(
        file = paste(rootname, ".pdf", sep = ""),
        width = round(width / resolution),
        height = round(height / resolution)
      )
      pheatmap(
        mat = mat_data,
        main = "Top 30 genes from VST normalized read counts",
        annotation_col = column_data,
        cluster_rows = FALSE,
        show_rownames = TRUE,
        cluster_cols = FALSE
      )
      dev.off()

      cat(paste(
        "\nExport expression heatmap to ",
        rootname,
        ".(png/pdf)",
        "\n",
        sep = ""
      ))
    },
    error = function(e) {
      dev.off()
      cat(
        paste(
          "\nFailed to export expression heatmap to ",
          rootname,
          ".(png/pdf)",
          "\n",
          sep = ""
        )
      )
    }
  )
}


# Function to generate main effect contrasts with different reference levels
generate_main_effect_contrasts <- function(dds, factors, factor_levels) {
  contrasts <- list()

  for (factor in factors) {
    other_factors <- setdiff(factors, factor)

    for (other_factor in other_factors) {
      other_levels <- factor_levels[[other_factor]]

      for (other_level in other_levels) {
        # dds_subset <- dds[colData(dds)[[other_factor]] == other_level, ]
        dds_subset <- dds

        for (ref_level in factor_levels[[factor]]) {
          colData(dds_subset)[[factor]] <- relevel(colData(dds_subset)[[factor]], ref = ref_level)

          levels <- factor_levels[[factor]]

          for (level in levels) {
            if (level != ref_level) {
              specificity_group <- paste(other_factor, other_level, sep = "_")
              contrast <- paste(sort(c(level, ref_level)), collapse = "_vs_")
              contrasts <- append(contrasts, list(list(effect_type = "main",
                                                       specificity_group = specificity_group,
                                                       numerator = level,
                                                       denominator = ref_level,
                                                       contrast = paste(factor, level, "vs", ref_level, sep = "_"),
                                                       subset = dds_subset)))
            }
          }
        }
      }
    }
  }

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

            contrasts <- append(contrasts, list(list(effect_type = "interaction",
                                                     specificity_group = specificity_group,
                                                     numerator = numerator,
                                                     denominator = denominator,
                                                     contrast = interaction,
                                                     subset = dds_subset)))
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
            dds_subset <- dds[colData(dds)[[factor1]] == ref_level1,]
            colData(dds_subset)[[factor2]] <- relevel(colData(dds_subset)[[factor2]], ref = ref_level2)

            contrasts <- append(contrasts, list(list(effect_type = "interaction",
                                                     specificity_group = specificity_group,
                                                     numerator = numerator,
                                                     denominator = denominator,
                                                     contrast = interaction,
                                                     subset = dds_subset)))
          }
        }
      }
    }
  }

  return(contrasts)
}

# Function to get number of significant genes
get_contrast_res <- function(contrast) {
  dds_subset <- contrast$subset
  dds_subset <- DESeq(dds_subset, test = "Wald")
  res <- results(dds_subset,
                 name = contrast$contrast,
                 alpha = args$fdr,
                 lfcThreshold = ifelse(args$use_lfc_thresh, args$lfcthreshold, 0),
                 independentFiltering = TRUE)
  return(res)
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
      print("Test mode is ON, each sample will be limited to 100 rows")
      isoforms <- isoforms %>%
        dplyr::slice_head(n = 100)
    } else {
      print("Test mode is OFF, processing all rows")
    }

    print(paste("Load ", nrow(isoforms), " rows from ", filenames[i], sep = ""))
    colnames(isoforms)[colnames(isoforms) == read_colname] <- paste(prefixes[i], read_colname, sep = " ")
    colnames(isoforms)[colnames(isoforms) == rpkm_colname] <- paste(prefixes[i], rpkm_colname, sep = " ")
    if (is.null(collected_isoforms)) {
      collected_isoforms <- isoforms
    } else {
      collected_isoforms <- merge(collected_isoforms,
                                  isoforms,
                                  by = intersect_by,
                                  sort = FALSE)
    }
  }
  print(paste(
    "Number of rows common for all loaded files ",
    nrow(collected_isoforms),
    sep = ""
  ))

  return(collected_isoforms)
}

assert_args <- function(args) {
  print("Check input parameters")

  if (length(args$treated) == 1 || length(args$untreated) == 1) {
    args$batchfile <- NULL # reset batchfile to NULL. We don't need it for DESeq even if it was provided
  }

  if (!is.null(args$batchfile)) {
    batch_metadata <- read.table(
      args$batchfile,
      sep = get_file_type(args$batchfile),
      row.names = 1,
      col.names = c("name", "batch"),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    rownames(batch_metadata) <- gsub("'|\"| ", "_", rownames(batch_metadata))
    if (all(is.element(c(args$ualias, args$talias), rownames(batch_metadata)))) {
      args$batchfile <- batch_metadata # dataframe
    } else {
      cat("\nMissing values in batch metadata file. Skipping multi-factor analysis\n")
      args$batchfile <- NULL
    }
  }

  return(args)
}


get_args <- function() {
  parser <- ArgumentParser(description = "Run BioWardrobe DESeq/DESeq2 for untreated-vs-treated groups (condition-1-vs-condition-2)")
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
  parser$add_argument("-m",
                      "--meta",
                      help = "Metadata file to describe relation between samples, where first column corresponds to --name, formatted as CSV/TSV",
                      type = "character",
                      required = "True")
  parser$add_argument("-bf", "--batchfile",
                      help = "Metadata file for multi-factor analysis. Headerless TSV/CSV file. First column - names from --ualias and --talias, second column - batch group name. Default: None", type =
                        "character"
  )
  parser$add_argument(
    "--design_formula",
    help = "Design formula for DESeq2. Must be provided as a string, e.g., '~ condition + batch'.",
    type = "character",
    required = TRUE
  )
  parser$add_argument("-ci",
                      "--contrast_indices",
                      help = "Comma-separated list of integers representing contrast indices (e.g., 1,2,3)",
                      type = "character",
                      required = "True"
  )
  parser$add_argument("--test_mode",
                      help = "Run for test, only first 100 rows",
                      action = "store_true",
                      default = FALSE)
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
    "--cluster",
    help = paste(
      "Hopach clustering method to be run on normalized read counts for the",
      "exploratory visualization part of the analysis. Default: do not run",
      "clustering"
    ),
    type = "character",
    choices = c("row", "column", "both")
  )
  parser$add_argument(
    "--batchcorrection",
    help = paste(
      "Specifies the batch correction method to be applied.
      - 'combatseq' applies ComBat_seq at the beginning of the analysis, removing batch effects from the design formula before differential expression analysis.
      - 'limmaremovebatcheffect' applies removeBatchEffect from the limma package after differential expression analysis, incorporating batch effects into the model during DE analysis.
      - Default: none"
    ),
    type = "character",
    choices = c("none", "combatseq", "limmaremovebatcheffect"),
    default = "none"
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
  parser$add_argument("-o",
                      "--output",
                      help = "Output prefix. Default: deseq",
                      type = "character",
                      default = "./deseq"
  )
  parser$add_argument("-d",
                      "--digits",
                      help = "Precision, number of digits to print. Default: 3",
                      type = "integer",
                      default = 3
  )
  parser$add_argument(
    "-p",
    "--threads",
    help = "Threads",
    type = "integer",
    default = 1
  )
  args <- assert_args(parser$parse_args(gsub(
    "'|\"| ", "_", commandArgs(trailingOnly = TRUE)
  )))
  return(args)
}

# Function to extract design formula, factor levels, and generate specific contrasts based on contrast_vector
generate_contrasts <- function(dds, contrast_vector) {
  # Extract the design formula from the DESeq2 object
  design_formula <- design(dds)

  # Get the levels of each factor in the design
  factors <- all.vars(design_formula)
  factor_levels <- lapply(factors, function(f) levels(colData(dds)[[f]]))
  names(factor_levels) <- factors

  # Print the factor levels for debugging/verification
  print(factor_levels)

  # Generate main effect contrasts
  main_contrasts <- generate_main_effect_contrasts(dds, factors, factor_levels)

  # Generate interaction effect contrasts
  interaction_contrasts <- generate_interaction_effect_contrasts(dds)

  # Combine main and interaction contrasts
  all_contrasts <- c(main_contrasts, interaction_contrasts)

  # Select contrasts based on the contrast_vector
  selected_contrasts <- all_contrasts[contrast_vector]

  return(selected_contrasts)
}

# Function to export results (plots, reports, etc.)
export_contrast_results <- function(res, contrast_name, dds, column_data, output) {
  # Normalized counts for PCA and heatmap
  normCounts <- counts(dds, normalized = TRUE)

  # Generate and export plots for each contrast
  export_ma_plot(res, paste0(output, "_", contrast_name, "_ma_plot"))
  # Heatmap for the top 30 most expressed genes
  vst <- varianceStabilizingTransformation(dds, blind = FALSE)
  vsd <- assay(vst)
  rownames(vsd) <- rownames(normCounts)
  mat <- vsd[order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:30],]
  export_heatmap(mat, column_data, paste0(output, "_", contrast_name, "_heatmap"))

  # Export the DESeq2 results to file
  results_filename <- paste0(output, "_", contrast_name, "_results.tsv")
  write.table(as.data.frame(res), file = results_filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  print(paste("Results for", contrast_name, "exported to", results_filename))
}

# Function to clean DESeq2 results (handle NA values)
clean_deseq_results <- function(DESeqRes) {
  DESeqRes$log2FoldChange[is.na(DESeqRes$log2FoldChange)] <- 0
  DESeqRes$pvalue[is.na(DESeqRes$pvalue)] <- 1
  DESeqRes$padj[is.na(DESeqRes$padj)] <- 1
  return(DESeqRes)
}

# Function to add metadata columns and create the final results data frame
add_metadata_to_results <- function(collected_isoforms, DESeqRes, read_count_cols, digits) {
  collected_isoforms <- data.frame(
    cbind(collected_isoforms[, !colnames(collected_isoforms) %in% read_count_cols], DESeqRes),
    check.names = F,
    check.rows = F
  )
  collected_isoforms[, "'-LOG10(pval)'"] <- format(-log(as.numeric(collected_isoforms$pvalue), 10), digits = digits)
  collected_isoforms[, "'-LOG10(padj)'"] <- format(-log(as.numeric(collected_isoforms$padj), 10), digits = digits)
  return(collected_isoforms)
}

# Function to export DESeq2 results to a TSV file
export_deseq_report <- function(collected_isoforms, output_prefix) {
  collected_isoforms_filename <- paste0(output_prefix, "_gene_exp_table.tsv")
  write.table(
    collected_isoforms,
    file = collected_isoforms_filename,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  print(paste("Export DESeq report to", collected_isoforms_filename))
}

# Function to export normalized counts and filtered counts to GCT format
export_gct_data <- function(normCounts, collected_isoforms, column_data, output_prefix) {
  row_metadata <- collected_isoforms %>%
    dplyr::mutate_at("GeneId", toupper) %>%
    dplyr::distinct(GeneId, .keep_all = TRUE) %>%
    remove_rownames() %>%
    column_to_rownames("GeneId") %>%
    dplyr::select(log2FoldChange, pvalue, padj) %>%
    arrange(desc(log2FoldChange))

  col_metadata <- column_data %>%
    mutate_at(colnames(.), as.vector)

  rownames(normCounts) <- toupper(collected_isoforms[, c("GeneId")])

  print("Exporting GCT data")
  print("Row metadata:")
  print(head(row_metadata))
  print("Column metadata:")
  print(head(col_metadata))
  print("Normalized counts:")
  print(head(normCounts))
  # Export unfiltered data
  export_gct(
    counts_mat = normCounts,
    row_metadata = row_metadata,
    col_metadata = col_metadata,
    location = paste0(output_prefix, "_counts_all.gct")
  )

  # Filter rows by padj <= FDR and export
  row_metadata <- row_metadata %>%
    dplyr::filter(.$padj <= args$fdr)

  filtered_normCounts <- normCounts[as.vector(rownames(row_metadata)),]

  export_gct(
    counts_mat = filtered_normCounts,
    row_metadata = row_metadata,
    col_metadata = col_metadata,
    location = paste0(output_prefix, "_counts_filtered.gct")
  )
}

# Function to cluster data and re-order based on clustering results
cluster_and_reorder <- function(normCounts, col_metadata, row_metadata, args) {
  if (!is.null(args$cluster)) {
    if (args$cluster == "column" || args$cluster == "both") {
      print("Clustering filtered read counts by columns")
      clustered_data <- get_clustered_data(expression_data = normCounts, dist = args$columndist, transpose = TRUE)
      col_metadata <- cbind(col_metadata, clustered_data$clusters) # Add cluster labels
      col_metadata <- col_metadata[clustered_data$order,] # Reorder based on clustering results
      print("Reordered samples")
      print(col_metadata)
    }
    if (args$cluster == "row" || args$cluster == "both") {
      print("Clustering filtered normalized read counts by rows")
      clustered_data <- get_clustered_data(expression_data = normCounts, dist = args$rowdist, transpose = FALSE)
      normCounts <- clustered_data$expression # Adjust for row centering
      row_metadata <- cbind(row_metadata, clustered_data$clusters) # Add cluster labels
      row_metadata <- row_metadata[clustered_data$order,] # Reorder based on clustering results
      print("Reordered features")
      print(head(row_metadata))
    }
  }
  return(list(normCounts = normCounts, col_metadata = col_metadata, row_metadata = row_metadata))
}

# Main DESeq2 workflow function
run_deseq2_analysis <- function(countData, column_data, design_formula, contrast_vector, output, batchfile = NULL, batchcorrection = NULL, threads = 1) {
  # Register multiple threads for parallel processing
  register(MulticoreParam(threads))

  # Define design and apply batch correction if applicable
  if (!is.null(batchfile) && batchcorrection == "combatseq") {
    # ComBat-seq adjusts the raw count matrix to remove batch effects before DESeq2
    countData <- ComBat_seq(as.matrix(countData), batch = batchfile, covar_mod = model.matrix(as.formula(design_formula), data = column_data))
    design <- as.formula(design_formula)

  } else if (!is.null(batchfile) && batchcorrection == "limmaremovebatcheffect") {
    # Limma's removeBatchEffect function is used on normalized data post DESeq2
    design <- as.formula(paste(design_formula, "+ batch"))

  } else {
    # No batch correction applied, keep original design formula
    design <- as.formula(design_formula)
  }

  # Create DESeq2 dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = column_data, design = design)
  dds <- DESeq2::DESeq(dds, test = "Wald", quiet = FALSE, parallel = TRUE)

  print("DESeq2 object created.")
  print(dds)

  # Generate all contrasts
  all_contrasts <- generate_contrasts(dds, contrast_vector)

  print("All contrasts list: ")
  print(all_contrasts)
  print(contrast_vector)
  print(all_contrasts[[1]])

  # Iterate through each index in contrast_vector and process the corresponding contrast
  for (contrast_index in 1:length(contrast_vector)) {
    contrast_name <- paste("contrast_", contrast_index, sep = "")

    print(paste("Processing contrast:", contrast_name))

    contrast <- all_contrasts[[contrast_index]]  # Get the contrast from the list
    # Get DESeq2 results for the specific contrast
    res <- get_contrast_res(contrast)

    print("DESeq2 results obtained.")
    print(summary(res))

    # Clean results
    DESeqRes <- clean_deseq_results(as.data.frame(res[, c(1, 2, 5, 6)]))

    print("DESeq2 results cleaned.")
    print(DESeqRes)

    # Add metadata to results
    final_isoforms <- add_metadata_to_results(expression_data_df, DESeqRes, read_counts_data_df, digits = 4)

    print("Metadata added to results.")

    # Export DESeq2 report
    export_deseq_report(final_isoforms, paste0(output, "_", contrast_name))
    print("Exported DESeq2 report.")
    print(str(final_isoforms))

    # Export normalized read counts to GCT format
    normCounts <- counts(dds, normalized = TRUE)
    export_gct_data(normCounts, final_isoforms, column_data, paste0(output, "_", contrast_name))
    print("Exported normalized counts to GCT format.")

    # Optionally perform clustering and export
    # if (!is.null(args$cluster)) {
    #   clustered_data <- cluster_and_reorder(normCounts, col_metadata, row_metadata, args)
    #   normCounts <- clustered_data$normCounts
    #   col_metadata <- column_data
    #   row_metadata <- expression_data_df
    # }
  }

  print("DESeq2 analysis complete.")
}


# Parse the arguments
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
design_formula <- as.formula(args$design)
print("Load design formula")
print(design_formula)

# Load expression data
expression_data_df <- load_expression_data(args$input, args$name, READ_COL, RPKM_COL, INTERSECT_BY)
print("Expression data")
print(head(expression_data_df))


# Select all columns with read counts data, reorder them based on the row names from metadata_df
read_counts_columns <- grep(
  paste(READ_COL, sep = ""),
  colnames(expression_data_df),
  value = TRUE,
  ignore.case = TRUE
)
read_counts_data_df <- expression_data_df[read_counts_columns]
colnames(read_counts_data_df) <- lapply(colnames(read_counts_data_df), function(s) {
  paste(head(unlist(strsplit(s, " ", fixed = TRUE)), -1), collapse = " ")
})
print("Read counts data")
print(head(read_counts_data_df))

# Convert both rownames of metadata_df and colnames of read_counts_data_df to lowercase
print("Check if metadata file rows are present in read counts data columns")
print(setdiff(rownames(metadata_df), colnames(read_counts_data_df)))
rownames(metadata_df) <- tolower(rownames(metadata_df))
colnames(read_counts_data_df) <- tolower(colnames(read_counts_data_df))
print("Read counts data columns")
print(colnames(read_counts_data_df))
print("Metadata file rows")
print(rownames(metadata_df))
# Trim whitespace from column names and row names
colnames(read_counts_data_df) <- trimws(colnames(read_counts_data_df))
rownames(metadata_df) <- trimws(rownames(metadata_df))

# Remove any non-visible characters
colnames(read_counts_data_df) <- gsub("[^[:alnum:]_]", "", colnames(read_counts_data_df))
rownames(metadata_df) <- gsub("[^[:alnum:]_]", "", rownames(metadata_df))

# Convert the comma-separated string into a vector of integers
contrast_vector <- as.integer(unlist(strsplit(args$contrast_indices, ",")))

# Ensure colnames of countData and rownames of column_data match exactly
metadata_df <- metadata_df[match(colnames(read_counts_data_df), rownames(metadata_df)),]

# Check if any mismatch after reordering
stopifnot(all(rownames(metadata_df) == colnames(read_counts_data_df)))

print("countData")
print(read_counts_data_df)
print("metadata_df")
print(metadata_df)

run_deseq2_analysis(
  countData = read_counts_data_df,
  column_data = metadata_df,
  design_formula = args$design_formula,
  contrast_vector = contrast_vector,
  output = args$output,
  fdr = args$fdr,
  use_lfc_thresh = args$use_lfc_thresh,
  lfcthreshold = args$lfcthreshold,
  batchfile = args$batchfile,
  batchcorrection = args$batchcorrection,
  threads = args$threads
)