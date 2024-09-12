#!/usr/bin/env Rscript
options(warn = -1)
options("width" = 400)

suppressMessages(library(argparse))
suppressMessages(library(BiocParallel))
suppressMessages(library(pheatmap))
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))

################################################################################################
# v0.0.1
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
################################################################################################


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
  if (length(args$input) != length(args$name)) {
    print("Exiting: --input and --name have different number of values")
    quit(save = "no",
         status = 1,
         runLast = FALSE)
  }
  tryCatch(
    expr = {
      # Try to load design formula
      design_formula <- as.formula(args$design)
    },
    error = function(e) {
      print(paste(
        "Exiting: failed to load --design ",
        args$design,
        " as formula",
        sep = ""
      ))
      quit(save = "no",
           status = 1,
           runLast = FALSE)
    }
  )
  tryCatch(
    expr = {
      # Try to load reduced formula
      reduced_formula <- as.formula(args$reduced)
    },
    error = function(e) {
      print(paste(
        "Exiting: failed to load --reduced ",
        args$reduced,
        " as formula",
        sep = ""
      ))
      quit(save = "no",
           status = 1,
           runLast = FALSE)
    }
  )
  return(args)
}


get_args <- function() {
  parser <- ArgumentParser(description = "Run DeSeq2 for multi-factor analysis using LRT (likelihood ratio or chi-squared test)")
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
  parser$add_argument(
    "-d",
    "--design",
    help = "Design formula. Should start with ~. See DeSeq2 manual for details",
    type = "character",
    required = "True"
  )
  parser$add_argument(
    "-r",
    "--reduced",
    help = "Reduced formula to compare against with the term(s) of interest removed. Should start with ~. See DeSeq2 manual for details",
    type = "character",
    required = "True"
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
  parser$add_argument("-o",
                      "--output",
                      help = "Output prefix for generated files",
                      type = "character",
                      default = "./deseq")
  parser$add_argument(
    "-p",
    "--threads",
    help = "Threads number",
    type = "integer",
    default = 1
  )
  parser$add_argument("--test_mode",
                      help = "Run for test, only first 100 rows",
                      action = "store_true",
                      default = FALSE)
  args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
  return(args)
}

generate_lrt_md <- function(deseq_results, full_formula, reduced_formula, output_file, alpha = 0.1) {
  # Initialize the markdown content
  md_content <- ""

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

    outliers <- gsub(".*: ", "", summary_output[6])  # Outliers
    low_counts <- gsub(".*: ", "", summary_output[7])  # Low counts (independent filtering)
    mean_count <- gsub("[^0-9]", "", summary_output[8])

    # Add a summary of significant genes
    lrt_summary <- paste0(
      "### Results Summary\n\n",
      "From this LRT analysis, **", significant_genes, " genes** (out of ", total_genes, " tested) are identified as significant with a padj value < ", alpha, ".\n\n"
    )

    # Add the information about outliers and low counts
    lrt_summary <- paste0(lrt_summary,
      "**Outliers**<sup>1</sup>: ", outliers, " of genes were detected as outliers and excluded from analysis.\n\n",
      "**Low counts**<sup>2</sup>: ", low_counts, " of genes were removed due to low counts (mean <", mean_count,  ") and independent filtering.\n\n",
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
            dds_subset <- dds[colData(dds)[[factor1]] == ref_level1, ]
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
get_num_significant_genes <- function(contrast) {
  dds_subset <- contrast$subset
  dds_subset <- DESeq(dds_subset, test = "Wald")
  res <- results(dds_subset, 
                 name = contrast$contrast,
                 alpha = args$fdr,
                 lfcThreshold = ifelse(args$use_lfc_thresh, args$lfcthreshold, 0),
                 independentFiltering = TRUE)
  return(sum(res$padj < args$fdr, na.rm = TRUE))
}

# Main function to generate the dataframe with all possible contrasts
generate_contrasts <- function(dds) {
  # Extract the design formula and model matrix
  design_formula <- design(dds)
  
  # Get the levels of each factor in the design
  factors <- all.vars(design_formula)
  factor_levels <- lapply(factors, function(f) levels(colData(dds)[[f]]))
  names(factor_levels) <- factors
  print(factor_levels)
  
  # stop("Stop here")
  
  # Generate main and interaction effect contrasts
  main_contrasts <- generate_main_effect_contrasts(dds, factors, factor_levels)
  interaction_contrasts <- generate_interaction_effect_contrasts(dds)
  all_contrasts <- c(main_contrasts, interaction_contrasts)
  
  # Create a dataframe to store the contrasts and number of significant genes
  contrast_df <- data.frame(
    effect = character(),
    specificity_group = character(),
    contrast = character(),
    numerator = character(),
    denominator = character(),
    num_significant_genes = integer(),
    stringsAsFactors = FALSE
  )
  
  # Run DESeq2 for each contrast and store the results
  count = 1
  for (contrast in all_contrasts) {
    print(count)
    num_significant_genes <- get_num_significant_genes(contrast)
    contrast_df <- rbind(contrast_df, data.frame(
      effect = contrast$effect_type,
      contrast = contrast$contrast,
      specificity_group = contrast$specificity_group,
      numerator = contrast$numerator,
      denominator = contrast$denominator,
      num_significant_genes = num_significant_genes,
      stringsAsFactors = FALSE
    ))
    count = count + 1
  }
  
  # Remove duplicate contrasts
  contrast_df <- contrast_df %>% 
    group_by(specificity_group) %>% 
    distinct(contrast, .keep_all = TRUE) %>% 
    arrange(specificity_group) %>%
    mutate(contrast_number = row_number()) %>%
    select(contrast_number, everything())
  
  # Sort the dataframe by the number of significant genes
  return(contrast_df)
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
design_formula <- as.formula(args$design)
print("Load design formula")
print(design_formula)

# Load reduced formula
reduced_formula <- as.formula(args$reduced)
print("Load reduced formula")
print(reduced_formula)

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

tryCatch(
  expr = {
    # Try reorder columns in read_counts_data_df based on metadata_df rownames
    print("Reorder read count columns based on the row names from the provided metadata file")
    read_counts_data_df <- read_counts_data_df[, rownames(metadata_df)]
    print(head(read_counts_data_df))
  },
  error = function(e) {
    print(
      "Exiting: failed to reorder read count columns based on the row names from the provided metadata file"
    )
    print(paste("Count data columns -", paste(
      colnames(read_counts_data_df), collapse = " "
    ), sep = " "))
    print(paste("Metadata file rows -", paste(rownames(metadata_df), collapse = " "), sep = " "))
    quit(save = "no",
         status = 1,
         runLast = FALSE)
  }
)

dse <- DESeqDataSetFromMatrix(countData = read_counts_data_df,
                              colData = metadata_df,
                              design = design_formula)

print("Run DESeq2 using Wald")
dsq_wald <- DESeq(
  dse,
  test = "Wald",
  quiet = FALSE,
  parallel = TRUE
)

print("Run DESeq2 using LRT")
dsq_lrt <- DESeq(
  dse,
  test = "LRT",
  reduced = reduced_formula,
  quiet = FALSE,
  parallel = TRUE
)

print("Generate contrasts")
all_contrasts <- generate_contrasts(dsq_wald)

dsq_lrt_res <- results(dsq_lrt,
               alpha = args$fdr,
               independentFiltering = TRUE)

print("LRT Results description")
print(mcols(dsq_lrt_res))

# Filter DESeq2 output
res_filtered <- as.data.frame(dsq_lrt_res[, c(2, 5, 6)])
res_filtered$log2FoldChange[is.na(res_filtered$log2FoldChange)] <- 0
res_filtered[is.na(res_filtered)] <- 1
# Export results to TSV file
expression_data_df <- data.frame(
  cbind(expression_data_df[, ], res_filtered),
  check.names = F,
  check.rows = F
)
expression_data_df[, "-LOG10(pval)"] <- -log(as.numeric(expression_data_df$pval), 10)
expression_data_df[, "-LOG10(padj)"] <- -log(as.numeric(expression_data_df$padj), 10)

contrasts_filename <- paste(args$output, "_contrasts_table.tsv", sep = "")
print(all_contrasts)
print(contrasts_filename)
write.table(
  all_contrasts,
  file = contrasts_filename,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

print(paste("Export contrasts to", contrasts_filename, sep = " "))

lrt_report_filename <- paste(args$output, "_lrt_result.md", sep = "")
summary(dsq_lrt_res)
generate_lrt_md(dsq_lrt_res, args$design, args$reduced, lrt_report_filename)
print(paste("Export LRT markdown report to", lrt_report_filename, sep = " "))

results_filename <- paste(args$output, "_gene_exp_table.tsv", sep = "")
write.table(
  expression_data_df,
  file = results_filename,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
print(paste("Export results to", results_filename, sep = " "))
graphics.off()
