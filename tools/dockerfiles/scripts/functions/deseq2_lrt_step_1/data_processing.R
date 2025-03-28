#!/usr/bin/env Rscript

# --- Data processing functions ---

# Constants for column name patterns
READ_COL <- "Read"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- "GeneId"

# Function to check the delimiter used in a file
check_file_delimiter <- function(file_path) {
  # Default to tab delimiter
  delimiter <- "\t"
  
  # Check file extension
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    delimiter <- ","
  } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    delimiter <- "\t"
  } else {
    # Try to detect delimiter from file content
    first_line <- readLines(file_path, n = 1)
    
    if (grepl(",", first_line)) {
      delimiter <- ","
    } else if (grepl("\t", first_line)) {
      delimiter <- "\t"
    } else if (grepl(";", first_line)) {
      delimiter <- ";"
    }
  }
  
  message(paste("Detected delimiter for file", file_path, ":", delimiter))
  return(delimiter)
}

# Function to load expression data from multiple files
load_expression_data <- function(input_files, sample_names, read_col="Read", rpkm_col="Rpkm", intersect_by="GeneId") {
  # Process each input file
  expression_data_list <- list()
  
  for (i in 1:length(input_files)) {
    input_file <- input_files[i]
    sample_name <- sample_names[i]
    
    # Determine file type
    delimiter <- check_file_delimiter(input_file)
    
    # Load data
    message(paste("Loading file:", input_file))
    
    # Use tryCatch to handle errors in file reading
    data <- tryCatch({
      read.table(
        input_file, 
        sep = delimiter,
        header = TRUE, 
        quote = "",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        comment.char = "#"
      )
    }, error = function(e) {
      stop(paste("Error reading file", input_file, ":", e$message))
    })
    
    # Rename columns to include sample name
    data_cols <- colnames(data)
    for (col in data_cols) {
      if (grepl(read_col, col, ignore.case = TRUE) || grepl(rpkm_col, col, ignore.case = TRUE)) {
        # Rename to include sample name
        new_col <- paste(sample_name, col, sep = " ")
        colnames(data)[colnames(data) == col] <- new_col
      }
    }
    
    # Check for duplicate column names
    if (any(duplicated(colnames(data)))) {
      dup_cols <- colnames(data)[duplicated(colnames(data))]
      warning(paste("Duplicate column names found in file", input_file, ":", paste(dup_cols, collapse=", ")))
      
      # Make column names unique
      colnames(data) <- make.unique(colnames(data))
      message("Column names made unique using make.unique()")
    }
    
    expression_data_list[[i]] <- data
  }
  
  # Merge data by gene ID using a safer approach with explicit base:: namespace
  tryCatch({
    # Check for required column in each data frame
    for (i in seq_along(expression_data_list)) {
      if (!intersect_by %in% colnames(expression_data_list[[i]])) {
        stop(paste("Column", intersect_by, "not found in file", input_files[i]))
      }
    }
    
    # Merge data frames one by one
    merged_data <- expression_data_list[[1]]
    
    if (length(expression_data_list) > 1) {
      for (i in 2:length(expression_data_list)) {
        # Check for duplicate column names before merging
        common_cols <- base::intersect(
          base::setdiff(colnames(merged_data), intersect_by),
          base::setdiff(colnames(expression_data_list[[i]]), intersect_by)
        )
        
        if (length(common_cols) > 0) {
          warning(paste("Common columns found while merging file", input_files[i], ":", paste(common_cols, collapse=", ")))
          
          # Make names unique in current data frame before merge
          rename_cols <- base::setdiff(colnames(expression_data_list[[i]]), intersect_by)
          new_names <- paste0(rename_cols, "_", i)
          names(new_names) <- rename_cols
          
          # Rename columns manually without dplyr
          for (old_name in names(new_names)) {
            pos <- which(colnames(expression_data_list[[i]]) == old_name)
            if (length(pos) > 0) {
              colnames(expression_data_list[[i]])[pos] <- new_names[old_name]
            }
          }
          
          message("Made column names unique by adding suffix before merging")
        }
        
        merged_data <- merge(
          merged_data, 
          expression_data_list[[i]], 
          by = intersect_by, 
          all = TRUE, 
          sort = FALSE,
          suffixes = c("", paste0("_", i))
        )
      }
    }
    
    # Final check for duplicate column names in merged data
    if (any(duplicated(colnames(merged_data)))) {
      dup_cols <- colnames(merged_data)[duplicated(colnames(merged_data))]
      warning(paste("Duplicate column names in final merged data:", paste(dup_cols, collapse=", ")))
      
      # Make all column names unique
      colnames(merged_data) <- make.unique(colnames(merged_data))
      message("Final merged data column names made unique")
    }
    
  }, error = function(e) {
    stop(paste("Error merging expression data:", e$message))
  })
  
  # Replace NA values with 0
  merged_data[is.na(merged_data)] <- 0
  
  return(merged_data)
}

# Function to filter expression data based on RPKM values
filter_rpkm <- function(expression_data, rpkm_cutoff) {
  if (is.null(rpkm_cutoff) || rpkm_cutoff <= 0) {
    message("No RPKM filtering applied")
    return(expression_data)
  }
  
  message(paste("Filtering data with RPKM cutoff:", rpkm_cutoff))
  
  # Get RPKM columns
  rpkm_cols <- grep(RPKM_COL, colnames(expression_data), ignore.case = TRUE, value = TRUE)
  
  if (length(rpkm_cols) == 0) {
    warning("No RPKM columns found, skipping RPKM filtering")
    return(expression_data)
  }
  
  # Keep rows where any RPKM value exceeds cutoff
  keep_rows <- apply(expression_data[, rpkm_cols, drop = FALSE], 1, function(x) {
    any(x > rpkm_cutoff, na.rm = TRUE)
  })
  
  filtered_data <- expression_data[keep_rows, , drop = FALSE]
  
  message(paste("RPKM filtering: removed", sum(!keep_rows), "rows, kept", sum(keep_rows), "rows"))
  
  return(filtered_data)
}

# Function to process count data and apply batch correction
process_count_data <- function(args, counts, metadata, design_formula) {
  message("Processing count data...")
  
  # Set up variables for batch correction result
  batch_warning <- NULL
  countData <- counts
  
  # Apply batch correction if requested
  if (args$batchcorrection != "none") {
    message(paste("Applying batch correction method:", args$batchcorrection))
    
    if (!"batch" %in% colnames(metadata)) {
      batch_warning <- "Batch correction requested but no 'batch' column found in metadata. Batch correction skipped."
      warning(batch_warning)
    } else {
      # Check batch correction method
      if (args$batchcorrection == "combatseq") {
        # Apply ComBat-Seq batch correction
        countData <- apply_combatseq_correction(counts, metadata, design_formula)
      } else if (args$batchcorrection == "model") {
        # Include batch in model (no data modification needed)
        design_formula <- update(design_formula, ~ . + batch)
        message("Added batch to design formula: ", deparse(design_formula))
      } else {
        batch_warning <- paste("Unknown batch correction method:", args$batchcorrection, ". Batch correction skipped.")
        warning(batch_warning)
      }
    }
  }
  
  # Return processed data
  return(list(
    countData = countData,
    design_formula = design_formula,
    batch_warning = batch_warning
  ))
}

# Helper function for ComBat-Seq batch correction
apply_combatseq_correction <- function(counts, metadata, design_formula) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    warning("Package 'sva' is required for ComBat-Seq batch correction but not available. Skipping batch correction.")
    return(counts)
  }
  
  message("Applying ComBat-Seq batch correction...")
  
  # Extract batch information
  batch <- metadata$batch
  
  # Get model matrix from design formula
  mod <- model.matrix(design_formula, data = metadata)
  
  # Apply ComBat-Seq
  corrected_counts <- sva::ComBat_seq(
    counts = as.matrix(counts),
    batch = batch,
    group = mod
  )
  
  # Convert back to data frame
  corrected_counts <- as.data.frame(corrected_counts)
  
  # Ensure row and column names are preserved
  rownames(corrected_counts) <- rownames(counts)
  colnames(corrected_counts) <- colnames(counts)
  
  message("Batch correction applied successfully")
  
  return(corrected_counts)
}

# Function to scale expression data based on scaling_type
scale_expression_data <- function(expression_data, scaling_type) {
  # Ensure scaling_type is a character
  scaling_type <- as.character(scaling_type)
  
  if (scaling_type == "none" || is.null(scaling_type)) {
    message("No scaling applied to expression data")
    return(expression_data)
  }
  
  message(paste("Applying", scaling_type, "scaling to expression data"))
  
  # Apply scaling
  if (scaling_type == "minmax") {
    # Min-max scaling to range [-2, 2]
    scaled_data <- apply(expression_data, 1, function(x) {
      min_val <- min(x, na.rm = TRUE)
      max_val <- max(x, na.rm = TRUE)
      if (max_val == min_val) {
        return(rep(0, length(x)))
      }
      return(4 * (x - min_val) / (max_val - min_val) - 2)
    })
    # Transpose back to original orientation
    scaled_data <- t(scaled_data)
    
  } else if (scaling_type == "zscore") {
    # Z-score scaling (mean=0, sd=1)
    scaled_data <- t(apply(expression_data, 1, function(x) {
      mean_val <- mean(x, na.rm = TRUE)
      sd_val <- sd(x, na.rm = TRUE)
      if (sd_val == 0) {
        return(rep(0, length(x)))
      }
      return((x - mean_val) / sd_val)
    }))
  } else {
    stop(paste("Unsupported scaling type:", scaling_type))
  }
  
  return(scaled_data)
}

# Function to perform clustering based on cluster method
perform_clustering <- function(expression_data, cluster_method, row_distance, column_distance, k, kmax) {
  if (cluster_method == "none" || is.null(cluster_method)) {
    message("No clustering performed")
    return(list(expression_data = expression_data, row_order = 1:nrow(expression_data), col_order = 1:ncol(expression_data)))
  }
  
  message(paste("Performing", cluster_method, "clustering"))
  
  # Load required library
  require(hopach)
  
  # Initialize clustering results
  row_clustering <- NULL
  col_clustering <- NULL
  
  # Check what parameters are accepted by hopach in this version
  hopach_args <- names(formals(hopach::hopach))
  message("Available hopach parameters: ", paste(hopach_args, collapse=", "))
  
  # Perform row clustering if requested
  if (cluster_method %in% c("row", "both")) {
    message(paste("Row clustering with distance:", row_distance))
    
    tryCatch({
      # Try modern parameter names first
      if ("distmethod" %in% hopach_args) {
        row_clustering <- hopach::hopach(
          data = expression_data,
          diss = NULL,
          distmethod = as.character(row_distance),
          K = as.numeric(k),
          kmax = as.numeric(kmax)
        )
      } 
      # Fall back to old parameter names
      else if ("dist.method" %in% hopach_args) {
        row_clustering <- hopach::hopach(
          data = expression_data,
          dmat = NULL,
          dist.method = as.character(row_distance),
          k = as.numeric(k),
          kmax = as.numeric(kmax),
          khigh = NULL
        )
      }
      # If neither works, use minimal parameters
      else {
        row_clustering <- hopach::hopach(
          data = expression_data
        )
      }
    }, error = function(e) {
      # Attempt with default parameters if specific ones fail
      message(paste("Error in row clustering:", e$message))
      message("Attempting row clustering with default parameters")
      tryCatch({
        row_clustering <<- hopach::hopach(data = expression_data)
      }, error = function(e2) {
        message(paste("Error in basic row clustering:", e2$message))
        message("Skipping row clustering")
      })
    })
  }
  
  # Perform column clustering if requested
  if (cluster_method %in% c("column", "both")) {
    message(paste("Column clustering with distance:", column_distance))
    
    tryCatch({
      # Try modern parameter names first
      if ("distmethod" %in% hopach_args) {
        col_clustering <- hopach::hopach(
          data = t(expression_data),  # Transpose for column clustering
          diss = NULL,
          distmethod = as.character(column_distance),
          K = as.numeric(k),
          kmax = as.numeric(kmax)
        )
      } 
      # Fall back to old parameter names
      else if ("dist.method" %in% hopach_args) {
        col_clustering <- hopach::hopach(
          data = t(expression_data),  # Transpose for column clustering
          dmat = NULL,
          dist.method = as.character(column_distance),
          k = as.numeric(k),
          kmax = as.numeric(kmax),
          khigh = NULL
        )
      }
      # If neither works, use minimal parameters
      else {
        col_clustering <- hopach::hopach(
          data = t(expression_data)  # Transpose for column clustering
        )
      }
    }, error = function(e) {
      # Attempt with default parameters if specific ones fail
      message(paste("Error in column clustering:", e$message))
      message("Attempting column clustering with default parameters")
      tryCatch({
        col_clustering <<- hopach::hopach(data = t(expression_data))
      }, error = function(e2) {
        message(paste("Error in basic column clustering:", e2$message))
        message("Skipping column clustering")
      })
    })
  }
  
  # Handle the case where clustering fails completely
  if (is.null(row_clustering) && cluster_method %in% c("row", "both")) {
    message("Row clustering failed, using original row order")
    row_order <- 1:nrow(expression_data)
  } else if (!is.null(row_clustering)) {
    row_order <- row_clustering$final$order
  } else {
    row_order <- 1:nrow(expression_data)
  }
  
  if (is.null(col_clustering) && cluster_method %in% c("column", "both")) {
    message("Column clustering failed, using original column order")
    col_order <- 1:ncol(expression_data)
  } else if (!is.null(col_clustering)) {
    col_order <- col_clustering$final$order
  } else {
    col_order <- 1:ncol(expression_data)
  }
  
  # Reorder data based on clustering
  ordered_data <- expression_data[row_order, col_order]
  
  return(list(
    expression_data = ordered_data,
    row_order = row_order,
    col_order = col_order,
    row_clustering = row_clustering,
    col_clustering = col_clustering
  ))
}

# Function to filter expression data based on RPKM cutoff
filter_by_rpkm <- function(expression_data, sample_metadata, rpkm_cutoff) {
  if (is.null(rpkm_cutoff)) {
    message("No RPKM-based filtering applied")
    return(expression_data)
  }
  
  message(paste("Filtering data with RPKM cutoff:", rpkm_cutoff))
  
  # Find columns containing "Rpkm" in their names
  rpkm_cols <- grep("Rpkm", colnames(expression_data), ignore.case = TRUE)
  
  if (length(rpkm_cols) == 0) {
    warning("No columns with 'Rpkm' in their names found, no filtering applied")
    return(expression_data)
  }
  
  # Keep rows where any RPKM column exceeds the cutoff
  keep_rows <- apply(expression_data[, rpkm_cols, drop = FALSE], 1, function(x) {
    any(x > rpkm_cutoff, na.rm = TRUE)
  })
  
  message(paste("Filtered out", sum(!keep_rows), "of", nrow(expression_data), "rows"))
  
  return(expression_data[keep_rows, ])
}

# Function to load and validate metadata file
load_metadata <- function(metadata_file) {
  if (!file.exists(metadata_file)) {
    stop(paste("Metadata file not found:", metadata_file))
  }
  
  # Determine file type
  delimiter <- check_file_delimiter(metadata_file)
  
  # Load metadata
  metadata <- tryCatch({
    read.table(
      metadata_file, 
      sep = delimiter,
      header = TRUE, 
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = "#"
    )
  }, error = function(e) {
    stop(paste("Error reading metadata file:", e$message))
  })
  
  # Check if metadata has required columns
  required_cols <- c("sample", "treatment", "cond")
  missing_cols <- required_cols[!required_cols %in% colnames(metadata)]
  
  if (length(missing_cols) > 0) {
    logger::warn(paste("Metadata is missing these expected columns:", paste(missing_cols, collapse=", ")))
  }
  
  # Ensure sample column exists and is unique
  if ("sample" %in% colnames(metadata)) {
    if (length(unique(metadata$sample)) != nrow(metadata)) {
      stop("Metadata sample column contains duplicate values")
    }
  } else {
    stop("Metadata must contain a 'sample' column")
  }
  
  # Set row names to sample names for easier indexing
  row.names(metadata) <- metadata$sample
  
  return(metadata)
}

# Function to prepare metadata for DESeq2
prepare_metadata_for_deseq <- function(metadata, count_colnames) {
  # Extract sample names from count data column names
  sample_pattern <- "^([^\\s]+)\\s.*$"
  sample_names <- gsub(sample_pattern, "\\1", count_colnames)
  
  # Ensure all samples in count data have metadata
  missing_samples <- sample_names[!sample_names %in% metadata$sample]
  
  if (length(missing_samples) > 0) {
    stop(paste("Metadata missing for samples:", paste(missing_samples, collapse=", ")))
  }
  
  # Subset metadata to only include samples in count data
  metadata_subset <- metadata[metadata$sample %in% sample_names, , drop = FALSE]
  
  # Ensure order matches count data
  metadata_subset <- metadata_subset[match(sample_names, metadata_subset$sample), , drop = FALSE]
  
  # Ensure column order matches count data
  if (!identical(metadata_subset$sample, sample_names)) {
    logger::warn("Metadata sample order doesn't match count data. Reordering.")
    metadata_subset <- metadata_subset[match(sample_names, metadata_subset$sample), , drop = FALSE]
  }
  
  # Convert character columns to factors for DESeq2
  for (col in colnames(metadata_subset)) {
    if (is.character(metadata_subset[[col]])) {
      metadata_subset[[col]] <- as.factor(metadata_subset[[col]])
    }
  }
  
  return(metadata_subset)
} 