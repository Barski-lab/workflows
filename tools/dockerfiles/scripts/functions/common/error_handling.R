#!/usr/bin/env Rscript
#
# Common error handling functions for DESeq2 analysis
#

#' Handle errors with traceback and a user-friendly message
#'
#' @param expr The expression to evaluate
#' @param message Optional custom error message
#' @param exit_on_error Whether to exit the script on error
#' @return The result of the expression if successful
#' @export
with_error_handling <- function(expr, message = NULL, exit_on_error = TRUE) {
  tryCatch(
    expr,
    error = function(e) {
      error_msg <- ifelse(is.null(message), e$message, paste0(message, ": ", e$message))
      log_error(error_msg)
      print(rlang::last_trace())
      if (exit_on_error) {
        quit(save = "no", status = 1, runLast = FALSE)
      }
      return(NULL)
    }
  )
}

#' Check if required packages are installed
#'
#' @param packages Vector of package names to check
#' @param exit_on_missing Whether to exit if packages are missing
#' @return Boolean indicating if all packages are installed
#' @export
check_required_packages <- function(packages, exit_on_missing = TRUE) {
  missing_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    log_error(paste0("Required packages are missing: ", paste(missing_packages, collapse = ", ")))
    
    if (exit_on_missing) {
      log_error("Please install missing packages and try again.")
      quit(save = "no", status = 1, runLast = FALSE)
    }
    
    return(FALSE)
  }
  
  return(TRUE)
}

# --- Error handling functions ---

# Detailed error reporting function
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