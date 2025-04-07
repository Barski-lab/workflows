#!/usr/bin/env Rscript
#
# Main entry point for DESeq2 LRT Step 2 Analysis
#

# Display startup message
message("Starting DESeq2 LRT Step 2 Analysis")
message("Working directory:", getwd())

# Docker paths
workflow_file <- "/usr/local/bin/functions/deseq2_lrt_step_2/workflow.R"
message("Looking for workflow file at:", workflow_file)

# Source the workflow file
if (file.exists(workflow_file)) {
  message("Found workflow file. Sourcing:", workflow_file)
  source(workflow_file)
  
  # Initialize the environment
  initialize_environment()
  
  # Run the workflow with memory management
  results <- main_with_memory_management()
  
  # Verify outputs
  verify_file <- "/usr/local/bin/verify_outputs.R"
  if (file.exists(verify_file)) {
    message("Verifying outputs with:", verify_file)
    source(verify_file)
    
    # Get the output prefix from command-line arguments
    args <- get_args()
    output_prefix <- if (!is.null(args$output_prefix)) args$output_prefix else 
                     if (!is.null(args$output)) args$output else "deseq-lrt-step-2"
    
    # Verify all required outputs were created
    message("Verifying all required outputs were created...")
    verify_workflow_outputs("lrt_step2", output_prefix, fail_on_missing = FALSE)
  } else {
    message("Warning: Verification file not found, skipping output verification")
  }
} else {
  message("ERROR: Workflow file not found at", workflow_file)
  message("Checking alternative locations...")
  
  # Try relative path
  relative_path <- "functions/deseq2_lrt_step_2/workflow.R"
  if (file.exists(relative_path)) {
    message("Found workflow file at relative path:", relative_path)
    source(relative_path)
    initialize_environment()
    results <- main_with_memory_management()
    
    # Verify outputs
    verify_path <- "verify_outputs.R"
    if (file.exists(verify_path)) {
      source(verify_path)
      
      # Get the output prefix from command-line arguments
      args <- get_args()
      output_prefix <- if (!is.null(args$output_prefix)) args$output_prefix else 
                       if (!is.null(args$output)) args$output else "deseq-lrt-step-2"
      
      # Verify all required outputs were created
      message("Verifying all required outputs were created...")
      verify_workflow_outputs("lrt_step2", output_prefix, fail_on_missing = FALSE)
    }
  } else {
    # Last resort - try to find it
    message("Attempting to locate workflow.R file...")
    system("find /usr/local -name workflow.R | grep deseq2_lrt_step_2/", intern = FALSE)
    stop("Could not find workflow.R file. Please verify your installation.")
  }
}

message("DESeq2 LRT Step 2 analysis completed.")
