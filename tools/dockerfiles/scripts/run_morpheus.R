#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(morpheus))
suppressMessages(library(argparse))
suppressMessages(library(htmlwidgets))


get_args <- function(){
    parser <- ArgumentParser(description="Morpheus heatmap from GCT file")
    parser$add_argument(
        "--gct",
        help=paste(
            "Path to the input GCT file."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--output",
        help=paste(
            "Output prefix for generated files"
        ),
        type="character", default="./heatmap"
    )
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}

# Parse arguments
args <- get_args()

print("Used parameters")
print(args)

print(paste("Loading GCT data from", args$gct))
gct_data <- read.gct(args$gct)

is_all_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}
morpheus_html <- morpheus(
    x=gct_data$data,
    rowAnnotations=if(nrow(gct_data$rowAnnotations) == 0) NULL else gct_data$rowAnnotations %>% dplyr::mutate_if(is_all_numeric, as.numeric),
    columnAnnotations=if(nrow(gct_data$columnAnnotations) == 0) NULL else gct_data$columnAnnotations
)

html_location <- paste0(args$output, ".html")
print(paste("Saving heatmap to", html_location))
htmlwidgets::saveWidget(
    morpheus_html,
    file=html_location
)