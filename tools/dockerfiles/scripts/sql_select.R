#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)


suppressMessages(library(argparse))
suppressMessages(library(sqldf))


###########################################################################################
# v0.0.1
# Input CSV/TSV file should have header
# Value provided through --sql parameter will be appended to "SELECT * FROM raw_data WHERE"
# Column names in --sql query are case insensitive
# Output is always TSV file with header
###########################################################################################


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = ","
    if (ext == "tsv"){
        separator = "\t"
    }
    return (separator)
}

load_and_run_query <- function(location, sql, columns) {
    print(paste("Loading raw data from", location))
    raw_data <- read.table(location, sep=get_file_type(location), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
    print(head(raw_data))
    full_sql_query <- paste("SELECT", columns, "FROM raw_data WHERE", sql)
    print(paste("Applying filter", full_sql_query))
    filtered_data <- sqldf(full_sql_query)
    return (filtered_data)
}   


assert_args <- function(args){
    print("Check input parameters")
    if(is.null(args$output)){
        args$output <- paste(head(unlist(strsplit(basename(args$raw), ".", fixed = TRUE)), 1), "filtered.tsv", sep="_")
    }
    return (args)
}


parser <- ArgumentParser(description='Filter CSV/TSV file based on the provided SQL query')
parser$add_argument("--raw",     help='Input CSV/TSV file', type="character", required="True")
parser$add_argument("--sql",     help='SQL query parameters that will be appended to `SELECT *` statement', type="character", required="True")
parser$add_argument("--columns", help="Comma-separated list of column to be print in the output. Default: all", type="character", default="*")
parser$add_argument("--header",  help='Print header in the output. Default: false', action='store_true')
parser$add_argument("--output",  help='Filename for filtered output TSV file', type="character")
args <- assert_args(parser$parse_args(gsub("'", "\"", commandArgs(trailingOnly = TRUE))))


filtered_data <- load_and_run_query(args$raw, args$sql, args$columns)
print(paste("Exporting filtered data to", args$output))
print(head(filtered_data))
write.table(filtered_data,
            file= args$output,
            sep="\t",
            row.names=FALSE,
            col.names=args$header,
            quote=FALSE
)