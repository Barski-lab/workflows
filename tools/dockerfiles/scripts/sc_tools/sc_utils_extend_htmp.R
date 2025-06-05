#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(tools))
suppressMessages(library(modules))
suppressMessages(library(tidyverse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))


args = commandArgs(trailingOnly=TRUE)
if(length(args) != 2){
    print("Not enough input parameters. Exiting without throwing the error.")
    quit(save="no", status=0, runLast=FALSE)
}
gct_location <- args[1]
metadata_location <- args[2]
gct_data <- morpheus::read.gct(gct_location)
metadata <- read.table(
                metadata_location,
                sep="\t",
                header=TRUE,
                check.names=FALSE,
                stringsAsFactors=FALSE,
                quote=""
            ) %>%
            mutate(label=paste(chr, paste(start, end, sep="-"), sep=":"))    # we need to overwrite the labels to include only chr:start-end 
row_metadata <- gct_data$rowAnnotations %>%
                select("regions_group") %>%                                  # this is the only column we don't have in the metadata
                rownames_to_column("label") %>%
                left_join(metadata, by="label") %>%
                remove_rownames() %>%
                column_to_rownames("label") %>%
                select(c("gene_id", "region", "log2FoldChange", "padj", "regions_group")) %>%
                rename(category = regions_group)
io$export_gct(
    counts_mat=gct_data$data,
    row_metadata=row_metadata,
    col_metadata=gct_data$columnAnnotations,
    location=basename(gct_location)
)
score_limits <- stats::quantile(
    gct_data$data,
    c(0.01, 0.98),
    na.rm=TRUE, names=FALSE
)
graphics$morpheus_html_heatmap(
    gct_location=basename(gct_location),
    rootname=file_path_sans_ext(basename(gct_location)),
    color_scheme=list(
        scalingMode="fixed",
        stepped=FALSE,
        values=as.list(score_limits),
        colors=c("white", "darkblue")
    )
)
io$export_data(
    row_metadata %>% tibble::rownames_to_column(var="feature"),
    paste0(file_path_sans_ext(basename(gct_location)), ".tsv")
)