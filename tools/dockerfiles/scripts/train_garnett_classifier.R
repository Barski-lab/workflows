#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)

suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(garnett))
suppressMessages(library(argparse))
suppressMessages(library(tidyverse))


load_training_data <- function(args) {
    tryCatch(
        expr = {
            training_data <<- UpdateSeuratObject(readRDS(args$rds))
        },
        error = function(e){
            load(args$rds, verbose=TRUE)
            training_data <<- UpdateSeuratObject(eval(as.name(args$objname)))
        }
    )
    return (training_data)
}


export_formatted_cell_markers_data <- function(data, location){
    file.create(location)
    for(i in 1:nrow(data)) {
        row <- data[i,]
        cat(paste0(">", row[1], "\n"), file=location, append=TRUE)
        cat(paste0("expressed: ", row[2], "\n\n"), file=location, append=TRUE)
    }
}


export_rds <- function(data, location){
    tryCatch(
        expr = {
            saveRDS(data, location)
            print(paste("Export data as RDS to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data as RDS to", location, sep=" "))
        }
    )
}


load_raw_cell_markers_data <- function (location, species, gene_sep=",") {
    raw_cell_markers_data <- read.table(
        location,
        sep="\t",
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    exp_columns <- paste0(
        "exp_",
        seq_len(
            max(length(strsplit(raw_cell_markers_data$geneSymbol, gene_sep)))
        )
    )
    formatted_cell_markers_data <- raw_cell_markers_data %>%
                                   filter(.$speciesType==species) %>%
                                   separate("geneSymbol", all_of(exp_columns), gene_sep) %>%
                                   gather("not_used", "Marker", all_of(exp_columns), na.rm=TRUE) %>%
                                   dplyr::select("cellName", "Marker") %>%
                                   dplyr::rename("Type"="cellName") %>%
                                   mutate("Marker"=str_trim(.$Marker)) %>%
                                   filter(.$Marker!="NA") %>%
                                   filter(.$Marker!="") %>%
                                   distinct() %>%
                                   group_by(Type) %>%
                                   mutate(Expressed=paste0(Marker, collapse=", ")) %>%
                                   dplyr::select("Type", "Expressed") %>%
                                   distinct()
    return (formatted_cell_markers_data)
}


get_monocle_data <- function (seurat_data, features, assay="RNA", matrix_slot="data") {
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    gene_metadata <- data.frame(
        "gene_short_name"=features,
        row.names=features,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    cell_metadata <- data.frame(
        seurat_data@meta.data,
        row.names=colnames(seurat_data),
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    expression_matrix <- GetAssayData(seurat_data, slot=matrix_slot)[features, ]
    monocle_data <- new_cell_data_set(
        expression_matrix,
        cell_metadata=cell_metadata,
        gene_metadata=gene_metadata
    )
    DefaultAssay(seurat_data) <- backup_assay
    return (monocle_data)
}


export_markers_plot <- function(data, rootname, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot_markers(data)))
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot_markers(data)))
                dev.off()
            }
            print(paste("Export markers plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export markers plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description="Trains Garnett classifier")
    parser$add_argument("--markers", help="Path to the Garnett marker genes file", type="character", required="True")
    parser$add_argument("--rds",     help="Path to the Seurat training data file in rds or Robj format", type="character", required="True")
    parser$add_argument("--objname", help="Object name to load Seurat training data from when provided through --rds file is in Robj format", type="character", default="tiss")
    parser$add_argument("--pdf",     help="Export plots in PDF. Default: false", action="store_true")
    parser$add_argument("--output",  help="Output prefix. Default: ./garnett", type="character", default="./garnett")
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}


args <- get_args()

print("Loading Seurat training data")
training_data <- load_training_data(args)

print(training_data)
print(head(training_data@meta.data))

print("Convereting Seurat training data to Monocle3")
monocle_data <- get_monocle_data(training_data, rownames(training_data))

# print("Checking markers")
# marker_check <- check_markers(
#     monocle_data,
#     args$markers,
#     db="none"
# )

# export_markers_plot(marker_check, paste(args$output, "markers", sep="_"), args$pdf)

# classifier <- train_cell_classifier(
#     cds=monocle_data,
#     marker_file=args$markers,
#     db="none"
# )

# export_rds(classifier,  paste(args$output, "classifier.rds", sep="_"))