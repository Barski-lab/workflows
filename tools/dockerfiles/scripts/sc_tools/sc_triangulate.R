#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(knitr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(stringr))
suppressMessages(library(modules))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))

## ----
export_all_plots <- function(seurat_data, args){
    Idents(seurat_data) <- "new.ident"                                                               # safety measure
    for (reduction in c("rnaumap", "atacumap", "wnnumap")){
        if (!(reduction %in% names(seurat_data@reductions))) {next}                                  # skip missing reductions
        graphics$dim_plot(
            data=seurat_data,
            reduction=reduction,
            plot_title=paste0(
                "Cells UMAP with integrated labels (",
                reduction, " dim. reduction)"
            ),
            legend_title="Integrated label",
            group_by=paste("custom", args$target, "pruned", sep="_"),
            label=FALSE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_tril_rd", reduction, sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction=reduction,
            plot_title=paste0(
                "Cells UMAP with winning annotations (",
                reduction, " dim. reduction)"
            ),
            legend_title="Winning annotation",
            group_by=paste("custom", args$target, "final_annotation", sep="_"),
            label=FALSE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_tria_rd", reduction, sep="_"),
            pdf=args$pdf
        )
        graphics$feature_plot(
            data=seurat_data,
            features=paste("custom", args$target, "confidence", sep="_"),,
            labels=NULL,                                                                    # we already have dataset names on the right side of each plot
            from_meta=TRUE,
            reduction=reduction,
            plot_title=paste0(
                "Cells UMAP with integration confidence scores (",
                reduction, " dim. reduction)"
            ),
            label=FALSE,
            alpha=1,
            order=TRUE,
            gradient_colors=c("yellow", "#662404"),
            color_scales=c(0, 1),
            color_limits=c(0, 1),
            theme=args$theme,
            rootname=paste(args$output, "umap_tric_rd", reduction, sep="_"),
            pdf=args$pdf
        )
    }
    gc(verbose=FALSE)
}

## ----
get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell Label Integration Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "genes expression and/or chromatin accessibility information stored in the RNA",
            "and/or ATAC assays correspondingly. Additionally, 'rnaumap', and/or 'atacumap',",
            "and/or 'wnnumap' dimensionality reductions should be present."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--barcodes",
        help=paste(
            "Path to the TSV/CSV file to optionally prefilter and extend Seurat object",
            "metadata be selected barcodes. First column should be named as 'barcode'.",
            "If file includes any other columns they will be added to the Seurat object",
            "metadata ovewriting the existing ones if those are present.",
            "Default: all cells used, no extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--source",
        help=paste(
            "Columns from the metadata of the loaded Seurat object to select",
            "conflicting cells annotations."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--target",
        help=paste(
            "Suffix to be used as part of the columns names to save label",
            "integration result. Default: sctri"
        ),
        type="character", default="sctri"
    )
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--verbose",
        help="Print debug information. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--h5seurat",
        help="Save Seurat data to h5seurat file. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--h5ad",
        help="Save raw counts from the RNA and/or ATAC assay(s) to h5ad file(s). Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--loupe",
        help=paste(
            "Save raw counts from the RNA assay to Loupe file.",
            "By enabling this feature you accept the End-User",
            "License Agreement available at https://10xgen.com/EULA.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--cbbuild",
        help="Export results to UCSC Cell Browser. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./sc",
        type="character", default="./sc"
    )
    parser$add_argument(
        "--theme",
        help=paste(
            "Color theme for all generated plots.",
            "Default: classic"
        ),
        type="character", default="classic",
        choices=c("gray", "bw", "linedraw", "light", "dark", "minimal", "classic", "void")
    )
    parser$add_argument(
        "--cpus",
        help="Number of cores/cpus to use. Default: 1",
        type="integer", default=1
    )
    parser$add_argument(
        "--memory",
        help=paste(
            "Maximum memory in GB allowed to be shared between the workers",
            "when using multiple --cpus.",
            "Default: 32"
        ),
        type="integer", default=32
    )
    parser$add_argument(
        "--seed",
        help="Seed number for random values. Default: 42",
        type="integer", default=42
    )
    args <- parser$parse_args(str_subset(commandArgs(trailingOnly=TRUE), "\\.R$", negate=TRUE))  # to exclude itself when executed from the sc_report_wrapper.R
    print(args)
    return (args)
}

## ----
args <- get_args()
prod$parallel(args)

## ----
print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
debug$print_info(seurat_data, args)

## ----
if (!any(c("rnaumap", "atacumap", "wnnumap") %in% names(seurat_data@reductions))){
    print(
        paste(
            "Loaded Seurat object includes neither of the required reductions:",
            "'rnaumap', and/or 'atacumap', and/or 'wnnumap'.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (!any(c("RNA", "ATAC") %in% names(seurat_data@assays))){
    print(
        paste(
            "Loaded Seurat object includes neither of the required assays:",
            "'RNA' and/or 'ATAC'.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)                  # sets identities to new.ident
    debug$print_info(seurat_data, args)
}

## ----
print("Adjusting --source parameters to include only columns present in Seurat object metadata")    # need to do it after filtering by --barcodes
args$source <- unique(args$source)
args$source <- args$source[args$source %in% colnames(seurat_data@meta.data)]

## ----
seurat_data <- analyses$integrate_labels(
    seurat_data=seurat_data,
    source_columns=args$source,
    args=args
)
print(head(seurat_data))

## ----
export_all_plots(seurat_data=seurat_data, args=args)

## ----
if(args$cbbuild){
    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/rna", sep=""),
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/atac", sep=""),
        )
    } else if ("RNA" %in% names(seurat_data@assays)){
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            rootname=paste(args$output, "_cellbrowser", sep=""),
        )
    } else {
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            rootname=paste(args$output, "_cellbrowser", sep=""),
        )
    }
}

## ----
io$export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))

## ----
if(args$h5seurat){
    io$export_h5seurat(seurat_data, paste(args$output, "_data.h5seurat", sep=""))
}

## ----
if(args$h5ad && ("RNA" %in% names(seurat_data@assays))){
    io$export_h5ad(
        data=seurat_data,
        location=paste(args$output, "_rna_counts.h5ad", sep=""),
        assay="RNA",
        slot="counts"
    )
}

## ----
if(args$h5ad && ("ATAC" %in% names(seurat_data@assays))){
    io$export_h5ad(
        data=seurat_data,
        location=paste(args$output, "_atac_counts.h5ad", sep=""),
        assay="ATAC",
        slot="counts"
    )
}

## ----
if(args$loupe && ("RNA" %in% names(seurat_data@assays))){
    ucsc$export_loupe(
        seurat_data=seurat_data,
        assay="RNA",
        rootname=paste0(args$output, "_rna_counts")
    )
}