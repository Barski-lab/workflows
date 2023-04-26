#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(argparse))


HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(filter <- modules::use(file.path(HERE, "modules/filter.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


export_plots <- function(seurat_data, da_cells, da_thresholds, args) {
    Idents(seurat_data) <- "new.ident"                                                                       # safety measure

    graphics$daseq_permutations(
        data=da_cells,
        plot_title=paste(
            "DA scores random permutations plot for",
            args$second, "vs", args$first, "comparison"
        ),
        x_label="DA score",
        y_label="Ranked by DA score cells",
        y_intercepts=if(is.null(args$ranges)) round(da_thresholds, 2) else round(c(da_thresholds, args$ranges), 2),
        palette_colors=if(is.null(args$ranges)) c("black", "black") else c("black", "black", "red", "red"),
        theme=args$theme,
        rootname=paste(args$output, "da_perm", sep="_"),
        pdf=args$pdf
    )

    for (reduction in c("rnaumap", "atacumap", "wnnumap")){
        if (!(reduction %in% names(seurat_data@reductions))) {next}                                          # skip missing reductions
        for (i in 1:length(args$resolution)) {
            current_resolution <- args$resolution[i]
            current_cluster <- paste0("da_", args$second, "_vs_", args$first, "_res.", current_resolution)
            graphics$dim_plot(
                data=seurat_data,
                reduction=reduction,
                cells=rownames(seurat_data@meta.data[seurat_data@meta.data[[current_cluster]] != 0, ]),      # remove cluster 0 as it includes non-DA cells
                plot_title=paste0(
                    "Clustered DA cells subpopulations UMAP (", reduction,
                    " dim. reduction). Resolution ", current_resolution
                ),
                legend_title="DA cluster",
                group_by=current_cluster,
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_rd", reduction, "res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$dim_plot(
                data=seurat_data,
                reduction=reduction,
                plot_title=paste0(
                    "Split by grouping condition clustered DA cells subpopulations UMAP (",
                    reduction, " dim. reduction). Resolution ", current_resolution
                ),
                legend_title="DA cluster",
                group_by=current_cluster,
                split_by=args$splitby,
                label=TRUE,
                label_color="black",
                palette_colors=c("lightgrey", graphics$D40_COLORS),         # adding grey for cluster 0 as we can't use cells and split_by together
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_cnd_rd", reduction, "res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
        graphics$feature_plot(
            data=seurat_data,
            features=paste("custom", "da_score", args$second, "vs", args$first, sep="_"),
            labels=NULL,                                                                    # we already have dataset names on the right side of each plot
            from_meta=TRUE,
            reduction=reduction,
            split_by="new.ident",
            plot_title=paste(
                "Split by dataset cells UMAP with DA scores for",
                args$second, "vs", args$first, "comparison"
            ),
            label=FALSE,
            alpha=0.5,
            order=FALSE,                                                                    # otherwise white will be on top of blue
            gradient_colors=c("blue", "white", "red"),
            color_scales=c(-1, 1),
            color_limits=c(-1, 1),
            combine_guides="collect",
            theme=args$theme,
            rootname=paste(args$output, "umap_spl_idnt_rd", reduction, "da_scr", sep="_"),
            pdf=args$pdf
        )
    }
}


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell Differential Abundance Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include genes",
            "expression information stored in the RNA assay and selected with the --reduction",
            "parameter dimensionality reduction. Additionally, 'rnaumap', and/or 'atacumap',",
            "and/or 'wnnumap' dimensionality reductions should be present."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--reduction",
        help=paste(
            "Dimensionality reduction to be used for DA analysis.",
            "Default: pca"
        ),
        type="character", default="pca"
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use when running DA analysis (from 1 to 50).",
            "If single value N is provided, use from 1 to N PCs. If multiple",
            "values are provided, subset to only selected PCs.",
            "Default: from 1 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--knn",
        help=paste(
            "Array of k values for kNN graph construction when calculating the",
            "score vector for each cell to represent the DA behavior in the",
            "neighborhood.",
            "Default: calculated based on the cells number"
        ),
        type="integer", nargs="*"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to optionally extend Seurat object metadata with",
            "categorical values using samples identities. First column - 'library_id'",
            "should correspond to all unique values from the 'new.ident' column of the",
            "loaded Seurat object. If any of the provided in this file columns are already",
            "present in the Seurat object metadata, they will be overwritten.",
            "Default: no extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to split cells into two groups",
            "to run --second vs --first DA analysis. May include columns from the",
            "extra metadata added with --metadata parameter."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--first",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby to define",
            "the first group of cells for DA analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--second",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby to define",
            "the second group of cells for DA analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--resolution",
        help=paste(
            "Clustering resolution applied to DA cells to identify DA cells populations.",
            "Can be set as an array.",
            "Default: 0.01, 0.03, 0.05"
        ),
        type="double", default=c(0.01, 0.03, 0.05), nargs="*"
    )
    parser$add_argument(
        "--ranges",
        help=paste(
            "DA scores ranges for to filter out not significant cells.",
            "Default: calculated based on the permutation test"
        ),
        type="double", nargs=2
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
        help="Save Seurat data to h5ad file. Default: false",
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
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}

args <- get_args()
print("Input parameters")
print(args)
if (length(args$dimensions) == 1) {
    print("Adjusting --dimensions parameter as only a single value was provided")
    args$dimensions <- c(1:args$dimensions[1])
    print(paste("--dimensions was adjusted to", paste(args$dimensions, collapse=", ")))
}
if(!is.null(args$ranges)){
    args$ranges <- sort(args$ranges, decreasing=TRUE)   # for consistency with other parts of the program
}

print(
    paste(
        "Setting parallelization to", args$cpus, "cores, and", args$memory,
        "GB of memory allowed to be shared between the processes"
    )
)
prod$parallel(args)

print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)

if (!(args$reduction %in% names(seurat_data@reductions))){
    print(
        paste(
            "Loaded Seurat object doesn't include selected",
            "reduction", args$reduction, "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

if (!("RNA" %in% names(seurat_data@assays))){
    print(
        paste(
            "Loaded Seurat object doesn't include required RNA assay.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

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

print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"
debug$print_info(seurat_data, args)

if (!is.null(args$metadata)){
    print("Extending Seurat object with the extra metadata fields")
    seurat_data <- io$extend_metadata(
        seurat_data=seurat_data,
        location=args$metadata,
        seurat_ref_column="new.ident",
        meta_ref_column="library_id"
    )
    debug$print_info(seurat_data, args)
}

print("Filtering Seurat object to include only selected groups of cells")
seurat_data <- io$apply_metadata_filters(seurat_data, args$splitby, c(args$first, args$second))
debug$print_info(seurat_data, args)

print(
    paste0(
        "Running ", args$second, " vs ", args$first, " differential abundance ",
        "analysis for datasets split by ", args$splitby, " using ",
        paste(args$dimensions, collapse=","), " dimensions from the ",
        args$reduction, " dimensionality reduction.",
        ifelse(
            !is.null(args$knn),
            paste0(" User provided array of k values for kNN graph construction: ", paste(args$knn, collapse=","), "."),
            ""
        ),
        ifelse(
            !is.null(args$ranges),
            paste0(" User provided thresholds for DA scores: ", paste(args$ranges, collapse=",")),
            ""
        )
    )
)

da_results <- analyses$da_analyze(seurat_data, args)                         # will add new metadata column with DA predictions

seurat_data <- da_results$seurat_data                                        # for easy access
da_cells <- da_results$da_cells
da_thresholds <- da_results$thresholds                                       # thresholds identified by DASeq by permutation test
rm(da_results)                                                               # remove unused data
gc(verbose=FALSE)

# need to adjust args$ranges (if they were provided) the same way as DASeq does it internally
if(!is.null(args$ranges)){
    if(da_thresholds[1] > args$ranges[1]){
        print(paste("Adjusting maximum value of user provided DA scores ranges to", da_thresholds[1]))
        args$ranges[1] <- da_thresholds[1]
    }
    if(da_thresholds[2] < args$ranges[2]){
        paste(paste("Adjusting minimum value of user provided DA scores ranges to", da_thresholds[2]))
        args$ranges[2] <- da_thresholds[2]
    }
}

export_plots(seurat_data, da_cells, da_thresholds, args)

if(args$cbbuild){
    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        print("Exporting RNA and ATAC assays to UCSC Cellbrowser jointly")
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            features=args$genes,                                             # can be NULL
            is_nested=TRUE,
            label_field=paste0("Clustering (da_", args$second, "_vs_", args$first, " ", args$resolution[1], ")"),                   # always use only the first resolution
            rootname=paste(args$output, "_cellbrowser/rna", sep=""),
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            features=nearest_peaks,                                          # use nearest to the genes if interest peaks
            is_nested=TRUE,
            label_field=paste0("Clustering (da_", args$second, "_vs_", args$first, " ", args$resolution[1], ")"),                   # always use only the first resolution
            rootname=paste(args$output, "_cellbrowser/atac", sep=""),
        )
    } else {
        print("Exporting RNA assay to UCSC Cellbrowser")
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            features=args$genes,                                             # can be NULL
            label_field=paste0("Clustering (da_", args$second, "_vs_", args$first, " ", args$resolution[1], ")"),                   # always use only the first resolution
            rootname=paste(args$output, "_cellbrowser", sep=""),
        )
    }
}

print("Exporting results to RDS file")
io$export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))
if(args$h5seurat){
    print("Exporting results to h5seurat file")
    io$export_h5seurat(seurat_data, paste(args$output, "_data.h5seurat", sep=""))
}

if(args$h5ad){
    print("Exporting results to h5ad file")
    io$export_h5ad(seurat_data, paste(args$output, "_data.h5ad", sep=""))
}