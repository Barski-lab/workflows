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


export_all_dimensionality_plots <- function(seurat_data, args) {
    Idents(seurat_data) <- "new.ident"                                                                                         # safety measure
    selected_features=c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal", "frip", "blacklist_fraction")
    selected_labels=c("UMI", "Peaks", "TSS enrichment score", "Nucleosome signal", "FRiP", "Bl. regions")

    graphics$corr_plot(
        data=seurat_data,
        reduction="qclsi",
        highlight_dims=args$dimensions,
        qc_columns=selected_features,
        qc_labels=selected_labels,
        plot_title="Correlation plots between QC metrics and cells LSI dimensions",
        combine_guides="collect",
        theme=args$theme,
        rootname=paste(args$output, "qc_dim_corr", sep="_"),
        pdf=args$pdf
    )
    graphics$feature_plot(
        data=seurat_data,
        features=selected_features,
        labels=selected_labels,
        from_meta=TRUE,
        reduction="atacumap",
        plot_title="QC metrics on cells UMAP",
        label=FALSE,
        alpha=0.4,
        max_cutoff="q99",                                                                   # to prevent outlier cells to distort coloring
        combine_guides="keep",
        theme=args$theme,
        rootname=paste(args$output, "umap_qc_mtrcs", sep="_"),
        pdf=args$pdf
    )

    graphics$dim_plot(
        data=seurat_data,
        reduction="atacumap",
        plot_title="Cells UMAP",
        legend_title="Dataset",
        group_by="new.ident",
        label=FALSE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap", sep="_"),
        pdf=args$pdf
    )

    if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
        graphics$dim_plot(
            data=seurat_data,
            reduction="atacumap",
            plot_title="Split by dataset cells UMAP",
            legend_title="Dataset",
            group_by="new.ident",
            split_by="new.ident",
            label=FALSE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_spl_idnt", sep="_"),
            pdf=args$pdf
        )
    }

    if (
        all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
        length(unique(as.vector(as.character(seurat_data@meta.data$condition)))) > 1
    ){
        graphics$dim_plot(
            data=seurat_data,
            reduction="atacumap",
            plot_title="Split by grouping condition cells UMAP",
            legend_title="Dataset",
            group_by="new.ident",
            split_by="condition",
            label=FALSE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_spl_cnd", sep="_"),
            pdf=args$pdf
        )
    }

    graphics$dim_plot(
        data=seurat_data,
        reduction="atacumap",
        plot_title="Split by the UMI per cell counts cells UMAP",
        legend_title="Dataset",
        group_by="new.ident",
        split_by="quartile_nCount_ATAC",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_umi", sep="_"),
        pdf=args$pdf
    )

    graphics$dim_plot(
        data=seurat_data,
        reduction="atacumap",
        plot_title="Split by the peaks per cell counts cells UMAP",
        legend_title="Dataset",
        group_by="new.ident",
        split_by="quartile_nFeature_ATAC",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_peak", sep="_"),
        pdf=args$pdf
    )

    graphics$dim_plot(
        data=seurat_data,
        reduction="atacumap",
        plot_title="Split by the TSS enrichment score cells UMAP",
        legend_title="Dataset",
        group_by="new.ident",
        split_by="quartile_TSS.enrichment",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_tss", sep="_"),
        pdf=args$pdf
    )

    graphics$dim_plot(
        data=seurat_data,
        reduction="atacumap",
        plot_title="Split by the nucleosome signal cells UMAP",
        legend_title="Dataset",
        group_by="new.ident",
        split_by="quartile_nucleosome_signal",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_ncls", sep="_"),
        pdf=args$pdf
    )

    graphics$dim_plot(
        data=seurat_data,
        reduction="atacumap",
        plot_title="Split by the FRiP cells UMAP",
        legend_title="Dataset",
        group_by="new.ident",
        split_by="quartile_frip",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_frip", sep="_"),
        pdf=args$pdf
    )

    graphics$dim_plot(
        data=seurat_data,
        reduction="atacumap",
        plot_title="Split by the genomic blacklist regions fraction cells UMAP",
        legend_title="Dataset",
        group_by="new.ident",
        split_by="quartile_blacklist_fraction",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_blck", sep="_"),
        pdf=args$pdf
    )

}


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell ATAC-Seq Dimensionality Reduction Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "chromatin accessibility information stored in the ATAC assay."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to optionally extend Seurat object metadata with",
            "categorical values using samples identities. First column - 'library_id'",
            "should correspond to all unique values from the 'new.ident' column of the",
            "loaded Seurat object. If any of the provided in this file columns are already",
            "present in the Seurat object metadata, they will be overwritten. When combined",
            "with --barcodes parameter, first the metadata will be extended, then barcode",
            "filtering will be applied.",
            "Default: no extra metadata is added"
        ),
        type="character"
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
        "--norm",
        help=paste(
            "TF-IDF normalization method applied to chromatin accessibility counts.",
            "log-tfidf - Stuart & Butler et al. 2019,",
            "tf-logidf - Cusanovich & Hill et al. 2018,",
            "logtf-logidf - Andrew Hill,",
            "idf - 10x Genomics,",
            "Default: log-tfidf"
        ),
        type="character",
        default="log-tfidf",
        choices=c("log-tfidf", "tf-logidf", "logtf-logidf", "idf")
    )
    parser$add_argument(
        "--ntgr",
        help=paste(
            "Integration method used for joint analysis of multiple datasets. Automatically",
            "set to 'none' if loaded Suerat object includes only one dataset.",
            "Default: signac"
        ),
        type="character",
        default="signac",
        choices=c("signac", "harmony", "none")
    )
    parser$add_argument(
        "--ntgrby",
        help=paste(
            "Column(s) from the Seurat object metadata to define the variable(s) that should",
            "be integrated out when running multiple datasets integration with harmony. May",
            "include columns from the extra metadata added with --metadata parameter. Ignored",
            "if --ntgr is not set to harmony.",
            "Default: new.ident"
        ),
        type="character", default=c("new.ident"), nargs="*"
    )
    parser$add_argument(
        "--minvarpeaks",
        help=paste(
            "Minimum percentile for identifying the top most common peaks as highly variable.",
            "For example, setting to 5 will use the the top 95 percent most common among all cells",
            "peaks as highly variable. These peaks are used for datasets integration, scaling",
            "and dimensionality reduction.",
            "Default: 0 (use all available peaks)"
        ),
        type="integer", default=0
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use for datasets integration and UMAP projection (from 2 to 50).",
            "If single value N is provided, use from 2 to N LSI components. If multiple values are",
            "provided, subset to only selected LSI components. In combination with --ntgr set to",
            "harmony, multiple values will result in using all dimensions starting from 1(!) to",
            "the max of the provided values. Default: from 2 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--uspread",
        help=paste(
            "The effective scale of embedded points on UMAP. In combination with '--mindist'",
            "it determines how clustered/clumped the embedded points are.",
            "Default: 1"
        ),
        type="double", default=1
    )
    parser$add_argument(
        "--umindist",
        help=paste(
            "Controls how tightly the embedding is allowed compress points together on UMAP.",
            "Larger values ensure embedded points are moreevenly distributed, while smaller",
            "values allow the algorithm to optimise more accurately with regard to local structure.",
            "Sensible values are in the range 0.001 to 0.5.",
            "Default:  0.3"
        ),
        type="double", default=0.3
    )
    parser$add_argument(
        "--uneighbors",
        help=paste(
            "Determines the number of neighboring points used in UMAP. Larger values will result",
            "in more global structure being preserved at the loss of detailed local structure.",
            "In general this parameter should often be in the range 5 to 50.",
            "Default: 30"
        ),
        type="integer", default=30
    )
    parser$add_argument(
        "--umetric",
        help=paste(
            "The metric to use to compute distances in high dimensional space for UMAP.",
            "Default: cosine"
        ),
        type="character", default="cosine",
        choices=c(
            "euclidean", "manhattan", "chebyshev", "minkowski", "canberra", "braycurtis",
            "mahalanobis", "wminkowski", "seuclidean", "cosine", "correlation", "haversine",
            "hamming", "jaccard", "dice", "russelrao", "kulsinski", "ll_dirichlet", "hellinger",
            "rogerstanimoto", "sokalmichener", "sokalsneath", "yule"
        )
    )
    # The default method for RunUMAP has changed from calling Python UMAP via reticulate to
    # the R-native UWOT using the cosine metric. To use Python UMAP via reticulate, set
    # umap.method to 'umap-learn' and metric to 'correlation'
    parser$add_argument(
        "--umethod",
        help=paste(
            "UMAP implementation to run. If set to 'umap-learn' use --umetric 'correlation'",
            "Default: uwot"
        ),
        type="character", default="uwot",
        choices=c("uwot", "uwot-learn", "umap-learn")
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
    args$dimensions <- c(2:args$dimensions[1])                                              # skipping the first LSI component
    print(paste("--dimensions was adjusted to", paste(args$dimensions, collapse=", ")))
}
args$minvarpeaks <- paste0("q", args$minvarpeaks)                                         # need to have it in a form of "qN", for example "q0"

print(
    paste(
        "Setting parallelization to", args$cpus, "cores, and", args$memory,
        "GB of memory allowed to be shared between the processes"
    )
)
prod$parallel(args)

print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
print("Setting default assay to ATAC")
DefaultAssay(seurat_data) <- "ATAC"
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

if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)
}
debug$print_info(seurat_data, args)

print("Running ATAC analysis")
seurat_data <- analyses$atac_analyze(seurat_data, args)                   # adds "atac_lsi" and "atacumap" reductions
seurat_data <- filter$collapse_fragments_list(seurat_data)                # collapse repetitive fragments if ATAC assay was splitted when running integration
debug$print_info(seurat_data, args)

print("Quantifying QC metrics")
seurat_data <- qc$quartile_qc_metrics(
    seurat_data=seurat_data,
    features=c(
        "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment",
        "nucleosome_signal", "frip", "blacklist_fraction"
    ),
    prefix="quartile"                                                       # we will use this prefix in ucsc$export_cellbrowser function
)
debug$print_info(seurat_data, args)

export_all_dimensionality_plots(
    seurat_data=seurat_data,
    args=args
)

if ("qclsi" %in% names(seurat_data@reductions)){                            # we only needed it for qc correlation plots
    print("Removing qclsi reduction")
    seurat_data[["qclsi"]] <- NULL
    debug$print_info(seurat_data, args)
}

if(args$cbbuild){
    print("Exporting ATAC assay to UCSC Cellbrowser")
    print("Reordering reductions to have atacumap on the first place")                      # will be shown first in UCSC Cellbrowser
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c("atacumap", reduc_names[reduc_names!="atacumap"])              # atacumap will be added by this time
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]
    debug$print_info(seurat_data, args)
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="ATAC",
        slot="counts",
        short_label="ATAC",
        palette_colors=graphics$D40_COLORS,                                                 # to have colors correspond to the plots
        rootname=paste(args$output, "_cellbrowser", sep=""),
    )
}

DefaultAssay(seurat_data) <- "ATAC"
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