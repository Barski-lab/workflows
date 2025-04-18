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
suppressMessages(filter <- modules::use(file.path(HERE, "modules/filter.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))
suppressMessages(logger <- modules::use(file.path(HERE, "modules/logger.R")))

## ----
export_all_dimensionality_plots <- function(seurat_data, args) {
    Idents(seurat_data) <- "new.ident"                                                                                         # safety measure
    selected_features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "S.Score", "G2M.Score")
    selected_labels=c("RNA reads", "Genes", "Mitochondrial %", "Novelty score", "S score", "G2M score")

    graphics$elbow_plot(
        data=seurat_data,
        x_intercept=max(args$dimensions),                  # need to take the max because dimensions was already adjusted to array of values
        reduction="pca",                                   # this PCA is always after intergration
        plot_title="Elbow plot",
        theme=args$theme,
        rootname=paste(args$output, "elbow", sep="_"),
        pdf=args$pdf
    )
    graphics$corr_plot(
        data=seurat_data,
        reduction="pca",                                   # this PCA is always after intergration
        highlight_dims=args$dimensions,
        qc_columns=selected_features,
        qc_labels=selected_labels,
        plot_title="Correlation between QC metrics and principal components",
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
        reduction="rnaumap",
        plot_title="UMAP, QC metrics",
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
        reduction="rnaumap",
        plot_title="UMAP, colored by dataset",
        legend_title="Dataset",
        group_by="new.ident",
        label=FALSE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap", sep="_"),
        pdf=args$pdf
    )
    if ("Phase" %in% colnames(seurat_data@meta.data)){
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="UMAP, colored by dataset, split by cell cycle phase",
            legend_title="Dataset",
            group_by="new.ident",
            split_by="Phase",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_spl_ph", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="ccpca",
            plot_title="PCA, colored by cell cycle phase",
            legend_title="Phase",
            group_by="Phase",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "ccpca", sep="_"),
            pdf=args$pdf
        )
    }
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title="UMAP, colored by dataset, split by mitochondrial percentage",
        legend_title="Dataset",
        group_by="new.ident",
        split_by="quartile_mito_percentage",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_mito", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title="UMAP, colored by dataset, split by RNA reads per cell",
        legend_title="Dataset",
        group_by="new.ident",
        split_by="quartile_nCount_RNA",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_umi", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title="UMAP, colored by dataset, split by genes per cell",
        legend_title="Dataset",
        group_by="new.ident",
        split_by="quartile_nFeature_RNA",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_gene", sep="_"),
        pdf=args$pdf
    )

    if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="UMAP, split by dataset",
            legend_title="Dataset",
            group_by="new.ident",
            split_by="new.ident",
            label=FALSE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_spl_idnt", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="ccpca",
            plot_title="PCA, colored by cell cycle phase, split by dataset",
            legend_title="Phase",
            group_by="Phase",
            split_by="new.ident",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "ccpca_spl_idnt", sep="_"),
            pdf=args$pdf
        )
    }

    if (
        all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
        length(unique(as.vector(as.character(seurat_data@meta.data$condition)))) > 1
    ){
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="UMAP, colored by dataset, split by grouping condition",
            legend_title="Dataset",
            group_by="new.ident",
            split_by="condition",
            label=FALSE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        if ("Phase" %in% colnames(seurat_data@meta.data)){
            graphics$dim_plot(
                data=seurat_data,
                reduction="rnaumap",
                plot_title="UMAP, colored by grouping condition, split by cell cycle phase",
                legend_title="Condition",
                group_by="condition",
                split_by="Phase",
                label=FALSE,
                alpha=0.5,
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_gr_cnd_spl_ph", sep="_"),
                pdf=args$pdf
            )
            graphics$dim_plot(
                data=seurat_data,
                reduction="ccpca",
                plot_title="PCA, colored by cell cycle phase, split by grouping condition",
                legend_title="Phase",
                group_by="Phase",
                split_by="condition",
                label=FALSE,
                alpha=0.5,
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "ccpca_spl_cnd", sep="_"),
                pdf=args$pdf
            )
        }
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="UMAP, colored by grouping condition, split by mitochondrial percentage",
            legend_title="Condition",
            group_by="condition",
            split_by="quartile_mito_percentage",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_gr_cnd_spl_mito", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="UMAP, colored by grouping condition, split by RNA reads per cell",
            legend_title="Condition",
            group_by="condition",
            split_by="quartile_nCount_RNA",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_gr_cnd_spl_umi", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="UMAP, colored by grouping condition, split by genes per cell",
            legend_title="Condition",
            group_by="condition",
            split_by="quartile_nFeature_RNA",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_gr_cnd_spl_gene", sep="_"),
            pdf=args$pdf
        )
    } 
}

## ----
get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell RNA-Seq Dimensionality Reduction Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include genes",
            "expression information stored in the RNA assay."
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
        "--cellcycle",
        help=paste(
            "Path to the TSV/CSV file with the information for cell cycle score assignment.",
            "First column - 'phase', second column 'gene_id'. If loaded Seurat object already",
            "includes cell cycle scores in 'S.Score', 'G2M.Score', and 'CC.Difference' metatada",
            "columns they will be overwritten.",
            "Default: skip cell cycle score assignment."
        ),
        type="character"
    )
    parser$add_argument(
        "--norm",
        help=paste(
            "Normalization method applied to genes expression counts. If loaded Seurat object",
            "includes multiple datasets, normalization will be run independently for each of",
            "them, unless integration is disabled with 'none' or set to 'harmony'",
            "Default: sct"
        ),
        type="character",
        default="sctglm",
        choices=c("sct", "log", "sctglm")
    )
    parser$add_argument(
        "--ntgr",
        help=paste(
            "Integration method used for joint analysis of multiple datasets. Automatically",
            "set to 'none' if loaded Seurat object includes only one dataset.",
            "Default: seurat"
        ),
        type="character",
        default="seurat",
        choices=c("seurat", "harmony", "none")
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
        "--highvargenes",
        help=paste(
            "Number of highly variable genes used in datasets integration, scaling and",
            "dimensionality reduction.",
            "Default: 3000"
        ),
        type="integer", default=3000
    )
    parser$add_argument(
        "--regressmt",
        help=paste(
            "Regress the percentage of RNA reads mapped to mitochondrial genes as a",
            "confounding source of variation.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--removegenes",
        help=paste(
            "Regex pattern to identify the genes which should",
            "be removed from the list of the most variable genes,",
            "so they don't impact neither integration nor PCA.",
            "Default: none"
        ),
        type="character"
    )
    cell_cycle_group <- parser$add_mutually_exclusive_group()
    cell_cycle_group$add_argument(
        "--regressccfull",
        help=paste(
            "Regress all signals associated with cell cycle phase.",
            "Ignored if --cellcycle is not provided. Mutually exclusive",
            "with --regressccdiff parameter.",
            "Default: false"
        ),
        action="store_true"
    )
    cell_cycle_group$add_argument(
        "--regressccdiff",
        help=paste(
            "Regress only differences in cell cycle phase among proliferating",
            "cells. Signals separating non-cycling and cycling cells will be",
            "maintained. Ignored if --cellcycle is not provided. Mutually",
            "exclusive with --regressccfull",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use for datasets integration (if provided RDS",
            "file includes multiple datasets and --ntgr is not set to 'harmony')",
            "and UMAP projection (from 1 to 50). Set to 0 to use auto-estimated",
            "dimensionality. Default: 10"
        ),
        type="integer", default=10
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
        help="Save raw counts from the RNA assay to h5ad file. Default: false",
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
        "--scope",
        help=paste(
            "Save Seurat data to SCope compatible loom file.",
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
        "--lowmem",
        help=paste(
            "Attempts to minimize RAM usage when integrating multiple datasets",
            "with SCTransform algorithm (slows down the computation). Ignored if",
            "'--ntgr' is not set to 'seurat' or if '--norm' is not set to either",
            "'sct' or 'sctglm'.",
            "Default: false"
        ),
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
    logger$setup(
        file.path(dirname(ifelse(args$output == "", "./", args$output)), "error_report.txt"),
        header="Single-Cell RNA-Seq Dimensionality Reduction Analysis (sc_rna_reduce.R)"
    )
    print(args)
    return (args)
}

## ----
args <- get_args()
if (args$dimensions != 0){                       # we extend it to array only when we don't want to auto-estimate it
    print("Adjusting --dimensions parameter")
    args$dimensions <- c(1:args$dimensions)
    print(paste("--dimensions was adjusted to", paste(args$dimensions, collapse=", ")))
}
prod$parallel(args)

## ----
print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
debug$print_info(seurat_data, args)

## ----
if (!("RNA" %in% names(seurat_data@assays))){
    logger$info(
        paste(
            "Loaded Seurat object doesn't",
            "include the required RNA assay.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"

## ----
print(paste("Loading cell cycle data from", args$cellcycle))
cell_cycle_data <- io$load_cell_cycle_data(seurat_data, args$cellcycle)

## ----
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

## ----
if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)    # sets identities to new.ident
}
debug$print_info(seurat_data, args)

## ----
if (!is.null(args$removegenes)){
    genes_to_remove <- grep(                                              # might return an empty list if nothing found
        args$removegenes,
        as.vector(as.character(rownames(seurat_data))),                   # with RNA assay set as default the rownames should be genes
        value=TRUE,
        ignore.case=TRUE
    )
    if (length(genes_to_remove) > 0){                                     # found some genes based on pattern
        print(
            paste(
                "Based on the", args$removegenes, "pattern",
                "the following genes will be removed from",
                "the list of the most variable genes:",
                paste(genes_to_remove, collapse=", ")
            )
        )
        args$removegenes <- genes_to_remove                               # for easy access
    } else {
        args$removegenes <- NULL                                          # no genes found, reset it to NULL
    }
}

## ----
print("Running RNA analysis")
seurat_data <- analyses$rna_analyze(seurat_data, args, cell_cycle_data)   # adds "pca" and "rnaumap" reductions
seurat_data <- filter$collapse_fragments_list(seurat_data)                # collapse repetitive fragments if ATAC assay was present in the Seurat object and was splitted
debug$print_info(seurat_data, args)

## ----
print("Quantifying QC metrics")
seurat_data <- qc$quartile_qc_metrics(
    seurat_data=seurat_data,
    features=c("nCount_RNA", "nFeature_RNA", "mito_percentage"),
    prefix="quartile"                                                     # we use this prefix in DEFAULT_META_FIELDS from ucsc.R
)
debug$print_info(seurat_data, args)

## ----
export_all_dimensionality_plots(
    seurat_data=seurat_data,
    args=args
)

## ----
if(args$cbbuild){
    print("Reordering reductions to have rnaumap on the first place")                      # will be shown first in UCSC Cellbrowser
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c("rnaumap", reduc_names[reduc_names!="rnaumap"])               # rnaumap will be added by this time
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        rootname=paste(args$output, "_cellbrowser", sep=""),
    )
}

## ----
DefaultAssay(seurat_data) <- "RNA"                                                         # better to stick to RNA assay by default https://www.biostars.org/p/395951/#395954 
io$export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))

## ----
if(args$h5seurat){
    io$export_h5seurat(seurat_data, paste(args$output, "_data.h5seurat", sep=""))
}

## ----
if(args$h5ad){
    io$export_h5ad(
        data=seurat_data,
        location=paste(args$output, "_counts.h5ad", sep=""),
        assay="RNA",
        slot="counts"
    )
}

## ----
if(args$loupe){
    ucsc$export_loupe(
        seurat_data=seurat_data,
        assay="RNA",
        rootname=paste0(args$output, "_counts")
    )
}

## ----
if(args$scope){
    io$export_scope_loom(                                                                  # we save only counts slot from the RNA assay 
        seurat_data,
        paste(args$output, "_data.loom", sep="")
    )
}