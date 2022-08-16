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
    selected_features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "S.Score", "G2M.Score")
    selected_labels=c("UMI", "Genes", "Mitochondrial %", "Novelty score", "S score", "G to M score")

    graphics$elbow_plot(
        data=seurat_data,
        ndims=50,
        reduction="pca",
        plot_title="Elbow plot (from cells PCA)",
        theme=args$theme,
        rootname=paste(args$output, "elbow", sep="_"),
        pdf=args$pdf
    )
    graphics$corr_plot(
        data=seurat_data,
        reduction="pca",
        highlight_dims=args$dimensions,
        qc_columns=selected_features,
        qc_labels=selected_labels,
        plot_title="Correlation plots between QC metrics and cells PCA components",
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
        reduction="rnaumap",
        plot_title="Cells UMAP",
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
            plot_title="Split by cell cycle phase cells UMAP",
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
    }
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title="Split by the percentage of transcripts mapped to mitochondrial genes cells UMAP",
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
        plot_title="Split by the UMI per cell counts cells UMAP",
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
        plot_title="Split by the genes per cell counts cells UMAP",
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

    if (seurat_data@meta.data$new.ident != seurat_data@meta.data$condition){
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
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
        if ("Phase" %in% colnames(seurat_data@meta.data)){
            graphics$dim_plot(
                data=seurat_data,
                reduction="rnaumap",
                plot_title="Grouped by condition split by cell cycle cells UMAP",
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
        }
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="Grouped by condition split by the percentage of transcripts mapped to mitochondrial genes cells UMAP",
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
            plot_title="Grouped by condition split by the UMI per cell counts cells UMAP",
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
            plot_title="Grouped by condition split by the genes per cell counts cells UMAP",
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


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell RNA-Seq Dimensionality Reduction Analysis")
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
            "includes cell cycle scores in 'S.Score' and 'G2M.Score' metatada columns they will",
            "be removed.",
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
        default="sct",
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
            "Regress the percentage of transcripts mapped to mitochondrial genes as a",
            "confounding source of variation.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--regressgenes",
        help=paste(
            "Genes which expression should be regressed as a confounding source of variation.",
            "Default: None"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--regresscellcycle",
        help=paste(
            "Regress cell cycle scores as a confounding source of variation.",
            "Ignored if --cellcycle is not provided.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use in UMAP projection (from 1 to 50). If single value N",
            "is provided, use from 1 to N PCs. If multiple values are provided, subset to",
            "only selected PCs. In combination with --ntgr set to harmony, selected principle",
            "components will be used in Harmony integration.",
            "Default: from 1 to 10"
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

print(
    paste(
        "Setting parallelization to", args$cpus, "cores, and", args$memory,
        "GB of memory allowed to be shared between the processes"
    )
)
prod$parallel(args)

print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"
debug$print_info(seurat_data, args)

print(paste("Loading cell cycle data from", args$cellcycle))
cell_cycle_data <- io$load_cell_cycle_data(args$cellcycle)
if(!is.null(cell_cycle_data) && any(c("S.Score", "G2M.Score") %in% colnames(seurat_data@meta.data))){
    print(
        paste(
            "Removing S.Score and G2M.Score columns from",
            "the metadata of the loaded Seurat object",
            "as cell cycle scores will be re-assigned."
        )
    )
    seurat_data[["S.Score"]] <- NULL
    seurat_data[["G2M.Score"]] <- NULL
}

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

if (!is.null(args$regressgenes)){
    print("Adjusting genes to be regressed to include only those that are present in the loaded Seurat object")
    args$regressgenes <- unique(args$regressgenes)
    args$regressgenes <- args$regressgenes[args$regressgenes %in% as.vector(as.character(rownames(seurat_data)))]     # with RNA assay set as default the rownames should be genes
    print(args$regressgenes)
    if (!is.null(args$regressgenes) && length(args$regressgenes) > 0){
        print("Calculating the percentage of transcripts mapped to the genes that should be regressed out")
        seurat_data <- qc$add_gene_expr_percentage(
            seurat_data=seurat_data,
            target_genes=args$regressgenes
        )
        debug$print_info(seurat_data, args)
    }
}

print("Running RNA analysis")
seurat_data <- analyses$rna_analyze(seurat_data, args, cell_cycle_data)   # adds "pca" and "rnaumap" reductions
seurat_data <- filter$collapse_fragments_list(seurat_data)                # collapse repetitive fragments if ATAC assay was present in the Seurat object and was splitted
debug$print_info(seurat_data, args)

print("Quantifying QC metrics to evaluate unwanted sources of variation")
seurat_data <- qc$quartile_qc_metrics(
    seurat_data=seurat_data,
    features=c("nCount_RNA", "nFeature_RNA", "mito_percentage"),          # use only those QC metrics that we allow to regress out
    prefix="quartile"                                                     # we will use this prefix in ucsc$export_cellbrowser function
)
debug$print_info(seurat_data, args)

export_all_dimensionality_plots(
    seurat_data=seurat_data,
    args=args
)

if(args$cbbuild){
    print("Exporting RNA assay to UCSC Cellbrowser")
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        rootname=paste(args$output, "_cellbrowser", sep=""),
    )
}

DefaultAssay(seurat_data) <- "RNA"                                                         # better to stick to RNA assay by default https://www.biostars.org/p/395951/#395954 
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