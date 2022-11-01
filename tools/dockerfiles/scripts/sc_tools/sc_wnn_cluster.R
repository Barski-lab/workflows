#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(argparse))
suppressMessages(library(GenomicRanges))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


export_all_clustering_plots <- function(seurat_data, args){
    Idents(seurat_data) <- "new.ident"                                                               # safety measure
    downsampled_to <- analyses$get_min_ident_size(SplitObject(seurat_data, split.by="new.ident"))    # need to split it for consistency
    downsampled_data <- subset(seurat_data, downsample=downsampled_to)
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        graphics$dim_plot(
            data=seurat_data,
            reduction="wnnumap",
            plot_title=paste("Clustered cells UMAP. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste("wsnn_res", current_resolution, sep="."),
            label=TRUE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
            graphics$dim_plot(
                data=seurat_data,
                reduction="wnnumap",
                plot_title=paste("Split by dataset clustered cells UMAP. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="new.ident",
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cluster split by dataset cells composition plot.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="new.ident",
                x_label="Dataset",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_clst_spl_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by dataset split by cluster cells composition plot.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Dataset",
                group_by="new.ident",
                split_by=paste("wsnn_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_idnt_spl_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
        if (all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition)))){
            graphics$dim_plot(
                data=seurat_data,
                reduction="wnnumap",
                plot_title=paste("Split by grouping condition clustered cells UMAP. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="condition",
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_cnd_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cluster split by condition cells composition plot.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="condition",
                x_label="Condition",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_clst_spl_cnd_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by condition split by cluster cells composition plot.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution
                ),
                legend_title="Condition",
                group_by="condition",
                split_by=paste("wsnn_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_cnd_spl_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
        if ("Phase" %in% colnames(seurat_data@meta.data)){
            graphics$dim_plot(
                data=seurat_data,
                reduction="wnnumap",
                plot_title=paste("Split by cell cycle phase clustered cells UMAP. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="Phase",
                label=TRUE,
                label_color="black",
                alpha=0.5,
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_ph_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
                graphics$composition_plot(
                    data=downsampled_data,
                    plot_title=paste(
                        "Grouped by cell cycle phase split by dataset cells composition plot.",
                        "Downsampled to", downsampled_to, "cells per dataset.",
                        "Resolution", current_resolution),
                    legend_title="Phase",
                    group_by="Phase",
                    split_by="new.ident",
                    x_label="Dataset",
                    y_label="Cells percentage",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cmp_gr_ph_spl_idnt_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
            }
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cell cycle phase split by cluster cells composition plot.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Phase",
                group_by="Phase",
                split_by=paste("wsnn_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_ph_spl_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
    }
    rm(downsampled_data)
    gc(verbose=FALSE)
}


export_all_expression_plots <- function(seurat_data, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        Idents(seurat_data) <- paste("wsnn_res", current_resolution, sep=".")
        graphics$dot_plot(
            data=seurat_data,
            features=args$genes,
            plot_title=paste("Log normalized scaled average gene expression per cluster. Resolution", current_resolution),
            x_label="Genes",
            y_label="Clusters",
            cluster_idents=FALSE,
            theme=args$theme,
            rootname=paste(args$output, "xpr_avg_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (length(args$genes) > 0){
            for (i in 1:length(args$genes)) {
                current_gene <- args$genes[i]
                graphics$feature_plot(
                    data=seurat_data,
                    features=current_gene,
                    labels=current_gene,
                    reduction="wnnumap",
                    plot_title=paste("Log normalized gene expression on cells UMAP. Resolution", current_resolution),
                    label=TRUE,
                    order=TRUE,
                    max_cutoff="q99",  # to prevent cells with overexpressed gene from distoring the color bar
                    combine_guides="keep",
                    width=800,
                    height=800,
                    theme=args$theme,
                    rootname=paste(args$output, "xpr_per_cell_res", current_resolution, current_gene, sep="_"),
                    pdf=args$pdf
                )
                graphics$expression_density_plot(
                    data=seurat_data,
                    features=current_gene,
                    reduction="wnnumap",
                    plot_title=paste("Log normalized gene expression density on cells UMAP. Resolution", current_resolution),
                    joint=FALSE,
                    width=800,
                    height=800,
                    theme=args$theme,
                    rootname=paste(args$output, "xpr_per_cell_sgnl_res", current_resolution, current_gene, sep="_"),
                    pdf=args$pdf
                )
                graphics$vln_plot(
                    data=seurat_data,
                    features=current_gene,
                    labels=current_gene,
                    plot_title=paste("Log normalized gene expression density per cluster. Resolution", current_resolution),
                    legend_title="Cluster",
                    log=TRUE,
                    pt_size=0,
                    combine_guides="collect",
                    width=800,
                    height=600,
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "xpr_dnst_res", current_resolution, current_gene, sep="_"),
                    pdf=args$pdf
                )
            }
        }
    }
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
}

export_all_coverage_plots <- function(seurat_data, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"                                          # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                           # safety measure

    genome_annotation <- Annotation(seurat_data)                                               # safety measure to build the coverage plot
    if( !("gene_biotype" %in% base::colnames(GenomicRanges::mcols(genome_annotation))) ){
        print("Genome annotation doesn't have 'gene_biotype' column. Adding NA")
        genome_annotation$gene_biotype <- NA
        Annotation(seurat_data) <- genome_annotation
    }

    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        for (i in 1:length(args$genes)) {
            current_gene <- args$genes[i]
            graphics$coverage_plot(
                data=seurat_data,
                assay="ATAC",
                region=current_gene,
                group_by=paste("wsnn_res", current_resolution, sep="."),
                plot_title=paste(
                    "Tn5 insertion frequency plot around", current_gene, "gene.",
                    "Resolution", current_resolution
                ),
                idents=NULL,                                                   # to include all values from the default "new.ident" column
                cells=colnames(seurat_data),                                   # limit to only those cells that are in out seurat_data
                features=current_gene,
                expression_assay="RNA",
                expression_slot="data",                                        # use scaled counts
                extend_upstream=2500,
                extend_downstream=2500,
                show_annotation=TRUE,
                show_peaks=TRUE,
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cvrg_res", current_resolution, current_gene, sep="_"),
                pdf=args$pdf
            )
        }
    }
}


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell WNN Cluster Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "genes expression and chromatin accessibility information stored in the RNA",
            "and ATAC assays correspondingly. Additionally, 'pca', 'rnaumap', 'atac_lsi'",
            "and 'atacumap' dimensionality reductions should be present."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--rnadimensions",
        help=paste(
            "Dimensionality from the 'pca' reduction to use when constructing weighted",
            "nearest-neighbor graph before clustering (from 1 to 50). If single value N",
            "is provided, use from 1 to N dimensions. If multiple values are provided,",
            "subset to only selected dimensions.",
            "Default: from 1 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--atacdimensions",
        help=paste(
            "Dimensionality from the 'atac_lsi' reduction to use when constructing weighted",
            "nearest-neighbor graph before clustering (from 1 to 50). If single value N",
            "is provided, use from 2 to N dimensions. If multiple values are provided,",
            "subset to only selected dimensions.",
            "Default: from 2 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--algorithm",
        help=paste(
            "Algorithm for modularity optimization when running clustering.",
            "Default: louvain"
        ),
        type="character", default="slm",
        choices=c(
            "louvain", "mult-louvain", "slm", "leiden"
        )
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
        "--resolution",
        help=paste(
            "Clustering resolution applied to the constructed weighted nearest-neighbor",
            "graph. Can be set as an array.",
            "Default: 0.3, 0.5, 1.0"
        ),
        type="double", default=c(0.3, 0.5, 1.0), nargs="*"
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment used in the loaded Seurat",
            "object. File should be saved in TSV format with tbi-index file."
        ),
        type="character"
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to build gene expression and Tn5 insertion frequency plots",
            "for the nearest peaks. If '--fragments' is not provided only gene expression",
            "plots will be built.",
            "Default: None"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--diffgenes",
        help=paste(
            "Identify differentially expressed genes (putative gene markers) between each",
            "pair of clusters for all resolutions.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--diffpeaks",
        help=paste(
            "Identify differentially accessible peaks between each pair of clusters for all resolutions.",
            "Default: false"
        ),
        action="store_true"
    )

    parser$add_argument(
        "--rnalogfc",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "on average have log fold change difference in expression between every",
            "tested pair of clusters not lower than this value. Ignored if '--diffgenes'",
            "is not set.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--rnaminpct",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "are detected in not lower than this fraction of cells in either of the",
            "two tested clusters. Ignored if '--diffgenes' is not set.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--rnaonlypos",
        help=paste(
            "For putative gene markers identification return only positive markers.",
            "Ignored if '--diffgenes' is not set.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--rnatestuse",
        help=paste(
            "Statistical test to use for putative gene markers identification.",
            "Ignored if '--diffgenes' is not set.",
            "Default: wilcox"
        ),
        type="character", default="wilcox",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
    )
    parser$add_argument(
        "--ataclogfc",
        help=paste(
            "For differentially accessible peaks identification include only those peaks that",
            "on average have log fold change difference in the chromatin accessibility between",
            "every tested pair of clusters not lower than this value. Ignored if '--diffpeaks'",
            "is not set.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--atacminpct",
        help=paste(
            "For differentially accessible peaks identification include only those peaks that",
            "are detected in not lower than this fraction of cells in either of the two tested",
            "clusters. Ignored if '--diffpeaks' is not set.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--atactestuse",
        help=paste(
            "Statistical test to use for differentially accessible peaks identification.",
            "Ignored if '--diffpeaks' is not set.",
            "Default: LR"
        ),
        type="character", default="LR",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
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
if (length(args$rnadimensions) == 1) {
    print("Adjusting --rnadimensions parameter as only a single value was provided")
    args$rnadimensions <- c(1:args$rnadimensions[1])
    print(paste("--rnadimensions was adjusted to", paste(args$rnadimensions, collapse=", ")))
}
if (length(args$atacdimensions) == 1) {
    print("Adjusting --atacdimensions parameter as only a single value was provided")
    args$atacdimensions <- c(2:args$atacdimensions[1])                                              # skipping the first LSI component
    print(paste("--atacdimensions was adjusted to", paste(args$atacdimensions, collapse=", ")))
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
DefaultAssay(seurat_data) <- "RNA"                                                               # for consistency with other scripts
debug$print_info(seurat_data, args)

if (!all(c("pca", "atac_lsi") %in% names(seurat_data@reductions))){
    print("Loaded Seurat object doesn't have 'pca' and/or 'atac_lsi' reduction(s). Exiting.")
    quit(save="no", status=1, runLast=FALSE)
}

if (!is.null(args$fragments)){
    print(paste("Loading fragments data from", args$fragments))
    seurat_data <- io$replace_fragments(args$fragments, seurat_data)                             # will change the default assay to ATAC
    debug$print_info(seurat_data, args)
}

print(
    paste(
        "Running weighted nearest-neighbor analysis using", paste(args$rnadimensions, collapse=", "),
        "dimensions from 'pca' and", paste(args$atacdimensions, collapse=", "), "dimensions from 'atac_lsi'",
        "dimensionality reductions."
    )
)
seurat_data <- analyses$add_wnn_clusters(                       # will add 'wnnumap' reduction
    seurat_data=seurat_data,
    graph_name="wsnn",                                          # will be used in all the plot generating functions
    reductions=list("pca", "atac_lsi"),
    dimensions=list(args$rnadimensions, args$atacdimensions),   # should be the same order as reductions
    args=args
)
debug$print_info(seurat_data, args)

export_all_clustering_plots(seurat_data=seurat_data, args=args)

nearest_peaks <- NULL
if (!is.null(args$genes)){
    print("Adjusting genes of interest to include only those that are present in the loaded Seurat object")
    args$genes <- unique(args$genes)
    DefaultAssay(seurat_data) <- "RNA"                                                                       # need it for rownames to return genes
    args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]                 # with RNA assay set as default the rownames should be genes
    DefaultAssay(seurat_data) <- "ATAC"                                                                      # Annotation needs the default assay to be ATAC
    args$genes <- args$genes[args$genes %in% as.vector(as.character(Annotation(seurat_data)$gene_name))]     # just in case check if the same genes are present in Annotation
    all_peaks <- StringToGRanges(rownames(seurat_data), sep=c("-", "-"))                                     # rownames are peaks when default assay is ATAC
    nearest_peaks <- sapply(
        args$genes,
        function(gene)
        GRangesToString(
            all_peaks[
                subjectHits(
                    distanceToNearest(
                        LookupGeneCoords(seurat_data, gene=gene, assay="ATAC"),
                        all_peaks
                    )
                )
            ]
        )
    )
    rm(all_peaks)
    print(nearest_peaks)
}

if(args$cbbuild){
    print("Exporting RNA assay to UCSC Cellbrowser")
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        features=args$genes,                                   # can be NULL
        is_nested=TRUE,
        rootname=paste(args$output, "_cellbrowser/rna", sep=""),
    )
    print("Exporting ATAC assay to UCSC Cellbrowser")
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="ATAC",
        slot="counts",
        short_label="ATAC",
        features=nearest_peaks,                               # use nearest to the genes if interest peaks
        is_nested=TRUE,
        rootname=paste(args$output, "_cellbrowser/atac", sep=""),
    )
}

if (!is.null(args$genes) || args$diffgenes) {
    print("Normalizing counts in RNA assay")
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
    if (!is.null(args$genes)){
        print("Generating genes expression plots")
        export_all_expression_plots(seurat_data=seurat_data, args=args)
    }
    if (!is.null(args$genes) && !is.null(args$fragments)){
        print("Generating coverage plots")
        export_all_coverage_plots(seurat_data=seurat_data, args=args)
    }
    if(args$diffgenes){
        print("Identifying differentially expressed genes between each pair of clusters for all resolutions")
        args$logfc <- args$rnalogfc                                  # need the proper names for get_putative_markers
        args$minpct <- args$rnaminpct
        args$onlypos <- args$rnaonlypos
        args$testuse <- args$rnatestuse
        all_rna_putative_markers <- analyses$get_putative_markers(   # will change default assay to RNA
            seurat_data=seurat_data,
            assay="RNA",
            resolution_prefix="wsnn_res",
            args=args
        )
        io$export_data(
            all_rna_putative_markers,
            paste(args$output, "_gene_markers.tsv", sep="")
        )
        rm(all_rna_putative_markers)
    }
}

if (args$diffpeaks){
    print("Identifying differentially accessible peaks between each pair of clusters for all resolutions")
    args$logfc <- args$ataclogfc                                     # need the proper names for get_putative_markers
    args$minpct <- args$atacminpct
    args$testuse <- args$atactestuse
    all_atac_putative_markers <- analyses$get_putative_markers(      # will change default assay to ATAC
        seurat_data=seurat_data,
        assay="ATAC",
        resolution_prefix="wsnn_res",
        latent_vars="nCount_ATAC",                                   # to remove the influence of sequencing depth
        args=args
    )
    io$export_data(
        all_atac_putative_markers,
        paste(args$output, "_peak_markers.tsv", sep="")
    )
    rm(all_atac_putative_markers)
}

DefaultAssay(seurat_data) <- "RNA"
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