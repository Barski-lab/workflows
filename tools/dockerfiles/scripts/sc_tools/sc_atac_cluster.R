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
            reduction="atacumap",
            plot_title=paste("Clustered cells UMAP. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste("atac_res", current_resolution, sep="."),
            label=TRUE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        graphics$silhouette_plot(
            data=seurat_data,
            reduction="atac_lsi",
            dims=args$dimensions,
            downsample=500,
            plot_title=paste("Silhouette scores. Downsampled to max 500 cells per cluster. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste("atac_res", current_resolution, sep="."),
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "slh_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
            graphics$dim_plot(
                data=seurat_data,
                reduction="atacumap",
                plot_title=paste("Split by dataset clustered cells UMAP. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("atac_res", current_resolution, sep="."),
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
                group_by=paste("atac_res", current_resolution, sep="."),
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
                split_by=paste("atac_res", current_resolution, sep="."),
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
                reduction="atacumap",
                plot_title=paste("Split by grouping condition clustered cells UMAP. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("atac_res", current_resolution, sep="."),
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
                group_by=paste("atac_res", current_resolution, sep="."),
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
                split_by=paste("atac_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_cnd_spl_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
    }
    rm(downsampled_data)
    gc(verbose=FALSE)
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
                group_by=paste("atac_res", current_resolution, sep="."),
                plot_title=paste(
                    "Tn5 insertion frequency plot around", current_gene, "gene.",
                    "Resolution", current_resolution
                ),
                idents=NULL,                                                               # to include all values from the default "new.ident" column
                cells=colnames(seurat_data),                                               # limit to only those cells that are in out seurat_data
                features=if("RNA" %in% names(seurat_data@assays)) current_gene else NULL,  # will fail if features are provided without "RNA" assay
                expression_assay="RNA",
                expression_slot="data",                                                    # use scaled counts
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
    parser <- ArgumentParser(description="Single-cell ATAC-Seq Cluster Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "chromatin accessibility information stored in the ATAC assay, as well as",
            "'atac_lsi' and 'atacumap' dimensionality reductions applied to that assay."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use when constructing nearest-neighbor graph before clustering",
            "(from 1 to 50). If single value N is provided, use from 2 to N dimensions. If",
            "multiple values are provided, subset to only selected dimensions.",
            "Default: from 2 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--ametric",
        help=paste(
            "Distance metric used when constructing nearest-neighbor graph before clustering.",
            "Default: euclidean"
        ),
        type="character", default="euclidean",
        choices=c(
            "euclidean", "cosine", "manhattan", "hamming"
        )
    )
    parser$add_argument(
        "--algorithm",
        help=paste(
            "Algorithm for modularity optimization when running clustering.",
            "Default: slm"
        ),
        type="character", default="slm",
        choices=c(
            "louvain", "mult-louvain", "slm", "leiden"
        )
    )
    parser$add_argument(
        "--resolution",
        help=paste(
            "Clustering resolution applied to the constructed nearest-neighbor graph.",
            "Can be set as an array.",
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
            "Genes of interest to build Tn5 insertion frequency plots for the nearest peaks.",
            "If loaded Seurat object includes genes expression information in the RNA assay",
            "it will be additionally shown on the right side of the plots.",
            "Ignored if '--fragments' is not provided.",
            "Default: None"
        ),
        type="character", nargs="*"
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
        "--logfc",
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
        "--minpct",
        help=paste(
            "For differentially accessible peaks identification include only those peaks that",
            "are detected in not lower than this fraction of cells in either of the two tested",
            "clusters. Ignored if '--diffpeaks' is not set.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--testuse",
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
if (length(args$dimensions) == 1) {
    print("Adjusting --dimensions parameter as only a single value was provided")
    args$dimensions <- c(2:args$dimensions[1])                                              # skipping the first LSI component
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
print("Setting default assay to ATAC")
DefaultAssay(seurat_data) <- "ATAC"
debug$print_info(seurat_data, args)

if (!all(c("atac_lsi", "atacumap") %in% names(seurat_data@reductions))){
    print("Loaded Seurat object doesn't have 'atac_lsi' and/or 'atacumap' reduction(s). Exiting.")
    quit(save="no", status=1, runLast=FALSE)
}

if (!is.null(args$fragments)){
    print(paste("Loading fragments data from", args$fragments))
    seurat_data <- io$replace_fragments(args$fragments, seurat_data)
    debug$print_info(seurat_data, args)
}

print(paste("Clustering ATAC data using", paste(args$dimensions, collapse=", "), "dimensions"))
seurat_data <- analyses$add_clusters(
    seurat_data=seurat_data,
    assay="ATAC",
    graph_name="atac",                          # will be used in all the plot generating functions
    reduction="atac_lsi",
    args=args
)
debug$print_info(seurat_data, args)

export_all_clustering_plots(seurat_data=seurat_data, args=args)

nearest_peaks <- NULL
if (!is.null(args$genes)){
    DefaultAssay(seurat_data) <- "ATAC"         # safety measure
    print("Adjusting genes of interest to include only those that are present in the loaded Seurat object")
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(Annotation(seurat_data)$gene_name))]
    all_peaks <- StringToGRanges(rownames(seurat_data), sep=c("-", "-"))
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
    print("Exporting ATAC assay to UCSC Cellbrowser")
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="ATAC",
        slot="counts",
        short_label="ATAC",
        features=nearest_peaks,                               # use nearest to the genes if interest peaks
        rootname=paste(args$output, "_cellbrowser", sep=""),
    )
}

if (!is.null(args$genes) && !is.null(args$fragments)){
    if ("RNA" %in% names(seurat_data@assays)){
        print("Normalizing counts in RNA assay to show average gene expression alongside the coverage plots")
        DefaultAssay(seurat_data) <- "RNA"                         # for now we will always use RNA even if SCT can be present
        seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
        DefaultAssay(seurat_data) <- "ATAC"                        # safety measure
    }
    export_all_coverage_plots(
        seurat_data=seurat_data,
        args=args
    )
}

if (args$diffpeaks){
    print("Identifying differentially accessible peaks between each pair of clusters for all resolutions")
    all_putative_markers <- analyses$get_putative_markers(
        seurat_data=seurat_data,
        assay="ATAC",
        resolution_prefix="atac_res",
        latent_vars="nCount_ATAC",                                  # to remove the influence of sequencing depth
        args=args
    )
    io$export_data(
        all_putative_markers,
        paste(args$output, "_peak_markers.tsv", sep="")
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