#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(forcats))
suppressMessages(library(argparse))
suppressMessages(library(tidyselect))

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
    print(paste("Downsampling datasets to", downsampled_to, "cells per sample"))
    downsampled_data <- subset(seurat_data, downsample=downsampled_to)
    if ("Phase" %in% colnames(seurat_data@meta.data)){
        if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
                graphics$composition_plot(
                    data=downsampled_data,
                    plot_title=paste(
                        "Grouped by cell cycle phase split by dataset cells composition plot.",
                        "Downsampled to", downsampled_to, "cells per dataset"),
                    legend_title="Phase",
                    group_by="Phase",
                    split_by="new.ident",
                    x_label="Dataset",
                    y_label="Cells percentage",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cmp_gr_ph_spl_idnt", sep="_"),
                    pdf=args$pdf
                )
        }
    }
    for (i in 1:length(args$resolution)){
        current_resolution <- args$resolution[i]
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title=paste("Clustered cells UMAP. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste("rna_res", current_resolution, sep="."),
            label=TRUE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        graphics$silhouette_plot(
            data=seurat_data,
            reduction="pca",
            dims=args$dimensions,
            downsample=500,
            plot_title=paste("Silhouette scores. Downsampled to max 500 cells per cluster. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste("rna_res", current_resolution, sep="."),
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "slh_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
            graphics$dim_plot(
                data=seurat_data,
                reduction="rnaumap",
                plot_title=paste("Split by dataset clustered cells UMAP. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("rna_res", current_resolution, sep="."),
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
                group_by=paste("rna_res", current_resolution, sep="."),
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
                split_by=paste("rna_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_idnt_spl_clst_res", current_resolution, sep="_"),
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
                plot_title=paste("Split by grouping condition clustered cells UMAP. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("rna_res", current_resolution, sep="."),
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
                group_by=paste("rna_res", current_resolution, sep="."),
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
                split_by=paste("rna_res", current_resolution, sep="."),
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
                reduction="rnaumap",
                plot_title=paste("Split by cell cycle phase clustered cells UMAP. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("rna_res", current_resolution, sep="."),
                split_by="Phase",
                label=TRUE,
                label_color="black",
                alpha=0.5,
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_ph_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cell cycle phase split by cluster cells composition plot.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Phase",
                group_by="Phase",
                split_by=paste("rna_res", current_resolution, sep="."),
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
    DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure
    if (!is.null(args$genes) && length(args$genes) > 0){
        for (i in 1:length(args$genes)){
            current_gene <- args$genes[i]
            graphics$feature_plot(
                data=seurat_data,
                features=current_gene,
                labels=current_gene,
                reduction="rnaumap",
                plot_title="Log normalized gene expression on cells UMAP",
                label=FALSE,
                order=TRUE,
                max_cutoff="q99",  # to prevent cells with overexpressed gene from distoring the color bar
                combine_guides="keep",
                width=800,
                height=800,
                theme=args$theme,
                rootname=paste(args$output, "xpr_per_cell", current_gene, sep="_"),
                pdf=args$pdf
            )
            graphics$expression_density_plot(
                data=seurat_data,
                features=current_gene,
                reduction="rnaumap",
                plot_title="Log normalized gene expression density on cells UMAP",
                joint=FALSE,
                width=800,
                height=800,
                theme=args$theme,
                rootname=paste(args$output, "xpr_per_cell_sgnl", current_gene, sep="_"),
                pdf=args$pdf
            )
        }
    }
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        Idents(seurat_data) <- paste("rna_res", current_resolution, sep=".")
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
        if (!is.null(args$genes) && length(args$genes) > 0){
            for (i in 1:length(args$genes)){
                current_gene <- args$genes[i]
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
    Idents(seurat_data) <- "new.ident"                            # safety measure
}


export_heatmaps <- function(seurat_data, markers, args){
    DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        grouped_markers <- markers %>%
                           dplyr::filter(.$resolution==current_resolution) %>%                     # shouldn't fail even if resolution is not present
                           dplyr::select(-c("resolution")) %>%
                           dplyr::group_by(cluster) %>%
                           dplyr::top_n(
                               n=tidyselect::all_of(floor(60/length(unique(markers$cluster)))),
                               wt=avg_log2FC
                           )
        if (nrow(grouped_markers) > 0){                                                            # in case we don't have any markers for specific resolution
            column_annotations <- c(paste("rna_res", current_resolution, sep="."))
            if (length(unique(as.vector(as.character(seurat_data@meta.data$new.ident)))) > 1){
                column_annotations <- c(column_annotations, "new.ident")                           # several datasets found
            }
            if (
                all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
                length(unique(as.vector(as.character(seurat_data@meta.data$condition)))) > 1
            ){
                column_annotations <- c(column_annotations, "condition")                           # several conditions found
            }
            custom_fields <- grep("^custom_", colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            if (length(custom_fields) > 0){
                column_annotations <- c(column_annotations, custom_fields)                         # adding all custom fields
            }
            graphics$feature_heatmap(                                                              # install.packages("magick") for better rasterization
                data=seurat_data,
                assay="RNA",
                slot="data",
                features=grouped_markers$feature,
                split_rows=forcats::fct_inorder(as.character(grouped_markers$cluster)),            # fct_inorder fails with numeric
                show_rownames=TRUE,
                scale_to_max=TRUE,
                group_by=column_annotations,
                palette_colors=graphics$D40_COLORS,
                heatmap_colors=c("black", "yellow"),
                plot_title=paste(
                    "Normalized gene expression heatmap. Resolution", current_resolution
                ),
                rootname=paste(args$output, "xpr_htmp_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
    }
    Idents(seurat_data) <- "new.ident"                            # safety measure
}


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell RNA-Seq Cluster Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include genes",
            "expression information stored in the RNA assay, as well as 'pca' and 'rnaumap'",
            "dimensionality reductions applied to that assay."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use when constructing nearest-neighbor graph before clustering",
            "(from 1 to 50). If single value N is provided, use from 1 to N dimensions. If",
            "multiple values are provided, subset to only selected dimensions.",
            "Default: from 1 to 10"
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
            "Default: louvain"
        ),
        type="character", default="louvain",
        choices=c(
            "louvain", "mult-louvain", "slm", "leiden"
        )
    )
    parser$add_argument(
        "--resolution",
        help=paste(
            "Clustering resolution applied to the constructed nearest-neighbor graph.",
            "Can be set as an array but only the first item from the list will be used",
            "for cluster labels and gene markers in the UCSC Cell Browser when running",
            "with --cbbuild and --diffgenes parameters.",
            "Default: 0.3, 0.5, 1.0"
        ),
        type="double", default=c(0.3, 0.5, 1.0), nargs="*"
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to build genes expression plots.",
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
        "--logfc",
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
        "--minpct",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "are detected in not lower than this fraction of cells in either of the",
            "two tested clusters. Ignored if '--diffgenes' is not set.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--onlypos",
        help=paste(
            "For putative gene markers identification return only positive markers.",
            "Ignored if '--diffgenes' is not set.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--testuse",
        help=paste(
            "Statistical test to use for putative gene markers identification.",
            "Ignored if '--diffgenes' is not set.",
            "Default: wilcox"
        ),
        type="character", default="wilcox",
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

if (!all(c("pca", "rnaumap") %in% names(seurat_data@reductions))){
    print("Loaded Seurat object doesn't have 'pca' and/or 'rnaumap' reduction(s). Exiting.")
    quit(save="no", status=1, runLast=FALSE)
}

print(paste("Clustering RNA data using", paste(args$dimensions, collapse=", "), "principal components"))
seurat_data <- analyses$add_clusters(
    seurat_data=seurat_data,
    assay="RNA",
    graph_name="rna",                          # will be used in all the plot generating functions
    reduction="pca",
    args=args
)
debug$print_info(seurat_data, args)

export_all_clustering_plots(
    seurat_data=seurat_data,
    args=args
)

if (!is.null(args$genes)){
    print("Adjusting genes of interest to include only those that are present in the loaded Seurat object")
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]     # with RNA assay set as default the rownames should be genes
    print(args$genes)
}

all_markers <- NULL
if (!is.null(args$genes) || args$diffgenes) {
    print("Normalizing counts in RNA assay before evaluating genes expression or identifying putative gene markers")
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
    if (!is.null(args$genes)){
        print("Generating genes expression plots")
        export_all_expression_plots(seurat_data=seurat_data, args=args)
    }
    if(args$diffgenes){
        print("Identifying differentially expressed genes between each pair of clusters for all resolutions")
        all_markers <- analyses$get_markers_by_res(                                        # either NULL or definitely not empty
            seurat_data=seurat_data,
            assay="RNA",
            resolution_prefix="rna_res",
            args=args
        )
        if (!is.null(all_markers)){
            io$export_data(
                all_markers,
                paste(args$output, "_gene_markers.tsv", sep="")
            )
            export_heatmaps(
                seurat_data=seurat_data,
                markers=all_markers,
                args=args
            )
        }
    }
}

if(args$cbbuild){
    print("Exporting RNA assay to UCSC Cellbrowser")
    if(!is.null(all_markers)){
        all_markers <- all_markers %>%
                       dplyr::filter(.$resolution==args$resolution[1]) %>%                 # won't fail even if resolution is not present
                       dplyr::select(-c("resolution"))
    }
    print("Reordering reductions to have rnaumap on the first place")                      # will be shown first in UCSC Cellbrowser
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c("rnaumap", reduc_names[reduc_names!="rnaumap"])               # we checked before that rnaumap is present
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]
    debug$print_info(seurat_data, args)
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        markers=all_markers,                                                               # can be NULL
        palette_colors=graphics$D40_COLORS,                                                # to have colors correspond to the plots
        label_field=paste0("Clustering (rna ", args$resolution[1], ")"),                   # always use only the first resolution
        rootname=paste(args$output, "_cellbrowser", sep="")
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