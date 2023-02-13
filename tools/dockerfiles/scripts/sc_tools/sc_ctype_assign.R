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
    print(paste("Downsampling datasets to", downsampled_to, "cells per sample"))
    downsampled_data <- subset(seurat_data, downsample=downsampled_to)
    for (reduction in c("rnaumap", "atacumap", "wnnumap")){
        if (!(reduction %in% names(seurat_data@reductions))) {next}                                  # skip missing reductions
        graphics$dim_plot(
            data=seurat_data,
            reduction=reduction,
            plot_title=paste0(
                "Cells UMAP with assigned cell types (",
                reduction, " dim. reduction)"
            ),
            legend_title="Cell type",
            group_by=args$target,
            label=FALSE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_rd", reduction, sep="_"),
            pdf=args$pdf
        )
        if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
            graphics$dim_plot(
                data=seurat_data,
                reduction=reduction,
                plot_title=paste0(
                    "Split by dataset cells UMAP with assigned cell types (",
                    reduction, " dim. reduction)"
                ),
                legend_title="Cell type",
                group_by=args$target,
                split_by="new.ident",
                label=FALSE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_idnt_rd", reduction, sep="_"),
                pdf=args$pdf
            )
        }
        if (
            all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
            length(unique(as.vector(as.character(seurat_data@meta.data$condition)))) > 1
        ){
            graphics$dim_plot(
                data=seurat_data,
                reduction=reduction,
                plot_title=paste0(
                    "Split by grouping condition cells UMAP with assigned cell types (",
                    reduction, " dim. reduction)"
                ),
                legend_title="Cell type",
                group_by=args$target,
                split_by="condition",
                label=FALSE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_cnd_rd", reduction, sep="_"),
                pdf=args$pdf
            )
        }
        if ("Phase" %in% colnames(seurat_data@meta.data)){
            graphics$dim_plot(
                data=seurat_data,
                reduction=reduction,
                plot_title=paste0(
                    "Split by cell cycle phase cells UMAP with assigned cell types (",
                    reduction, " dim. reduction)"
                ),
                legend_title="Cell type",
                group_by=args$target,
                split_by="Phase",
                label=FALSE,
                label_color="black",
                alpha=0.5,
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_ph_rd", reduction, sep="_"),
                pdf=args$pdf
            )
        }
    }

    if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
        graphics$composition_plot(
            data=downsampled_data,
            plot_title=paste(
                "Grouped by cell type split by dataset cells composition plot.",
                "Downsampled to", downsampled_to, "cells per dataset."
            ),
            legend_title="Cell type",
            group_by=args$target,
            split_by="new.ident",
            x_label="Dataset",
            y_label="Cells counts",
            palette_colors=graphics$D40_COLORS,
            bar_position="dodge",
            theme=args$theme,
            rootname=paste(args$output, "cmp_gr_ctyp_spl_idnt", sep="_"),
            pdf=args$pdf
        )
        graphics$composition_plot(
            data=downsampled_data,
            plot_title=paste(
                "Grouped by dataset split by cell type cells composition plot.",
                "Downsampled to", downsampled_to, "cells per dataset."
            ),
            legend_title="Dataset",
            group_by="new.ident",
            split_by=args$target,
            x_label="Cell type",
            y_label="Cells counts",
            palette_colors=graphics$D40_COLORS,
            bar_position="dodge",
            theme=args$theme,
            rootname=paste(args$output, "cmp_gr_idnt_spl_ctyp", sep="_"),
            pdf=args$pdf
        )
        if ("Phase" %in% colnames(seurat_data@meta.data)){
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cell cycle phase split by dataset cells composition plot.",
                    "Downsampled to", downsampled_to, "cells per dataset."
                ),
                legend_title="Phase",
                group_by="Phase",
                split_by="new.ident",
                x_label="Dataset",
                y_label="Cells counts",
                palette_colors=graphics$D40_COLORS,
                bar_position="dodge",
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_ph_spl_idnt", sep="_"),
                pdf=args$pdf
            )
        }
    }
    if (
        all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
        length(unique(as.vector(as.character(seurat_data@meta.data$condition)))) > 1
    ){
        graphics$composition_plot(
            data=downsampled_data,
            plot_title=paste(
                "Grouped by cell type split by condition cells composition plot.",
                "Downsampled to", downsampled_to, "cells per dataset."
            ),
            legend_title="Cell type",
            group_by=args$target,
            split_by="condition",
            x_label="Condition",
            y_label="Cells counts",
            palette_colors=graphics$D40_COLORS,
            bar_position="dodge",
            theme=args$theme,
            rootname=paste(args$output, "cmp_gr_ctyp_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$composition_plot(
            data=downsampled_data,
            plot_title=paste(
                "Grouped by condition split by cell type cells composition plot.",
                "Downsampled to", downsampled_to, "cells per dataset."
            ),
            legend_title="Condition",
            group_by="condition",
            split_by=args$target,
            x_label="Cell type",
            y_label="Cells counts",
            palette_colors=graphics$D40_COLORS,
            bar_position="dodge",
            theme=args$theme,
            rootname=paste(args$output, "cmp_gr_cnd_spl_ctyp", sep="_"),
            pdf=args$pdf
        )
    }

    if ("Phase" %in% colnames(seurat_data@meta.data)){
        graphics$composition_plot(
            data=downsampled_data,
            plot_title=paste(
                "Grouped by cell cycle phase split by cell type cells composition plot.",
                "Downsampled to", downsampled_to, "cells per dataset."
            ),
            legend_title="Phase",
            group_by="Phase",
            split_by=args$target,
            x_label="Cell type",
            y_label="Cells counts",
            palette_colors=graphics$D40_COLORS,
            bar_position="dodge",
            theme=args$theme,
            rootname=paste(args$output, "cmp_gr_ph_spl_ctyp", sep="_"),
            pdf=args$pdf
        )
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
    }
    if( !("tx_id" %in% base::colnames(GenomicRanges::mcols(genome_annotation))) ){             # https://github.com/stuart-lab/signac/issues/1159
        print("Genome annotation doesn't have 'tx_id' column. Adding from 'transcript_id'")
        genome_annotation$tx_id <- genome_annotation$transcript_id
    }
    Annotation(seurat_data) <- genome_annotation

    if (!is.null(args$genes) && length(args$genes) > 0){
        for (i in 1:length(args$genes)){
            current_gene <- args$genes[i]
            graphics$coverage_plot(
                data=seurat_data,
                assay="ATAC",
                region=current_gene,
                group_by=args$target,
                plot_title=paste(
                    "Tn5 insertion frequency plot around", current_gene, "gene."
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
                rootname=paste(args$output, "cvrg", current_gene, sep="_"),
                pdf=args$pdf
            )
        }
    }
}


export_all_expression_plots <- function(seurat_data, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    Idents(seurat_data) <- args$target
    graphics$dot_plot(
        data=seurat_data,
        features=args$genes,
        plot_title=paste("Log normalized scaled average gene expression per cell type."),
        x_label="Genes",
        y_label="Cell type",
        cluster_idents=FALSE,
        theme=args$theme,
        rootname=paste(args$output, "xpr_avg", sep="_"),
        pdf=args$pdf
    )
    if (!is.null(args$genes) && length(args$genes) > 0){
        for (i in 1:length(args$genes)){
            current_gene <- args$genes[i]
            graphics$vln_plot(
                data=seurat_data,
                features=current_gene,
                labels=current_gene,
                plot_title=paste("Log normalized gene expression density per cell type"),
                legend_title="Cell type",
                log=TRUE,
                pt_size=0,
                combine_guides="collect",
                width=800,
                height=600,
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "xpr_dnst", current_gene, sep="_"),
                pdf=args$pdf
            )
            for (reduction in c("rnaumap", "atacumap", "wnnumap")){
                if (!(reduction %in% names(seurat_data@reductions))) {next}                                  # skip missing reductions
                graphics$feature_plot(
                    data=seurat_data,
                    features=current_gene,
                    labels=current_gene,
                    reduction=reduction,
                    plot_title=paste0("Log normalized gene expression on cells UMAP with assigned cell types (", reduction, " dim. reduction)"),
                    label=FALSE,
                    order=TRUE,
                    max_cutoff="q99",  # to prevent cells with overexpressed gene from distoring the color bar
                    combine_guides="keep",
                    width=800,
                    height=800,
                    theme=args$theme,
                    rootname=paste(args$output, "xpr_per_cell_rd", reduction, current_gene, sep="_"),
                    pdf=args$pdf
                )
                graphics$expression_density_plot(
                    data=seurat_data,
                    features=current_gene,
                    reduction=reduction,
                    plot_title=paste0("Log normalized gene expression density on cells UMAP with assigned cell types (", reduction, " dim. reduction)"),
                    joint=FALSE,
                    width=800,
                    height=800,
                    theme=args$theme,
                    rootname=paste(args$output, "xpr_per_cell_sgnl_rd", reduction, current_gene, sep="_"),
                    pdf=args$pdf
                )
            }
        }
    }
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
}


export_heatmaps <- function(seurat_data, markers, args){
    DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure
    grouped_markers <- markers %>%
                       dplyr::group_by(cluster) %>%
                       dplyr::top_n(
                           n=tidyselect::all_of(floor(60/length(unique(markers$cluster)))),
                           wt=avg_log2FC
                       )
    column_annotations <- c(args$target)
    if (length(unique(as.vector(as.character(seurat_data@meta.data$new.ident)))) > 1){
        column_annotations <- c(column_annotations, "new.ident")                           # several datasets found
    }
    if (
        all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
        length(unique(as.vector(as.character(seurat_data@meta.data$condition)))) > 1
    ){
        column_annotations <- c(column_annotations, "condition")                           # several conditions found
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
        plot_title="Normalized gene expression heatmap",
        rootname=paste(args$output, "xpr_htmp", sep="_"),
        pdf=args$pdf
    )
    Idents(seurat_data) <- "new.ident"                            # safety measure
}


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell Manual Cell Type Assignment")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "genes expression and/or chromatin accessibility information stored in the RNA",
            "and ATAC assays correspondingly. Additionally, 'rnaumap', and/or 'atacumap',",
            "and/or 'wnnumap' dimensionality reductions should be present."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--celltypes",
        help=paste(
            "Path to the TSV/CSV file for manual cell type assignment for each of the clusters.",
            "First column - 'cluster', second column may have arbitrary name."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--source",
        help=paste(
            "Column from the metadata of the loaded Seurat object to select clusters from."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--target",
        help=paste(
            "Column from the metadata of the loaded Seurat object to save manually",
            "assigned cell types. Should start with 'custom_', otherwise, it won't",
            "be shown in UCSC Cell Browser."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--diffgenes",
        help=paste(
            "Identify differentially expressed genes (putative gene markers) for",
            "assigned cell types. Ignored if loaded Seurat object doesn't include",
            "genes expression information stored in the RNA assay.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--diffpeaks",
        help=paste(
            "Identify differentially accessible peaks for assigned cell types. Ignored",
            "if loaded Seurat object doesn't include chromatin accessibility information",
            "stored in the ATAC assay.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--rnalogfc",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "on average have log fold change difference in expression between every",
            "tested pair of cell types not lower than this value. Ignored if '--diffgenes'",
            "is not set or RNA assay is not present.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--rnaminpct",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "are detected in not lower than this fraction of cells in either of the",
            "two tested cell types. Ignored if '--diffgenes' is not set or RNA assay",
            "is not present.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--rnaonlypos",
        help=paste(
            "For putative gene markers identification return only positive markers.",
            "Ignored if '--diffgenes' is not set or RNA assay is not present.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--rnatestuse",
        help=paste(
            "Statistical test to use for putative gene markers identification.",
            "Ignored if '--diffgenes' is not set or RNA assay is not present.",
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
            "every tested pair of cell types not lower than this value. Ignored if '--diffpeaks'",
            "is not set or ATAC assay is not present.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--atacminpct",
        help=paste(
            "For differentially accessible peaks identification include only those peaks that",
            "are detected in not lower than this fraction of cells in either of the two tested",
            "cell types. Ignored if '--diffpeaks' is not set or ATAC assay is not present.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--atactestuse",
        help=paste(
            "Statistical test to use for differentially accessible peaks identification.",
            "Ignored if '--diffpeaks' is not set or ATAC assay is not present.",
            "Default: LR"
        ),
        type="character", default="LR",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment used in the loaded Seurat",
            "object. File should be saved in TSV format with tbi-index file. Ignored if the",
            "loaded Seurat object doesn't include ATAC assay."
        ),
        type="character"
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to build gene expression and/or Tn5 insertion frequency plots",
            "for the nearest peaks. To build gene expression plots the loaded Seurat object",
            "should include RNA assay. To build Tn5 insertion frequency plots for the nearest",
            "peaks the loaded Seurat object should include ATAC assay as well as the --fragments",
            "file should be provided.",
            "Default: None"
        ),
        type="character", nargs="*"
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

print(
    paste(
        "Setting parallelization to", args$cpus, "cores, and", args$memory,
        "GB of memory allowed to be shared between the processes"
    )
)
prod$parallel(args)

print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
debug$print_info(seurat_data, args)

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

seurat_data <- io$extend_metadata(
    seurat_data=seurat_data,
    location=args$celltypes,
    seurat_ref_column=args$source,
    meta_ref_column="cluster",
    seurat_target_columns=args$target
)
debug$print_info(seurat_data, args)

if ( (!is.null(args$fragments)) && ("ATAC" %in% names(seurat_data@assays)) ){
    print(paste("Loading fragments data from", args$fragments))
    seurat_data <- io$replace_fragments(args$fragments, seurat_data)                                             # will change the default assay to ATAC
    debug$print_info(seurat_data, args)
}

export_all_clustering_plots(seurat_data=seurat_data, args=args)

nearest_peaks <- NULL                                                                                            # will be used only in UCSC Cell Browser build for ATAC assay
if (!is.null(args$genes)){
    print("Adjusting genes of interest to include only those that are present in the loaded Seurat object")
    args$genes <- unique(args$genes)
    if("RNA" %in% names(seurat_data@assays)){
        DefaultAssay(seurat_data) <- "RNA"                                                                       # need it for rownames to return genes
        args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]                 # with RNA assay set as default the rownames should be genes
    }
    if("ATAC" %in% names(seurat_data@assays)){
        DefaultAssay(seurat_data) <- "ATAC"                                                                      # Annotation needs the default assay to be ATAC
        args$genes <- args$genes[args$genes %in% as.vector(as.character(Annotation(seurat_data)$gene_name))]     # check if genes of interest are present in Annotation
    }
    if("RNA" %in% names(seurat_data@assays) && length(args$genes) > 0){                                          # check the length in case all the genes have been removed
        DefaultAssay(seurat_data) <- "RNA"                                                                       # need it for rownames to return genes
        seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
        print("Generating genes expression plots")
        export_all_expression_plots(seurat_data=seurat_data, args=args)                                          # changes default assay to RNA
    }
    if("ATAC" %in% names(seurat_data@assays)){
        DefaultAssay(seurat_data) <- "ATAC"                                                                      # Annotation needs the default assay to be ATAC
        all_peaks <- StringToGRanges(rownames(seurat_data), sep=c("-", "-"))                                     # rownames are peaks when default assay is ATAC
        nearest_peaks <- sapply(                                                                                 # we still might need them for UCSC even when fragments are not provided
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
        if (!is.null(args$fragments) && length(args$genes) > 0){                                                 # check the length in case all the genes have been removed
            print("Generating coverage plots")
            export_all_coverage_plots(seurat_data=seurat_data, args=args)                                        # changes default assay to ATAC
        }
        rm(all_peaks)
        print(nearest_peaks)
    }
}

all_rna_markers <- NULL
if (args$diffgenes && ("RNA" %in% names(seurat_data@assays))){
    print("Normalizing counts in RNA assay before identifying putative gene markers")
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data <- NormalizeData(seurat_data, verbose=FALSE)                              # this might be the second time we normalize data, but it should be ok to do it
    print("Identifying differentially expressed genes between each pair of cell types")
    args$logfc <- args$rnalogfc                                                           # need the proper names for get_markers
    args$minpct <- args$rnaminpct
    args$onlypos <- args$rnaonlypos
    args$testuse <- args$rnatestuse
    all_rna_markers <- analyses$get_markers(                                              # will change default assay to RNA
        seurat_data=seurat_data,
        assay="RNA",
        group_by=args$target,
        args=args
    )
    args <- args[names(args) %in% c("logfc", "minpct", "onlypos", "testuse") == FALSE]    # to remove temporary added items
    if (!is.null(all_rna_markers)){
        io$export_data(
            all_rna_markers,
            paste(args$output, "_gene_markers.tsv", sep="")
        )
        export_heatmaps(                                                                  # will change default assay to RNA
            seurat_data=seurat_data,
            markers=all_rna_markers,
            args=args
        )
    }
}

all_atac_markers <- NULL
if (args$diffpeaks && ("ATAC" %in% names(seurat_data@assays))){
    print("Identifying differentially accessible peaks between each pair of cell types")
    DefaultAssay(seurat_data) <- "ATAC"                                                   # safety measure
    args$logfc <- args$ataclogfc                                                          # need the proper names for get_markers
    args$minpct <- args$atacminpct
    args$onlypos <- FALSE                                                                 # need to overwrite what was set for RNA
    args$testuse <- args$atactestuse
    all_atac_markers <- analyses$get_markers(                                             # will change default assay to ATAC
        seurat_data=seurat_data,
        assay="ATAC",
        group_by=args$target,
        latent_vars="nCount_ATAC",                                                        # to remove the influence of sequencing depth
        args=args
    )
    args <- args[names(args) %in% c("logfc", "minpct", "onlypos", "testuse") == FALSE]    # to remove temporary added items
    if (!is.null(all_atac_markers)){
        io$export_data(
            all_atac_markers,
            paste(args$output, "_peak_markers.tsv", sep="")
        )
    }
}

if(args$cbbuild){
    print("Attemting to reorder reductions")
    reduc_names <- names(seurat_data@reductions)
    possible_suffixes <- c("wsnn", "rna", "atac")
    possible_names <- c("wnnumap", "rnaumap", "atacumap")
    first_reduc_name <- NULL
    for (i in 1:length(possible_suffixes)){
        if (grepl(possible_suffixes[i], args$source) && possible_names[i] %in% reduc_names){
            first_reduc_name <- possible_names[i]
            break
        }
    }
    if (!is.null(first_reduc_name)){
        print(paste("Moving", first_reduc_name, "reduction on top"))
        ordered_reduc_names <- c(first_reduc_name, reduc_names[reduc_names!=first_reduc_name])
        seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]
        debug$print_info(seurat_data, args)
    }

    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        print("Exporting RNA and ATAC assays to UCSC Cellbrowser jointly")
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            markers=all_rna_markers,                                         # can be NULL
            label_field <- base::gsub("custom_", "Custom ", args$target),
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/rna", sep="")
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            markers=all_atac_markers,                                        # can be NULL
            label_field <- base::gsub("custom_", "Custom ", args$target),
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/atac", sep="")
        )
    } else if ("RNA" %in% names(seurat_data@assays)){
        print("Exporting RNA assay to UCSC Cellbrowser")
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            markers=all_rna_markers,                                         # can be NULL
            label_field <- base::gsub("custom_", "Custom ", args$target),
            rootname=paste(args$output, "_cellbrowser", sep="")
        )
    } else {
        print("Exporting ATAC assay to UCSC Cellbrowser")
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            markers=all_atac_markers,                                        # can be NULL
            label_field <- base::gsub("custom_", "Custom ", args$target),
            rootname=paste(args$output, "_cellbrowser", sep="")
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