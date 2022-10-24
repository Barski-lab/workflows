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
            label=TRUE,
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
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_idnt_rd", reduction, sep="_"),
                pdf=args$pdf
            )
        }
        if (all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition)))){
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
                label=TRUE,
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
                label=TRUE,
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
    
    if (all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition)))){
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
        Annotation(seurat_data) <- genome_annotation
    }
    for (i in 1:length(args$genes)) {
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
                label=TRUE,
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
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
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
            "assigned cell types."
        ),
        type="character", required="True"
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

if(args$cbbuild){
    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        print("Exporting RNA and ATAC assays to UCSC Cellbrowser jointly")
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            features=args$genes,                                   # can be NULL
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/rna", sep=""),
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            features=nearest_peaks,                               # use nearest to the genes if interest peaks
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/atac", sep=""),
        )
    } else if ("RNA" %in% names(seurat_data@assays)){
        print("Exporting RNA assay to UCSC Cellbrowser")
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            features=args$genes,                                   # can be NULL
            rootname=paste(args$output, "_cellbrowser", sep=""),
        )
    } else {
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