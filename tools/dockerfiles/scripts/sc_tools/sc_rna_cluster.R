#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(knitr))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(stringr))
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

## ----
export_all_qc_plots <- function(seurat_data, args){

    for (i in 1:length(args$resolution)){
        current_resolution <- args$resolution[i]
        current_cluster_column <- paste("rna_res", current_resolution, sep=".")
        Idents(seurat_data) <- current_cluster_column            # othervise it will be split by dataset in vln_plot

        graphics$geom_bar_plot(
            data=seurat_data@meta.data,
            x_axis=current_cluster_column,
            color_by=current_cluster_column,
            x_label="Cluster",
            y_label="Cells",
            legend_title="Cluster",
            plot_title="Number of cells per cluster",
            plot_subtitle=paste(
                "All cells;",
                "resolution", current_resolution
            ),
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "cell_cnts_gr_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )

        graphics$geom_point_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_RNA",
            y_axis="nFeature_RNA",
            split_by=current_cluster_column,
            color_by="mito_percentage",
            highlight_rows=which(seurat_data@meta.data$rna_doublets == "doublet"),
            gradient_colors=c("lightslateblue", "orange", "red"),                       # to have the same scale, even if only blue values are present
            color_limits=c(0, 100),
            color_break=ceiling(max(seurat_data@meta.data$mito_percentage)),            # we take the round max, because data is already filtered
            x_label="RNA reads per cell",
            y_label="Genes per cell",
            legend_title="Mitochondrial %",
            plot_title="Genes vs RNA reads per cell",
            plot_subtitle=paste(
                "Split by cluster;",
                "all cells;",
                "resolution", current_resolution
            ),
            scale_x_log10=TRUE,
            scale_y_log10=TRUE,
            show_lm=TRUE,
            show_density=TRUE,
            density_bins=4,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "gene_umi_spl_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )

        graphics$geom_point_plot(
            data=seurat_data@meta.data,
            x_axis="mito_percentage",
            y_axis="nCount_RNA",
            split_by=current_cluster_column,
            color_by="mito_percentage",
            highlight_rows=which(seurat_data@meta.data$rna_doublets == "doublet"),
            gradient_colors=c("lightslateblue", "orange", "red"),
            color_limits=c(0, 100),
            color_break=ceiling(max(seurat_data@meta.data$mito_percentage)),            # we take the round max, because data is already filtered
            x_label="Mitochondrial %",
            y_label="RNA reads per cell",
            legend_title="Mitochondrial %",
            plot_title="RNA reads vs mitochondrial % per cell",
            plot_subtitle=paste(
                "Split by cluster;",
                "all cells;",
                "resolution", current_resolution
            ),
            scale_x_log10=TRUE,
            scale_y_log10=TRUE,
            show_lm=FALSE,
            show_density=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umi_mito_spl_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )

        if (nrow(seurat_data@meta.data[seurat_data@meta.data$rna_doublets == "doublet", ]) > 0){      # show plot only if we still have doublets
            graphics$composition_plot(
                data=seurat_data,
                plot_title="Percentage of RNA doublets per cluster",
                plot_subtitle=paste(
                    "All cells;",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by="rna_doublets",
                split_by=current_cluster_column,
                x_label="Cluster",
                y_label="Cell percentage",
                palette_colors=c("#00AEAE", "#0BFFFF"),
                theme=args$theme,
                rootname=paste(args$output, "rnadbl_gr_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }

        graphics$vln_plot(
            data=seurat_data,
            features=c("nCount_RNA", "nFeature_RNA", "mito_percentage"),
            labels=c("RNA reads", "Genes", "Mitochondrial %"),
            scale_y_log10=c(TRUE, TRUE, FALSE),
            from_meta=TRUE,
            show_box_plots=TRUE,
            plot_title="Distribution of QC metrics per cell colored by cluster",
            plot_subtitle=paste(
                "All cells;",
                "resolution", current_resolution
            ),
            legend_title="Cluster",
            pt_size=0,
            combine_guides="collect",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "qc_mtrcs_dnst_gr_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
    }
}

## ----
export_all_clustering_plots <- function(seurat_data, args){
    Idents(seurat_data) <- "new.ident"                                                               # safety measure
    datasets_count <- length(unique(as.vector(as.character(seurat_data@meta.data$new.ident))))
    conditions_count <- length(unique(as.vector(as.character(seurat_data@meta.data$condition))))
    not_default_conditions <- all(
        as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))
    )

    if (datasets_count > 1){
        Idents(seurat_data) <- "new.ident"
        min_dataset_size <- analyses$get_min_ident_size(SplitObject(seurat_data, split.by="new.ident"))
        print(paste("Downsampling to", min_dataset_size, "cells per datasets"))
        downsampled_per_dataset <- subset(seurat_data, downsample=min_dataset_size)    # downsample per "new.ident"
        if (conditions_count > 1 && not_default_conditions){
            Idents(downsampled_per_dataset) <- "condition"
            min_condition_size <- analyses$get_min_ident_size(SplitObject(downsampled_per_dataset, split.by="condition"))
            print(paste("Additionally downsampling to", min_condition_size, "cells per grouping condition"))
            downsampled_per_condition <- subset(downsampled_per_dataset, downsample=min_condition_size)
            Idents(downsampled_per_dataset) <- "new.ident"
        }
    }

    if ("Phase" %in% colnames(seurat_data@meta.data)){                             # plots not related to clustering resolution but more useful in
        if (datasets_count > 1){                                                   # the clustering pipeline rather than in dim. reduction step
            graphics$dim_plot(
                data=downsampled_per_dataset,
                reduction="rnaumap",
                plot_title="UMAP colored by cell cycle phase",
                plot_subtitle=paste(
                    "Split by dataset;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset"
                ),
                legend_title="Phase",
                group_by="Phase",
                split_by="new.ident",
                label=FALSE,
                label_color="black",
                palette_colors=graphics$CC_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_gr_ph_spl_idnt", sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_per_dataset,
                plot_title="Composition plot colored by cell cycle phase",
                plot_subtitle=paste(
                    "Split by dataset;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset"
                ),
                legend_title="Phase",
                group_by="Phase",
                split_by="new.ident",
                x_label="Dataset",
                y_label="Cell percentage",
                palette_colors=graphics$CC_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_ph_spl_idnt", sep="_"),
                pdf=args$pdf
            )
            if (conditions_count > 1 && not_default_conditions){
                graphics$dim_plot(
                    data=downsampled_per_condition,
                    reduction="rnaumap",
                    plot_title="UMAP colored by cell cycle phase",
                    plot_subtitle=paste(
                        "Split by grouping condition;",
                        "first downsampled to", min_dataset_size,
                        "cells per dataset,",
                        "then downsampled to", min_condition_size,
                        "cells per grouping condition"
                    ),
                    legend_title="Phase",
                    group_by="Phase",
                    split_by="condition",
                    label=FALSE,
                    label_color="black",
                    palette_colors=graphics$CC_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "umap_gr_ph_spl_cnd", sep="_"),
                    pdf=args$pdf
                )
                graphics$composition_plot(
                    data=downsampled_per_condition,
                    plot_title="Composition plot colored by cell cycle phase",
                    plot_subtitle=paste(
                        "Split by grouping condition;",
                        "first downsampled to", min_dataset_size,
                        "cells per dataset,",
                        "then downsampled to", min_condition_size,
                        "cells per grouping condition"
                    ),
                    legend_title="Phase",
                    group_by="Phase",
                    split_by="condition",
                    x_label="Condition",
                    y_label="Cell percentage",
                    palette_colors=graphics$CC_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cmp_gr_ph_spl_cnd", sep="_"),
                    pdf=args$pdf
                )
            }
        }
    }

    for (i in 1:length(args$resolution)){
        current_resolution <- args$resolution[i]
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="UMAP colored by cluster",
            plot_subtitle=paste(
                "All cells;",
                "resolution", current_resolution
            ),
            legend_title="Cluster",
            group_by=paste("rna_res", current_resolution, sep="."),
            label=FALSE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_gr_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        graphics$silhouette_plot(
            data=seurat_data,
            reduction="pca",
            dims=args$dimensions,
            downsample=500,
            plot_title="Silhouette scores",
            plot_subtitle=paste(
                "All cells;",
                "resolution", current_resolution
            ),
            legend_title="Cluster",
            group_by=paste("rna_res", current_resolution, sep="."),
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "slh_gr_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (datasets_count > 1){
            graphics$dim_plot(
                data=downsampled_per_dataset,
                reduction="rnaumap",
                plot_title="UMAP colored by cluster",
                plot_subtitle=paste(
                    "Split by dataset;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset;",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("rna_res", current_resolution, sep="."),
                split_by="new.ident",
                show_density=TRUE,
                label=FALSE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_gr_clst_spl_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_per_dataset,
                plot_title="Composition plot colored by cluster",
                plot_subtitle=paste(
                    "Split by dataset;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset;",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("rna_res", current_resolution, sep="."),
                split_by="new.ident",
                x_label="Dataset",
                y_label="Cell percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_clst_spl_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_per_dataset,
                plot_title="Composition plot colored by dataset",
                plot_subtitle=paste(
                    "Split by cluster;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset;",
                    "resolution", current_resolution
                ),
                legend_title="Dataset",
                group_by="new.ident",
                split_by=paste("rna_res", current_resolution, sep="."),
                bar_position="dodge",
                x_label="Cluster",
                y_label="Cell counts",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_idnt_spl_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            if ("Phase" %in% colnames(seurat_data@meta.data)){
                graphics$dim_plot(
                    data=downsampled_per_dataset,
                    reduction="rnaumap",
                    plot_title="UMAP colored by cluster",
                    plot_subtitle=paste(
                        "Split by cell cycle phase;",
                        "downsampled to", min_dataset_size,
                        "cells per dataset;",
                        "resolution", current_resolution
                    ),
                    legend_title="Cluster",
                    group_by=paste("rna_res", current_resolution, sep="."),
                    split_by="Phase",
                    label=FALSE,
                    label_color="black",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "umap_gr_clst_spl_ph_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
                graphics$composition_plot(
                    data=downsampled_per_dataset,
                    plot_title="Composition plot colored by cell cycle phase",
                    plot_subtitle=paste(
                        "Split by cluster;",
                        "downsampled to", min_dataset_size,
                        "cells per dataset;",
                        "resolution", current_resolution
                    ),
                    legend_title="Phase",
                    group_by="Phase",
                    split_by=paste("rna_res", current_resolution, sep="."),
                    bar_position="dodge",
                    x_label="Cluster",
                    y_label="Cell counts",
                    palette_colors=graphics$CC_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cmp_gr_ph_spl_clst_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
            }
            if (conditions_count > 1 && not_default_conditions){
                graphics$dim_plot(
                    data=downsampled_per_condition,
                    reduction="rnaumap",
                    plot_title="UMAP colored by cluster",
                    plot_subtitle=paste(
                        "Split by grouping condition;",
                        "first downsampled to", min_dataset_size,
                        "cells per dataset,",
                        "then downsampled to", min_condition_size,
                        "cells per grouping condition;",
                        "resolution", current_resolution
                    ),
                    legend_title="Cluster",
                    group_by=paste("rna_res", current_resolution, sep="."),
                    split_by="condition",
                    show_density=TRUE,
                    label=FALSE,
                    label_color="black",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "umap_gr_clst_spl_cnd_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
                graphics$composition_plot(
                    data=downsampled_per_condition,
                    plot_title="Composition plot colored by cluster",
                    plot_subtitle=paste(
                        "Split by grouping condition;",
                        "first downsampled to", min_dataset_size,
                        "cells per dataset,",
                        "then downsampled to", min_condition_size,
                        "cells per grouping condition;",
                        "resolution", current_resolution
                    ),
                    legend_title="Cluster",
                    group_by=paste("rna_res", current_resolution, sep="."),
                    split_by="condition",
                    x_label="Condition",
                    y_label="Cell percentage",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cmp_gr_clst_spl_cnd_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
                graphics$composition_plot(
                    data=downsampled_per_condition,
                    plot_title="Composition plot colored by grouping condition",
                    plot_subtitle=paste(
                        "Split by cluster;",
                        "first downsampled to", min_dataset_size,
                        "cells per dataset,",
                        "then downsampled to", min_condition_size,
                        "cells per grouping condition;",
                        "resolution", current_resolution
                    ),
                    legend_title="Condition",
                    group_by="condition",
                    split_by=paste("rna_res", current_resolution, sep="."),
                    bar_position="dodge",
                    x_label="Cluster",
                    y_label="Cell counts",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cmp_gr_cnd_spl_clst_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
            }
        } else {
            if ("Phase" %in% colnames(seurat_data@meta.data)){                # the same plots but with not downsampled data
                graphics$dim_plot(
                    data=seurat_data,
                    reduction="rnaumap",
                    plot_title="UMAP colored by cluster",
                    plot_subtitle=paste(
                        "Split by cell cycle phase;",
                        "resolution", current_resolution
                    ),
                    legend_title="Cluster",
                    group_by=paste("rna_res", current_resolution, sep="."),
                    split_by="Phase",
                    label=FALSE,
                    label_color="black",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "umap_gr_clst_spl_ph_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
                graphics$composition_plot(
                    data=seurat_data,
                    plot_title="Composition plot colored by cell cycle phase",
                    plot_subtitle=paste(
                        "Split by cluster;",
                        "resolution", current_resolution
                    ),
                    legend_title="Phase",
                    group_by="Phase",
                    split_by=paste("rna_res", current_resolution, sep="."),
                    bar_position="dodge",
                    x_label="Cluster",
                    y_label="Cell counts",
                    palette_colors=graphics$CC_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cmp_gr_ph_spl_clst_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
            }
        }
    }
    rm(downsampled_per_dataset, downsampled_per_condition)
    gc(verbose=FALSE)
}

## ----
export_all_expression_plots <- function(seurat_data, args) {
    DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure

    for (i in 1:length(args$genes)){
        current_gene <- args$genes[i]
        graphics$feature_plot(
            data=seurat_data,
            features=current_gene,
            labels=current_gene,
            reduction="rnaumap",
            plot_title="UMAP colored by gene expression",
            legend_title="Expression",
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
            plot_title="UMAP colored by gene expression density",
            joint=FALSE,
            width=800,
            height=800,
            theme=args$theme,
            rootname=paste(args$output, "xpr_per_cell_sgnl", current_gene, sep="_"),
            pdf=args$pdf
        )
    }

    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        Idents(seurat_data) <- paste("rna_res", current_resolution, sep=".")
        graphics$dot_plot(
            data=seurat_data,
            features=args$genes,
            plot_title="Average gene expression",
            plot_subtitle=paste(
                "Resolution", current_resolution
            ),
            x_label="Genes",
            y_label="Clusters",
            cluster_idents=FALSE,
            theme=args$theme,
            rootname=paste(args$output, "xpr_avg_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        graphics$vln_plot(
            data=seurat_data,
            features=args$genes,
            labels=args$genes,
            plot_title="Gene expression density",
            plot_subtitle=paste(
                "Resolution", current_resolution
            ),
            legend_title="Cluster",
            rotate_labels=TRUE,
            pt_size=0,
            combine_guides="collect",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "xpr_dnst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
    }

    Idents(seurat_data) <- "new.ident"                            # safety measure
}

## ----
export_heatmaps <- function(seurat_data, args){
    DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure
    datasets_count <- length(unique(as.vector(as.character(seurat_data@meta.data$new.ident))))
    conditions_count <- length(unique(as.vector(as.character(seurat_data@meta.data$condition))))
    not_default_conditions <- all(
        as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))
    )

    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        cluster_column <- paste("rna_res", current_resolution, sep=".")
        clusters_order <- levels(seurat_data@meta.data[[cluster_column]])
        if (is.null(clusters_order)){
            clusters_order <- unique(seurat_data@meta.data[[cluster_column]])
        }

        base::print(
            base::paste(
                "Checking if RNA markers for", cluster_column,
                "metadata column have been calculated"
            )
        )
        filtered_markers <- NULL
        base::tryCatch(
            expr = {
                filtered_markers <- seurat_data@misc$markers$RNA[[cluster_column]] %>%                  # should be ok even if selected markers have 0 rows
                                    dplyr::mutate(
                                        cluster=base::factor(
                                            cluster,
                                            levels=clusters_order
                                        )
                                    ) %>%
                                    dplyr::arrange(cluster) %>%                                         # to have upregulated genes on the diagonal
                                    dplyr::filter(.$p_val_adj <= 0.05) %>%                              # to have only significant gene markers
                                    dplyr::filter(.$pct.1 >= 0.1) %>%                                   # to have at least 10% of cells expressing this gene
                                    dplyr::group_by(feature) %>%
                                    dplyr::arrange(desc(pct.1), .by_group=TRUE) %>%                     # sort all duplicated features by desc(pct.1)
                                    dplyr::slice_head(n=1) %>%                                          # choose the feature with the highest pct.1
                                    dplyr::ungroup() %>%
                                    dplyr::group_by(cluster) %>%
                                    dplyr::arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>%
                                    dplyr::group_modify(~ .x %>%
                                        dplyr::slice_head(n=analyses$get_fraction(.x, 0.25))            # take 25% of the features
                                    ) %>%
                                    dplyr::arrange(desc(avg_log2FC), .by_group=TRUE) %>%
                                    dplyr::ungroup()
            },
            error = function(e){
                base::print(
                    base::paste(
                        "Failed to find RNA markers for", cluster_column,
                        "metadata column due to", e
                    )
                )
            }
        )

        if (!is.null(filtered_markers) && nrow(filtered_markers) > 0){                             # in case we don't have any markers for specific resolution
            column_annotations <- c("Cluster")
            seurat_data@meta.data$Cluster <- seurat_data@meta.data[[cluster_column]]               # will be overwritten for each resolution
            if (conditions_count > 1 && not_default_conditions){
                column_annotations <- c(column_annotations, "Condition")                           # several conditions found
                seurat_data@meta.data$Condition <- seurat_data@meta.data$condition
            }
            if (datasets_count > 1){
                column_annotations <- c(column_annotations, "Dataset")                             # several datasets found
                seurat_data@meta.data$Dataset <- seurat_data@meta.data$new.ident
            }
            graphics$feature_heatmap(                                                              # install.packages("magick") for better rasterization
                data=seurat_data,
                assay="RNA",
                slot="data",
                features=filtered_markers$feature,
                highlight_features=args$genes,                                                     # can be NULL
                show_rownames=FALSE,
                scale_to_max=FALSE,
                scale="row",                                                                       # will calculate z-score
                heatmap_colors=c("darkblue", "black", "yellow"),
                group_by=column_annotations,
                palette_colors=graphics$D40_COLORS,
                plot_title=paste0(
                    "Gene expression heatmap (resolution ", current_resolution, ")"
                ),
                rootname=paste(args$output, "xpr_htmp_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            io$export_data(
                filtered_markers,
                paste(args$output, "_xpr_htmp_res_", current_resolution, ".tsv", sep="")
            )
        }
    }
    seurat_data@meta.data$Cluster <- NULL                         # safety measure
    seurat_data@meta.data$Condition <- NULL                       # safety measure
    seurat_data@meta.data$Dataset <- NULL                         # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure
}

## ----
get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell RNA-Seq Cluster Analysis")
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
            "Dimensionality to use when constructing nearest-neighbor",
            "graph before clustering (from 1 to 50). Set to 0 to use",
            "auto-estimated dimensionality. Default: 10"
        ),
        type="integer", default=10
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
        "--cbbuild",
        help="Export results to UCSC Cell Browser. Default: false",
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
print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"
debug$print_info(seurat_data, args)

## ----
if (!all(c("pca", "rnaumap") %in% names(seurat_data@reductions))){
    print("Loaded Seurat object doesn't have 'pca' and/or 'rnaumap' reduction(s). Exiting.")
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if ("Phase" %in% colnames(seurat_data@meta.data)){                                 # safety measure to have the proper order on the plots
    seurat_data@meta.data[["Phase"]] <- factor(
        seurat_data@meta.data[["Phase"]],
        levels=c("G1", "S", "G2M")
    )
}

## ----
print("Adjusting --dimensions parameter")
if (length(args$dimensions) == 1 && args$dimensions == 0){
    print("Estimating datasets dimensionality")
    estimated_dimensionality <- analyses$get_dimensionality(seurat_data, "pca")    # won't be bigger than the max pca dimensions
    print(paste("Estimated dimensionality:", estimated_dimensionality))
    args$dimensions <- c(1:estimated_dimensionality)
} else {
    args$dimensions <- c(
        1:min(args$dimensions, ncol(Embeddings(seurat_data, "pca")))               # use min to make sure it's within the max pca dimensions
    )
}
print(paste("--dimensions was adjusted to", paste(args$dimensions, collapse=", ")))

## ----
print(paste("Clustering RNA data using", paste(args$dimensions, collapse=", "), "principal components"))
seurat_data <- analyses$add_clusters(
    seurat_data=seurat_data,
    assay="RNA",
    graph_name="rna",                          # will be used in all the plot generating functions
    reduction="pca",
    args=args
)
debug$print_info(seurat_data, args)

## ----
export_all_qc_plots(
    seurat_data=seurat_data,
    args=args
)

## ----
export_all_clustering_plots(
    seurat_data=seurat_data,
    args=args
)

## ----
if (!is.null(args$genes)){
    print("Adjusting genes of interest to include only those that are present in the loaded Seurat object")
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]     # with RNA assay set as default the rownames should be genes
    if (length(args$genes) == 0){
        print("Neither of the provided genes of interest are present in the loaded Seurat object. Setting back to NULL")
        args$genes <- NULL                                                                       # all genes were filtered out, so return to NULL
    } else {
        print(paste("--genes was adjusted to", paste(args$genes, collapse=", ")))
    }
}

## ----
if (!is.null(args$genes) || args$diffgenes) {                                                    # we check it so we don't normalize it without reason
    print("Normalizing counts in RNA assay before evaluating genes expression or identifying putative gene markers")
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
    if (!is.null(args$genes)) {
        print("Generating genes expression plots")
        export_all_expression_plots(seurat_data=seurat_data, args=args)
    }
    if(args$diffgenes){
        print("Identifying differentially expressed genes between each pair of clusters for all resolutions")
        seurat_data <- analyses$get_markers_by_res(
            seurat_data=seurat_data,
            assay="RNA",
            resolution_prefix="rna_res",
            args=args
        )
        debug$print_info(seurat_data, args)
        io$export_markers(
            data=seurat_data,
            assay="RNA",
            markers_regex="^rna_res",
            location=paste0(args$output, "_gene_markers.tsv")
        )
        export_heatmaps(
            seurat_data=seurat_data,
            args=args
        )
    }
}

## ----
if(args$cbbuild){
    print("Reordering reductions to have rnaumap on the first place")                      # will be shown first in UCSC Cellbrowser
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c("rnaumap", reduc_names[reduc_names!="rnaumap"])               # we checked before that rnaumap is present
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        label_field=paste0("Clustering (rna ", args$resolution[1], ")"),                   # always use only the first resolution
        rootname=paste(args$output, "_cellbrowser", sep="")
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
        active_cluster=paste("rna_res", args$resolution[1], sep="."),
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