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
suppressMessages(library(GenomicRanges))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))

## ----
export_all_qc_plots <- function(seurat_data, args){
    Idents(seurat_data) <- args$target                                # othervise it will be split by dataset in vln_plot

    graphics$geom_bar_plot(
        data=seurat_data@meta.data,
        x_axis=args$target,
        color_by=args$target,
        x_label="Cell type",
        y_label="Cells",
        legend_title="Cell type",
        plot_title="Number of cells per cell type",
        plot_subtitle="All cells",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "cell_cnts_gr_ctyp", sep="_"),
        pdf=args$pdf
    )

    # need to select certain parameters depending on what assays are present
    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        color_by <- "mito_percentage"
        highlight_rows <- which(seurat_data@meta.data$rna_doublets == "doublet" | seurat_data@meta.data$atac_doublets == "doublet")
        gradient_colors <- c("lightslateblue", "orange", "red")
        color_limits <- c(0, 100)
        color_break <- ceiling(max(seurat_data@meta.data$mito_percentage))
        legend_title <- "Mitochondrial %"
        ncol <- 2
        width <- 2400
        height <- 1200
    } else if ("RNA" %in% names(seurat_data@assays)){
        color_by <- "mito_percentage"
        highlight_rows <- which(seurat_data@meta.data$rna_doublets == "doublet")
        gradient_colors <- c("lightslateblue", "orange", "red")
        color_limits <- c(0, 100)
        color_break <- ceiling(max(seurat_data@meta.data$mito_percentage))
        legend_title <- "Mitochondrial %"
        ncol <- 1
        width <- 1200
        height <- 800
    } else {
        color_by <- "frip"
        highlight_rows <- which(seurat_data@meta.data$atac_doublets == "doublet")
        gradient_colors <- c("orange", "lightslateblue", "lightslateblue")
        color_limits <- c(0, 1)
        color_break <- min(seurat_data@meta.data$frip)
        legend_title <- "FRiP"
        ncol <- 1
        width <- 1200
        height <- 1600
    }

    if ("RNA" %in% names(seurat_data@assays)){
        graphics$geom_point_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_RNA",
            y_axis="nFeature_RNA",
            split_by=args$target,
            color_by=color_by,
            highlight_rows=highlight_rows,
            gradient_colors=gradient_colors,
            color_limits=color_limits,
            color_break=color_break,
            legend_title=legend_title,
            x_label="RNA reads per cell",
            y_label="Genes per cell",
            plot_title="Genes vs RNA reads per cell",
            plot_subtitle=paste(
                "Split by cell type;",
                "all cells"
            ),
            scale_x_log10=TRUE,
            scale_y_log10=TRUE,
            show_lm=TRUE,
            show_density=TRUE,
            density_bins=4,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "gene_umi_spl_ctyp", sep="_"),
            pdf=args$pdf
        )

        graphics$geom_point_plot(
            data=seurat_data@meta.data,
            x_axis="mito_percentage",
            y_axis="nCount_RNA",
            split_by=args$target,
            color_by=color_by,
            highlight_rows=highlight_rows,
            gradient_colors=gradient_colors,
            color_limits=color_limits,
            color_break=color_break,
            legend_title=legend_title,
            x_label="Mitochondrial %",
            y_label="RNA reads per cell",
            plot_title="RNA reads vs mitochondrial % per cell",
            plot_subtitle=paste(
                "Split by cell type;",
                "all cells"
            ),
            scale_x_log10=TRUE,
            scale_y_log10=TRUE,
            show_lm=FALSE,
            show_density=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umi_mito_spl_ctyp", sep="_"),
            pdf=args$pdf
        )

        if (nrow(seurat_data@meta.data[seurat_data@meta.data$rna_doublets == "doublet", ]) > 0){
            graphics$composition_plot(
                data=seurat_data,
                plot_title="Percentage of RNA doublets per cell type",
                plot_subtitle="All cells",
                legend_title="Cell type",
                group_by="rna_doublets",
                split_by=args$target,
                x_label="Cell type",
                y_label="Cell percentage",
                palette_colors=c("#00AEAE", "#0BFFFF"),
                theme=args$theme,
                rootname=paste(args$output, "rnadbl_gr_ctyp", sep="_"),
                pdf=args$pdf
            )
        }
    }

    if ("ATAC" %in% names(seurat_data@assays)){
        graphics$geom_point_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_ATAC",
            y_axis="TSS.enrichment",
            split_by=args$target,
            color_by=color_by,
            highlight_rows=highlight_rows,
            gradient_colors=gradient_colors,
            color_limits=color_limits,
            color_break=color_break,
            legend_title=legend_title,
            x_label="ATAC fragments in peaks per cell",
            y_label="TSS enrichment score",
            plot_title="TSS enrichment score vs ATAC fragments in peaks per cell",
            plot_subtitle=paste(
                "Split by cell type;",
                "all cells"
            ),
            scale_x_log10=TRUE,
            scale_y_log10=FALSE,
            show_density=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "tss_frgm_spl_ctyp", sep="_"),
            pdf=args$pdf
        )

        if (nrow(seurat_data@meta.data[seurat_data@meta.data$atac_doublets == "doublet", ]) > 0){
            graphics$composition_plot(
                data=seurat_data,
                plot_title="Percentage of ATAC doublets per cell type",
                plot_subtitle="All cells",
                legend_title="Cell type",
                group_by="atac_doublets",
                split_by=args$target,
                x_label="Cell type",
                y_label="Cell percentage",
                palette_colors=c("#00DCDC", "#0BFFFF"),
                theme=args$theme,
                rootname=paste(args$output, "atacdbl_gr_ctyp", sep="_"),
                pdf=args$pdf
            )
        }
    }

    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        graphics$geom_point_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_ATAC",
            y_axis="nCount_RNA",
            split_by=args$target,
            color_by=color_by,
            highlight_rows=highlight_rows,
            gradient_colors=gradient_colors,
            color_limits=color_limits,
            color_break=color_break,
            legend_title=legend_title,
            x_label="ATAC fragments in peaks per cell",
            y_label="RNA reads per cell",
            plot_title="RNA reads vs ATAC fragments in peaks per cell",
            plot_subtitle=paste(
                "Split by cell type;",
                "all cells"
            ),
            scale_x_log10=TRUE,
            scale_y_log10=TRUE,
            show_density=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "rna_atac_cnts_spl_ctyp", sep="_"),
            pdf=args$pdf
        )

        if (
            nrow(seurat_data@meta.data[seurat_data@meta.data$rna_doublets == "doublet", ]) > 0 &&
            nrow(seurat_data@meta.data[seurat_data@meta.data$atac_doublets == "doublet", ]) > 0
        ){
            seurat_data@meta.data <- seurat_data@meta.data %>%
                                        dplyr::mutate(
                                            doublets_overlap = dplyr::case_when(
                                                rna_doublets == "singlet" & atac_doublets == "singlet" ~ "Singlet",
                                                rna_doublets == "doublet" & atac_doublets == "singlet" ~ "Only RNA",
                                                rna_doublets == "singlet" & atac_doublets == "doublet" ~ "Only ATAC",
                                                rna_doublets == "doublet" & atac_doublets == "doublet" ~ "RNA and ATAC"
                                            )
                                        ) %>%
                                        dplyr::mutate(
                                            doublets_overlap=base::factor(
                                                doublets_overlap,
                                                levels=c("Only RNA", "RNA and ATAC", "Only ATAC", "Singlet")
                                            )
                                        )
            graphics$composition_plot(
                data=seurat_data,
                plot_title="Percentage of RNA and ATAC doublets per cell type",
                plot_subtitle="All cells",
                legend_title="Cell type",
                group_by="doublets_overlap",
                split_by=args$target,
                x_label="Cell type",
                y_label="Cell percentage",
                palette_colors=c("#00AEAE", "#008080", "#00DCDC", "#0BFFFF"),
                theme=args$theme,
                rootname=paste(args$output, "vrlpdbl_gr_ctyp", sep="_"),
                pdf=args$pdf
            )
        }
    }

    graphics$vln_plot(                                                           # won't fail even if some columns are not present
        data=seurat_data,
        features=c(
            "nCount_RNA", "nFeature_RNA", "mito_percentage",
            "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "nFeature_ATAC", "frip", "blacklist_fraction"
        ),
        labels=c(
            "RNA reads", "Genes", "Mitochondrial %",
            "ATAC fragments in peaks", "TSS enrichment score", "Nucl. signal", "Peaks", "FRiP", "Bl. regions"
        ),
        scale_y_log10=c(
            TRUE, TRUE, FALSE,
            TRUE, FALSE, FALSE, TRUE, FALSE, FALSE
        ),
        from_meta=TRUE,
        show_box_plots=TRUE,
        plot_title="Distribution of QC metrics per cell colored by cell type",
        plot_subtitle="All cells",
        legend_title="Cell type",
        pt_size=0,
        combine_guides="collect",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "qc_mtrcs_dnst_gr_ctyp", sep="_"),
        pdf=args$pdf
    )
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

    graphics$dim_plot(
        data=seurat_data,
        reduction=args$reduction,
        plot_title="UMAP colored by cell type",
        plot_subtitle="All cells",
        legend_title="Cell type",
        group_by=args$target,
        label=FALSE,
        label_color="black",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "umap_gr_ctyp", sep="_"),
        pdf=args$pdf
    )
    if (datasets_count > 1){
        graphics$dim_plot(
            data=downsampled_per_dataset,
            reduction=args$reduction,
            plot_title="UMAP colored by cell type",
            plot_subtitle=paste(
                "Split by dataset;",
                "downsampled to", min_dataset_size,
                "cells per dataset"
            ),
            legend_title="Cell type",
            group_by=args$target,
            split_by="new.ident",
            show_density=TRUE,
            label=FALSE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_gr_ctyp_spl_idnt", sep="_"),
            pdf=args$pdf
        )
        graphics$composition_plot(
            data=downsampled_per_dataset,
            plot_title="Composition plot colored by cell type",
            plot_subtitle=paste(
                "Split by dataset;",
                "downsampled to", min_dataset_size,
                "cells per dataset"
            ),
            legend_title="Cell type",
            group_by=args$target,
            split_by="new.ident",
            x_label="Dataset",
            y_label="Cell percentage",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "cmp_gr_ctyp_spl_idnt", sep="_"),
            pdf=args$pdf
        )
        graphics$composition_plot(
            data=downsampled_per_dataset,
            plot_title="Composition plot colored by dataset",
            plot_subtitle=paste(
                "Split by cell type;",
                "downsampled to", min_dataset_size,
                "cells per dataset"
            ),
            legend_title="Dataset",
            group_by="new.ident",
            split_by=args$target,
            bar_position="dodge",
            x_label="Cell type",
            y_label="Cell counts",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "cmp_gr_idnt_spl_ctyp", sep="_"),
            pdf=args$pdf
        )
        if ("Phase" %in% colnames(seurat_data@meta.data)){
            graphics$dim_plot(
                data=downsampled_per_dataset,
                reduction=args$reduction,
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
            graphics$dim_plot(
                data=downsampled_per_dataset,
                reduction=args$reduction,
                plot_title="UMAP colored by cell type",
                plot_subtitle=paste(
                    "Split by cell cycle phase;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset"
                ),
                legend_title="Cell type",
                group_by=args$target,
                split_by="Phase",
                label=FALSE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_gr_ctyp_spl_ph", sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_per_dataset,
                plot_title="Composition plot colored by cell cycle phase",
                plot_subtitle=paste(
                    "Split by cell type;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset"
                ),
                legend_title="Phase",
                group_by="Phase",
                split_by=args$target,
                bar_position="dodge",
                x_label="Cell type",
                y_label="Cell counts",
                palette_colors=graphics$CC_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_ph_spl_ctyp", sep="_"),
                pdf=args$pdf
            )
        }
        if (conditions_count > 1 && not_default_conditions){
            graphics$dim_plot(
                data=downsampled_per_condition,
                reduction=args$reduction,
                plot_title="UMAP colored by cell type",
                plot_subtitle=paste(
                    "Split by grouping condition;",
                    "first downsampled to", min_dataset_size,
                    "cells per dataset,",
                    "then downsampled to", min_condition_size,
                    "cells per grouping condition"
                ),
                legend_title="Cell type",
                group_by=args$target,
                split_by="condition",
                show_density=TRUE,
                label=FALSE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_gr_ctyp_spl_cnd", sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_per_condition,
                plot_title="Composition plot colored by cell type",
                plot_subtitle=paste(
                    "Split by grouping condition;",
                    "first downsampled to", min_dataset_size,
                    "cells per dataset,",
                    "then downsampled to", min_condition_size,
                    "cells per grouping condition"
                ),
                legend_title="Cell type",
                group_by=args$target,
                split_by="condition",
                x_label="Condition",
                y_label="Cell percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_ctyp_spl_cnd", sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_per_condition,
                plot_title="Composition plot colored by grouping condition",
                plot_subtitle=paste(
                    "Split by cell type;",
                    "first downsampled to", min_dataset_size,
                    "cells per dataset,",
                    "then downsampled to", min_condition_size,
                    "cells per grouping condition"
                ),
                legend_title="Condition",
                group_by="condition",
                split_by=args$target,
                bar_position="dodge",
                x_label="Cell type",
                y_label="Cell counts",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_cnd_spl_ctyp", sep="_"),
                pdf=args$pdf
            )
            if ("Phase" %in% colnames(seurat_data@meta.data)){
                graphics$dim_plot(
                    data=downsampled_per_condition,
                    reduction=args$reduction,
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
    } else {
        if ("Phase" %in% colnames(seurat_data@meta.data)){                # the same plots but with not downsampled data
            graphics$dim_plot(
                data=seurat_data,
                reduction=args$reduction,
                plot_title="UMAP colored by cell type",
                plot_subtitle="Split by cell cycle phase",
                legend_title="Cell type",
                group_by=args$target,
                split_by="Phase",
                label=FALSE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_gr_ctyp_spl_ph", sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=seurat_data,
                plot_title="Composition plot colored by cell cycle phase",
                plot_subtitle="Split by cell type",
                legend_title="Phase",
                group_by="Phase",
                split_by=args$target,
                bar_position="dodge",
                x_label="Cell type",
                y_label="Cell counts",
                palette_colors=graphics$CC_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_ph_spl_ctyp", sep="_"),
                pdf=args$pdf
            )
        }
    }
    rm(downsampled_per_dataset, downsampled_per_condition)
    gc(verbose=FALSE)
}

## ----
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

    for (i in 1:length(args$genes)){
        current_gene <- args$genes[i]
        graphics$coverage_plot(
            data=seurat_data,
            assay="ATAC",
            region=current_gene,
            group_by=args$target,
            plot_title="ATAC fragment coverage",
            plot_subtitle=current_gene,
            idents=NULL,                                                               # to include all values from the default "new.ident" column
            cells=colnames(seurat_data),                                               # limit to only those cells that are in out seurat_data
            features=if("RNA" %in% names(seurat_data@assays)) current_gene else NULL,  # will fail if features are provided without "RNA" assay
            expression_assay="RNA",
            expression_slot="data",                                                    # use scaled counts
            extend_upstream=args$upstream,
            extend_downstream=args$downstream,
            show_annotation=TRUE,
            show_peaks=TRUE,
            show_tile=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "cvrg", current_gene, sep="_"),
            pdf=args$pdf
        )
    }
}

## ----
export_all_expression_plots <- function(seurat_data, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- args$target

    groups_count <- 1                                    # used only to adjust the height of the plot
    if (!is.null(args$splitby)){
        groups_count <- length(unique(as.vector(as.character(seurat_data@meta.data[[args$splitby]]))))
    }

    graphics$dot_plot(
        data=seurat_data,
        features=args$genes,
        plot_title="Average gene expression",
        x_label="Genes",
        y_label="Cell type",
        cluster_idents=FALSE,
        theme=args$theme,
        rootname=paste(args$output, "xpr_avg", sep="_"),
        pdf=args$pdf
    )

    graphics$vln_plot(
        data=seurat_data,
        features=args$genes,
        labels=args$genes,
        plot_title="Gene expression density",
        legend_title="Cell type",
        rotate_labels=TRUE,
        pt_size=0,
        combine_guides="collect",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "xpr_dnst", sep="_"),
        pdf=args$pdf
    )

    for (i in 1:length(args$genes)){
        current_gene <- args$genes[i]
        graphics$feature_plot(
            data=seurat_data,
            features=current_gene,
            labels=current_gene,
            reduction=args$reduction,
            plot_title="UMAP colored by gene expression",
            plot_subtitle=if (!is.null(args$splitby)) paste("Split by", args$splitby) else NULL,
            legend_title="Expression",
            label=FALSE,
            order=TRUE,
            split_by=args$splitby,
            pt_size=1,                                # need to have it fixed, otherwise dots are different on the splitted plot
            max_cutoff="q99",                         # to prevent cells with overexpressed gene from distoring the color bar
            combine_guides="collect",
            width=800,
            height=ifelse(groups_count == 2, 400, 800),
            theme=args$theme,
            rootname=paste(args$output, "xpr_per_cell", current_gene, sep="_"),
            pdf=args$pdf
        )
        graphics$expression_density_plot(
            data=seurat_data,
            features=current_gene,
            reduction=args$reduction,
            plot_title="UMAP colored by gene expression density",
            joint=FALSE,
            width=800,
            height=800,
            theme=args$theme,
            rootname=paste(args$output, "xpr_per_cell_sgnl", current_gene, sep="_"),
            pdf=args$pdf
        )
    }
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
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

    clusters_order <- levels(seurat_data@meta.data[[args$target]])
    if (is.null(clusters_order)){
        clusters_order <- unique(seurat_data@meta.data[[args$target]])
    }

    base::print(
        base::paste(
            "Checking if RNA markers for", args$target,
            "metadata column have been calculated"
        )
    )
    filtered_markers <- NULL
    base::tryCatch(
        expr = {
            filtered_markers <- seurat_data@misc$markers$RNA[[args$target]] %>%
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
                    "Failed to find RNA markers for", args$target,
                    "metadata column due to", e
                )
            )
        }
    )

    if (!is.null(filtered_markers) && (nrow(filtered_markers) > 0)){
        column_annotations <- c("Cell type")
        colnames(seurat_data@meta.data)[colnames(seurat_data@meta.data) == args$target] <- "Cell type"
        if (conditions_count > 1 && not_default_conditions){
            column_annotations <- c(column_annotations, "Condition")                           # several conditions found
            colnames(seurat_data@meta.data)[colnames(seurat_data@meta.data) == "condition"] <- "Condition"
        }
        if (datasets_count > 1){
            column_annotations <- c(column_annotations, "Dataset")                             # several datasets found
            colnames(seurat_data@meta.data)[colnames(seurat_data@meta.data) == "new.ident"] <- "Dataset"
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
            plot_title="Gene expression heatmap",
            rootname=paste(args$output, "xpr_htmp", sep="_"),
            pdf=args$pdf
        )
        io$export_data(
            filtered_markers,
            paste(args$output, "xpr_htmp.tsv", sep="_")
        )
    }
}

## ----
get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell Manual Cell Type Assignment")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "genes expression and/or chromatin accessibility information stored in the RNA",
            "and ATAC assays correspondingly."
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
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to additionally split",
            "every cluster selected with --source into smaller groups.",
            "Default: do not split"
        ),
        type="character"
    )
    parser$add_argument(
        "--reduction",
        help=paste(
            "Dimensionality reduction to be used in the generated plots. If not",
            "provided it will be automatically defined on the basis of the --source",
            "parameter as follows: rna_res.* - rnaumap, atac_res.* - atacumap,",
            "wsnn_res.* - wnnumap.",
            "Default: defined automatically"
        ),
        type="character"
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
        "--upstream",
        help=paste(
            "Number of bases to extend the genome coverage region for",
            "a specific gene upstream. Ignored if --genes or --fragments",
            "parameters are not provided. Default: 2500"
        ),
        type="integer", default=2500
    )
    parser$add_argument(
        "--downstream",
        help=paste(
            "Number of bases to extend the genome coverage region for",
            "a specific gene downstream. Ignored if --genes or --fragments",
            "parameters are not provided. Default: 2500"
        ),
        type="integer", default=2500
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
        help="Save raw counts from the RNA and/or ATAC assay(s) to h5ad file(s). Default: false",
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
            "Save Seurat data to SCope compatible loom file. Only",
            "not normalized raw counts from the RNA assay will be",
            "saved. If loaded Seurat object doesn't have RNA assay",
            "this parameter will be ignored. Default: false"
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
debug$print_info(seurat_data, args)

## ----
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

## ----
if (!is.null(args$reduction) && !(args$reduction %in% names(seurat_data@reductions))){
    print(
        paste(
            "Loaded Seurat object doesn't include the",
            "selected", args$reduction, "reduction.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (is.null(args$reduction)){
    print(
        paste(
            "Attempting to automatically define the",
            "--reduction parameter on the basis of",
            args$source
        )
    )
    reduc_names <- names(seurat_data@reductions)
    possible_suffixes <- c("wsnn", "rna", "atac")
    possible_names <- c("wnnumap", "rnaumap", "atacumap")
    for (i in 1:length(possible_suffixes)){
        if (grepl(possible_suffixes[i], args$source) && possible_names[i] %in% reduc_names){
            args$reduction <- possible_names[i]
            print(paste("--reduction was set to", args$reduction))
            break
        }
    }
    if (is.null(args$reduction)){
        print(
            paste(
                "Failed to automatically define the",
                "--reduction parameter. Exiting."
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
}

## ----
if ("Phase" %in% colnames(seurat_data@meta.data)){                                 # safety measure to have the proper order on the plots
    seurat_data@meta.data[["Phase"]] <- factor(
        seurat_data@meta.data[["Phase"]],
        levels=c("G1", "S", "G2M")
    )
}

## ----
seurat_data <- io$extend_metadata(
    seurat_data=seurat_data,
    location=args$celltypes,
    seurat_ref_column=args$source,
    meta_ref_column="cluster",
    split_by=args$splitby,                 # can be NULL if --splitby not provided
    seurat_target_columns=args$target
)
debug$print_info(seurat_data, args)

## ----
if ( (!is.null(args$fragments)) && ("ATAC" %in% names(seurat_data@assays)) ){
    print(paste("Loading ATAC fragments data from", args$fragments))
    seurat_data <- io$replace_fragments(args$fragments, seurat_data)                                             # will change the default assay to ATAC
    debug$print_info(seurat_data, args)
}

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
    if("RNA" %in% names(seurat_data@assays)){
        DefaultAssay(seurat_data) <- "RNA"                                                                       # need it for rownames to return genes
        args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]                 # with RNA assay set as default the rownames should be genes
    }
    if("ATAC" %in% names(seurat_data@assays)){
        DefaultAssay(seurat_data) <- "ATAC"                                                                      # Annotation needs the default assay to be ATAC
        args$genes <- args$genes[args$genes %in% as.vector(as.character(Annotation(seurat_data)$gene_name))]     # check if genes of interest are present in Annotation
    }
    if (length(args$genes) == 0){
        print("Neither of the provided genes of interest are present in the loaded Seurat object. Setting back to NULL")
        args$genes <- NULL                                                                                       # all genes were filtered out, so return to NULL
    } else {
        print(paste("--genes was adjusted to", paste(args$genes, collapse=", ")))
    }
}

## ----
if ( ("RNA" %in% names(seurat_data@assays)) && (!is.null(args$genes) || args$diffgenes) ){                       # we check genes or diffgenes so we don't normalize it without reason
    print("Normalizing counts in RNA assay before evaluating genes expression or identifying putative gene markers")
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
    if (!is.null(args$genes)) {
        print("Generating genes expression plots")
        export_all_expression_plots(seurat_data=seurat_data, args=args)                                          # changes default assay to RNA
    }
    if(args$diffgenes){
        print("Identifying differentially expressed genes between each pair of cell types")
        args$logfc <- args$rnalogfc                                                                              # need the proper names for get_markers
        args$minpct <- args$rnaminpct
        args$onlypos <- args$rnaonlypos
        args$testuse <- args$rnatestuse
        seurat_data <- analyses$get_markers(
            seurat_data=seurat_data,
            assay="RNA",
            group_by=args$target,
            args=args
        )
        debug$print_info(seurat_data, args)
        args <- args[names(args) %in% c("logfc", "minpct", "onlypos", "testuse") == FALSE]                       # to remove temporary added items
        io$export_markers(
            data=seurat_data,
            assay="RNA",
            markers_regex=args$target,
            location=paste0(args$output, "_gene_markers.tsv")
        )
        export_heatmaps(                                                                                         # will change default assay to RNA
            seurat_data=seurat_data,
            args=args
        )
    }
}

## ----
if ("ATAC" %in% names(seurat_data@assays)){                                                                      # no reason to check for genes and fragments or diffpeaks, because we only change the assay
    DefaultAssay(seurat_data) <- "ATAC"                                                                          # safety measure
    if(!is.null(args$genes) && !is.null(args$fragments)){
        print("Generating coverage plots")                                                                       # by now "RNA" assay, if present, will be already normalized
        export_all_coverage_plots(seurat_data=seurat_data, args=args)                                            # changes default assay to ATAC
    }
    if(args$diffpeaks){
        print("Identifying differentially accessible peaks between each pair of cell types")
        args$logfc <- args$ataclogfc                                                                             # need the proper names for get_markers
        args$minpct <- args$atacminpct
        args$onlypos <- FALSE                                                                                    # need to overwrite what was set for RNA
        args$testuse <- args$atactestuse
        seurat_data <- analyses$get_markers(                                                                # will change default assay to ATAC
            seurat_data=seurat_data,
            assay="ATAC",
            group_by=args$target,
            latent_vars="nCount_ATAC",                                                                           # to remove the influence of sequencing depth
            args=args
        )
        debug$print_info(seurat_data, args)
        args <- args[names(args) %in% c("logfc", "minpct", "onlypos", "testuse") == FALSE]                       # to remove temporary added items
        io$export_markers(
            data=seurat_data,
            assay="ATAC",
            markers_regex=args$target,
            location=paste0(args$output, "_peak_markers.tsv")
        )
    }
}

## ----
if(args$cbbuild){
    print(
        paste(
            "Reordering reductions to have", args$reduction,                               # will be shown first in UCSC Cellbrowser
            "on the first place"
        )
    )
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c(args$reduction, reduc_names[reduc_names!=args$reduction])     # we checked before that args$reduction is present
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]

    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            label_field <- base::gsub("custom_", "Custom ", args$target),
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/rna", sep="")
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            label_field <- base::gsub("custom_", "Custom ", args$target),
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/atac", sep="")
        )
    } else if ("RNA" %in% names(seurat_data@assays)){
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            label_field <- base::gsub("custom_", "Custom ", args$target),
            rootname=paste(args$output, "_cellbrowser", sep="")
        )
    } else {
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            label_field <- base::gsub("custom_", "Custom ", args$target),
            rootname=paste(args$output, "_cellbrowser", sep="")
        )
    }
}

## ----
io$export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))

## ----
if(args$h5seurat){
    io$export_h5seurat(seurat_data, paste(args$output, "_data.h5seurat", sep=""))
}

## ----
if(args$h5ad && ("RNA" %in% names(seurat_data@assays))){
    io$export_h5ad(
        data=seurat_data,
        location=paste(args$output, "_rna_counts.h5ad", sep=""),
        assay="RNA",
        slot="counts"
    )
}

## ----
if(args$h5ad && ("ATAC" %in% names(seurat_data@assays))){
    io$export_h5ad(
        data=seurat_data,
        location=paste(args$output, "_atac_counts.h5ad", sep=""),
        assay="ATAC",
        slot="counts"
    )
}

## ----
if(args$loupe && ("RNA" %in% names(seurat_data@assays))){
    ucsc$export_loupe(
        seurat_data=seurat_data,
        assay="RNA",
        active_cluster=args$target,
        rootname=paste0(args$output, "_rna_counts")
    )
}

## ----
if(args$scope && ("RNA" %in% names(seurat_data@assays))){
    io$export_scope_loom(                                                                  # we save only counts slot from the RNA assay 
        seurat_data,
        paste(args$output, "_data.loom", sep="")
    )
}