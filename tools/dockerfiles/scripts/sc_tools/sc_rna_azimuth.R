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

    graphics$feature_plot(
        data=seurat_data,
        features=c("prediction_mapping_score", "prediction_confidence_score"),
        labels=c("Prediction mapping score", "Prediction confidence score"),
        from_meta=TRUE,
        reduction="refumap",
        plot_title="UMAP, QC metrics",
        plot_subtitle="All cells",
        legend_title="Score",
        label=FALSE,
        alpha=0.5,
        pt_size=1,
        max_cutoff="q99",                                                                   # to prevent outlier cells to distort coloring
        combine_guides="collect",
        order=TRUE,
        gradient_colors=c("black", "blue", "lightblue"),
        color_scales=c(0, 1),
        color_limits=c(0, 1),
        theme=args$theme,
        rootname=paste(args$output, "umap_qc_mtrcs", sep="_"),
        pdf=args$pdf
    )

    graphics$vln_plot(                                                           # won't fail even if some columns are not present
        data=seurat_data,
        features=c(
            "nCount_RNA", "nFeature_RNA", "mito_percentage",
            "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "nFeature_ATAC", "frip", "blacklist_fraction",
            "prediction_mapping_score", "prediction_confidence_score"
        ),
        labels=c(
            "RNA reads", "Genes", "Mitochondrial %",
            "ATAC fragments in peaks", "TSS enrichment score", "Nucl. signal", "Peaks", "FRiP", "Bl. regions",
            "Prediction mapping score", "Prediction confidence score"
        ),
        scale_y_log10=c(
            TRUE, TRUE, FALSE,
            TRUE, FALSE, FALSE, TRUE, FALSE, FALSE,
            FALSE, FALSE
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

    highlight_rows <- which(seurat_data@meta.data$rna_doublets == "doublet")
    if ("ATAC" %in% names(seurat_data@assays)){
        highlight_rows <- which(seurat_data@meta.data$rna_doublets == "doublet" | seurat_data@meta.data$atac_doublets == "doublet")
    }

    graphics$geom_point_plot(
        data=seurat_data@meta.data,
        x_axis="nCount_RNA",
        y_axis="nFeature_RNA",
        split_by=args$target,
        color_by="mito_percentage",
        highlight_rows=highlight_rows,
        gradient_colors=c("lightslateblue", "orange", "red"),
        color_limits=c(0, 100),
        color_break=ceiling(max(seurat_data@meta.data$mito_percentage)),
        legend_title="Mitochondrial %",
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
        color_by="mito_percentage",
        highlight_rows=highlight_rows,
        gradient_colors=c("lightslateblue", "orange", "red"),
        color_limits=c(0, 100),
        color_break=ceiling(max(seurat_data@meta.data$mito_percentage)),
        legend_title="Mitochondrial %",
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

    if ("ATAC" %in% names(seurat_data@assays)){
        graphics$geom_point_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_ATAC",
            y_axis="TSS.enrichment",
            split_by=args$target,
            color_by="mito_percentage",
            highlight_rows=highlight_rows,
            gradient_colors=c("lightslateblue", "orange", "red"),
            color_limits=c(0, 100),
            color_break=ceiling(max(seurat_data@meta.data$mito_percentage)),
            legend_title="Mitochondrial %",
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

        graphics$geom_point_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_ATAC",
            y_axis="nCount_RNA",
            split_by=args$target,
            color_by="mito_percentage",
            highlight_rows=highlight_rows,
            gradient_colors=c("lightslateblue", "orange", "red"),
            color_limits=c(0, 100),
            color_break=ceiling(max(seurat_data@meta.data$mito_percentage)),
            legend_title="Mitochondrial %",
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

            if (nrow(seurat_data@meta.data[seurat_data@meta.data$rna_doublets == "doublet", ]) > 0){
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

    graphics$dim_plot(
        data=seurat_data,
        reduction="refumap",
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
            reduction="refumap",
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
                reduction="refumap",
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
                reduction="refumap",
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
                reduction="refumap",
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
                    reduction="refumap",
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
                reduction="refumap",
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
            features=current_gene,
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
            reduction="refumap",
            plot_title="UMAP colored by gene expression",
            legend_title="Expression",
            label=FALSE,
            order=TRUE,
            pt_size=1,                                # need to have it fixed, otherwise dots are different on the splitted plot
            max_cutoff="q99",                         # to prevent cells with overexpressed gene from distoring the color bar
            combine_guides="collect",
            width=800,
            height=800,
            theme=args$theme,
            rootname=paste(args$output, "xpr_per_cell", current_gene, sep="_"),
            pdf=args$pdf
        )
        graphics$expression_density_plot(
            data=seurat_data,
            features=current_gene,
            reduction="refumap",
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
    parser <- ArgumentParser(description="Single-Cell RNA-Seq Reference Mapping")     # it can be used after any pipeline that produces Seurat object with RNA assay
    parser$add_argument(                                                              # if ATAC assay present it remains unchanged, but correspondent plots are produced
        "--query",                                                                    # the pipeline adds refumap redaction and several metadata columns
        help=paste(
            "Path to the RDS file to load the query Seurat object from. This file",
            "should include genes expression information stored in the RNA assay",
            "and, optionally, chromatin accessibility information stored in the",
            "ATAC assay. The later is used only for plots."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--reference",
        help=paste(
            "Path to the RDS file to load the reference Seurat object",
            "from. This file can be downloaded as ref.Rds from the",
            "https://azimuth.hubmapconsortium.org/references/"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--annoyidx",
        help=paste(
            "Path to the annoy index file generated for the provided reference",
            "Seurat object. This file can be downloaded as idx.annoy from the",
            "https://azimuth.hubmapconsortium.org/references/"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--source",
        help=paste(
            "Column from the metadata of the reference Seurat object",
            "to select the reference annotations."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--diffgenes",
        help=paste(
            "Identify differentially expressed genes (putative",
            "gene markers) for the predicted cell types.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--diffpeaks",
        help=paste(
            "Identify differentially accessible peaks for the predicted",
            "cell types. Ignored if the query Seurat object doesn't include",
            "chromatin accessibility information stored in the ATAC assay.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--rnalogfc",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "on average have a log fold change difference in the expression between",
            "every tested pair of the predicted cell types not lower than this value.",
            "Ignored if '--diffgenes is not set.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--rnaminpct",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "are detected in not lower than this fraction of cells in either of the",
            "two tested predicted cell types. Ignored if '--diffgenes' is not set.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--rnaonlypos",
        help=paste(
            "For putative gene markers identification return only",
            "upregulated markers. Ignored if '--diffgenes' is not set.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--rnatestuse",
        help=paste(
            "Statistical test to use for putative gene markers",
            "identification. Ignored if '--diffgenes' is not set.",
            "Default: wilcox"
        ),
        type="character", default="wilcox",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
    )
    parser$add_argument(
        "--ataclogfc",
        help=paste(
            "For differentially accessible peaks identification include",
            "only those peaks that on average have a log fold change",
            "difference in the chromatin accessibility between every",
            "tested pair of the predicted cell types not lower than this",
            "value. Ignored if '--diffpeaks is not set or if the query",
            "Seurat object doesn't include ATAC assay.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--atacminpct",
        help=paste(
            "For differentially accessible peaks identification include",
            "only those peaks that are detected in not lower than this",
            "fraction of cells in either of the two tested predicted",
            "cell types. Ignored if '--diffpeaks' is not set or if the",
            "query Seurat object doesn't include ATAC assay.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--atactestuse",
        help=paste(
            "Statistical test to use for differentially accessible peaks",
            "identification. Ignored if '--diffpeaks' is not set or if",
            "the query Seurat object doesn't include ATAC assay.",
            "Default: LR"
        ),
        type="character", default="LR",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment",
            "used in the query Seurat object. File should be saved",
            "in TSV format with tbi-index file. Ignored if the",
            "query Seurat object doesn't include ATAC assay."
        ),
        type="character"
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to build gene expression and Tn5",
            "insertion frequency plots for the nearest peaks. To",
            "build Tn5 insertion frequency plots for the nearest",
            "peaks the query Seurat object should include ATAC assay",
            "as well as the --fragments file should be provided.",
            "Default: None"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--upstream",
        help=paste(
            "Number of bases to extend the genome coverage region for",
            "a specific gene upstream. Ignored if --genes or --fragments",
            "parameters are not provided or when the query Seurat object",
            "doesn't include ATAC assay.",
            "Default: 2500"
        ),
        type="integer", default=2500
    )
    parser$add_argument(
        "--downstream",
        help=paste(
            "Number of bases to extend the genome coverage region for",
            "a specific gene downstream. Ignored if --genes or --fragments",
            "parameters are not provided or when the query Seurat object",
            "doesn't include ATAC assay.",
            "Default: 2500"
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
        help=paste(
            "Save raw counts from the RNA and ATAC (if present)",
            "assays to h5ad file(s).",
            "Default: false"
        ),
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
            "Save Seurat data to SCope compatible loom file. Only not",
            "normalized raw counts from the RNA assay will be saved.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--tmpdir",
        help=paste(
            "Directory to keep temporary files. Default: either /tmp",
            "or defined by the environment variables TMPDIR, TMP, TEMP."
        ),
        type="character", default=tempdir()
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
print(paste("Loading query Seurat data from", args$query))
seurat_data <- readRDS(args$query)
debug$print_info(seurat_data, args)

## ----
if (!("RNA" %in% names(seurat_data@assays))){
    print(
        paste(
            "Loaded query Seurat object doesn't include",
            "the required RNA assay. Exiting."
        )
    )
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
ref_dir <- base::file.path(args$tmpdir, "reference")
if (!base::dir.exists(ref_dir)) {
    base::dir.create(ref_dir, recursive=TRUE)
}
print(paste("Copying reference RDS and annoy index files to the", ref_dir))
base::file.copy(args$reference, base::file.path(ref_dir, "ref.Rds"))               # only these names are recognized by Azimuth
base::file.copy(args$annoyidx, base::file.path(ref_dir, "idx.annoy"))

## ----
seurat_data <- analyses$rna_reference_map(                                         # sets RNA as a default assay, and new.ident as the current identities
    seurat_data=seurat_data,                                                       # adds prediction_confidence_score, prediction_mapping_score, and prediction_cell_type
    reference_dir=ref_dir,                                                         # metadata column, as well as the refumap reduction
    args=args
)
args$target <- "prediction_cell_type"                                               # for easy access and consistency in code between the scripts
debug$print_info(seurat_data, args)

## ----
if ( (!is.null(args$fragments)) && ("ATAC" %in% names(seurat_data@assays)) ){
    print(paste("Loading ATAC fragments data from", args$fragments))
    seurat_data <- io$replace_fragments(args$fragments, seurat_data)               # will change the default assay to ATAC
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
    DefaultAssay(seurat_data) <- "RNA"                                                                           # need it for rownames to return genes
    args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]                     # with RNA assay set as default the rownames should be genes
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
if (!is.null(args$genes) || args$diffgenes){                                                                     # we check genes or diffgenes so we don't normalize it without reason
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
            "Reordering reductions to have refumap",
            "on the first place"
        )
    )
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c("refumap", reduc_names[reduc_names!="refumap"])
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]

    if ("ATAC" %in% names(seurat_data@assays)){
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            label_field <- "Prediction cell type",
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/rna", sep="")
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            label_field <- "Prediction cell type",
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/atac", sep="")
        )
    } else {
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            label_field <- "Prediction cell type",
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
if(args$h5ad){
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
if(args$loupe){
    ucsc$export_loupe(
        seurat_data=seurat_data,
        assay="RNA",
        active_cluster=args$target,
        rootname=paste0(args$output, "_rna_counts")
    )
}

## ----
if(args$scope){
    io$export_scope_loom(                                                                  # we save only counts slot from the RNA assay 
        seurat_data,
        paste(args$output, "_data.loom", sep="")
    )
}