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
suppressMessages(library(argparse))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(filter <- modules::use(file.path(HERE, "modules/filter.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))

## ----
export_all_qc_plots <- function(seurat_data, suffix, args, macs2_peaks=FALSE){
    Idents(seurat_data) <- "new.ident"                                                                # safety measure
    peak_type <- ifelse(macs2_peaks, "- MACS2", "- 10x")
    selected_features <- c(
        "nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi",
        "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "nFeature_ATAC", "frip", "blacklist_fraction"
    )
    selected_labels <- c(
        "RNA reads", "Genes", "Mitochondrial %", "Novelty score",
        paste(c("ATAC fragments\nin peaks", "TSS enrichment\nscore", "Nucl. signal", "Peaks", "FRiP", "Bl. regions"), peak_type)
    )
    selected_scales <- c(
        TRUE, TRUE, FALSE, FALSE,
        TRUE, FALSE, FALSE, TRUE, FALSE, FALSE
    )
    datasets_count <- length(unique(as.vector(as.character(seurat_data@meta.data$new.ident))))
    conditions_count <- length(unique(as.vector(as.character(seurat_data@meta.data$condition))))
    not_default_conditions <- all(
        as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))
    )

    qc_metrics_pca <- qc$qc_metrics_pca(
        seurat_data=seurat_data,
        qc_columns=selected_features,
        qc_labels=selected_labels,
        orq_transform=TRUE
    )

    graphics$pca_plot(
        pca_data=qc_metrics_pca,
        pcs=c(1, 2),
        plot_title="QC metrics PCA",
        plot_subtitle=paste0(
            graphics$expand_qc_suffix(suffix),
            "; ", "PC1/PC2"
        ),
        legend_title="QC metrics",
        color_by="labels",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, suffix, paste(c(1, 2) ,collapse="_"), "qc_mtrcs_pca", sep="_"),
        pdf=args$pdf
    )
    graphics$pca_plot(
        pca_data=qc_metrics_pca,
        pcs=c(2, 3),
        plot_title="QC metrics PCA",
        plot_subtitle=paste0(
            graphics$expand_qc_suffix(suffix),
            "; ", "PC2/PC3"
        ),
        legend_title="QC metrics",
        color_by="labels",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, suffix, paste(c(2, 3) ,collapse="_"), "qc_mtrcs_pca", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_bar_plot(
        data=seurat_data@meta.data,
        x_axis="new.ident",
        color_by="new.ident",
        x_label="Dataset",
        y_label="Cells",
        legend_title="Dataset",
        plot_title="Number of cells per dataset",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        width=ifelse(datasets_count > 1, 1200, 400),
        rootname=paste(args$output, suffix, "cell_cnts", sep="_"),
        pdf=args$pdf
    )

    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nCount_RNA",
        group_by="new.ident",
        split_by="new.ident",
        x_left_intercept=args$minumis,
        x_label="RNA reads per cell",
        y_label="Distribution",
        legend_title="Dataset",
        plot_title="Distribution of RNA reads per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        scale_x_log10=TRUE,
        show_zoomed=FALSE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "umi_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nFeature_RNA",
        group_by="new.ident",
        split_by="new.ident",
        x_left_intercept=args$mingenes,
        x_right_intercept=args$maxgenes,
        x_label="Genes per cell",
        y_label="Distribution",
        legend_title="Dataset",
        plot_title="Distribution of genes per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        scale_x_log10=TRUE,
        show_zoomed=FALSE,
        show_ranked=TRUE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, suffix, "gene_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_point_plot(
        data=seurat_data@meta.data,
        x_axis="nCount_RNA",
        y_axis="nFeature_RNA",
        split_by="new.ident",
        x_left_intercept=args$minumis,
        y_low_intercept=args$mingenes,
        y_high_intercept=args$maxgenes,
        color_by="mito_percentage",
        highlight_rows=which(seurat_data@meta.data$rna_doublets == "doublet" | seurat_data@meta.data$atac_doublets == "doublet"),    # need a single |
        gradient_colors=c("lightslateblue", "orange", "red"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        x_label="RNA reads per cell",
        y_label="Genes per cell",
        legend_title="Mitochondrial %",
        plot_title="Genes vs RNA reads per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        show_lm=TRUE,
        show_density=TRUE,
        density_bins=4,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        width=ifelse(datasets_count > 1, 1200, 600),
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "gene_umi", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_point_plot(
        data=seurat_data@meta.data,
        x_axis="mito_percentage",
        y_axis="nCount_RNA",
        split_by="new.ident",
        x_left_intercept=args$maxmt,
        y_low_intercept=args$minumis,
        color_by="mito_percentage",
        highlight_rows=which(seurat_data@meta.data$rna_doublets == "doublet" | seurat_data@meta.data$atac_doublets == "doublet"),    # need a single |
        gradient_colors=c("lightslateblue", "orange", "red"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        x_label="Mitochondrial %",
        y_label="RNA reads per cell",
        legend_title="Mitochondrial %",
        plot_title="RNA reads vs mitochondrial % per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        show_lm=FALSE,
        show_density=TRUE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        width=ifelse(datasets_count > 1, 1200, 600),
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "umi_mito", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="mito_percentage",
        group_by="new.ident",
        split_by="new.ident",
        x_right_intercept=args$maxmt,
        x_label="Percentage of RNA reads mapped to mitochondrial genes per cell",
        y_label="Distribution",
        legend_title="Dataset",
        plot_title="Distribution of RNA reads mapped to mitochondrial genes per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        show_zoomed=FALSE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "mito_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="log10_gene_per_log10_umi",
        group_by="new.ident",
        split_by="new.ident",
        x_left_intercept=args$minnovelty,
        x_label="log10 Genes / log10 UMI per cell",
        y_label="Distribution",
        legend_title="Dataset",
        plot_title="Distribution of novelty score per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        show_zoomed=FALSE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "nvlt_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nCount_ATAC",
        group_by="new.ident",
        split_by="new.ident",
        x_left_intercept=args$minfragments,
        x_label="ATAC fragments in peaks per cell",
        y_label="Distribution",
        legend_title="Dataset",
        plot_title="Distribution of ATAC fragments in peaks per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        scale_x_log10=TRUE,
        show_zoomed=FALSE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "frgm_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nFeature_ATAC",
        group_by="new.ident",
        split_by="new.ident",
        x_label="Peaks per cell",
        y_label="Distribution",
        legend_title="Dataset",
        plot_title="Distribution of peaks per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        scale_x_log10=FALSE,
        show_zoomed=FALSE,
        show_ranked=TRUE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, suffix, "peak_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="blacklist_fraction",
        group_by="new.ident",
        split_by="new.ident",
        x_right_intercept=args$maxblacklist,
        x_label="Fraction of ATAC fragments within genomic blacklist regions per cell",
        y_label="Distribution",
        legend_title="Dataset",
        plot_title="Distribution of ATAC fragments within genomic blacklist regions per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        scale_x_log10=FALSE,
        show_zoomed=FALSE,
        show_ranked=TRUE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, suffix, "blck_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_point_plot(
        data=seurat_data@meta.data,
        split_by="new.ident",
        x_axis="nCount_ATAC",
        x_label="ATAC fragments in peaks per cell",
        y_axis="nCount_RNA",
        y_label="RNA reads per cell",
        x_left_intercept=args$minfragments,
        y_low_intercept=args$minumis,
        alpha_intercept=1,
        color_by="mito_percentage",
        highlight_rows=which(seurat_data@meta.data$rna_doublets == "doublet" | seurat_data@meta.data$atac_doublets == "doublet"),    # need a single |
        gradient_colors=c("lightslateblue", "orange", "red"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        legend_title="Mitochondrial %",
        plot_title="RNA reads vs ATAC fragments in peaks per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        show_density=TRUE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        width=ifelse(datasets_count > 1, 1200, 600),
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "rna_atac_cnts", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_point_plot(
        data=seurat_data@meta.data,
        split_by="new.ident",
        x_axis="nCount_ATAC",
        x_label="ATAC fragments in peaks per cell",
        y_axis="TSS.enrichment",
        y_label="TSS enrichment score",
        x_left_intercept=args$minfragments,
        y_low_intercept=args$mintssenrich,
        alpha_intercept=1,
        color_by="mito_percentage",
        highlight_rows=which(seurat_data@meta.data$rna_doublets == "doublet" | seurat_data@meta.data$atac_doublets == "doublet"),    # need a single |
        gradient_colors=c("lightslateblue", "orange", "red"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        legend_title="Mitochondrial %",
        plot_title="TSS enrichment score vs ATAC fragments in peaks per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        scale_x_log10=TRUE,
        scale_y_log10=FALSE,
        show_density=TRUE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        width=ifelse(datasets_count > 1, 1200, 600),
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "tss_frgm", sep="_"),
        pdf=args$pdf
    )
    graphics$vln_plot(
        data=seurat_data,
        features=selected_features,
        labels=selected_labels,
        scale_y_log10=selected_scales,
        from_meta=TRUE,
        show_box_plots=TRUE,
        plot_title="Distribution of QC metrics per cell",
        plot_subtitle=graphics$expand_qc_suffix(suffix),
        legend_title="Dataset",
        hide_x_text=TRUE,
        pt_size=0,
        combine_guides="collect",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, suffix, "qc_mtrcs_dnst", sep="_"),
        pdf=args$pdf
    )

    if (nrow(seurat_data@meta.data[seurat_data@meta.data$rna_doublets == "doublet", ]) > 0){      # show plot only if not all rna doublets have been removed
        graphics$composition_plot(
            data=seurat_data,
            plot_title="Percentage of RNA doublets",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                sep="; "
            ),
            legend_title="Dataset",
            group_by="rna_doublets",
            split_by="new.ident",
            x_label="Dataset",
            y_label="Cell percentage",
            palette_colors=c("#00AEAE", "#0BFFFF"),
            theme=args$theme,
            width=ifelse(datasets_count > 1, 1200, 400),
            rootname=paste(args$output, suffix, "rnadbl", sep="_"),
            pdf=args$pdf
        )
    }

    if (nrow(seurat_data@meta.data[seurat_data@meta.data$atac_doublets == "doublet", ]) > 0){      # show plot only if not all atac doublets have been removed
        graphics$composition_plot(
            data=seurat_data,
            plot_title="Percentage of ATAC doublets",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                sep="; "
            ),
            legend_title="Dataset",
            group_by="atac_doublets",
            split_by="new.ident",
            x_label="Dataset",
            y_label="Cell percentage",
            palette_colors=c("#00DCDC", "#0BFFFF"),
            theme=args$theme,
            width=ifelse(datasets_count > 1, 1200, 400),
            rootname=paste(args$output, suffix, "atacdbl", sep="_"),
            pdf=args$pdf
        )
    }

    if (
        nrow(seurat_data@meta.data[seurat_data@meta.data$rna_doublets == "doublet", ]) > 0 &&             # no reason to show this plot if not both doublets are present
        nrow(seurat_data@meta.data[seurat_data@meta.data$atac_doublets == "doublet", ]) > 0
    ){
        seurat_data@meta.data <- seurat_data@meta.data %>%                                                # temporary column to see doublets overlap between RNA and ATAC
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
                                         levels=c("Only RNA", "RNA and ATAC", "Only ATAC", "Singlet")     # it sage to set all 4 levels even if some are removed from data
                                     )
                                 )
        graphics$composition_plot(
            data=seurat_data,
            plot_title="Percentage of RNA and ATAC doublets",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                sep="; "
            ),
            legend_title="Dataset",
            group_by="doublets_overlap",
            split_by="new.ident",
            x_label="Dataset",
            y_label="Cell percentage",
            palette_colors=c("#00AEAE", "#008080", "#00DCDC", "#0BFFFF"),
            theme=args$theme,
            width=ifelse(datasets_count > 1, 1200, 400),
            rootname=paste(args$output, suffix, "vrlpdbl", sep="_"),
            pdf=args$pdf
        )
    }

    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "ATAC"
    graphics$tss_plot(
        data=seurat_data,
        split_by="new.ident",
        group_by_value=args$mintssenrich,
        combine_guides="collect",
        plot_title="Signal enrichment around TSS",
        plot_subtitle=paste(
            graphics$expand_qc_suffix(suffix),
            "split by the minimum TSS enrichment score threshold",
            sep="; "
        ),
        theme=args$theme,
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "tss_nrch", sep="_"),
        pdf=args$pdf
    )
    graphics$fragments_hist(
        data=seurat_data,
        split_by="new.ident",
        group_by_value=args$maxnuclsignal,
        combine_guides="collect",
        plot_title="Histogram of ATAC fragment length",
        plot_subtitle=paste(
            graphics$expand_qc_suffix(suffix),
            "split by the maximum nucleosome signal threshold",
            sep="; "
        ),
        theme=args$theme,
        height=ifelse(datasets_count > 1, 800, 400),
        rootname=paste(args$output, suffix, "frgm_hist", sep="_"),
        pdf=args$pdf
    )
    DefaultAssay(seurat_data) <- backup_assay

    if (conditions_count > 1 && not_default_conditions){
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_RNA",
            group_by="new.ident",
            split_by="condition",
            x_left_intercept=args$minumis,
            x_label="RNA reads per cell",
            y_label="Distribution",
            legend_title="Dataset",
            plot_title="Distribution of RNA reads per cell",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                "split by grouping condition",
                sep="; "
            ),
            scale_x_log10=TRUE,
            show_zoomed=FALSE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, suffix, "umi_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nFeature_RNA",
            group_by="new.ident",
            split_by="condition",
            x_left_intercept=args$mingenes,
            x_right_intercept=args$maxgenes,
            x_label="Genes per cell",
            y_label="Distribution",
            legend_title="Dataset",
            plot_title="Distribution of genes per cell",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                "split by grouping condition",
                sep="; "
            ),
            scale_x_log10=TRUE,
            show_zoomed=FALSE,
            show_ranked=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, suffix, "gene_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="mito_percentage",
            group_by="new.ident",
            split_by="condition",
            x_right_intercept=args$maxmt,
            x_label="Percentage of RNA reads mapped to mitochondrial genes per cell",
            y_label="Distribution",
            legend_title="Dataset",
            plot_title="Distribution of RNA reads mapped to mitochondrial genes per cell",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                "split by grouping condition",
                sep="; "
            ),
            show_zoomed=FALSE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, suffix, "mito_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="log10_gene_per_log10_umi",
            group_by="new.ident",
            split_by="condition",
            x_left_intercept=args$minnovelty,
            x_label="log10 Genes / log10 UMI per cell",
            y_label="Distribution",
            legend_title="Dataset",
            plot_title="Distribution of novelty score per cell",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                "split by grouping condition",
                sep="; "
            ),
            show_zoomed=FALSE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, suffix, "nvlt_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_ATAC",
            group_by="new.ident",
            split_by="condition",
            x_left_intercept=args$minfragments,
            x_label="ATAC fragments in peaks per cell",
            y_label="Distribution",
            legend_title="Dataset",
            plot_title="Distribution of ATAC fragments in peaks per cell",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                "split by grouping condition",
                sep="; "
            ),
            scale_x_log10=TRUE,
            show_zoomed=FALSE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, suffix, "frgm_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nFeature_ATAC",
            group_by="new.ident",
            split_by="condition",
            x_label="Peaks per cell",
            y_label="Distribution",
            legend_title="Dataset",
            plot_title="Distribution of peaks per cell",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                "split by grouping condition",
                sep="; "
            ),
            scale_x_log10=FALSE,
            show_zoomed=FALSE,
            show_ranked=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, suffix, "peak_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="blacklist_fraction",
            group_by="new.ident",
            split_by="condition",
            x_right_intercept=args$maxblacklist,
            x_label="Fraction of ATAC fragments within genomic blacklist regions per cell",
            y_label="Distribution",
            legend_title="Dataset",
            plot_title="Distribution of ATAC fragments within genomic blacklist regions per cell",
            plot_subtitle=paste(
                graphics$expand_qc_suffix(suffix),
                "split by grouping condition",
                sep="; "
            ),
            scale_x_log10=FALSE,
            show_zoomed=FALSE,
            show_ranked=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, suffix, "blck_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
    }
}

## ----
get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell Multiome ATAC-Seq and RNA-Seq Filtering Analysis")
    parser$add_argument(
        "--mex",
        help=paste(
            "Path to the folder with feature-barcode matrix from Cell Ranger ARC Count/Aggregate",
            "experiment in MEX format. The rows consist of all the genes and peaks concatenated",
            "together and the columns are restricted to those barcodes that are identified as cells."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--identity",
        help=paste(
            "Path to the metadata TSV/CSV file to set the datasets identities, if '--mex' points to",
            "the Cell Ranger ARC Aggregate outputs. The aggr.csv file can be used."
        ),
        type="character"
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment observed in the experiment in TSV",
            "format. Tbi-index file is required."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--annotations",
        help="Path to the genome annotation file in GTF format",
        type="character", required="True"
    )
    parser$add_argument(
        "--seqinfo",
        help="Path to the headerless chromosome length file in TSV format",
        type="character", required="True"
    )
    parser$add_argument(
        "--grouping",
        help=paste(
            "Path to the TSV/CSV file to define datasets grouping. First column - 'library_id'",
            "with the values and order that correspond to the 'library_id' column from the",
            "'--identity' file, second column 'condition'.",
            "Default: each dataset is assigned to its own group."
        ),
        type="character"
    )
    parser$add_argument(
        "--blacklist",
        help="Path to the optional BED file with the genomic blacklist regions.",
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
        "--rnamincells",
        help=paste(
            "Include only genes detected in at least this many cells.",
            "Default: 5 (applied to all datasets)"
        ),
        type="integer", default=5
    )
    parser$add_argument(
        "--mingenes",
        help=paste(
            "Include cells where at least this many genes are detected. If multiple values",
            "provided, each of them will be applied to the correspondent dataset from the",
            "'--mex' input based on the '--identity' file. Any 0 will be replaced with the",
            "auto-estimated threshold (median - 2.5 * MAD) calculated per dataset.",
            "Default: 250 (applied to all datasets)"
        ),
        type="integer", default=250, nargs="*"
    )
    parser$add_argument(
        "--maxgenes",
        help=paste(
            "Include cells with the number of genes not bigger than this value. If multiple",
            "values provided, each of them will be applied to the correspondent dataset from",
            "the '--mex' input based on the '--identity' file. Any 0 will be replaced with the",
            "auto-estimated threshold (median + 5 * MAD) calculated per dataset.",
            "Default: 5000 (applied to all datasets)"
        ),
        type="integer", default=5000, nargs="*"
    )
    parser$add_argument(
        "--minumis",
        help=paste(
            "Include cells where at least this many RNA reads are detected. If multiple",
            "values provided, each of them will be applied to the correspondent dataset from",
            "the '--mex' input based on the '--identity' file. Any 0 will be replaced with the",
            "auto-estimated threshold (median - 2.5 * MAD) calculated per dataset.",
            "Default: 500 (applied to all datasets)"
        ),
        type="integer", default=500, nargs="*"
    )
    parser$add_argument(
        "--mitopattern",
        help=paste(
            "Regex pattern to identify mitochondrial genes.",
            "Default: '^mt-|^MT-'"
        ),
        type="character", default="^mt-|^MT-"
    )
    parser$add_argument(
        "--maxmt",
        help=paste(
            "Include cells with the percentage of RNA reads mapped",
            "to mitochondrial genes not bigger than this value.",
            "Set to 0 for using an auto-estimated threshold equal to",
            "the maximum among (median + 2 * MAD) values calculated per",
            "dataset. Default: 5 (applied to all datasets)"
        ),
        type="double", default=5
    )
    parser$add_argument(
        "--minnovelty",
        help=paste(
            "Include cells with the novelty score not lower than this value, calculated for",
            "as log10(genes)/log10(UMI) for RNA assay. If multiple values provided, each of them will",
            "be applied to the correspondent dataset from the '--mex' input based on the",
            "'--identity' file.",
            "Default: 0.8 (applied to all datasets)"
        ),
        type="double", default=0.8, nargs="*"
    )
    parser$add_argument(
        "--atacmincells",
        help=paste(
            "Include only peaks detected in at least this many cells.",
            "Default: 5 (applied to all datasets)"
        ),
        type="integer", default=5
    )
    parser$add_argument(                                                                    # when loading Cell Ranger data we have cut sites
        "--minfragments",                                                                   # per peak so we recalculate it to have fragments
        help=paste(
            "Include cells where at least this many ATAC fragments in peaks",
            "are detected. If multiple values provided, each of them will be",
            "applied to the correspondent dataset from the '--mex' input",
            "based on the '--identity' file. Any 0 will be replaced with the",
            "auto-estimated threshold (median - 2.5 * MAD) calculated per dataset.",
            "Default: 1000 (applied to all datasets)"
        ),
        type="integer", default=1000, nargs="*"
    )
    parser$add_argument(
        "--maxnuclsignal",
        help=paste(
            "Include cells with the nucleosome signal not bigger than this value.",
            "Nucleosome signal quantifies the approximate ratio of mononucleosomal",
            "to nucleosome-free ATAC fragments. If multiple values provided, each of",
            "them will be applied to the correspondent dataset from the '--mex' input",
            "based on the '--identity' file.",
            "Default: 4 (applied to all datasets)"
        ),
        type="double", default=4, nargs="*"
    )
    parser$add_argument(
        "--mintssenrich",
        help=paste(
            "Include cells with the TSS enrichment score not lower than this value.",
            "Score is calculated based on the ratio of ATAC fragments centered at the TSS",
            "to ATAC fragments in TSS-flanking regions. If multiple values provided, each",
            "of them will be applied to the correspondent dataset from the '--mex' input",
            "based on the '--identity' file.",
            "Default: 2 (applied to all datasets)"
        ),
        type="double", default=2, nargs="*"
    )
    parser$add_argument(
        "--minfrip",
        help=paste(
            "Include cells with the FRiP not lower than this",
            "value. FRiP is calculated for ATAC fragments.",
            "Default: 0.15 (applied to all datasets)"
        ),
        type="double", default=0.15
    )
    parser$add_argument(
        "--maxblacklist",
        help=paste(
            "Include cells with the fraction of ATAC fragments in genomic blacklist regions",
            "not bigger than this value. If multiple values provided, each of them",
            "will be applied to the correspondent dataset from the '--mex' input based",
            "on the '--identity' file.",
            "Default: 0.05 (applied to all datasets)"
        ),
        type="double", default=0.05, nargs="*"
    )
    parser$add_argument(
        "--callby",
        help=paste(
            "Replace Cell Ranger ARC peaks with MACS2 peaks called for cells grouped by",
            "the column from the optionally provided --barcodes file. If --barcodes file",
            "was not provided MACS2 peaks can be still called per dataset by setting --callby",
            "to new.ident. Peaks are called only after applying all RNA related thresholds,",
            "maximum nucleosome signal, and minimum TSS enrichment scores filters.",
            "Default: do not call peaks"
        ),
        type="character"
    )
    parser$add_argument(
        "--qvalue",
        help=paste(
            "Minimum FDR (q-value) cutoff for MACS2 peak detection.",
            "Ignored if --callby is not provided. Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--removedoublets",
        help=paste(
            "Remove cells that were identified as doublets. For RNA assay cells",
            "with UMI < 200 will not be evaluated. Default: do not remove doublets"
        ),
        type="character",
        choices=c("union", "onlyrna", "onlyatac", "intersect")
    )
    parser$add_argument(
        "--rnadbr",
        help=paste(
            "Expected RNA doublet rate. Default: 1 percent per thousand",
            "cells captured with 10x genomics"
        ),
        type="double"
    )
    parser$add_argument(
        "--rnadbrsd",
        help=paste(
            "Uncertainty range in the RNA doublet rate, interpreted as",
            "a +/- around the value provided in --rnadbr. Set to 0 to",
            "disable. Set to 1 to make the threshold depend entirely",
            "on the misclassification rate. Default: 40 percents of the",
            "value provided in --rnadbr"
        ),
        type="double"
    )
    parser$add_argument(
        "--atacdbr",
        help=paste(
            "Expected ATAC doublet rate. Default: 1 percent per thousand",
            "cells captured with 10x genomics"
        ),
        type="double"
    )
    parser$add_argument(
        "--atacdbrsd",
        help=paste(
            "Uncertainty range in the ATAC doublet rate, interpreted as",
            "a +/- around the value provided in --atacdbr. Set to 0 to",
            "disable. Set to 1 to make the threshold depend entirely",
            "on the misclassification rate. Default: 40 percents of the",
            "value provided in --atacdbr"
        ),
        type="double"
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
        help="Save raw counts from the RNA and ATAC assays to h5ad files. Default: false",
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
        "--tmpdir",
        help=paste(
            "Directory to keep temporary files. Default: either /tmp",
            "or defined by environment variables TMPDIR, TMP, TEMP."
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
            "when using multiple '--cpus'.",
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
print(paste("Loading datasets identities from", args$identity))
cell_identity_data <- io$load_cell_identity_data(args$identity)                          # identities are always prepended with letters to keep the order
io$export_data(
    cell_identity_data[, "library_id", drop=FALSE],                                      # we need only the first column with the prepended letters
    paste(args$output, "_meta.tsv", sep="")
)

## ----
print(paste("Loading datasets grouping from", args$grouping))
grouping_data <- io$load_grouping_data(args$grouping, cell_identity_data)

## ----
print(paste("Loading genomic blacklist regions from", args$blacklist))
blacklist_data <- io$load_blacklist_data(args$blacklist)

## ----
print(paste("Loading chromosome length information from", args$seqinfo))
seqinfo_data <- io$load_seqinfo_data(args$seqinfo)

## ----
print(paste("Loading genome annotation data from", args$annotations))
annotation_data <- io$load_annotation_data(args$annotations)

## ----
print(paste("Loading gene/peak-barcode matrices from", args$mex))
print(paste("Loading ATAC fragments from", args$fragments))
seurat_data <- io$load_10x_multiome_data(                                                # identities are set to the "new.ident" column
    args=args,
    cell_identity_data=cell_identity_data,
    grouping_data=grouping_data,
    seqinfo_data=seqinfo_data,
    annotation_data=annotation_data
)
debug$print_info(seurat_data, args)

## ----
seurat_data <- qc$estimate_doublets(                                                    # we want to search for doublets before we applied any filters
    seurat_data=seurat_data,
    assay="RNA",
    target_column="rna_doublets",
    dbl_rate=args$rnadbr,
    dbl_rate_sd=args$rnadbrsd
)

## ----
seurat_data <- qc$estimate_doublets(                                                    # we want to search for doublets before we applied any filters
    seurat_data=seurat_data,
    assay="ATAC",
    target_column="atac_doublets",
    dbl_rate=args$atacdbr,
    dbl_rate_sd=args$atacdbrsd
)
debug$print_info(seurat_data, args)

## ----
idents_before_filtering <- sort(unique(as.vector(as.character(Idents(seurat_data)))))    # A->Z sorted identities
if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)       # sets identities to new.ident
}
debug$print_info(seurat_data, args)
idents_after_filtering <- sort(unique(as.vector(as.character(Idents(seurat_data)))))     # A->Z sorted identities

## ----
print("Adjusting input parameters")
for (key in names(args)){
    if (key %in% c("mingenes", "maxgenes", "minumis", "minnovelty", "minfragments", "maxnuclsignal", "mintssenrich", "maxblacklist", "maxmt")){
        if (length(args[[key]]) == 1){
            print(paste("Extending filtering parameter", key, "to have a proper size"))
            args[[key]] <- rep(args[[key]][1], length(idents_after_filtering))           # we use number of identities after filtering as we may have potentially removed some of them
        } else {
            if (length(args[[key]]) != length(idents_before_filtering)){                 # we use the original number of identities to make sure
                print(                                                                   # that a user provided the the correct number of filtering
                    paste(                                                               # parameter from the very beginning (a.k.a the same as in --identity)
                        "The size of the filtering parameter", key, "is not",
                        "equal to the number of originally provided datasets.",
                        "Exiting"
                    )
                )
                quit(save="no", status=1, runLast=FALSE)
            }
            if (length(idents_after_filtering) != length(idents_before_filtering)){
                filtered_params <- c()
                for (i in 1:length(idents_before_filtering)){
                    if (idents_before_filtering[i] %in% idents_after_filtering){
                        filtered_params <- append(filtered_params, args[[key]][i])
                    } else {
                        print(
                            paste(
                                "Excluding value", args[[key]][i], "from the", key, "parameter as",
                                "identity", idents_before_filtering[i], "is not present anymore"
                            )
                        )
                    }
                }
                args[[key]] <- filtered_params
            }
        }
    }
}
print("Adjusted parameters")
print(args)

## ----
print("Adding RNA QC metrics for not filtered datasets")
seurat_data <- qc$add_rna_qc_metrics(seurat_data, args)

## ----
print("Adding ATAC QC metrics for not filtered datasets")
seurat_data <- qc$add_atac_qc_metrics(seurat_data, args)

## ----
print("Adding peak QC metrics for peaks called by Cell Ranger ARC")
seurat_data <- qc$add_peak_qc_metrics(seurat_data, blacklist_data, args)
debug$print_info(seurat_data, args)

## ----
args <- qc$update_qc_thresholds(
    seurat_data=seurat_data,
    args=args,
    qc_keys=c("minumis", "mingenes", "maxgenes", "maxmt"),
    qc_columns=c("nCount_RNA", "nFeature_RNA", "nFeature_RNA", "mito_percentage"),
    qc_coef=c(-2.5, -2.5, 5, 2),
    remove_rna_doublets=TRUE
)

## ----
args <- qc$update_qc_thresholds(
    seurat_data=seurat_data,
    args=args,
    qc_keys=c("minfragments"),
    qc_columns=c("nCount_ATAC"),
    qc_coef=c(-2.5),
    remove_atac_doublets=TRUE
)
print("Filtering thresholds adjusted based on the MAD values")
args$maxmt <- max(args$maxmt)                                      # need to convert it back to a single value (take max when it was autoestimated)
print(args)

## ----
export_all_qc_plots(
    seurat_data=seurat_data,
    suffix="raw",
    args=args,
    macs2_peaks=FALSE
)

## ----
print("Applying filters based on RNA QC metrics")
seurat_data <- filter$apply_rna_qc_filters(seurat_data, args)                            # cleans up all reductions
debug$print_info(seurat_data, args)

## ----
print("Applying filters based on ATAC QC metrics")
seurat_data <- filter$apply_atac_qc_filters(seurat_data, args)                           # cleans up all reductions
debug$print_info(seurat_data, args)

## ----
if (
    !is.null(args$removedoublets) &&
    ( args$removedoublets %in% c("union", "onlyrna", "onlyatac", "intersect") )
){
    print("Filtering by estimated doublets")
    seurat_data <- filter$remove_doublets(
        seurat_data,
        what_to_remove=args$removedoublets
    )
    debug$print_info(seurat_data, args)
}

## ----
if (!is.null(args$callby)){
    print("Forced to replace Cell Ranger ARC peaks with MACS2 peaks")
    seurat_data <- analyses$call_macs2_peaks(
        seurat_data=seurat_data,
        seqinfo_data=seqinfo_data,
        annotation_data=annotation_data,
        args=args
    )
    debug$print_info(seurat_data, args)
    print("Updating ATAC QC metrics after calling MACS2 peaks")
    seurat_data <- qc$add_atac_qc_metrics(seurat_data, args)                             # with the new peaks, we have different number of ATAC UMI counted per cell, so all the metrics should be updated
    print("Updating peak QC metrics after calling MACS2 peaks")
    seurat_data <- qc$add_peak_qc_metrics(seurat_data, blacklist_data, args)             # recalculate peak QC metrics
    debug$print_info(seurat_data, args)
    export_all_qc_plots(                                                                 # after RNA and ATAC filters have been applied
        seurat_data=seurat_data,
        suffix="mid_fltr",
        args=args,
        macs2_peaks=TRUE
    )
    print("Applying filters based on updated ATAC QC metrics after calling MACS2 peaks")
    seurat_data <- filter$apply_atac_qc_filters(seurat_data, args)     # cleans up all reductions
    debug$print_info(seurat_data, args)
}

## ----
print("Applying filters based on peaks QC metrics")
seurat_data <- filter$apply_peak_qc_filters(seurat_data, args)                           # cleans up all reductions

## ----
print("Updating ATAC QC metrics after all filtering thresholds applied")
seurat_data <- qc$add_atac_qc_metrics(seurat_data, args)                                 # recalculate ATAC QC metrics

## ----
print("Updating peak QC metrics after all filtering thresholds applied")
seurat_data <- qc$add_peak_qc_metrics(seurat_data, blacklist_data, args)                 # recalculate peak QC metrics
debug$print_info(seurat_data, args)

## ----
export_all_qc_plots(                                                                     # after all filters have been applied
    seurat_data=seurat_data,
    suffix="fltr",
    args=args,
    macs2_peaks=!is.null(args$callby)                                                    # can be both from 10x or MACS2
)

## ----
print("Adding genes vs RNA reads per cell as gene_rnaumi dimensionality reduction")
seurat_data@reductions[["gene_rnaumi"]] <- CreateDimReducObject(
    embeddings=as.matrix(
        seurat_data@meta.data[, c("nCount_RNA", "nFeature_RNA"), drop=FALSE] %>%
        dplyr::rename("GRU_1"="nCount_RNA", "GRU_2"="nFeature_RNA") %>%
        dplyr::mutate(GRU_1=log10(GRU_1), GRU_2=log10(GRU_2))
    ),
    key="GRU_",
    assay="RNA"
)

## ----
print("Adding RNA reads vs ATAC fragments in peaks per cell as rnaumi_atacfrgm dimensionality reduction")
seurat_data@reductions[["rnaumi_atacfrgm"]] <- CreateDimReducObject(
    embeddings=as.matrix(
        seurat_data@meta.data[, c("nCount_ATAC", "nCount_RNA"), drop=FALSE] %>%
        dplyr::rename("RAU_1"="nCount_ATAC", "RAU_2"="nCount_RNA") %>%
        dplyr::mutate(RAU_1=log10(RAU_1), RAU_2=log10(RAU_2))
    ),
    key="RAU_",
    assay="RNA"
)

## ----
print("Adding TSS enrichment score vs ATAC fragments in peaks per cell as tss_atacfrgm dimensionality reduction")
seurat_data@reductions[["tss_atacfrgm"]] <- CreateDimReducObject(
    embeddings=as.matrix(
        seurat_data@meta.data[, c("nCount_ATAC", "TSS.enrichment"), drop=FALSE] %>%
        dplyr::rename("TAU_1"="nCount_ATAC", "TAU_2"="TSS.enrichment") %>%
        dplyr::mutate(TAU_1=log10(TAU_1))
    ),
    key="TAU_",
    assay="ATAC"
)

## ----
if(args$cbbuild){
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        is_nested=TRUE,
        rootname=paste(args$output, "_cellbrowser/rna", sep=""),
    )
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="ATAC",
        slot="counts",
        short_label="ATAC",
        is_nested=TRUE,
        rootname=paste(args$output, "_cellbrowser/atac", sep=""),
    )
}

## ----
DefaultAssay(seurat_data) <- "RNA"                                       # better to stick to RNA assay by default https://www.biostars.org/p/395951/#395954 
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
        rootname=paste0(args$output, "_rna_counts")
    )
}