#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(knitr))
suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(modules))
suppressMessages(library(argparse))
suppressMessages(library(tidyverse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


# Book chapter with the related materials
# https://www.ncbi.nlm.nih.gov/books/NBK27130/#:~:text=The%20antigen%20receptors%20on%20B,cell%20receptor%E2%80%94that%20are%20associated
# https://www.sc-best-practices.org/air_repertoire/ir_profiling.html#raw-data
# https://www.sciencedirect.com/science/article/pii/S0958166920301051

## ----
export_all_plots <- function(seurat_data, args){
    Idents(seurat_data) <- "new.ident"                                                               # safety measure
    datasets_count <- length(unique(as.vector(as.character(seurat_data@meta.data$new.ident))))
    conditions_count <- length(unique(as.vector(as.character(seurat_data@meta.data$condition))))
    donor_count <- length(unique(as.vector(as.character(seurat_data@meta.data$donor))))              # by this time it should be already present in the seurat object
    not_default_conditions <- all(
        as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))
    )
    not_default_donor <- all(
        as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$donor))
    )
    selected_features <- c("clonalFrequency_TRA", "clonalFrequency_TRB", "clonalFrequency_IGH", "clonalFrequency_IGL", "clonalFrequency_both")
    selected_labels <- c("TRA", "TRB", "IGH", "IGL", "Both")
    max_frequency <- max(
        seurat_data@meta.data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            max_frequency = max(
                dplyr::c_across(
                    tidyselect::any_of(selected_features)      # some columns might be missing
                )
            )
        ) %>%
        dplyr::pull(max_frequency)
    )

    graphics$clonotype_quant_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains=seurat_data@misc$vdj$chains,
        x_label="Dataset",
        y_label="Unique clonotypes %",
        legend_title="Dataset",
        plot_title="Percentage of unique clonotypes per dataset",
        plot_subtitle=paste(
            "Split by chain;",
            "filtered by clonotype frequency per donor >=",
            args$minfrequency
        ),
        group_by="new.ident",
        min_frequency=args$minfrequency,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        combine_guides="collect",
        width=ifelse(datasets_count > 1, 1200, 400),
        rootname=paste(args$output, "cl_qnt_gr_idnt_spl_ch", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_abundance_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains=seurat_data@misc$vdj$chains,
        x_label="Clonotype frequency",
        y_label="Distribution",
        legend_title="Dataset",
        plot_title="Distribution of clonotype frequencies per dataset",
        plot_subtitle=paste(
            "Split by chain;",
            "filtered by clonotype frequency per donor >=",
            args$minfrequency
        ),
        group_by="new.ident",
        min_frequency=args$minfrequency,
        scale=TRUE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        combine_guides="collect",
        height=400,
        rootname=paste(args$output, "cl_dnst_gr_idnt_spl_ch", sep="_"),
        pdf=args$pdf
    )

    top_clones_per_dataset <- floor(                                       # we don't want to exceed the number of available colors
        length(graphics$D40_COLORS)/
        length(seurat_data@misc$vdj$chains)/                               # each chain should have a separate set of colors
        datasets_count                                                     # for the worst scenario when neither of datasets have shared clonotypes
    )
    graphics$clonotype_alluvial_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains=seurat_data@misc$vdj$chains,
        x_label="Dataset",
        y_label="Proportion",
        legend_title="Clonotype",
        plot_title="Proportion of top shared clonotypes between datasets",
        plot_subtitle=paste0(
            "Split by chain; ",
            "filtered by clonotype frequency per donor >= ",
            args$minfrequency, "; ",
            "top ", top_clones_per_dataset, " clonotypes ",
            "selected from each dataset"
        ),
        group_by="new.ident",
        min_frequency=args$minfrequency,
        top_clones=top_clones_per_dataset,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        combine_guides="collect",
        ncol=3,
        legend_position="bottom",
        rootname=paste(args$output, "allu_gr_idnt_spl_ch", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_homeostasis_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains=seurat_data@misc$vdj$chains,
        x_label="Dataset",
        y_label="Proportion",
        legend_title="Clonotype size",
        plot_title="Proportion of clonotype frequencies per dataset",
        plot_subtitle=paste(
            "Split by chain;",
            "not filtered by clonotype frequency"
        ),
        group_by="new.ident",
        min_frequency=0,
        palette_colors=c("#080808", "#571357", "#BD3265", "#F28D30", "#FCFF9A"),
        theme=args$theme,
        combine_guides="collect",
        ncol=3,
        rootname=paste(args$output, "hmst_gr_idnt_spl_ch", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_overlap_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains=seurat_data@misc$vdj$chains,
        x_label="Dataset",
        y_label="Dataset",
        plot_title="Overlap of clonotypes between datasets",
        plot_subtitle=paste(
            "Split by chain;",
            "filtered by clonotype frequency per donor >=",
            args$minfrequency
        ),
        group_by="new.ident",
        min_frequency=args$minfrequency,
        methods=c("raw", "jaccard"),
        gradient_colors=c("lightgrey", "blue", "black", "orange"),
        theme=args$theme,
        combine_guides="collect",
        ncol=3,
        nrow=2,
        rootname=paste(args$output, "vrlp_gr_idnt_spl_ch", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_diversity_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains=seurat_data@misc$vdj$chains,
        x_label=if(donor_count > 1 && not_default_donor) "Donor" else NULL,
        y_label="Score",
        legend_title="Dataset",
        plot_title="Diversity of clonotypes per dataset",
        plot_subtitle=paste(
            "Split by chain;",
            "not filtered by clonotype frequency"
        ),
        group_by="new.ident",
        split_by=if(donor_count > 1 && not_default_donor) "donor" else NULL,    # optionally split by donor (score values remain the same)
        min_frequency=0,
        theme=args$theme,
        ncol=1,
        combine_guides="collect",
        rootname=paste(args$output, "dvrs_gr_idnt_spl_ch", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_feature_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains=seurat_data@misc$vdj$chains,
        x_label="Gene",
        y_label="Proportion",
        plot_title="Distribution of gene usage per dataset",
        plot_subtitle=paste(
            "Split by chain;",
            "filtered by clonotype frequency per donor >=",
            args$minfrequency
        ),
        group_by="new.ident",
        order_by="variance",
        min_frequency=args$minfrequency,
        scale=TRUE,
        theme=args$theme,
        ncol=ifelse("TRA" %in% seurat_data@misc$vdj$chains, 5, 3),
        width=ifelse("TRA" %in% seurat_data@misc$vdj$chains, 5*500, 3*500),
        height=datasets_count * 400,
        combine_guides="collect",
        rootname=paste(args$output, "gene_gr_idnt_spl_ch", sep="_"),
        pdf=args$pdf
    )

    graphics$feature_plot(
        data=seurat_data,
        features=selected_features,
        labels=selected_labels,
        from_meta=TRUE,
        reduction="rnaumap",
        legend_title="Frequency",
        plot_title="UMAP colored by clonotype frequency",
        plot_subtitle=paste(
            "Split by chain;",
            "filtered by clonotype frequency per donor >=",
            args$minfrequency
        ),
        label=FALSE,
        order=TRUE,
        pt_size=0.5,                                                  # need to have it fixed, otherwise dots are different on the splitted plot
        alpha=0.5,
        gradient_colors=c("lightgrey", "lightgrey", "darkred", "orange"),
        color_limits=c(0, max_frequency),
        color_scales=c(0, args$minfrequency-0.001*args$minfrequency, args$minfrequency, max_frequency),
        color_breaks=c(0, args$minfrequency, max_frequency),
        ncol=3,
        height=400,
        combine_guides="collect",
        theme=args$theme,
        rootname=paste(args$output, "umap_cl_freq_spl_ch", sep="_"),
        pdf=args$pdf
    )

    if (donor_count > 1 && not_default_donor){

        graphics$clonotype_quant_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Donor",
            y_label="Unique clonotypes %",
            legend_title="Donor",
            plot_title="Percentage of unique clonotypes per donor",
            plot_subtitle=paste(
                "Split by chain;",
                "filtered by clonotype frequency per donor >=",
                args$minfrequency
            ),
            group_by="donor",
            min_frequency=args$minfrequency,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            combine_guides="collect",
            width=ifelse(donor_count > 1, 1200, 400),
            rootname=paste(args$output, "cl_qnt_gr_dnr_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_abundance_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Clonotype frequency",
            y_label="Distribution",
            legend_title="Donor",
            plot_title="Distribution of clonotype frequencies per donor",
            plot_subtitle=paste(
                "Split by chain;",
                "filtered by clonotype frequency per donor >=",
                args$minfrequency
            ),
            group_by="donor",
            min_frequency=args$minfrequency,
            scale=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            combine_guides="collect",
            height=400,
            rootname=paste(args$output, "cl_dnst_gr_dnr_spl_ch", sep="_"),
            pdf=args$pdf
        )

        top_clones_per_donor <- floor(
            length(graphics$D40_COLORS)/
            length(seurat_data@misc$vdj$chains)/
            donor_count
        )
        graphics$clonotype_alluvial_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Donor",
            y_label="Proportion",
            legend_title="Clonotype",
            plot_title="Proportion of top shared clonotypes between donors",
            plot_subtitle=paste0(
                "Split by chain; ",
                "filtered by clonotype frequency per donor >= ",
                args$minfrequency, "; ",
                "top ", top_clones_per_donor, " clonotypes ",
                "selected from each donor"
            ),
            group_by="donor",
            min_frequency=args$minfrequency,
            top_clones=top_clones_per_donor,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            combine_guides="collect",
            ncol=3,
            legend_position="bottom",
            rootname=paste(args$output, "allu_gr_dnr_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_homeostasis_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Donor",
            y_label="Proportion",
            legend_title="Clonotype size",
            plot_title="Proportion of clonotype frequencies per donor",
            plot_subtitle=paste(
                "Split by chain;",
                "not filtered by clonotype frequency"
            ),
            group_by="donor",
            min_frequency=0,
            palette_colors=c("#080808", "#571357", "#BD3265", "#F28D30", "#FCFF9A"),
            theme=args$theme,
            combine_guides="collect",
            ncol=3,
            rootname=paste(args$output, "hmst_gr_dnr_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_overlap_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Donor",
            y_label="Donor",
            plot_title="Overlap of clonotypes between donors",
            plot_subtitle=paste(
                "Split by chain;",
                "filtered by clonotype frequency per donor >=",
                args$minfrequency
            ),
            group_by="donor",
            min_frequency=args$minfrequency,
            methods=c("raw", "jaccard"),
            gradient_colors=c("lightgrey", "blue", "black", "orange"),
            theme=args$theme,
            combine_guides="collect",
            ncol=3,
            nrow=2,
            rootname=paste(args$output, "vrlp_gr_dnr_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_diversity_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label=NULL,
            y_label="Score",
            legend_title="Donor",
            plot_title="Diversity of clonotypes per donor",
            plot_subtitle=paste(
                "Split by chain;",
                "not filtered by clonotype frequency"
            ),
            group_by="donor",
            split_by=NULL,
            min_frequency=0,
            theme=args$theme,
            ncol=1,
            combine_guides="collect",
            rootname=paste(args$output, "dvrs_gr_dnr_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_feature_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Gene",
            y_label="Proportion",
            plot_title="Distribution of gene usage per donor",
            plot_subtitle=paste(
                "Split by chain;",
                "filtered by clonotype frequency per donor >=",
                args$minfrequency
            ),
            group_by="donor",
            order_by="variance",
            min_frequency=args$minfrequency,
            scale=TRUE,
            theme=args$theme,
            ncol=ifelse("TRA" %in% seurat_data@misc$vdj$chains, 5, 3),
            width=ifelse("TRA" %in% seurat_data@misc$vdj$chains, 5*500, 3*500),
            height=donor_count * 400,
            combine_guides="collect",
            rootname=paste(args$output, "gene_gr_dnr_spl_ch", sep="_"),
            pdf=args$pdf
        )
    }

    if (conditions_count > 1 && not_default_conditions){

        graphics$clonotype_quant_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Condition",
            y_label="Unique clonotypes %",
            legend_title="Condition",
            plot_title="Percentage of unique clonotypes per grouping condition",
            plot_subtitle=paste(
                "Split by chain;",
                "filtered by clonotype frequency per donor >=",
                args$minfrequency
            ),
            group_by="condition",
            min_frequency=args$minfrequency,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            combine_guides="collect",
            width=ifelse(conditions_count > 1, 1200, 400),
            rootname=paste(args$output, "cl_qnt_gr_cnd_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_abundance_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Clonotype frequency",
            y_label="Distribution",
            legend_title="Condition",
            plot_title="Distribution of clonotype frequencies per grouping condition",
            plot_subtitle=paste(
                "Split by chain;",
                "filtered by clonotype frequency per donor >=",
                args$minfrequency
            ),
            group_by="condition",
            min_frequency=args$minfrequency,
            scale=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            combine_guides="collect",
            height=400,
            rootname=paste(args$output, "cl_dnst_gr_cnd_spl_ch", sep="_"),
            pdf=args$pdf
        )

        top_clones_per_condition <- floor(
            length(graphics$D40_COLORS)/
            length(seurat_data@misc$vdj$chains)/
            conditions_count
        )
        graphics$clonotype_alluvial_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Condition",
            y_label="Proportion",
            legend_title="Clonotype",
            plot_title="Proportion of top shared clonotypes between grouping conditions",
            plot_subtitle=paste0(
                "Split by chain; ",
                "filtered by clonotype frequency per donor >= ",
                args$minfrequency, "; ",
                "top ", top_clones_per_condition, " clonotypes ",
                "selected from each grouping condition"
            ),
            group_by="condition",
            min_frequency=args$minfrequency,
            top_clones=top_clones_per_condition,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            combine_guides="collect",
            ncol=3,
            legend_position="bottom",
            rootname=paste(args$output, "allu_gr_cnd_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_homeostasis_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Condition",
            y_label="Proportion",
            legend_title="Clonotype size",
            plot_title="Proportion of clonotype frequencies per grouping condition",
            plot_subtitle=paste(
                "Split by chain;",
                "not filtered by clonotype frequency"
            ),
            group_by="condition",
            min_frequency=0,
            palette_colors=c("#080808", "#571357", "#BD3265", "#F28D30", "#FCFF9A"),
            theme=args$theme,
            combine_guides="collect",
            ncol=3,
            rootname=paste(args$output, "hmst_gr_cnd_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_overlap_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Condition",
            y_label="Condition",
            plot_title="Overlap of clonotypes between grouping conditions",
            plot_subtitle=paste(
                "Split by chain;",
                "filtered by clonotype frequency per donor >=",
                args$minfrequency
            ),
            group_by="condition",
            min_frequency=args$minfrequency,
            methods=c("raw", "jaccard"),
            gradient_colors=c("lightgrey", "blue", "black", "orange"),
            theme=args$theme,
            combine_guides="collect",
            ncol=3,
            nrow=2,
            rootname=paste(args$output, "vrlp_gr_cnd_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_diversity_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label=NULL,
            y_label="Score",
            legend_title="Condition",
            plot_title="Diversity of clonotypes per grouping condition",
            plot_subtitle=paste(
                "Split by chain;",
                "not filtered by clonotype frequency"
            ),
            group_by="condition",
            split_by=NULL,
            min_frequency=0,
            theme=args$theme,
            ncol=1,
            combine_guides="collect",
            rootname=paste(args$output, "dvrs_gr_cnd_spl_ch", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_feature_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains=seurat_data@misc$vdj$chains,
            x_label="Gene",
            y_label="Proportion",
            plot_title="Distribution of gene usage per grouping condition",
            plot_subtitle=paste(
                "Split by chain;",
                "filtered by clonotype frequency per donor >=",
                args$minfrequency
            ),
            group_by="condition",
            order_by="variance",
            min_frequency=args$minfrequency,
            scale=TRUE,
            theme=args$theme,
            ncol=ifelse("TRA" %in% seurat_data@misc$vdj$chains, 5, 3),
            width=ifelse("TRA" %in% seurat_data@misc$vdj$chains, 5*500, 3*500),
            height=conditions_count * 400,
            combine_guides="collect",
            rootname=paste(args$output, "gene_gr_cnd_spl_ch", sep="_"),
            pdf=args$pdf
        )
    }

    gc(verbose=FALSE)
}

## ----
get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell Immune Profiling Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This",
            "file should include gene expression information stored",
            "in the RNA assay, as well as pca and rnaumap dimensionality",
            "reductions applied to that assay. If loaded Seurat object",
            "includes multiple datasets, it should have a donor column",       # donor is a required column when aggregating multiple Cell Ranger Multi samples
            "to define grouping for clonotype calling."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--contigs",
        help=paste(
            "Path to the file with high-level annotations of each",
            "high-confidence contig from cell-associated barcodes",
            "from the Cell Ranger Multi or Cell Ranger Aggregate",
            "experiments in TSV/CSV format."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to optionally extend Seurat",
            "object metadata with categorical values using samples",
            "identities. First column - 'library_id' should correspond",
            "to all unique values from the 'new.ident' column of the",
            "loaded Seurat object. If any of the provided in this file",
            "columns are already present in the Seurat object metadata,",
            "they will be overwritten. When combined with --barcodes",
            "parameter, first the metadata will be extended, then barcode",
            "filtering will be applied. Default: no extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--barcodes",
        help=paste(
            "Path to the TSV/CSV file to optionally prefilter and",
            "extend Seurat object metadata be selected barcodes.",
            "First column should be named as 'barcode'. If file",
            "includes any other columns they will be added to the",
            "Seurat object metadata ovewriting the existing ones",
            "if those are present. Default: all cells used, no",
            "extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--cloneby",
        help=paste(
            "Defines how to call the clonotype. gene: based on VDJC gene",
            "sequence. nt: based on the nucleotide sequence. aa: based on",
            "the amino acid sequence. strict: based on the combination of",
            "the nucleotide and gene sequences. Default: gene"
        ),
        type="character", default="gene",
        choices=c("gene", "nt", "aa", "strict")
    )
    parser$add_argument(
        "--minfrequency",
        help=paste(
            "Minimum frequency (number of cells) per clonotype to be reported.",
            "Default: 3"
        ),
        type="integer", default=3
    )
    parser$add_argument(
        "--filter",
        help=paste(
            "Stringency filters to be applied.",
            "cells: remove cells with more than 2 chains.",
            "chains: remove chains exceeding 2 (select the most expressed ones).",
            "Default: do not apply any filters."
        ),
        type="character",
        choices=c("cells", "chains")
    )
    parser$add_argument(
        "--removepartial",
        help=paste(
            "Remove cells with only one chain detected.",
            "Default: keep all cells if at least one chain detected"
        ),
        action="store_true"
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
debug$print_info(seurat_data, args)

## ----
if (!("RNA" %in% names(seurat_data@assays))){
    print(
        paste(
            "Loaded Seurat object doesn't",
            "include RNA assay. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (!("rnaumap" %in% names(seurat_data@reductions))){
    print(
        paste(
            "Loaded Seurat object doesn't",
            "include rnaumap reduction. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
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

## ----
if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)         # sets identities to new.ident
    debug$print_info(seurat_data, args)
}

## ----
print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"

print("Extending Seurat data with congtigs annotations")
seurat_data <- io$load_10x_vdj_data(seurat_data, args)
debug$print_info(seurat_data, args)

## ----
export_all_plots(
    seurat_data=seurat_data,
    args=args
)

## ----
io$export_clonotypes(
    data=seurat_data,
    location=paste0(args$output, "_clonotypes.tsv")
)

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
