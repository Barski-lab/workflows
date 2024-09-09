#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(knitr))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(modules))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(filter <- modules::use(file.path(HERE, "modules/filter.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))

## ----
export_all_qc_plots <- function(seurat_data, suffix, args){
    Idents(seurat_data) <- "new.ident"                                                                # safety measure
    selected_features <- c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi")
    selected_labels <- c("RNA reads", "Genes", "Mitochondrial %", "Novelty score")
    selected_scales <- c(TRUE, TRUE, FALSE, FALSE)
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
        highlight_rows=which(seurat_data@meta.data$rna_doublets == "doublet"),
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
        highlight_rows=which(seurat_data@meta.data$rna_doublets == "doublet"),
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

    if (nrow(seurat_data@meta.data[seurat_data@meta.data$rna_doublets == "doublet", ]) > 0){      # show plot only if we still have doublets
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
    }
}

## ----
get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell RNA-Seq Filtering Analysis")
    parser$add_argument(
        "--mex",
        help=paste(
            "Path to the folder with feature-barcode matrix from Cell Ranger Count/Aggregate",
            "experiment in MEX format. If multiple locations provided data is assumed to be not",
            "aggregated (outputs from the multiple Cell Ranger Count experiments) and will be",
            "merged before the analysis."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--identity",
        help=paste(
            "Path to the metadata TSV/CSV file to set the datasets identities, if '--mex' points to",
            "the Cell Ranger Aggregate outputs. The aggregation.csv file can be used. In case of",
            "using feature-barcode matrices the file with identities should include at least one",
            "column - 'library_id', and a row with aliases per each experiment from the '--mex'",
            "input. The order of rows should correspond to the order of feature-barcode matrices",
            "provided in the '--mex' parameter."
        ),
        type="character"
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
            "Include only genes detected in at least this many cells. Ignored when '--mex'",
            "points to the feature-barcode matrices from the multiple Cell Ranger Count",
            "experiments.",
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
        "--minnovelty",
        help=paste(
            "Include cells with the novelty score not lower than this value, calculated for",
            "as log10(genes)/log10(UMI). If multiple values provided, each of them will",
            "be applied to the correspondent dataset from the '--mex' input based on the",
            "'--identity' file.",
            "Default: 0.8 (applied to all datasets)"
        ),
        type="double", default=0.8, nargs="*"
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
        "--removedoublets",
        help=paste(
            "Remove cells that were identified as doublets. Cells with",            # from scDblFinder: it might be necessary to remove cells with
            "RNA UMI < 200 will not be evaluated. Default: do not remove",          # a very low coverage (e.g. <200 reads) to avoid errors
            "doublets"
        ),
        action="store_true"
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
cell_identity_data <- io$load_cell_identity_data(args$identity)                         # identities are always prepended with letters to keep the order

io$export_data(
    cell_identity_data[, "library_id", drop=FALSE],                                     # we need only the first column with the prepended letters
    paste(args$output, "_meta.tsv", sep="")
)

## ----
print(paste("Loading datasets grouping from", args$grouping))
grouping_data <- io$load_grouping_data(args$grouping, cell_identity_data)

## ----
print("Loading feature-barcode matrices from:")
for (location in args$mex){print(location)}

seurat_data <- io$load_10x_rna_data(                                                    # identities are set to the "new.ident" column
    args=args,
    cell_identity_data=cell_identity_data,
    grouping_data=grouping_data
)
debug$print_info(seurat_data, args)

## ----
seurat_data <- qc$estimate_doublets(                                                    # we want to search for doublets before we apply any filters
    seurat_data=seurat_data,
    assay="RNA",
    target_column="rna_doublets",
    dbl_rate=args$rnadbr,
    dbl_rate_sd=args$rnadbrsd
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
    if (key %in% c("mingenes", "maxgenes", "minumis", "minnovelty", "maxmt")){           # we convert maxmt to vector as well, because we want to use it in the update_qc_thresholds function
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
print("Adding QC metrics for not filtered datasets")
seurat_data <- qc$add_rna_qc_metrics(seurat_data, args)
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
print("Filtering thresholds adjusted based on the MAD values")
args$maxmt <- max(args$maxmt)                                      # need to convert it back to a single value (take max when it was autoestimated)
print(args)

## ----
export_all_qc_plots(
    seurat_data=seurat_data,
    suffix="raw",
    args=args
)

## ----
print("Applying filters based on QC metrics")
seurat_data <- filter$apply_rna_qc_filters(seurat_data, args)                              # cleans up all reductions
debug$print_info(seurat_data, args)

## ----
if (!is.null(args$removedoublets) && args$removedoublets){
    print("Filtering by estimated doublets")
    seurat_data <- filter$remove_doublets(
        seurat_data,
        what_to_remove="onlyrna"
    )
    debug$print_info(seurat_data, args)
}

## ----
export_all_qc_plots(                                                                       # after all filters have been applied
    seurat_data=seurat_data,
    suffix="fltr",
    args=args
)

## ----
print("Adding genes vs RNA reads per cell as gene_rnaumi dimensionality reduction")
seurat_data@reductions[["gene_rnaumi"]] <- CreateDimReducObject(
    embeddings=as.matrix(
        seurat_data@meta.data[, c("nCount_RNA", "nFeature_RNA"), drop=FALSE] %>%           # can't be included into the add_rna_qc_metrics function
        dplyr::rename("GRU_1"="nCount_RNA", "GRU_2"="nFeature_RNA") %>%                    # as apply_rna_qc_filters removes all reductions
        dplyr::mutate(GRU_1=log10(GRU_1), GRU_2=log10(GRU_2))
    ),
    key="GRU_",
    assay="RNA"
)

## ----
if(args$cbbuild){
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        rootname=paste(args$output, "_cellbrowser", sep=""),
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