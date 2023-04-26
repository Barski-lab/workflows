#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(forcats))
suppressMessages(library(modules))
suppressMessages(library(argparse))
suppressMessages(library(SummarizedExperiment))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(filter <- modules::use(file.path(HERE, "modules/filter.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))


export_raw_plots <- function(seurat_data, args){
    DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure
    for (reduction in c("rnaumap", "atacumap", "wnnumap")) {
        if (!(reduction %in% names(seurat_data@reductions))) {next}                                  # skip missing reductions
        graphics$dim_plot(
            data=seurat_data,
            reduction=reduction,
            plot_title=paste0(
                "Cells UMAP split by ", args$splitby, " ",
                ifelse(
                    (!is.null(args$groupby) && !is.null(args$subset)),
                    paste("subsetted to", paste(args$subset, collapse=", "), "values from the", args$groupby, "column "),
                    ""
                ),
                "(", reduction, " dim. reduction)"
            ),
            legend_title=ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                args$groupby,
                args$splitby
            ),
            split_by=args$splitby,
            group_by=ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                args$groupby,
                args$splitby
            ),
            highlight_group=if(!is.null(args$groupby) && !is.null(args$subset)) args$subset else NULL,
            label=FALSE,
            label_color="black",
            palette_colors=if(!is.null(args$groupby) && !is.null(args$subset)) c("#4CB6BA", "#FB1C0D") else c("darkred", "darkorchid4"),
            theme=args$theme,
            rootname=paste(args$output, "umap_rd", reduction, sep="_"),
            pdf=args$pdf
        )
    }
    gc(verbose=FALSE)
}


export_processed_plots <- function(seurat_data, de_results, args){
    DefaultAssay(seurat_data) <- "RNA"                                  # safety measure
    Idents(seurat_data) <- "new.ident"                                  # safety measure
    
    if (!is.null(de_results$bulk)){
        graphics$mds_html_plot(
            norm_counts_data=de_results$bulk$counts_data,               # not filtered pseudobulk aggregated normalized counts data in a form of SummarizedExperiment
            rootname=paste(args$output, "mds_plot", sep="_")
        )

        pca_data <- qc$counts_pca(                                      # adds 'group' column to identify the datasets
            SummarizedExperiment::assay(de_results$bulk$counts_data)
        )

        graphics$pca_plot(
            pca_data=pca_data,
            pcs=c(1, 2),
            plot_title="Normalized read counts PCA (1, 2). All genes",
            legend_title="Dataset",
            color_by="group",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "pca_1_2", sep="_"),
            pdf=args$pdf
        )

        graphics$pca_plot(
            pca_data=pca_data,
            pcs=c(2, 3),
            plot_title="Normalized read counts PCA (2, 3). All genes",
            legend_title="Dataset",
            color_by="group",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "pca_2_3", sep="_"),
            pdf=args$pdf
        )
    }

    genes_to_highlight <- get_genes_to_highlight(de_results, args)  # may return empty vector
    graphics$volcano_plot(
        data=de_results$de_genes,                                   # this is not filtered differentially expressed features
        x_axis="log2FoldChange",
        y_axis="padj",
        x_cutoff=0,
        y_cutoff=args$padj,
        x_label="log2 FC",
        y_label="-log10 Padj",
        plot_title=paste(
            "Differentially expressed",
            ifelse(
                (!is.null(args$genes) && length(args$genes) > 0),
                "user provided",
                paste("top", length(genes_to_highlight))
            ),
            "genes"
        ),
        plot_subtitle=paste0(
            args$second, " vs ", args$first, " for cells split by ", args$splitby, ". ",
            "Use ", args$test, " test. ",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste("Subsetted to", paste(args$subset, collapse=", "), "values from", args$groupby, "column. "),
                ""
            ),
            ifelse(
                (!is.null(args$batchby)),
                paste0("Model batch effect by ", args$batchby, ". "),
                ""
            ),
            "Displayed Padj threshold equals to ", args$padj
        ),
        caption=paste(nrow(de_results$de_genes), "genes"),
        features=genes_to_highlight,
        theme=args$theme,
        rootname=paste(args$output, "dxpr_vlcn", sep="_"),
        pdf=args$pdf
    )

    if (!is.null(genes_to_highlight) && length(genes_to_highlight) > 0){
        graphics$vln_plot(
            data=seurat_data,
            features=genes_to_highlight,
            labels=genes_to_highlight,
            plot_title=paste0(
                "Gene expression density for ", args$second, " vs ", args$first, " cells",
                ifelse(
                    (!is.null(args$groupby) && !is.null(args$subset)),
                    paste(" subsetted to", paste(args$subset, collapse=", "), "values from", args$groupby, "column"),
                    ""
                )
            ),
            legend_title=args$splitby,
            log=TRUE,
            pt_size=0,
            group_by=args$splitby,
            combine_guides="collect",
            palette_colors=c("darkred", "darkorchid4"),
            ncol=ceiling(sqrt(length(genes_to_highlight))),
            theme=args$theme,
            rootname=paste(args$output, "xpr_dnst", sep="_"),
            pdf=args$pdf
        )
        for (current_gene in genes_to_highlight) {
            for (reduction in c("rnaumap", "atacumap", "wnnumap")) {
                if (!(reduction %in% names(seurat_data@reductions))) {next}                                  # skip missing reductions
                graphics$feature_plot(
                    data=seurat_data,
                    features=current_gene,
                    labels=current_gene,
                    reduction=reduction,
                    plot_title=paste0(
                        "Gene expression on cells UMAP split by ", args$splitby, " ",
                        ifelse(
                            (!is.null(args$groupby) && !is.null(args$subset)),
                            paste("subsetted to", paste(args$subset, collapse=", "), "values from the", args$groupby, "column "),
                            ""
                        ),
                        "(", reduction, " dim. reduction)"
                    ),
                    label=FALSE,
                    order=TRUE,
                    split_by=args$splitby,
                    combine_guides="collect",
                    theme=args$theme,
                    width=1200,
                    height=600,
                    rootname=paste(args$output, "xpr_per_cell_rd", reduction, current_gene, sep="_"),
                    pdf=args$pdf
                )
            }
        }
    }

    print("Reordering heatmaps columns to correspond to the column metadata exported to GCT file")
    for (i in 1:length(colnames(de_results$cell$column_metadata))) {
        current_column <- colnames(de_results$cell$column_metadata)[i]
        seurat_data@meta.data[[current_column]] <- factor(
            seurat_data@meta.data[[current_column]],
            levels=unique(                                               # set levels to the order of appearance in column_metadata
                as.character(
                    de_results$cell$column_metadata[[current_column]]
                )
            )
        )
    }

    if ("HCL" %in% colnames(de_results$cell$row_metadata)){       # safe to take from cells, as if present in bulk it will be still the same
        split_rows <- forcats::fct_inorder(as.character(de_results$cell$row_metadata[["HCL"]]))
    } else if ("HCL.1" %in% colnames(de_results$cell$row_metadata)){
        split_rows <- forcats::fct_inorder(as.character(de_results$cell$row_metadata[["HCL.1"]]))
    } else {
        split_rows <- NULL
    }

    graphics$feature_heatmap(                                      # we build it only for cells level
        data=seurat_data,
        assay="RNA",                                               # for now we will always use RNA even if SCT may be present
        slot="data",
        features=rownames(de_results$cell$row_metadata),           # to make sure we use the same number and order of genes as in GCT file
        show_rownames=nrow(de_results$cell$row_metadata) <= 60,    # hide gene names if there are too many of them
        scale_to_max=TRUE,
        plot_title=paste0(
            "Normalized gene expression for cells",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste(
                    " subsetted to", paste(args$subset, collapse=", "),
                    "values from", args$groupby, "column. "
                ),
                ". "
            ),
            "padj <= ", args$padj
        ),
        group_by=colnames(de_results$cell$column_metadata),
        split_rows=split_rows,
        palette_colors=graphics$D40_COLORS,
        heatmap_colors=c("black", "yellow"),
        rootname=paste(args$output, "xpr_htmp", sep="_"),
        pdf=args$pdf
    )

}


get_genes_to_highlight <- function(de_results, args){
    if (!is.null(args$genes) && length(args$genes) > 0){
        print("Using user provided genes to highlight regardless of the significance score")
        genes_to_highlight <- args$genes[args$genes %in% as.vector(as.character(de_results$de_genes[, "gene"]))]
    } else {
        print(
            paste(
                "Identifying top 10 the most DE genes with P adjusted <= ", args$padj
            )
        )
        top_de_genes <- de_results$de_genes %>%
                        filter(.$padj<=args$padj) %>%
                        arrange(desc(log2FoldChange)) %>%
                        filter(row_number() > max(row_number()) - 10 | row_number() <= 10)
        genes_to_highlight <- as.vector(as.character(top_de_genes[, "gene"]))
    }
    print(paste("Genes to highlight", paste(genes_to_highlight, collapse=", ")))
    return (genes_to_highlight)
}


get_args <- function(){
    parser <- ArgumentParser(
        description="Single-cell Differential Expression Analysis"
    )
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include genes",
            "expression information stored in the RNA assay. Additionally, 'rnaumap', and/or",
            "'atacumap', and/or 'wnnumap' dimensionality reductions should be present."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to optionally extend Seurat object metadata with",
            "categorical values using samples identities. First column - 'library_id'",
            "should correspond to all unique values from the 'new.ident' column of the",
            "loaded Seurat object. If any of the provided in this file columns are already",
            "present in the Seurat object metadata, they will be overwritten. When combined",
            "with --barcodes parameter, first the metadata will be extended, then barcode",
            "filtering will be applied. Default: no extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--barcodes",
        help=paste(
            "Path to the TSV/CSV file to optionally prefilter and extend Seurat object",
            "metadata by selected barcodes. First column should be named as 'barcode'.",
            "If file includes any other columns they will be added to the Seurat object",
            "metadata ovewriting the existing ones if those are present.",
            "Default: all cells used, no extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Column from the Seurat object metadata to group cells for optional subsetting",
            "when combined with --subset parameter. May be one of the extra metadata columns",
            "added with --metadata or --barcodes parameters. Ignored if --subset is not set.",
            "Default: do not subset, include all cells into analysis."
        ),
        type="character"
    )
    parser$add_argument(
        "--subset",
        help=paste(
            "Values from the column set with --groupby parameter to subset cells",
            "before running differential expression analysis. Ignored if --groupby",
            "is not provided. Default: do not subset cells, include all of them."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to split cells into two groups to",
            "run --second vs --first differential expression analysis. May be one of",
            "the extra metadata columns added with --metadata or --barcodes parameters."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--first",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby parameter",
            "to define the first group of cells for differential expression analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--second",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby parameter",
            "to define the second group of cells for differential expression analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--test",
        help=paste(
            "Test type to use in differential expression analysis. If set to deseq or",
            "deseq-lrt, gene expression will be aggregated to the pseudobulk level per",
            "dataset. For deseq, the pair-wise Wald test will be used. For deseq-lrt, the",
            "reduced formula will look like ~1 if --batchby parameter is omitted or will",
            "be set to ~batchby to exclude the criteria if interest (defined by --splitby).",
            "For all other values of the --test parameter the FindMarkers function will",
            "be used (genes will be prefiltered by minimum percentage >= 0.1 and by minimum",
            "log2FoldChange >= 0.25 before running differential expression analysis).",
            "Default: use FindMarkers with Wilcoxon Rank Sum test."
        ),
        type="character", default="wilcoxon",
        choices=c(
            "wilcoxon",                   # (wilcox) Wilcoxon Rank Sum test
            "likelihood-ratio",           # (bimod) Likelihood-ratio test
            "t-test",                     # (t) Student's t-test
            "negative-binomial",          # (negbinom) Negative Binomial Generalized Linear Model (supports --batchby)
            "poisson",                    # (poisson) Poisson Generalized Linear Model (supports --batchby)
            "logistic-regression",        # (LR) Logistic Regression (supports --batchby)
            "mast",                       # (MAST) MAST package (supports --batchby)
            "deseq",                      # DESeq2 Wald test on pseudobulk aggregated gene expression
            "deseq-lrt"                   # DESeq2 LRT test on pseudobulk aggregated gene expression
        )
    )
    parser$add_argument(
        "--batchby",
        help=paste(
            "Column from the Seurat object metadata to group cells into batches.",
            "If --test is set to deseq or deseq-lrt the --batchby parameter will",
            "be used in the design formula in the following way ~splitby+batchby.",
            "If --test is set to negative-binomial, poisson, logistic-regression,",
            "or mast it will be used as a latent variable in the FindMarkers function.",
            "Not supported for --test values equal to wilcoxon, likelihood-ratio, or",
            "t-test. May be one of the extra metadata columns added with --metadata",
            "or --barcodes parameters. Default: do not model batch effect."
        ),
        type="character"
    )
    parser$add_argument(
        "--padj",
        help=paste(
            "In the exploratory visualization part of the analysis output only",
            "differentially expressed genes with adjusted P-value not bigger",
            "than this value. Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to label on the generated plots. Default: top 10",
            "genes with the highest and the lowest log2FoldChange values."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--exclude",
        help=paste(
            "Regex pattern to identify and exclude specific genes from the",
            "differential expression analysis (not case-sensitive). If any",
            "of such genes are provided in the --genes parameter, they will",
            "be excluded from there as well. Default: use all genes"
        ),
        type="character"
    )
    parser$add_argument(
        "--cluster",
        help=paste(
            "Hopach clustering method to be run on the normalized read counts for",
            "the exploratory visualization part of the analysis. Clustering by",
            "column is supported only when --test is set to deseq or deseq-lrt.",
            "Default: do not run clustering"
        ),
        type="character",
        choices=c("row", "column", "both")
    )
    parser$add_argument(
        "--rowdist",
        help=paste(
            "Distance metric for HOPACH row clustering. Ignored if --cluster is set",
            "to column or not provided. Default: cosangle"
        ),
        type="character", default="cosangle",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--columndist",
        help=paste(
            "Distance metric for HOPACH column clustering. Ignored if --cluster is set",
            "to row or not provided. Default: euclid"
        ),
        type="character", default="euclid",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--center",
        help=paste(
            "Apply mean centering for gene expression prior to running",
            "clustering by row. Ignored if --cluster is set to column or",
            "not provided. Default: do not centered"
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

print("Adjusting --test parameter")
args$test <- switch(
    args$test,
    "wilcoxon"            = "wilcox",
    "likelihood-ratio"    = "bimod",
    "t-test"              = "t",
    "negative-binomial"   = "negbinom",
    "poisson"             = "poisson",
    "logistic-regression" = "LR",
    "mast"                = "MAST",
    "deseq"               = "deseq",
    "deseq-lrt"           = "lrt"
)
print(paste("--test was adjusted to", args$test))

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
if (!("RNA" %in% names(seurat_data@assays))){
    print(
        paste(
            "Loaded Seurat object doesn't include required RNA assay.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}
if ( !is.null(args$cluster) && (args$cluster %in% c("column", "both")) && !(args$test %in% c("deseq", "lrt")) ){
    print(
        paste(
            "Clustering by column is not supported if gene expression data",
            "are not aggregated to pseudobulk form. Check --test parameter.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}
if (!is.null(args$batchby) && (args$test %in% c("wilcox", "bimod", "t"))){
    print(
        paste(
            "Modeling batch effect is not supported if --test parameter",
            "was set to wilcoxon, likelihood-ratio, or t-test. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"

excluded_genes <- c()
if (!is.null(args$exclude)){
    excluded_genes <- grep(
        args$exclude,
        as.vector(as.character(rownames(seurat_data))),    # with RNA assay set as default the rownames should be genes
        value=TRUE,
        ignore.case=TRUE
    )
    print(
        paste(
            "Based on the pattern", args$exclude, "the following genes will",
            "be excluded from the differential expression analysis:",
            paste(excluded_genes, collapse=", ")
        )
    )
}

if (!is.null(args$genes) && length(args$genes) > 0){
    print(
        paste(
            "Adjusting genes of interest to include only those genes that are",
            "present in the loaded Seurat object and have not been excluded."
        )
    )
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]     # with RNA assay set as default the rownames should be genes
    args$genes <- args$genes[!(args$genes %in% excluded_genes)]                                  # excluded_genes can be empty list
    print(args$genes)
}

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

if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)
    debug$print_info(seurat_data, args)
}

print("Subsetting Seurat object to include only cells from the tested conditions")
seurat_data <- io$apply_metadata_filters(seurat_data, args$splitby, c(args$first, args$second))
debug$print_info(seurat_data, args)

export_raw_plots(seurat_data, args)                                                        # will include only cells from --first and --second

if(!is.null(args$groupby) && !is.null(args$subset)){
    print("Subsetting Seurat object to include only selected groups of cells")
    seurat_data <- io$apply_metadata_filters(seurat_data, args$groupby, args$subset)
    debug$print_info(seurat_data, args)
}

if (!all(table(seurat_data@meta.data[[args$splitby]]) > 0)){                               # check if we accidentally removed cells we want to compare
    print(
        paste(
            "Not enough cells for comparison. Check --groupby",
            "and --subset parameters. Exiting."
        )
    )
    print(table(seurat_data@meta.data[[args$splitby]]))
    quit(save="no", status=1, runLast=FALSE)
}

print("Normalizing RNA counts after all filters applied")
seurat_data <- NormalizeData(seurat_data, verbose=FALSE)                                   # just in case rerun normalize as we use data slot in FindMarkers function

de_results <- analyses$rna_de_analyze(seurat_data, args, excluded_genes)
print(head(de_results$de_genes))

export_processed_plots(seurat_data, de_results, args)

print("Exporting differentially expressed genes")
io$export_data(
    de_results$de_genes,                                                                    # not filtered na-removed differentially expressed genes
    paste(args$output, "_de_genes.tsv", sep="")
)

if (!is.null(de_results$bulk)){                                                             # looks like we run pseudobulk analysis
    print("Exporting not filtered normalized pseudobulk read counts in GCT format")
    io$export_gct(                                                                          # will be used by GSEA
        counts_mat=SummarizedExperiment::assay(de_results$bulk$counts_data),
        row_metadata=NULL,                                                                  # should be NULL, because de_results$row_metadata is filtered by padj
        col_metadata=de_results$bulk$column_metadata,                                       # includes samples as row names
        location=paste(args$output, "_bulk_counts.gct", sep="")
    )
    print("Exporting CLS phenotype data")
    io$export_cls(
        categories=de_results$bulk$column_metadata[, args$splitby],                         # to have the same samples order as in GCT file
        paste(args$output, "_bulk_phntps.cls", sep="")
    )
}

print("Exporting filtered normalized cell read counts to GCT format")
io$export_gct(                                                                              # will be used by Morpheus
    counts_mat=de_results$cell$counts_mat,
    row_metadata=de_results$cell$row_metadata,                                              # includes genes as row names
    col_metadata=de_results$cell$column_metadata,                                           # includes cells as row names
    location=paste(args$output, "_cell_counts.gct", sep="")
)
