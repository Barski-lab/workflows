#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(filter <- modules::use(file.path(HERE, "modules/filter.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


export_raw_plots <- function(seurat_data, args){
    DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title=paste(
            "Cells UMAP split by", args$splitby, "column",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste("subsetted to", paste(args$subset, collapse=", "), "values from", args$groupby, "column"),
                ""
            )
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
        label=ifelse(
            (!is.null(args$groupby) && !is.null(args$subset)),
            TRUE,
            FALSE
        ),
        label_color="black",
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, "umap", sep="_"),
        pdf=args$pdf
    )
    gc(verbose=FALSE)
}


export_processed_plots <- function(seurat_data, de_results, args){
    DefaultAssay(seurat_data) <- "RNA"                                  # safety measure
    Idents(seurat_data) <- "new.ident"                                  # safety measure
    
    pca_data <- qc$counts_pca(de_results$vst_counts)                    # adds 'group' column to identify the datasets
    graphics$pca_plot(
        pca_data=pca_data,
        pcs=c(1, 2),
        plot_title=paste0(
            "VST normalized counts PCA subsetted to all DE genes regardless of Padj, PC1 and PC2",
            ifelse(
                (!is.null(args$batchby)),
                paste(", batch corrected by", args$batchby),
                ""
            )
        ),
        legend_title="Dataset",
        color_by="group",
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, "vst_pca_1_2", sep="_"),
        pdf=args$pdf
    )
    graphics$pca_plot(
        pca_data=pca_data,
        pcs=c(2, 3),
        plot_title=paste0(
            "VST normalized counts PCA subsetted to all DE genes regardless of Padj, PC2 and PC3",
            ifelse(
                (!is.null(args$batchby)),
                paste(", batch corrected by", args$batchby),
                ""
            )
        ),
        legend_title="Dataset",
        color_by="group",
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, "vst_pca_2_3", sep="_"),
        pdf=args$pdf
    )

    genes_to_highlight <- get_genes_to_highlight(de_results, args)
    graphics$volcano_plot(
        data=de_results$de_genes,                                   # this is not filtered differentially expressed features
        x_axis="log2FoldChange",
        y_axis="padj",
        x_cutoff=0,
        y_cutoff=args$alpha,
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
            args$second, " vs ", args$first, " pseudobulk DE analysis for samples split by ", args$splitby, ". ",
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
            "Displayed Padj threshold equals to ", args$alpha
        ),
        caption=paste(nrow(de_results$de_genes), "genes"),
        features=genes_to_highlight,
        rootname=paste(args$output, "dxpr_vlcn", sep="_"),
        pdf=args$pdf
    )

    if (!is.null(genes_to_highlight) && length(genes_to_highlight) > 0){
        seurat_data <- ScaleData(                                           # we need it only for gene expression heatmap
            seurat_data,
            features=genes_to_highlight,                                    # use those genes that we are interested in
            assay="RNA",
            verbose=FALSE
        )
        for (current_gene in genes_to_highlight) {
            graphics$vln_plot(
                data=seurat_data,
                features=current_gene,
                labels=current_gene,
                plot_title=paste(
                    "Log normalized gene expression density per dataset",
                    ifelse(
                        (!is.null(args$groupby) && !is.null(args$subset)),
                        paste("subsetted to", paste(args$subset, collapse=", "), "values from", args$groupby, "column"),
                        ""
                    )
                ),
                legend_title=args$splitby,
                log=TRUE,
                pt_size=0,
                split_by=args$splitby,
                combine_guides="collect",
                palette_colors=graphics$D40_COLORS,
                width=1000,
                height=400,
                rootname=paste(args$output, "xpr_dnst", current_gene, sep="_"),
                pdf=args$pdf
            )
            graphics$feature_plot(
                data=seurat_data,
                features=current_gene,
                labels=current_gene,
                reduction="rnaumap",
                plot_title=paste(
                    "Log normalized gene expression on cells UMAP per dataset",
                    ifelse(
                        (!is.null(args$groupby) && !is.null(args$subset)),
                        paste("subsetted to", paste(args$subset, collapse=", "), "values from", args$groupby, "column"),
                        ""
                    )
                ),
                label=FALSE,
                order=TRUE,
                split_by="new.ident",
                combine_guides="collect",
                rootname=paste(args$output, "xpr_per_cell", current_gene, sep="_"),
                pdf=args$pdf
            )
        }
        graphics$feature_heatmap(
            data=seurat_data,
            features=genes_to_highlight,
            plot_title=paste(
                "Log normalized gene expression heatmap per dataset",
                ifelse(
                    (!is.null(args$groupby) && !is.null(args$subset)),
                    paste("subsetted to", paste(args$subset, collapse=", "), "values from", args$groupby, "column"),
                    ""
                )
            ),
            assay="RNA",
            slot="scale.data",
            cells=NULL,
            group_by="new.ident",
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, "xpr_htmp", sep="_"),
            pdf=args$pdf
        )
    }
}


get_genes_to_highlight <- function(de_results, args){
    if (!is.null(args$genes) && length(args$genes) > 0){
        print("Using user provided genes to highlight regardless of the significance score")
        genes_to_highlight <- args$genes[args$genes %in% as.vector(as.character(de_results$de_genes[, "gene"]))]
    } else {
        print(
            paste(
                "Identifying top", args$topgenes, "the most DE genes with P adjusted <= ", args$alpha
            )
        )
        top_de_genes <- de_results$de_genes %>%
                        filter(.$padj<=args$alpha) %>%
                        arrange(desc(log2FoldChange)) %>%
                        filter(row_number() > max(row_number()) - all_of(args$topgenes) | row_number() <= all_of(args$topgenes))
        genes_to_highlight <- as.vector(as.character(top_de_genes[, "gene"]))
    }
    print(paste("Genes to highlight", paste(genes_to_highlight, collapse=", ")))
    return (genes_to_highlight)
}


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell Pseudobulk Differential Expression Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include genes",
            "expression information stored in the RNA assay. Certain plots require 'rnaumap'",
            "dimensionality reduction to be present."
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
            "present in the Seurat object metadata, they will be overwritten.",
            "Default: no extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to split cells into two groups",
            "to run --second vs --first pseudobulk DE analysis. May include columns",
            "from the extra metadata added with --metadata parameter."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--first",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby to define the",
            "first group of cells for pseudobulk DE analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--second",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby to define the",
            "second group of cells for pseudobulk DE analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--batchby",
        help=paste(
            "Column from the Seurat object metadata to define the variable that should",
            "be modeled as a batch effect when running pseudobulk DE analysis. May include",
            "columns from the extra metadata added with --metadata parameter.",
            "Default: do not model batch effect."
        ),
        type="character"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Column from the Seurat object metadata to group cells for optional subsetting",
            "when combined with --subset parameter. May include columns from the extra metadata",
            "added with --metadata parameter. Ignored if --subset is not provided.",
            "Default: do not subset, run pseudobulk DE analysis for all cells jointly"
        ),
        type="character"
    )
    parser$add_argument(
        "--subset",
        help=paste(
            "Value(s) from the column set with --groupby parameter to subset cells",
            "before running pseudobulk DE analysis. If multiple values are provided",
            "run analysis jointly for selected groups of cells. Ignored if --groupby",
            "is not set.",
            "Default: do not subset, run pseudobulk DE analysis for all cells jointly"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--lrt",
        help=paste(
            "Use LRT instead of the pair-wise Wald test. Determines if the increased",
            "likelihood of the data using the --splitby term is more than expected if",
            "that term is truly zero.",
            "Default: use Wald test"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--alpha",
        help=paste(
            "The maximum significance level alpha used by DESeq2 for independent filtering",
            "of results. The adjusted p values for the genes which do not pass this filter",
            "threshold will be set to NA and removed from the outputs. The same threshold",
            "will be used for identifying --topgenes the most DE genes.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to label on the generated plots.",
            "Default: --topgenes N genes with the highest and the lowest log2 fold change",
            "expression values."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--exclude",
        help=paste(
            "Regex pattern to identify and exclude non-coding RNA genes from the pseudobulk",
            "DE analysis (not case-sensitive). If any of such genes were provided in the --genes",
            "parameter, they will be excluded from there as well.",
            "Default: use all genes"
        ),
        type="character"
    )
    parser$add_argument(
        "--topgenes",
        help=paste(
            "Show N genes with the highest and N genes with the lowest log2 fold change",
            "expression values. Ignored if --genes are provided.",
            "Default: 10"
        ),
        type="integer", default=10
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
print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"
debug$print_info(seurat_data, args)

excluded_genes <- c()
if (!is.null(args$exclude)){
    excluded_genes <- grep(
        args$exclude,
        as.vector(as.character(rownames(seurat_data))),                                          # with RNA assay set as default the rownames should be genes
        value=TRUE,
        ignore.case=TRUE
    )
    print(
        paste(
            "Based on the pattern", args$exclude, "the following genes will",
            "be excluded from the pseudobulk DE analysis:",
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

print("Extending Seurat object with the extra metadata fields")
seurat_data <- io$extend_metadata(
    seurat_data=seurat_data,
    location=args$metadata,
    seurat_ref_column="new.ident",
    meta_ref_column="library_id"
)
debug$print_info(seurat_data, args)
export_raw_plots(seurat_data, args)

print("Filtering Seurat object to include only selected groups of cells")
seurat_data <- filter$apply_metadata_filters(seurat_data, args$splitby, c(args$first, args$second))
if(!is.null(args$groupby) && !is.null(args$subset)){
    seurat_data <- filter$apply_metadata_filters(seurat_data, args$groupby, args$subset)
}
debug$print_info(seurat_data, args)

print(
    paste(
        "Running", args$second, "vs", args$first, "differential expression",
        "analysis for pseudobulk RNA-Seq samples split by", args$splitby,
        ifelse(
            (!is.null(args$groupby) && !is.null(args$subset)),
            paste("subsetted to", paste(args$subset, collapse=", "), "values from", args$groupby, "column."),
            ""
        ),
        ifelse(
            (!is.null(args$batchby)),
            paste("Model batch effect by", args$batchby, "column."),
            ""
        )
    )
)
de_results <- analyses$rna_de_analyze(seurat_data, args, excluded_genes)
print(head(de_results$de_genes))

export_processed_plots(seurat_data, de_results, args)

print("Exporting differentially expressed genes")
io$export_data(
    de_results$de_genes,
    paste(args$output, "_de_genes.tsv", sep="")
)

print("Exporting VST normalized counts")
io$export_gct(
    de_results$vst_counts,
    paste(args$output, "_vst_counts.gct", sep="")
)

print("Exporting CLS phenotype data")
io$export_cls(
    de_results,
    paste(args$output, "_phenotypes.cls", sep=""),
    args
)
