#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(Seurat))
suppressMessages(library(modules))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))

set.seed(42)

# Book chapter with the related materials
# https://www.ncbi.nlm.nih.gov/books/NBK27130/#:~:text=The%20antigen%20receptors%20on%20B,cell%20receptor%E2%80%94that%20are%20associated
# https://www.sc-best-practices.org/air_repertoire/ir_profiling.html#raw-data


export_all_plots <- function(seurat_data, args){
    Idents(seurat_data) <- "new.ident"                                                               # safety measure
    datasets_count <- length(unique(as.vector(as.character(seurat_data@meta.data$new.ident))))
    conditions_count <- length(unique(as.vector(as.character(seurat_data@meta.data$condition))))

    graphics$clonotype_bar_plot(                                 # all NA are not taken into account
        data=seurat_data,
        clone_by=args$cloneby,
        chains=c(seurat_data@misc$vdj$chains, "both"),
        x_label="Dataset",
        y_label="Unique clonotypes %",
        legend_title="Dataset",
        plot_title="Percentage of unique clonotyps per dataset",
        split_by="new.ident",
        scale=TRUE,                                              # before calculating the percentage we exclude all the cells that don't have any clonotypes assigned
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        combine_guides="collect",
        width=1800,
        rootname=paste(args$output, "count_gr_idnt", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_bar_plot(                                 # all NA are not taken into account
        data=seurat_data,
        clone_by=args$cloneby,
        chains=c(seurat_data@misc$vdj$chains, "both"),
        x_label="Cluster",
        y_label="Unique clonotypes %",
        legend_title="Cluster",
        plot_title="Percentage of unique clonotyps per cluster",
        split_by=args$source,
        scale=TRUE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        combine_guides="collect",
        width=1800,
        rootname=paste(args$output, "count_gr_clst", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_homeostasis_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains="both",
        x_label="Dataset",
        y_label="Relative Abundance",
        legend_title="Clonotype group",
        plot_title=paste(
            "Clonal space homeostasis,",
            "split by dataset"
        ),
        split_by="new.ident",
        palette_colors=c("#2E86C1", "#5DADE2", "#A9DFBF", "#E74C3C", "#CB4335"),
        theme=args$theme,
        combine_guides="collect",
        rootname=paste(args$output, "hmst_spl_idnt", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_homeostasis_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains="both",
        x_label="Cluster",
        y_label="Relative Abundance",
        legend_title="Clonotype group",
        plot_title=paste(
            "Clonal space homeostasis,",
            "split by cluster"
        ),
        split_by=args$source,
        palette_colors=c("#2E86C1", "#5DADE2", "#A9DFBF", "#E74C3C", "#CB4335"),
        theme=args$theme,
        combine_guides="collect",
        rootname=paste(args$output, "hmst_spl_clst", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_overlap_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains="both",
        plot_title=paste(
            "Clonotypes similarity,",
            "split by cluster"
        ),
        split_by=args$source,
        x_label="Cluster",
        y_label="Cluster",
        method="morisita",
        theme=args$theme,
        combine_guides="collect",
        rootname=paste(args$output, "vrlp_spl_clst", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_overlap_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains="both",
        plot_title=paste(
            "Clonotypes similarity,",
            "split by dataset"
        ),
        split_by="new.ident",
        x_label="Dataset",
        y_label="Dataset",
        method="morisita",
        theme=args$theme,
        combine_guides="collect",
        rootname=paste(args$output, "vrlp_spl_idnt", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_network_plot(
        data=seurat_data,
        reduction="rnaumap",
        clone_by=args$cloneby,
        chains="both",
        plot_title=paste(
            "Clonotypes network,",
            "colored by cluster"
        ),
        legend_title="Cluster",
        group_by=args$source,
        theme=args$theme,
        pt_size=0.5,
        combine_guides="collect",
        rootname=paste(args$output, "ntwr_gr_clst", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_network_plot(
        data=seurat_data,
        reduction="rnaumap",
        clone_by=args$cloneby,
        chains="both",
        plot_title=paste(
            "Clonotypes network,",
            "colored by dataset"
        ),
        legend_title="Dataset",
        group_by="new.ident",
        theme=args$theme,
        pt_size=0.5,
        combine_guides="collect",
        rootname=paste(args$output, "ntwr_gr_idnt", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_diversity_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains="both",
        plot_title=paste(
            "Clonotypes diversity,",
            "colored by cluster,",
            "split by dataset"
        ),
        legend_title="Cluster",
        split_by="new.ident",
        group_by=args$source,
        x_label="Dataset",
        y_label="Index Score",
        theme=args$theme,
        width=1600,
        alpha=0.75,
        combine_guides="collect",
        rootname=paste(args$output, "dvrs_gr_clst_spl_idnt", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_diversity_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        chains="both",
        plot_title=paste(
            "Clonotypes diversity,",
            "colored by dataset,",
            "split by cluster"
        ),
        legend_title="Dataset",
        split_by=args$source,
        group_by="new.ident",
        x_label="Cluster",
        y_label="Index Score",
        theme=args$theme,
        width=1600,
        alpha=0.75,
        combine_guides="collect",
        rootname=paste(args$output, "dvrs_gr_idnt_spl_clst", sep="_"),
        pdf=args$pdf
    )

    features <- c("V", "J", "C")
    # if Idents are set to "new.ident" it will be used in
    # group_by(df[,ncol(df)], df[,y.axis], element.names)
    # because the last column in df will be identity


    graphics$clonotype_feature_plot(
        data=seurat_data,
        feature="D",
        chains="TRB",          # can't use "both" here
        plot_title=paste(
            "Relative usage of D",
            "genes, split by cluster"
        ),
        x_label="Gene",
        y_label="Mean gene usage",
        split_by=args$source,
        order_by="variance",
        scale=TRUE,
        theme=args$theme,
        combine_guides="collect",
        width=1600,
        height=600,
        rootname=paste(args$output, "gene_spl_clst", base::tolower("D"), sep="_"),
        pdf=args$pdf
    )

    for (i in 1:length(features)){
        current_feature <- features[i]
        graphics$clonotype_feature_plot(
            data=seurat_data,
            feature=current_feature,
            chains=seurat_data@misc$vdj$chains,          # can't use "both" here
            plot_title=paste(
                "Relative usage of", current_feature,
                "genes, split by cluster"
            ),
            x_label="Gene",
            y_label="Mean gene usage",
            split_by=args$source,
            order_by="variance",
            scale=TRUE,
            theme=args$theme,
            combine_guides="collect",
            width=1600,
            height=600,
            rootname=paste(args$output, "gene_spl_clst", base::tolower(current_feature), sep="_"),
            pdf=args$pdf
        )
        graphics$clonotype_feature_plot(
            data=seurat_data,
            feature=current_feature,
            chains=seurat_data@misc$vdj$chains,          # can't use "both" here
            plot_title=paste(
                "Relative usage of", current_feature,
                "genes, split by dataset"
            ),
            x_label="Gene",
            y_label="Mean gene usage",
            split_by="new.ident",
            order_by="variance",
            scale=TRUE,
            theme=args$theme,
            combine_guides="collect",
            width=1600,
            height=600,
            rootname=paste(args$output, "gene_spl_idnt", base::tolower(current_feature), sep="_"),
            pdf=args$pdf
        )
    }

    graphics$clonotype_chord_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        plot_title=paste(
            "Shared clonotype,",
            "colored by cluster"
        ),
        group_by=args$source,
        theme=args$theme,
        rootname=paste(args$output, "chrd_gr_clst", sep="_"),
        pdf=args$pdf
    )

    graphics$clonotype_chord_plot(
        data=seurat_data,
        clone_by=args$cloneby,
        plot_title=paste(
            "Shared clonotype,",
            "colored by dataset"
        ),
        group_by="new.ident",
        theme=args$theme,
        rootname=paste(args$output, "chrd_gr_idnt", sep="_"),
        pdf=args$pdf
    )

    if (
        all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
        conditions_count > 1
    ){
        graphics$clonotype_chord_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            plot_title=paste(
                "Shared clonotype,",
                "colored by grouping condition"
            ),
            group_by="condition",
            theme=args$theme,
            rootname=paste(args$output, "chrd_gr_cnd", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_bar_plot(                                 # all NA are not taken into account
            data=seurat_data,
            clone_by=args$cloneby,
            chains=c(seurat_data@misc$vdj$chains, "both"),
            x_label="Condition",
            y_label="Unique clonotypes %",
            legend_title="Condition",
            plot_title="Percentage of unique clonotyps per grouping condition",
            split_by="condition",
            scale=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            combine_guides="collect",
            width=1800,
            rootname=paste(args$output, "count_gr_cnd", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_homeostasis_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains="both",
            x_label="Condition",
            y_label="Relative Abundance",
            legend_title="Clonotype group",
            plot_title=paste(
                "Clonal space homeostasis,",
                "split by grouping condition"
            ),
            split_by="condition",
            palette_colors=c("#2E86C1", "#5DADE2", "#A9DFBF", "#E74C3C", "#CB4335"),
            theme=args$theme,
            combine_guides="collect",
            rootname=paste(args$output, "hmst_spl_cnd", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_overlap_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains="both",
            plot_title=paste(
                "Clonotypes similarity,",
                "split by grouping condition"
            ),
            split_by="condition",
            x_label="Condition",
            y_label="Condition",
            method="morisita",
            theme=args$theme,
            combine_guides="collect",
            rootname=paste(args$output, "vrlp_spl_cnd", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_network_plot(
            data=seurat_data,
            reduction="rnaumap",
            clone_by=args$cloneby,
            chains="both",
            plot_title=paste(
                "Clonotypes network,",
                "colored by grouping condition"
            ),
            legend_title="Condition",
            group_by="condition",
            theme=args$theme,
            pt_size=0.5,
            combine_guides="collect",
            rootname=paste(args$output, "ntwr_gr_cnd", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_diversity_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains="both",
            plot_title=paste(
                "Clonotypes diversity,",
                "colored by cluster,",
                "split by grouping condition"
            ),
            legend_title="Cluster",
            split_by="condition",
            group_by=args$source,
            x_label="Condition",
            y_label="Index Score",
            theme=args$theme,
            width=1600,
            alpha=0.75,
            combine_guides="collect",
            rootname=paste(args$output, "dvrs_gr_clst_spl_cnd", sep="_"),
            pdf=args$pdf
        )

        graphics$clonotype_diversity_plot(
            data=seurat_data,
            clone_by=args$cloneby,
            chains="both",
            plot_title=paste(
                "Clonotypes diversity,",
                "colored by grouping condition,",
                "split by cluster"
            ),
            legend_title="Dataset",
            split_by=args$source,
            group_by="condition",
            x_label="Condition",
            y_label="Index Score",
            theme=args$theme,
            width=1600,
            alpha=0.75,
            combine_guides="collect",
            rootname=paste(args$output, "dvrs_gr_cnd_spl_clst", sep="_"),
            pdf=args$pdf
        )
    }

    gc(verbose=FALSE)
}


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell Immune Profiling Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This",
            "file should include gene expression information stored",
            "in the RNA assay, as well as 'pca' and 'rnaumap'",
            "dimensionality reductions applied to that assay."
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
        "--source",
        help=paste(
            "Column from the metadata of the loaded Seurat",
            "object to select clusters from."
        ),
        type="character", required="True"
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
        "--groupby",
        help=paste(
            "Column from the metadata of the loaded Seurat object",
            "to group cells for clonotype frequency calculation.",
            "Default: group by dataset"
        ),
        type="character", default="new.ident"
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
        help="Save Seurat data to h5ad file. Default: false",
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

if (!("RNA" %in% names(seurat_data@assays))){
    print(
        paste(
            "Loaded Seurat object doesn't",
            "include RNA assay. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

if (!("rnaumap" %in% names(seurat_data@reductions))){
    print(
        paste(
            "Loaded Seurat object doesn't",
            "include rnaumap reduction. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
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
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)         # sets identities to new.ident
    debug$print_info(seurat_data, args)
}

print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"

print("Extending Seurat data with congtigs annotations")
seurat_data <- io$load_10x_vdj_data(seurat_data, args)
debug$print_info(seurat_data, args)

print("Quantifying clonotype counts")
seurat_data <- qc$quartile_qc_metrics(
    seurat_data=seurat_data,
    features=c("clntp_counts"),
    prefix="quartile"                                                                     # we will use this prefix in ucsc$export_cellbrowser function
)

export_all_plots(
    seurat_data=seurat_data,
    args=args
)

if(args$cbbuild){
    print("Exporting RNA assay to UCSC Cellbrowser")
    print("Reordering reductions to have rnaumap on the first place")                      # will be shown first in UCSC Cellbrowser
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c("rnaumap", reduc_names[reduc_names!="rnaumap"])               # we checked before that rnaumap is present
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]
    debug$print_info(seurat_data, args)
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        palette_colors=graphics$D40_COLORS,                                                # to have colors correspond to the plots
        rootname=paste(args$output, "_cellbrowser", sep="")
    )
}

DefaultAssay(seurat_data) <- "RNA"                                                         # better to stick to RNA assay by default https://www.biostars.org/p/395951/#395954 
print("Exporting results to RDS file")
# io$export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))
if(args$h5seurat){
    print("Exporting results to h5seurat file")
    io$export_h5seurat(seurat_data, paste(args$output, "_data.h5seurat", sep=""))
}

if(args$h5ad){
    print("Exporting results to h5ad file")
    io$export_h5ad(seurat_data, paste(args$output, "_data.h5ad", sep=""))
}

if(args$scope){
    print("Exporting results to SCope compatible loom file")
    io$export_scope_loom(                                                                  # we save only counts slot from the RNA assay 
        seurat_data,
        paste(args$output, "_data.loom", sep="")
    )
}
