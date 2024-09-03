#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(knitr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(stringr))
suppressMessages(library(modules))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


# https://bioconductor.org/packages/release/bioc/vignettes/TrajectoryUtils/inst/doc/overview.html
# https://github.com/dynverse/ti_slingshot/blob/master/package/R/ti_slingshot.R
# https://github.com/KirstLab/asc_seurat/blob/main/R/adapted_slingshot.R

## ----
export_all_plots <- function(seurat_data, args){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    lineages_count <- length(seurat_data@misc$trajectories[[args$reduction]]$slingshot@metadata$lineages)

    graphics$trajectory_plot(
        data=seurat_data,
        reduction=args$reduction,
        plot_title=paste(
            "Trajectory plot,",
            "colored by cluster,",
            "reduction", args$reduction
        ),
        legend_title="Cluster",
        color_cells="grouping",
        color_density="grouping",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "trjc_gr_clst", sep="_"),
        pdf=args$pdf
    )

    graphics$trajectory_plot(
        data=seurat_data,
        reduction=args$reduction,
        plot_title=paste(
            "Trajectory plot,",
            "colored by pseudotime,",
            "reduction", args$reduction
        ),
        color_cells="pseudotime",
        theme=args$theme,
        rootname=paste(args$output, "trjc_pstm", sep="_"),
        pdf=args$pdf
    )

    graphics$trajectory_graph(
        data=seurat_data,
        reduction=args$reduction,
        plot_title=paste(
            "Trajectory graph,",
            "colored by cluster"
        ),
        legend_title="Cluster",
        color_cells="grouping",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "grph_gr_clst", sep="_"),
        pdf=args$pdf
    )

    graphics$trajectory_graph(
        data=seurat_data,
        reduction=args$reduction,
        plot_title=paste(
            "Trajectory graph,",
            "colored by pseudotime"
        ),
        color_cells="pseudotime",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "grph_pstm", sep="_"),
        pdf=args$pdf
    )

    graphics$dendro_plot(
        data=seurat_data,
        reduction=args$reduction,
        plot_title=paste(
            "Dendrogram plot,",
            "colored by cluster"
        ),
        legend_title="Cluster",
        color_cells="grouping",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        height=ifelse(lineages_count > 1, 800, 400),
        rootname=paste(args$output, "dndr_gr_clst", sep="_"),
        pdf=args$pdf
    )

    graphics$dendro_plot(
        data=seurat_data,
        reduction=args$reduction,
        plot_title=paste(
            "Dendrogram plot,",
            "colored by pseudotime"
        ),
        color_cells="pseudotime",
        theme=args$theme,
        height=ifelse(lineages_count > 1, 800, 400),
        rootname=paste(args$output, "dndr_pstm", sep="_"),
        pdf=args$pdf
    )

    graphics$topology_plot(
        data=seurat_data,
        reduction=args$reduction,
        plot_title="Topology plot",
        theme=args$theme,
        rootname=paste(args$output, "tplg", sep="_"),
        pdf=args$pdf
    )

    graphics$trajectory_heatmap(
        data=seurat_data,
        reduction=args$reduction,
        assay="RNA",
        slot="data",
        plot_title="Gene expression heatmap",
        features=args$ngenes,
        height=16*args$ngenes,
        rootname=paste(args$output, "xpr_htmp", sep="_"),
        pdf=args$pdf
    )

    if (!is.null(args$genes) && length(args$genes) > 0){
        graphics$trajectory_expression(
            data=seurat_data,
            reduction=args$reduction,
            features=args$genes,
            plot_title="Gene expression along pseudotime",
            combine_guides="collect",
            theme=args$theme,
            rootname=paste(args$output, "xpr_pstm", sep="_"),
            pdf=args$pdf
        )
    }

    Idents(seurat_data) <- args$source
    ptime_column <- paste0("ptime_", args$reduction)
    ptime_ranges <- c(min(seurat_data@meta.data[[ptime_column]]), max(seurat_data@meta.data[[ptime_column]]))
    ptime_colors <- c("#140b34", "#38588c", "#25848e", "#2ab07f", "#83d44b", "#c0df25")
    datasets_count <- length(unique(as.vector(as.character(seurat_data@meta.data$new.ident))))
    conditions_count <- length(unique(as.vector(as.character(seurat_data@meta.data$condition))))

    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis=ptime_column,
        group_by="new.ident",
        split_by="new.ident",
        x_label="Pseudotime",
        y_label="Density",
        legend_title="Dataset",
        plot_title="Pseudotime density",
        plot_subtitle="Split by dataset",
        show_zoomed=TRUE,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "pstm_dnst_spl_idnt", sep="_"),
        pdf=args$pdf
    )

    graphics$trajectory_hist(
        data=seurat_data@meta.data,
        x_axis=ptime_column,
        group_by=args$source,
        split_by="new.ident",
        x_label="Pseudotime",
        y_label="Density",
        legend_title="Cluster",
            plot_title=paste(
                "Pseudotime histogram,",
                "colored by cluster,",
                "split by dataset"
            ),
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        rootname=paste(args$output, "pstm_hist_gr_clst_spl_idnt", sep="_"),
        pdf=args$pdf
    )

    if (
        all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
        conditions_count > 1
    ){
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis=ptime_column,
            group_by="condition",
            split_by="condition",
            x_label="Pseudotime",
            y_label="Density",
            legend_title="Condition",
            plot_title=paste(
                "Pseudotime density,",
                "split by grouping condition"
            ),
            show_zoomed=TRUE,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "pstm_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$trajectory_hist(
            data=seurat_data@meta.data,
            x_axis=ptime_column,
            group_by=args$source,
            split_by="condition",
            x_label="Pseudotime",
            y_label="Density",
            legend_title="Cluster",
                plot_title=paste(
                    "Pseudotime histogram,",
                    "colored by cluster,",
                    "split by grouping condition"
                ),
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "pstm_hist_gr_clst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
    }

    for (reduction in c("rnaumap", "atacumap", "wnnumap")){
        if (reduction == args$reduction || !(reduction %in% names(seurat_data@reductions))) {next}
        graphics$feature_plot(
            data=seurat_data,
            features=ptime_column,
            labels=NULL,
            from_meta=TRUE,
            reduction=reduction,
            plot_title=paste(
                "UMAP,",
                "colored by pseudotime,",
                "reduction", reduction
            ),
            label=TRUE,
            label_color="blue",
            label_size=6,
            alpha=0.5,
            pt_size=2,
            order=TRUE,
            gradient_colors=ptime_colors,
            color_scales=ptime_ranges,
            color_limits=ptime_ranges,
            theme=args$theme,
            rootname=paste(args$output, "umap_rd", reduction, sep="_"),
            pdf=args$pdf
        )
        if (datasets_count > 1){
            graphics$feature_plot(
                data=seurat_data,
                features=ptime_column,
                labels=NULL,
                from_meta=TRUE,
                reduction=reduction,
                split_by="new.ident",
                combine_guides="collect",
                plot_title=paste(
                    "UMAP,",
                    "colored by pseudotime,",
                    "split by dataset,",
                    "reduction", reduction
                ),
                label=TRUE,
                label_color="blue",
                label_size=6,
                alpha=0.5,
                pt_size=ifelse(datasets_count > 1, 1, 2),
                order=TRUE,
                gradient_colors=ptime_colors,
                color_scales=ptime_ranges,
                color_limits=ptime_ranges,
                height=ifelse(datasets_count == 2, 400, 800),
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_idnt_rd", reduction, sep="_"),
                pdf=args$pdf
            )
        }
        if (
            all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
            conditions_count > 1
        ){
            graphics$feature_plot(
                data=seurat_data,
                features=ptime_column,
                labels=NULL,
                from_meta=TRUE,
                reduction=reduction,
                split_by="condition",
                combine_guides="collect",
                plot_title=paste(
                    "UMAP,",
                    "colored by pseudotime,",
                    "split by grouping condition,",
                    "reduction", reduction
                ),
                label=TRUE,
                label_color="blue",
                label_size=6,
                alpha=0.5,
                pt_size=ifelse(conditions_count > 1, 1, 2),
                order=TRUE,
                gradient_colors=ptime_colors,
                color_scales=ptime_ranges,
                color_limits=ptime_ranges,
                height=ifelse(conditions_count == 2, 400, 800),
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_cnd_rd", reduction, sep="_"),
                pdf=args$pdf
            )
        }
    }
    SeuratObject::Idents(seurat_data) <- "new.ident"                                            # safety measure
}

## ----
get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell RNA-Seq Trajectory Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should",
            "include genes expression information stored in the RNA assay and",
            "dimensionality reduction specified in the --reduction parameter."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--reduction",
        help=paste(
            "Dimensionality reduction to be used in the trajectory analysis.",
            "Default: pca"
        ),
        type="character", default="pca"
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use (from 1 to 50).",
            "Default: use all available dimensions"
        ),
        type="integer"
    )
    parser$add_argument(
        "--source",
        help=paste(
            "Column from the metadata of the loaded",
            "Seurat object to select clusters from"
        ),
        type="character", required="True"
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
        "--start",
        help=paste(
            "Value from the metadata column defined with --source",
            "parameter to set the starting point for the trajectory.",
            "Default: defined automatically"
        ),
        type="character"
    )
    parser$add_argument(
        "--ngenes",
        help=paste(
            "Number of the most predictive genes to be shows",
            "on the gene expression heatmap. Default: 50"
        ),
        type="integer", default=50
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
if (!is.null(args$dimensions)){
    print("Adjusting --dimensions parameter")
    args$dimensions <- c(1:args$dimensions)
    print(paste("--dimensions was adjusted to", paste(args$dimensions, collapse=", ")))
}
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
if (!(args$reduction %in% names(seurat_data@reductions))){
    print(
        paste(
            "Loaded Seurat object doesn't include",
            args$reduction, "reduction. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)              # sets identities to new.ident
    debug$print_info(seurat_data, args)
}

## ----
print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"

print("Normalizing RNA counts")
seurat_data <- NormalizeData(seurat_data, verbose=FALSE)                                        # just in case rerun normalize, because we will show gene expression

## ----
if (!is.null(args$genes)){
    print(
        paste(
            "Adjusting genes of interest to include only those",
            "that are present in the loaded Seurat object"
        )
    )
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]    # with RNA assay set as default the rownames should be genes
    print(args$genes)
}

## ----
print("Running RNA trajectory analysis")
seurat_data <- analyses$add_trajectory(seurat_data, args)      # adds ptime_"reduction" column to metadata and misc object into trajectories list
debug$print_info(seurat_data, args)

## ----
export_all_plots(seurat_data, args)

## ----
if(args$cbbuild){
    print(paste("Moving", args$reduction, "reduction on top"))
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c(args$reduction, reduc_names[reduc_names!=args$reduction])
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]

    print("Attemting to identify label field")
    if (grepl("^custom_", args$source, ignore.case=TRUE)){
        label_field <- gsub("custom_", "Custom ", args$source)
    } else if (grepl("_res\\.", args$source, ignore.case=TRUE)) {
        split_line <- unlist(strsplit(args$source, split="_res\\."))
        label_field <- paste("Clustering (", split_line[1], " ", split_line[2], ")", sep="")
    } else {
        label_field <- NULL
    }
    print(paste("Label field", label_field))

    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            label_field=label_field,
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/rna", sep="")
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            # markers=all_atac_markers,                                        # can be NULL
            label_field=label_field,
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/atac", sep="")
        )
    } else {
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            label_field=label_field,
            rootname=paste(args$output, "_cellbrowser", sep="")
        )
    }
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