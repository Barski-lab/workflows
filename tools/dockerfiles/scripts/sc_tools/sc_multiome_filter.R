#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
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


call_peaks <- function(seurat_data, args) {
    print("Calling MACS2 peaks grouping cells by identity")
    group_by <- "new.ident"
    if (args$callpeaks == "cluster"){
        print(
            paste(
                "Forced to group cell by RNA cluster identified from PCA",
                "reduction using", paste(args$dimensions, collapse=", "),
                "dimensions and", args$resolution, "resolution"
            )
        )
        print("Running RNA analysis before calling peaks by cluster")
        seurat_data <- analyses$rna_analyze(seurat_data, args)
        seurat_data <- filter$collapse_fragments_list(seurat_data)                 # collapse repetitive fragments as we could split the datasets before integration
        debug$print_info(seurat_data, args)
        group_by <- paste0("peak_calling_res.", args$resolution)
        seurat_data <- FindNeighbors(
            seurat_data,
            reduction="pca",                                                # this is the same reduction we got after running RunPCA on our RNA data, it's bound to the specific assay
            dims=args$dimensions,
            graph.name=c("peak_calling_nn", "peak_calling"),
            verbose=FALSE
        )
        seurat_data <- FindClusters(
            seurat_data,
            graph.name="peak_calling",
            resolution=args$resolution,
            verbose=FALSE
        )
    }
    Idents(seurat_data) <- "new.ident"                                      # safety measure as FindClusters may change identity column
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "ATAC"
    macs2_peaks <- CallPeaks(
        seurat_data,
        group.by=group_by,
        verbose=args$verbose
    )
    macs2_counts <- FeatureMatrix(
        fragments=Fragments(seurat_data),
        sep=c("-", "-"),
        features=macs2_peaks,
        cells=colnames(seurat_data),
        verbose=FALSE
    )
    rm(macs2_peaks)                                                         # remove unused data
    genome_annotation <- rtracklayer::import(args$annotations, format="GFF")
    if( !("gene_biotype" %in% colnames(GenomicRanges::mcols(genome_annotation))) ){
        print("Loaded genome annotation doesn't have 'gene_biotype' column. Adding NA")
        genome_annotation$gene_biotype <- NA
    }
    atac_assay <- CreateChromatinAssay(
        counts=macs2_counts,
        sep=c("-", "-"),
        fragments=Fragments(seurat_data),
        min.cells=0,                                                        # setting something other than 0 will update nCount_ATAC, which bring some discrepancy to the QC plots
        min.features=-1,                                                    # as they check ncount.cell > min.features and by default it's 0, we will remove cells without peaks and won't be able to add new assay to our seurat_data
        annotation=genome_annotation
    )
    rm(macs2_counts, genome_annotation)                                            # remove unused data
    seurat_data[["ATAC"]] <- atac_assay
    rm(atac_assay)                                                          # remove unused data
    DefaultAssay(seurat_data) <- backup_assay
    Idents(seurat_data) <- "new.ident"                                      # safety measure
    gc()
    return (seurat_data)
}

export_all_dimensionality_plots <- function(seurat_data, suffix, args) {
    Idents(seurat_data) <- "new.ident"                                      # safety measure
    graphics$elbow_plot(
        data=seurat_data,
        ndims=50,
        plot_title=paste("Elbow plot (from cells PCA) for RNA assay (", suffix, ")", sep=""),
        rootname=paste(args$output, suffix, "elbow", sep="_"),
        pdf=args$pdf
    )
    graphics$corr_plot(
        data=seurat_data,
        reduction="pca",
        highlight_dims=args$dimensions,
        qc_columns=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi"),
        qc_labels=c("UMI(RNA)", "Genes", "Mitochondrial %", "Novelty score"),
        plot_title=paste(
            "Correlation plots between QC metrics and cells PCA components for RNA assay ",
            "(", suffix, ")", sep=""
        ),
        combine_guides="collect",
        rootname=paste(args$output, suffix, "qc_dim_corr", sep="_"),
        pdf=args$pdf
    )
}

export_all_clustering_plots <- function(seurat_data, suffix, args, cluster_prefix, macs2_peaks=FALSE){
    Idents(seurat_data) <- "new.ident"                                    # safety measure
    peak_type <- ifelse(macs2_peaks, "- MACS2", "- 10x")
    selected_features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "nFeature_ATAC", "frip", "blacklist_fraction")
    selected_labels=c(
        "UMI(RNA)", "Genes", "Mitochondrial %", "Novelty score",
        paste(
            c("UMI(ATAC)", "TSS enrichment score", "Nucl. signal", "Peaks", "FRiP", "Bl. regions"),
            peak_type
        )
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title=paste("Clustered cells UMAP for RNA assay (", suffix, "). Resolution ", args$resolution, sep=""),
        legend_title="Cluster",
        group_by=paste(paste(cluster_prefix, "res", sep="_"), args$resolution, sep="."),
        label=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "umap_res", args$resolution, sep="_"),
        pdf=args$pdf
    )
    Idents(seurat_data) <- paste(paste(cluster_prefix, "res", sep="_"), args$resolution, sep=".")
    graphics$feature_plot(
        data=seurat_data,
        features=selected_features,
        labels=selected_labels,
        from_meta=TRUE,
        reduction="rnaumap",
        plot_title=paste("QC metrics on cells UMAP for RNA assay (", suffix, "). Resolution ", args$resolution, sep=""),
        label=TRUE,
        alpha=0.4,
        max_cutoff="q99",                    # to prevent outlier cells to distort coloring
        combine_guides="keep",
        rootname=paste(args$output, suffix, "umap_qc_mtrcs_res", args$resolution, sep="_"),
        pdf=args$pdf
    )
    Idents(seurat_data) <- "new.ident"
}

export_all_qc_plots <- function(seurat_data, suffix, args, macs2_peaks=FALSE){
    Idents(seurat_data) <- "new.ident"                                                                # safety measure
    peak_type <- ifelse(macs2_peaks, "- MACS2", "- 10x")
    selected_features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "nFeature_ATAC", "frip", "blacklist_fraction")
    selected_labels=c(
        "UMI(RNA)", "Genes", "Mitochondrial %", "Novelty score",
        paste(c("UMI(ATAC)", "TSS enrichment score", "Nucl. signal", "Peaks", "FRiP", "Bl. regions"), peak_type)
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
        plot_title=paste(
            paste(
                paste0("PC", c(1, 2)),
                collapse=" and "
            ),
            " from the QC metrics PCA (", suffix, ")", sep=""
        ),
        legend_title="QC metrics",
        color_by="labels",
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, paste(c(1, 2) ,collapse="_"), "qc_mtrcs_pca", sep="_"),
        pdf=args$pdf
    )
    graphics$pca_plot(
        pca_data=qc_metrics_pca,
        pcs=c(2, 3),
        plot_title=paste(
            paste(
                paste0("PC", c(2, 3)),
                collapse=" and "
            ),
            " from the QC metrics PCA (", suffix, ")", sep=""
        ),
        legend_title="QC metrics",
        color_by="labels",
        palette_colors=graphics$D40_COLORS,
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
        plot_title=paste("Number of cells per dataset (", suffix, ")", sep=""),
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "cells_count", sep="_"),
        pdf=args$pdf
    )

    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nCount_RNA",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$rnaminumi,
        x_label="UMI per cell",
        y_label="Density",
        legend_title="Dataset",
        plot_title=paste("UMI per cell density for RNA assay (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "rna_umi_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nFeature_RNA",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$mingenes,
        x_right_intercept=args$maxgenes,
        x_label="Genes per cell",
        y_label="Density",
        legend_title="Dataset",
        plot_title=paste("Genes per cell density (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "gene_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_point_plot(
        data=seurat_data@meta.data,
        x_axis="nCount_RNA",
        y_axis="nFeature_RNA",
        facet_by="new.ident",
        x_left_intercept=args$rnaminumi,
        y_low_intercept=args$mingenes,
        y_high_intercept=args$maxgenes,
        color_by="mito_percentage",
        gradient_colors=c("lightslateblue", "red", "green"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        x_label="UMI per cell",
        y_label="Genes per cell",
        legend_title="Mitochondrial %",
        plot_title=paste("Genes vs UMI per cell correlation for RNA assay (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "gene_umi_corr", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="mito_percentage",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$maxmt,
        x_label="Percentage of transcripts mapped to mitochondrial genes per cell",
        y_label="Density",
        legend_title="Dataset",
        plot_title=paste("Percentage of transcripts mapped to mitochondrial genes per cell density (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "mito_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="log10_gene_per_log10_umi",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$minnovelty,
        x_label="log10 Genes / log10 UMI per cell",
        y_label="Density",
        legend_title="Dataset",
        plot_title=paste("Novelty score per cell density for RNA assay (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "nvlt_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nCount_ATAC",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$atacminumi,
        x_label="UMI per cell",
        y_label="Density",
        legend_title="Dataset",
        plot_title=paste("UMI per cell density for ATAC assay (", suffix, ") ", peak_type, sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "atac_umi_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nFeature_ATAC",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=0,
        x_label="Peaks per cell",
        y_label="Density",
        legend_title="Dataset",
        plot_title=paste("Peaks per cell density (", suffix, ") ", peak_type, sep=""),
        scale_x_log10=FALSE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "peak_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="blacklist_fraction",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$maxblacklist,
        x_label="Fraction of ATAC fragments within genomic blacklist regions per cell",
        y_label="Density",
        legend_title="Dataset",
        plot_title=paste("Fraction of ATAC fragments within genomic blacklist regions per cell density (", suffix, ") ", peak_type, sep=""),
        scale_x_log10=FALSE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "blck_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_point_plot(
        data=seurat_data@meta.data,
        facet_by="new.ident",
        x_axis="nCount_ATAC",
        x_label="UMI(ATAC) per cell",
        y_axis="nCount_RNA",
        y_label="UMI(RNA) per cell",
        x_left_intercept=args$atacminumi,
        y_low_intercept=args$rnaminumi,
        alpha_intercept=1,
        color_by="mito_percentage",
        gradient_colors=c("lightslateblue", "red", "green"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        legend_title="Mitochondrial %",
        plot_title=paste("UMI per cell correlation for RNA vs ATAC assays (", suffix, ") ", peak_type, sep=""),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "rna_atac_umi_corr", sep="_"),
        pdf=args$pdf
    )
    graphics$vln_plot(
        data=seurat_data,
        features=selected_features,
        labels=selected_labels,
        from_meta=TRUE,
        plot_title=paste("QC metrics per cell density (", suffix, ")", sep=""),
        legend_title="Dataset",
        hide_x_text=TRUE,
        pt_size=0,
        combine_guides="collect",
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "qc_mtrcs_dnst", sep="_"),
        pdf=args$pdf
    )

    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "ATAC"
    graphics$tss_plot(
        data=seurat_data,
        split_by="new.ident",
        group_by_value=args$mintssenrich,
        combine_guides="collect",
        plot_title=paste("TSS enrichment score (", suffix, ") ", peak_type, sep=""),
        rootname=paste(args$output, suffix, "tss_nrch", sep="_"),
        pdf=args$pdf
    )
    graphics$fragments_hist(
        data=seurat_data,
        split_by="new.ident",
        group_by_value=args$maxnuclsignal,
        combine_guides="collect",
        plot_title=paste("Fragments length histogram (", suffix, ") ", peak_type, sep=""),
        rootname=paste(args$output, suffix, "frgm_hist", sep="_"),
        pdf=args$pdf
    )
    DefaultAssay(seurat_data) <- backup_assay

    if (seurat_data@meta.data$new.ident != seurat_data@meta.data$condition){
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_RNA",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$rnaminumi,
            x_label="UMI per cell",
            y_label="Density",
            legend_title="Dataset",
            plot_title=paste("Split by grouping condition UMI per cell density for RNA assay (", suffix, ")", sep=""),
            scale_x_log10=TRUE,
            zoom_on_intercept=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "rna_umi_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nFeature_RNA",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$mingenes,
            x_right_intercept=args$maxgenes,
            x_label="Genes per cell",
            y_label="Density",
            legend_title="Dataset",
            plot_title=paste("Split by grouping condition genes per cell density (", suffix, ")", sep=""),
            scale_x_log10=TRUE,
            zoom_on_intercept=TRUE,
            show_ranked=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "gene_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="mito_percentage",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$maxmt,
            x_label="Percentage of transcripts mapped to mitochondrial genes per cell",
            y_label="Density",
            legend_title="Dataset",
            plot_title=paste("Split by grouping condition the percentage of transcripts mapped to mitochondrial genes per cell density (", suffix, ")", sep=""),
            zoom_on_intercept=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "mito_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="log10_gene_per_log10_umi",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$minnovelty,
            x_label="log10 Genes / log10 UMI per cell",
            y_label="Density",
            legend_title="Dataset",
            plot_title=paste("Split by grouping condition the novelty score per cell density for RNA assay (", suffix, ")", sep=""),
            zoom_on_intercept=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "nvlt_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_ATAC",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$atacminumi,
            x_label="UMI per cell",
            y_label="Density",
            legend_title="Dataset",
            plot_title=paste("Split by grouping condition UMI per cell density for ATAC assay (", suffix, ") ", peak_type, sep=""),
            scale_x_log10=TRUE,
            zoom_on_intercept=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "atac_umi_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nFeature_ATAC",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=0,
            x_label="Peaks per cell",
            y_label="Density",
            legend_title="Dataset",
            plot_title=paste("Split by grouping condition peaks per cell density (", suffix, ") ", peak_type, sep=""),
            scale_x_log10=FALSE,
            zoom_on_intercept=TRUE,
            show_ranked=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "peak_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="blacklist_fraction",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$maxblacklist,
            x_label="Fraction of ATAC fragments within genomic blacklist regions per cell",
            y_label="Density",
            legend_title="Dataset",
            plot_title=paste("Split by grouping condition the fraction of ATAC fragments within genomic blacklist regions per cell density (", suffix, ") ", peak_type, sep=""),
            scale_x_log10=FALSE,
            zoom_on_intercept=TRUE,
            show_ranked=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "blck_dnst_spl_cnd", sep="_"),
            pdf=args$pdf
        )
    }
}

get_args <- function(){
    parser <- ArgumentParser(description="Single-cell Multiome ATAC and RNA-Seq Filtering Analysis")
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
            "Path to the metadata TSV/CSV file to set the datasets identities. If '--mex' points to",
            "the Cell Ranger ARC Aggregate outputs, the aggr.csv file can be used. If Cell Ranger",
            "ARC Count outputs have been used in the '--mex' input, the file should include at least",
            "one column - 'library_id' and one row with the alias for Cell Ranger ARC Count experiment."
        ),
        type="character", required="True"
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
            "Path to the headerless TSV/CSV file with the list of barcodes to select",
            "cells of interest (one barcode per line). Prefilters input feature-barcode",
            "matrix to include only selected cells.",
            "Default: use all cells."
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
            "'--mex' input based on the '--identity' file.",
            "Default: 250 (applied to all datasets)"
        ),
        type="integer", default=250, nargs="*"
    )
    parser$add_argument(
        "--maxgenes",
        help=paste(
            "Include cells with the number of genes not bigger than this value. If multiple",
            "values provided, each of them will be applied to the correspondent dataset from",
            "the '--mex' input based on the '--identity' file.",
            "Default: 5000 (applied to all datasets)"
        ),
        type="integer", default=5000, nargs="*"
    )
    parser$add_argument(
        "--rnaminumi",
        help=paste(
            "Include cells where at least this many UMI (RNA transcripts) are detected.",
            "If multiple values provided, each of them will be applied to the correspondent",
            "dataset from the '--mex' input based on the '--identity' file.",
            "Default: 500 (applied to all datasets)"
        ),
        type="integer", default=500, nargs="*"
    )
    parser$add_argument(
        "--mitopattern",
        help=paste(
            "Regex pattern to identify mitochondrial genes.",
            "Default: '^Mt-'"
        ),
        type="character", default="^Mt-"
    )
    parser$add_argument(
        "--maxmt",
        help=paste(
            "Include cells with the percentage of transcripts mapped to mitochondrial",
            "genes not bigger than this value.",
            "Default: 5 (applied to all datasets)"
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
    parser$add_argument(
        "--atacminumi",
        help=paste(
            "Include cells where at least this many UMI (ATAC transcripts) are detected.",
            "If multiple values provided, each of them will be applied to the correspondent",
            "dataset from the '--mex' input based on the '--identity' file.",
            "Default: 1000 (applied to all datasets)"
        ),
        type="integer", default=1000, nargs="*"
    )
    parser$add_argument(
        "--maxnuclsignal",
        help=paste(
            "Include cells with the nucleosome signal not bigger than this value.",
            "Nucleosome signal quantifies the approximate ratio of mononucleosomal",
            "to nucleosome-free fragments. If multiple values provided, each of",
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
            "Score is calculated based on the ratio of fragments centered at the TSS",
            "to fragments in TSS-flanking regions. If multiple values provided, each",
            "of them will be applied to the correspondent dataset from the '--mex' input",
            "based on the '--identity' file.",
            "Default: 2 (applied to all datasets)"
        ),
        type="double", default=2, nargs="*"
    )
    parser$add_argument(
        "--minfrip",
        help=paste(
            "Include cells with the FRiP not lower than this value. If multiple values",
            "provided, each of them will be applied to the correspondent dataset from the",
            "'--mex' input based on the '--identity' file. FRiP is calculated for fragments.",
            "Default: 0.15 (applied to all datasets)"
        ),
        type="double", default=0.15, nargs="*"
    )
    parser$add_argument(
        "--maxblacklist",
        help=paste(
            "Include cells with the fraction of fragments in genomic blacklist regions",
            "not bigger than this value. If multiple values provided, each of them",
            "will be applied to the correspondent dataset from the '--mex' input based",
            "on the '--identity' file.",
            "Default: 0.05 (applied to all datasets)"
        ),
        type="double", default=0.05, nargs="*"
    )
    parser$add_argument(
        "--callpeaks",
        help=paste(
            "Call peaks with MACS2 instead of those that are provided by Cell Ranger ARC Count.",
            "Peaks are called per identity or per RNA cluster after applying all RNA related",
            "thresholds, maximum nucleosome signal, and minimum TSS enrichment scores filters.",
            "If set to 'cluster' RNA clusters are identified based on the parameters set with",
            "'--resolution', '--dimensions', '--highvargenes', '--norm', and '--ntgr'.",
            "Default: do not call peaks"
        ),
        type="character",
        choices=c("identity", "cluster")
    )
    parser$add_argument(
        "--norm",
        help=paste(
            "Normalization method applied to genes expression counts when identifying RNA based",
            "clusters before calling custom MACS2 peaks. Ignored if '--callpeaks' is not set to",
            "'cluster'.",
            "Default: sct"
        ),
        type="character",
        default="sct",
        choices=c("sct", "log", "sctglm")
    )
    parser$add_argument(
        "--highvargenes",
        help=paste(
            "Number of highly variable genes used in RNA datasets integration, scaling and",
            "dimensionality reduction when identifying RNA based clusters for calling",
            "custom MACS2 peaks. Ignored if '--callpeaks' is not set to 'cluster'.",
            "Default: 3000"
        ),
        type="integer", default=3000
    )
    parser$add_argument(
        "--ntgr",
        help=paste(
            "RNA datasets integration method used for identifying RNA based clusters",
            "before calling custom MACS2 peaks. Automatically set to 'none' if '--mex' points",
            "to the Cell Ranger ARC Count outputs (single, not aggregated dataset that",
            "doesn't require any integration). Ignored if '--callpeaks' is not set to",
            "'cluster'.",
            "Default: seurat"
        ),
        type="character",
        default="seurat",
        choices=c("seurat", "none")
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use in projection and clustering for RNA assay when identifying",
            "RNA based clusters for calling custom MACS2 peaks (from 1 to 50). If single",
            "number N is provided, use from 1 to N PCs. If multiple numbers are provided,",
            "subset to only selected PCs. Ignored if '--callpeaks' is not set to 'cluster'.",
            "Default: from 1 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--resolution",
        help=paste(
            "Resolution to be used when identifying RNA based clusters for calling",
            "custom MACS2 peaks. Ignored if '--callpeaks' is not set to 'cluster'.",
            "Default: 0.3"
        ),
        type="double", default=0.3
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
        "--lowmem",
        help=paste(
            "Attempts to minimize RAM usage when integrating multiple datasets",
            "with SCTransform algorithm (slows down the computation).",
            "Ignored if '--callpeaks' is not set to 'cluster', if '--ntgr' is not set",
            "to 'seurat', if '--norm' is not set to either 'sct' or 'sctglm'.",
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

print(paste("Loading datasets identities from", args$identity))
cell_identity_data <- io$load_cell_identity_data(args$identity)

print(paste("Loading datasets grouping from", args$grouping))
grouping_data <- io$load_grouping_data(args$grouping, cell_identity_data)

print(paste("Loading genomic blacklist regions from", args$blacklist))
blacklist_data <- io$load_blacklist_data(args$blacklist)

print(paste("Loading gene/peak-barcode matrices from", args$mex))
print(paste("Loading fragments from", args$fragments))
print(paste("Loading annotations from", args$annotations))
seurat_data <- io$load_10x_multiome_data(                           # identities are set to the "new.ident" column
    args=args,
    cell_identity_data=cell_identity_data,
    grouping_data=grouping_data
)
debug$print_info(seurat_data, args)

print("Adjusting input parameters")
idents_count <- length(unique(as.vector(as.character(Idents(seurat_data)))))
for (key in names(args)){
    if (key %in% c("mingenes", "maxgenes", "rnaminumi", "minnovelty", "atacminumi", "maxnuclsignal", "mintssenrich", "minfrip", "maxblacklist")){
        if (length(args[[key]]) != 1 && length(args[[key]]) != idents_count){
            print(paste("Filtering parameter", key, "has an ambiguous size. Exiting"))
            quit(save="no", status=1, runLast=FALSE)
        }
        if (length(args[[key]]) == 1){
            print(paste("Extending filtering parameter", key, "to have a proper size"))
            args[[key]] <- rep(args[[key]][1], idents_count)
        }
    }
}
if (length(args$dimensions) == 1) {                            # only one value was provided, so we need to inflate it to 1:N
    args$dimensions <- c(1:args$dimensions[1])
}
print("Adjusted parameters")
print(args)

print(paste("Loading barcodes of interest from", args$barcodes))
barcodes_data <- io$load_barcodes_data(args$barcodes, seurat_data)
print("Applying cell filters based on the loaded barcodes of interest")
seurat_data <- filter$apply_cell_filters(seurat_data, barcodes_data)
debug$print_info(seurat_data, args)

print("Adding RNA QC metrics for not filtered datasets")
seurat_data <- qc$add_rna_qc_metrics(seurat_data, args)
print("Adding ATAC QC metrics for not filtered datasets")
seurat_data <- qc$add_atac_qc_metrics(seurat_data, args)
print("Adding peak QC metrics for peaks called by Cell Ranger ARC")
seurat_data <- qc$add_peak_qc_metrics(seurat_data, blacklist_data, args)
debug$print_info(seurat_data, args)

export_all_qc_plots(
    seurat_data=seurat_data,
    suffix="raw",
    args=args,
    macs2_peaks=FALSE
)

print("Applying filters based on RNA QC metrics")
seurat_data <- filter$apply_rna_qc_filters(seurat_data, cell_identity_data, args)          # cleans up all reductions
debug$print_info(seurat_data, args)
print("Applying filters based on ATAC QC metrics")
seurat_data <- filter$apply_atac_qc_filters(seurat_data, cell_identity_data, args)         # cleans up all reductions
debug$print_info(seurat_data, args)

if (!is.null(args$callpeaks)){
    print("Forced to replace Cell Ranger ARC peaks with MACS2 peaks")
    seurat_data <- call_peaks(seurat_data, args)
    debug$print_info(seurat_data, args)
    print("Updating ATAC QC metrics after calling MACS2 peaks")
    seurat_data <- qc$add_atac_qc_metrics(seurat_data, args)                               # with the new peaks, we have different number of ATAC UMI counted per cell, so all the metrics should be updated
    print("Updating peak QC metrics after calling MACS2 peaks")
    seurat_data <- qc$add_peak_qc_metrics(seurat_data, blacklist_data, args)             # recalculate peak QC metrics
    debug$print_info(seurat_data, args)
    export_all_qc_plots(                                                                   # after RNA and ATAC filters have been applied
        seurat_data=seurat_data,
        suffix="mid_fltr",
        args=args,
        macs2_peaks=TRUE
    )
    if (args$callpeaks == "cluster"){                                                      # we have this plots only when we run RNA clustering for calling peaks per cluster
        export_all_dimensionality_plots(
            seurat_data=seurat_data,
            suffix="mid_fltr",
            args=args
        )
        export_all_clustering_plots(
            seurat_data=seurat_data,
            suffix="mid_fltr",
            args=args,
            cluster_prefix="peak_calling"
        )
    }
    print("Applying filters based on updated ATAC QC metrics after calling MACS2 peaks")
    seurat_data <- filter$apply_atac_qc_filters(seurat_data, cell_identity_data, args)     # cleans up all reductions
    debug$print_info(seurat_data, args)
}

print("Applying filters based on peaks QC metrics")
seurat_data <- filter$apply_peak_qc_filters(seurat_data, cell_identity_data, args)         # cleans up all reductions
print("Updating ATAC QC metrics after all filtering thresholds applied")
seurat_data <- qc$add_atac_qc_metrics(seurat_data, args)                                   # recalculate ATAC QC metrics
print("Updating peak QC metrics after all filtering thresholds applied")
seurat_data <- qc$add_peak_qc_metrics(seurat_data, blacklist_data, args)                 # recalculate peak QC metrics
debug$print_info(seurat_data, args)

export_all_qc_plots(                                                                       # after all filters have been applied
    seurat_data=seurat_data,
    suffix="fltr",
    args=args,
    macs2_peaks=!is.null(args$callpeaks)                                                   # can be both from 10x or MACS2
)

DefaultAssay(seurat_data) <- "RNA"                                                         # better to stick to RNA assay by default https://www.biostars.org/p/395951/#395954 
print("Exporting results to RDS file")
io$export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))
if(args$h5seurat){
    print("Exporting results to h5seurat file")
    io$export_h5seurat(seurat_data, paste(args$output, "_data.h5seurat", sep=""))
}

if(args$h5ad){
    print("Exporting results to h5ad file")
    io$export_h5ad(seurat_data, paste(args$output, "_data.h5ad", sep=""))
}