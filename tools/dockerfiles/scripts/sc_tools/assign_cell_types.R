#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(future))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(reticulate))
suppressMessages(library(sctransform))
suppressMessages(library(RColorBrewer))

# get not exported functions from Seurat package
IsSCT <- utils::getFromNamespace("IsSCT", "Seurat")        # returns TRUE if assay was SCT normalized


debug_seurat_data <- function(seurat_data, args) {
    if (args$verbose){
        print("Assays")
        print(seurat_data@assays)
        print("Reductions")
        print(seurat_data@reductions)
        print("Active assay")
        print(DefaultAssay(seurat_data))
        print("Metadata")
        print(head(seurat_data@meta.data))
        print(gc())
    }
}


setup_parallelization <- function (args) {
    print(
        paste(
            "Setting parallelization parameters to", args$cpus,
            "cores, and", args$memory, "GB of memory"
        )
    )
    invisible(capture.output(plan("multiprocess", workers=args$cpus)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(args$cpus)))
    options(future.globals.maxSize = args$memory * 1024^3)               # convert to bytes
}


export_rds <- function(data, location){
    tryCatch(
        expr = {
            saveRDS(data, location)
            print(paste("Export data as RDS to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data as RDS to", location, sep=" "))
        }
    )
}


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


load_ctype_data <- function (location) {
    print(paste("Loading cell types data from", location))
    ctype_data_data <- read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    return (ctype_data_data)
}


writeSparseMatrix <- function (inMat, outFname, sliceSize=1000){
    fnames <- c()
    mat <- inMat
    geneCount <- nrow(mat)
    startIdx <- 1
    while (startIdx < geneCount){
        endIdx <- min(startIdx + sliceSize - 1, geneCount)
        matSlice <- mat[startIdx:endIdx,]
        denseSlice <- as.matrix(matSlice)
        dt <- data.table(denseSlice)
        dt <- cbind(gene=rownames(matSlice), dt)
        writeHeader <- FALSE
        if (startIdx == 1) { 
            writeHeader <- TRUE
        }
        sliceFname <- paste0("temp", startIdx, ".txt")
        fwrite(dt, sep="\t", file=sliceFname, quote=FALSE, col.names=writeHeader)
        fnames <- append(fnames, sliceFname)
        startIdx <- startIdx + sliceSize
    }
    system(paste("cat", paste(fnames, collapse=" "), "| gzip >", outFname, sep=" "))
    unlink(fnames)
}


findMatrix <- function(object, matrix.slot){
    if (matrix.slot == "counts") {
        counts <- GetAssayData(object=object, slot="counts")
    } else if (matrix.slot == "scale.data") {
        counts <- GetAssayData(object=object, slot="scale.data")
    } else if (matrix.slot == "data") {
        counts <- GetAssayData(object=object)
    } else {
        error("matrix.slot can only be one of: counts, scale.data, data")
    }
}


ExportToCellbrowser <- function(
    object,
    matrix.slot,
    dir,
    cluster.field,
    features,
    meta.fields,
    meta.fields.names
) {
    if (!dir.exists(dir)) {
        dir.create(dir)
    }
    idents <- Idents(object)
    meta <- object@meta.data
    cellOrder <- colnames(object)
    counts <- findMatrix(object, matrix.slot)
    genes <- rownames(x=object)
    dr <- object@reductions

    gzPath <- file.path(dir, "exprMatrix.tsv.gz")
    if ((((ncol(counts)/1000)*(nrow(counts)/1000))>2000) && is(counts, 'sparseMatrix')){
        writeSparseMatrix(counts, gzPath);
    } else {
        mat <- as.matrix(counts)
        df <- as.data.frame(mat, check.names=FALSE)
        df <- data.frame(gene=genes, df, check.names=FALSE)
        z <- gzfile(gzPath, "w")
        write.table(x=df, sep="\t", file=z, quote=FALSE, row.names=FALSE)
        close(con=z)
    }
    embeddings = names(dr)
    embeddings.conf <- c()
    for (embedding in embeddings) {
        emb <- dr[[embedding]]
        df <- emb@cell.embeddings
        if (ncol(df) > 2){
            df <- df[, 1:2]
        }
        colnames(df) <- c("x", "y")
        df <- data.frame(cellId=rownames(df), df, check.names=FALSE)
        write.table(
            df[cellOrder,],
            sep="\t",
            file=file.path(dir, sprintf("%s.coords.tsv", embedding)),
            quote=FALSE,
            row.names=FALSE
        )
        embeddings.conf <- c(
            embeddings.conf,
            sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}', embedding)
        )
    }
    df <- data.frame(row.names=cellOrder, check.names=FALSE)
    enum.fields <- c()
    for (i in 1:length(meta.fields)){
        field <- meta.fields[i]
        name <- meta.fields.names[i]
        if (field %in% colnames(meta)){
            df[[name]] <- meta[[field]]
            if (!is.numeric(df[[name]])) {
                enum.fields <- c(enum.fields, name)
            }
        }
    }
    df <- data.frame(Cell=rownames(df), df, check.names=FALSE)
    write.table(
        as.matrix(df[cellOrder,]),
        sep="\t",
        file=file.path(dir, "meta.tsv"),
        quote=FALSE,
        row.names=FALSE
    )
    if (!is.null(features)){
        write.table(
            features,
            sep="\t",
            file=file.path(dir, "quickGenes.tsv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE
        )
    }
    config <- '
name="cellbrowser"
shortLabel="cellbrowser"
geneLabel="Gene"
exprMatrix="exprMatrix.tsv.gz"
meta="meta.tsv"
radius=3
alpha=0.5
geneIdType="auto"
clusterField="%s"
labelField="%s"
enumFields=%s
coords=%s'
    enum.string <- paste0("[", paste(paste0('"', enum.fields, '"'), collapse=", "), "]")
    coords.string <- paste0("[", paste(embeddings.conf, collapse = ",\n"), "]")
    config <- sprintf(
        config,
        cluster.field,
        cluster.field,
        enum.string,
        coords.string
    )
    if (!is.null(features)){
        config <- paste(config, 'quickGenesFile="quickGenes.tsv"', sep="\n")
    }
    confPath = file.path(dir, "cellbrowser.conf")
    cat(config, file=confPath)
    cb.dir <- paste(dir, "html_data", sep="/")
    cb <- import(module = "cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
}


export_cellbrowser_data <- function(seurat_data, suffix, args, assay="RNA", matrix_slot="data"){
    tryCatch(
        expr = {
            backup_assay <- DefaultAssay(seurat_data)
            DefaultAssay(seurat_data) <- assay
            meta_fields <- c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "Phase", "S.Score", "G2M.Score", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal", "frip", "blacklisted_fraction", "miQC.keep")
            meta_fields_names <- c("GEX UMIs", "Genes", "Mitochondrial %", "Novelty score", "Cell cycle phase", "S score", "G to M score", "ATAC UMIs", "Peaks", "TSS enrichment score", "Nucleosome signal", "FRiP", "Blacklisted fraction", "miQC prediction")

            if (!is.null(args$querysource)){                             # run in --ctype mode or used --downsample with --ref, so user knows which clustering resolution to use
                meta_fields <- c(meta_fields, args$querysource)
                meta_fields_names <- c(meta_fields_names, "Clustering")
            } else {                                                     # run in --ref mode, so user wants to compare automatically assigned cell types to all resolutions
                tryCatch(
                    expr = {
                        res_meta_fields <- grep("_res\\.", colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
                        res_fields_names <- paste(
                            "Clustering (", map(strsplit(res_meta_fields, "_res\\."), ~.x[2]) %>% unlist(), ")",
                            sep=""
                        )
                        meta_fields <- c(
                            meta_fields,
                            res_meta_fields
                        )
                        meta_fields_names <- c(
                            meta_fields_names,
                            res_fields_names
                        )
                    },
                    error = function(e){
                        print(paste("Failed to add clustering fields due to", e))
                    }
                )
            }
            meta_fields <- c(meta_fields, args$querytarget)
            meta_fields_names <- c(meta_fields_names, "Cell Type")

            Idents(seurat_data) <- "condition"
            if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
                meta_fields <- c("condition", meta_fields)
                meta_fields_names <- c("Condition", meta_fields_names)
            }

            Idents(seurat_data) <- "new.ident"
            if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
                meta_fields <- c("new.ident", meta_fields)
                meta_fields_names <- c("Identity", meta_fields_names)
            }

            # we do not want to lose custom columns added by run_seurat_wnn.R
            tryCatch(
                expr = {
                    custom_fields <- grep("^custom_", colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
                    custom_fields_names <- gsub("custom_", "", custom_fields)
                    meta_fields <- c(
                        meta_fields,
                        custom_fields
                    )
                    meta_fields_names <- c(
                        meta_fields_names,
                        custom_fields_names
                    )
                },
                error = function(e){
                    print(paste("Failed to add extra metadata fields due to", e))
                }
            )

            rootname <- paste(args$output, suffix, sep="_")
            ExportToCellbrowser(
                seurat_data,
                matrix.slot=matrix_slot,
                dir=rootname,
                cluster.field="Cell Type",
                features=args$features,
                meta.fields=meta_fields,
                meta.fields.names=meta_fields_names
            )
            print(paste("Export UCSC Cellbrowser data to", rootname, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export UCSC Cellbrowser data to", rootname, sep=" "))
        },
        finally = {
            DefaultAssay(seurat_data) <- backup_assay
        }
    )
}


export_dim_plot <- function(data, rootname, reduction, plot_title, legend_title, split_by=NULL, group_by=NULL, perc_split_by=NULL, perc_group_by=NULL, label=FALSE, palette=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- DimPlot(
                        data,
                        reduction=reduction,
                        split.by=split_by,
                        group.by=group_by,
                        label=label
                    ) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    guides(color=guide_legend(legend_title, override.aes=list(size=3)))

            if (!is.null(palette)){ plot <- plot + scale_color_brewer(palette=palette) }

            if(!is.null(perc_split_by) && !is.null(perc_group_by)){
                width <- 2 * width
                perc_data <- data@meta.data %>%
                             dplyr::group_by(across(all_of(perc_split_by)), across(all_of(perc_group_by))) %>%
                             summarise(counts=n(), .groups="drop_last") %>%               # drop the perc_group_by grouping level so we can get only groups defined by perc_split_by
                             mutate(freq=counts/sum(counts)*100) %>%                      # sum is taken for the group defined by perc_split_by
                             ungroup() %>%
                             complete(
                                 (!!as.symbol(perc_split_by)), (!!as.symbol(perc_group_by)),
                                 fill=list(counts=0, freq=0)
                             ) %>%
                             arrange(all_of(perc_split_by), all_of(perc_group_by))        # sort for consistency
                label_data <- data@meta.data %>%
                              dplyr::group_by(across(all_of(perc_split_by))) %>%
                              summarise(counts=n(), .groups="drop_last") %>%              # drops all grouping as we have only one level
                              arrange(all_of(perc_split_by))                              # sort for consistency
                perc_plot <- ggplot(perc_data, aes_string(x=perc_split_by, y="freq", fill=perc_group_by)) +
                             geom_col(position="dodge", width=0.9, linetype="solid", color="black", show.legend=FALSE) +
                             xlab("") +
                             ylab("Cells percentage") +
                             theme_gray() +
                             geom_label_repel(
                                label_data, mapping=aes(y=-Inf, label=counts),
                                color="black", fill="white", segment.colour=NA,
                                direction="y", size=3, show.legend=FALSE
                             ) +
                             RotatedAxis()
                plot <- plot + perc_plot
            }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export dim plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export dim plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_elbow_plot <- function(data, rootname, plot_title, ndims=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- ElbowPlot(
                        data,
                        ndims=ndims,
                        reduction="pca"
                    ) +
                    theme_gray() +
                    ggtitle(plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export elbow plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export elbow plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_depth_corr_plot <- function(data, rootname, plot_title, assay, reduction, ndims=NULL, pdf=FALSE, width=1200, height=400, resolution=100){
    tryCatch(
        expr = {
            plot <- DepthCor(
                        data,
                        assay=assay,
                        reduction=reduction,
                        n=ndims
                    ) +
                    theme_gray() +
                    ggtitle(plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export depth correlation plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export depth correlation plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_all_clustering_plots <- function(seurat_data, suffix, args) {
    default_umap <- ifelse("wnnumap" %in% names(seurat_data@reductions), "wnnumap", "umap")  # we are interested in WNN if the input data was generated by run_seurat_wnn.R
    export_dim_plot(
        data=seurat_data,
        reduction=default_umap,
        plot_title="Query UMAP grouped by cell type",
        legend_title="Cell type",
        group_by=args$querytarget,
        perc_split_by="new.ident",
        perc_group_by=args$querytarget,
        label=TRUE,
        rootname=paste(args$output, suffix, "umap_ctype", sep="_"),
        pdf=args$pdf
    )
    export_dim_plot(
        data=seurat_data,
        reduction=default_umap,
        plot_title="Query UMAP split by condition grouped by cell type",
        legend_title="Cell type",
        split_by="condition",
        group_by=args$querytarget,
        perc_split_by="condition",
        perc_group_by=args$querytarget,
        label=TRUE,
        rootname=paste(args$output, suffix, "umap_ctype_spl_by_cond", sep="_"),
        pdf=args$pdf
    )
}


export_vln_plot <- function(data, features, labels, rootname, plot_title, legend_title, log=FALSE, group_by=NULL, hide_x_text=FALSE, pt_size=NULL, palette=NULL, combine_guides=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots <- VlnPlot(
                         data,
                         features=features,
                         pt.size=pt_size,
                         group.by=group_by,
                         log=log,
                         combine=FALSE       # to return a list of gglots
                     )
            plots <- lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] +
                              ggtitle(labels[i]) +
                              theme_gray() +
                              theme(axis.title.x=element_blank()) +
                              guides(fill=guide_legend(legend_title)) +
                              stat_boxplot(width=0.15, geom="errorbar") +
                              geom_boxplot(width=0.15, outlier.alpha=0) +
                              RotatedAxis()
                if (!is.null(palette)){ plots[[i]] <- plots[[i]] + scale_fill_brewer(palette=palette) }
                if (hide_x_text){ plots[[i]] <- plots[[i]] + theme(axis.text.x=element_blank()) }
                return (plots[[i]])
            })
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export violin plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export violin plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_feature_plot <- function(data, features, labels, rootname, reduction, plot_title, split_by=NULL, label=FALSE, order=FALSE, min_cutoff=NA, max_cutoff=NA, pt_size=NULL, combine_guides=NULL, alpha=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots <- FeaturePlot(
                        data,
                        features=features,
                        pt.size=pt_size,
                        order=order,
                        min.cutoff=min_cutoff,
                        max.cutoff=max_cutoff,
                        reduction=reduction,
                        split.by=split_by,
                        label=label,
                        combine=FALSE       # to return a list of gglots
                    )
            plots <- lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] + ggtitle(labels[i]) + theme_gray()
                if (!is.null(alpha)) { plots[[i]]$layers[[1]]$aes_params$alpha <- alpha }
                return (plots[[i]])
            })
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export feature plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export feature plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_dot_plot <- function(data, features, rootname, plot_title, x_label, y_label, cluster_idents=FALSE, min_pct=0.01, col_min=-2.5, col_max=2.5, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- DotPlot(
                        data,
                        features=features,
                        cluster.idents=cluster_idents,
                        dot.min=min_pct,
                        col.min=col_min,
                        col.max=col_max,
                        scale=TRUE,
                        scale.by="size"  # for optimal perception
                    ) +
                    xlab(x_label) +
                    ylab(y_label) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    RotatedAxis()

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export dot plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export dot plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_all_dimensionality_plots <- function(seurat_data, suffix, assay, args) {
    export_elbow_plot(
        data=seurat_data,
        ndims=length(seurat_data@reductions[[args$refreduction]]),
        rootname=paste(args$output, suffix, "elbow", sep="_"),
        plot_title=paste(
            "Elbow plot for", args$refreduction, "reduction of the reference dataset"
        ),
        pdf=args$pdf
    )
    export_depth_corr_plot(
        data=seurat_data,
        assay=assay,
        reduction=args$refreduction,
        plot_title=paste(
            "Correlation plot between depth and", args$refreduction,
            "reduced dimension components of", assay, "assay"
        ),
        rootname=paste(args$output, suffix, "depth_corr", sep="_"),
        pdf=args$pdf
    )
}


export_all_expression_plots <- function(seurat_data, suffix, args, assay="RNA") {
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    Idents(seurat_data) <- args$querytarget
    export_dot_plot(
        data=seurat_data,
        features=args$features,
        plot_title="Scaled average log normalized gene expression per predicted cell type of query dataset",
        x_label="Genes",
        y_label="Cell types",
        cluster_idents=FALSE,  # no need to cluster cell types together
        rootname=paste(args$output, suffix, "avg_per_ctype", sep="_"),
        pdf=args$pdf
    )
    export_feature_plot(
        data=seurat_data,
        features=args$features,
        labels=args$features,
        reduction="umap",
        plot_title="Log normalized gene expression per cell of query dataset with predicted cell types",
        label=TRUE,
        order=TRUE,
        max_cutoff="q99",  # to prevent cells with overexpressed gene from distorting the color bar
        rootname=paste(args$output, suffix, "per_ctype_cell", sep="_"),
        combine_guides="keep",
        pdf=args$pdf
    )
    export_vln_plot(
        data=seurat_data,
        features=args$features,
        labels=args$features,
        plot_title="Log normalized gene expression densities per predicted cell type of query dataset",
        legend_title="Cell type",
        log=TRUE,
        pt_size=0,
        combine_guides="collect",
        rootname=paste(args$output, suffix, "dnst_per_ctype", sep="_"),
        pdf=args$pdf
    )
    Idents(seurat_data) <- "new.ident"
    DefaultAssay(seurat_data) <- backup_assay
}


get_args <- function(){
    parser <- ArgumentParser(description="Assigns cell types to clusters")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the query RDS file to load Seurat object from.",
            "RDS file can be produced by run_seurat.R or",
            "run_seurat_wnn.R scripts."
        ),
        type="character", required="True"
    )
    group <- parser$add_mutually_exclusive_group(required="True")
    group$add_argument(
        "--ref",
        help=paste(
            "Path to the reference RDS file for automatic cell type",
            "assignment. Either --ref or --ctype parameters should be",
            "provided."
        ),
        type="character"
    )
    group$add_argument(
        "--ctype",
        help=paste(
            "Path to the cell types metadata TSV/CSV file for manual",
            "cell type assignment. Either --ref or --ctype parameters",
            "should be provided."
        ),
        type="character"
    )
    parser$add_argument(
        "--refsource",
        help=paste(
            "Metadata column from the reference RDS file to select",
            "labels to transfer. Required if running with --ref"
        ),
        type="character"
    )
    parser$add_argument(
        "--refreduction",
        help=paste(
            "Reduction name from the reference RDS file to transfer",
            "embeddings from. Required if running with --ref"
        ),
        type="character"
    )
    parser$add_argument(
        "--querysource",
        help=paste(
            "Metadata column from the query RDS file to select clusters",
            "for manual cell type assignment. Required if running with",
            "--ctype. When running with --ref and --downsample parameter",
            "was provided, defines clusters each of which will be dowsampled",
            "to a fixed number of cells."
        ),
        type="character"
    )
    parser$add_argument(
        "--querytarget",
        help=paste(
            "Column to be added to metadata of the query RDS file for",
            "automatically or manually assigned cell types."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--highvarcount",
        help=paste(
            "Number of highly variable features to detect. Used for anchors",
            "detection between the reference and query datasets.",
            "Default: 3000"
        ),
        type="integer", default=3000
    )
    parser$add_argument(
        "--regressmt",
        help=paste(
            "Regress the percentage of transcripts mapped to mitochondrial genes",
            "as a confounding source of variation. Applied for query dataset only.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--ndim",
        help=paste(
            "Dimensions from the specified with --refreduction reference dataset's",
            "reduction that should be used in label transfer. If single value N is",
            "provided, use from 1 to N dimensions. If multiple values are provided,",
            "subset to only selected dimensions. Ignored when running with --ctype.",
            "Default: from 1 to 30"
        ),
        type="integer", default=30, nargs="*"
    )
    parser$add_argument(
        "--features",
        help=paste(
            "Features of interest to highlight.",
            "Default: None"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--verbose",
        help=paste(
            "Print debug information.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help=paste(
            "Output prefix for all generated files.",
            "Default: ./seurat"
        ),
        type="character", default="./seurat"
    )
    parser$add_argument(
        "--pdf",
        help=paste(
            "Export plots in PDF format in additions to PNG.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--cpus",
        help=paste(
            "Number of cores/cpus to use.",
            "Default: 1"
        ),
        type="integer", default=1
    )
    parser$add_argument(
        "--downsample",
        help=paste(
            "Downsample query dataset to include a fixed number of",
            "cells per cluster defined by --querysource parameter.",
            "Applied only when running automatic cell type assignment",
            "with both --ref and --querysource parameters provided.",
            "Default: do not downsample"
        ),
        type="integer"
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
setup_parallelization(args)

print(paste("Loading query Seurat data from", args$query))
query_seurat_data <- readRDS(args$query)
DefaultAssay(query_seurat_data) <- "RNA"                                  # safety measure in case we forgot to explicitly set RNA assay
debug_seurat_data(query_seurat_data, args)

if (!is.null(args$ctype)){
    print("Attempting to run manual cell types assignment")
    if(is.null(args$querysource)){
        print(paste("Missing required --querysource parameter"))
        quit(save="no", status=1, runLast=FALSE)
    }
    ctype_data <- load_ctype_data(args$ctype)
    print(
        paste(
            "Searching for clusters in", args$querysource, "column,",
            "adding cell types to", args$querytarget, "column"
        )
    )
    query_seurat_data[[args$querytarget]] <- ctype_data$type[
        match(
            as.vector(as.character(query_seurat_data@meta.data[, args$querysource])),
            as.vector(as.character(ctype_data$cluster))
        )
    ]
    debug_seurat_data(query_seurat_data, args)
} else {                                                                                        # no need to check for args$ref as we use required mutually exclusive group
    print("Attempting to run automatic cell types assignment")
    if(is.null(args$refsource) | is.null(args$refreduction)){
        print(paste("Missing --refsource and/or --refreduction requred parameters"))
        quit(save="no", status=1, runLast=FALSE)
    }
    if (length(args$ndim) == 1) {                                                               # only one value was provided, so we need to inflate it to 1:N
        args$ndim <- c(1:args$ndim[1])
    }
    if(!is.null(args$querysource) && !is.null(args$downsample)){
        Idents(query_seurat_data) <- args$querysource
        print(
            paste(
                "Downsampling query dataset to include only", args$downsample,
                "cells per cluster defined by", args$querysource, "column"
            )
        )
        query_seurat_data <- subset(query_seurat_data, downsample=args$downsample)
        debug_seurat_data(query_seurat_data, args)
        Idents(query_seurat_data) <- "new.ident"
    }
    print(paste("Loading reference Seurat data from", args$ref))
    ref_seurat_data <- readRDS(args$ref)
    debug_seurat_data(ref_seurat_data, args)
    ref_assay <- DefaultAssay(ref_seurat_data@reductions[[args$refreduction]])                  # get reference assay name based on the selected reduction
    ref_norm_method <- ifelse(IsSCT(ref_seurat_data[[ref_assay]]), "SCT", "LogNormalize")
    export_all_dimensionality_plots(ref_seurat_data, "ref", ref_assay, args)                    # for exploratory purposes to set --ndim
    print(
        paste0(
            "For ", args$refreduction, " reference reduction run on ", ref_assay, " assay ",
            "the normalizations method was ", ref_norm_method, ". Query dataset will be",
            "split by new.ident and re-normalized with ", ref_norm_method, "."
        )
    )
    Idents(query_seurat_data) <- "new.ident"
    query_identities <- unique(as.vector(as.character(Idents(query_seurat_data))))
    merged_seurat_data <- NULL                                                                  # to collect re-normalized datasets together
    vars_to_regress <- NULL
    if (args$regressmt) {
        print(
            paste(
                "Configuring to regress the percentage of transcripts",
                "mapped to mitochondrial genes as a confounding source",
                "of variation. Will be applied only for query dataset."
            )
        )
        vars_to_regress <- c("mito_percentage")
    }
    for (i in 1:length(query_identities)){
        current_identity <- query_identities[i]
        print(paste("Processing", current_identity))
        filtered_data <- subset(query_seurat_data, idents=current_identity)
        if (ref_norm_method == "LogNormalize"){
            print(
                paste(
                    "Running LogNormalize on RNA assay using", args$highvarcount,
                    "highly variable features"
                )
            )
            filtered_data <- NormalizeData(
                filtered_data,
                assay="RNA",                                                                    # Always use RNA assay
                verbose=FALSE
            )
            filtered_data <- FindVariableFeatures(
                filtered_data,
                nfeatures=args$highvarcount,
                verbose=FALSE,
                assay="RNA"
            )
            filtered_data <- ScaleData(
                filtered_data,
                vars.to.regress=vars_to_regress,
                verbose=FALSE,
                assay="RNA"
            )
        } else {
            print(
                paste(
                    "Running SCTransform on RNA assay using", args$highvarcount,
                    "highly variable features"
                )
            )
            filtered_data <- SCTransform(
                filtered_data,
                assay="RNA",                                                                    # still need to use RNA as SCTransform adds data to SCT assay based on RNA
                new.assay.name="SCT",                                                           # this is the default name, but we set it for better understanding of the logic
                variable.features.n=args$highvarcount,
                vars.to.regress=vars_to_regress,
                verbose=FALSE
            )
        }
        debug_seurat_data(filtered_data, args)
        print(
            paste(
                "Searching for transfer anchors using", paste(args$ndim, collapse=", "),
                "dimensions from the reference dataset"
            )
        )
        transfer_anchors <- FindTransferAnchors(
            reference=ref_seurat_data,
            query=filtered_data,
            normalization.method=ref_norm_method,
            reference.assay=ref_assay,
            features=as.vector(as.character(rownames(ref_seurat_data))),
            query.assay=ifelse(ref_norm_method == "LogNormalize", "RNA", "SCT"),
            reference.reduction=args$refreduction,
            dims=args$ndim
        )
        print(
            paste(
                "Transfering data using identified transfer anchors and",
                args$refsource, "column from the reference dataset.",
                "Adding predicted cell types to", args$querytarget, "column",
                "of the query dataset"
            )
        )
        predictions <- TransferData(
            anchorset=transfer_anchors,
            reference=ref_seurat_data,
            refdata=args$refsource
        )
        filtered_data[[args$querytarget]] <- predictions$predicted.id
        if (is.null(merged_seurat_data)){
            merged_seurat_data <- filtered_data
        } else {
            merged_seurat_data <- merge(merged_seurat_data, y=filtered_data)
        }
        rm(
            filtered_data,
            transfer_anchors,
            predictions,
            current_identity
        )
        gc()
    }
    debug_seurat_data(query_seurat_data, args)
    merged_seurat_data@reductions <- query_seurat_data@reductions
    query_seurat_data <- merged_seurat_data                                    # for consistency with the --ctype mode
    rm(ref_seurat_data, merged_seurat_data)
    gc()
}

debug_seurat_data(query_seurat_data, args)

if (!is.null(args$features)){
    print("Check genes of interest to include only those that are present in the datasets")
    backup_assay <- DefaultAssay(query_seurat_data)
    DefaultAssay(query_seurat_data) <- "RNA"
    args$features <- unique(args$features)
    args$features <- args$features[args$features %in% as.vector(as.character(rownames(query_seurat_data)))]
    print(args$features)
    DefaultAssay(query_seurat_data) <- backup_assay
}

print("Saving results")
export_rds(query_seurat_data, paste(args$output, "_ctype_data.rds", sep=""))
export_all_clustering_plots(query_seurat_data, "query", args)
if (!is.null(args$features)){
    export_all_expression_plots(query_seurat_data, "query_expr", args, assay="RNA")
}

print("Exporting UCSC Cellbrowser data")
export_cellbrowser_data(
    seurat_data=query_seurat_data,
    suffix="cellbrowser",
    args=args,
    assay="RNA",
    matrix_slot="counts"
)
