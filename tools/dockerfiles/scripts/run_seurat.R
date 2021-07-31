#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB should be good by default

suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(DESeq2))
suppressMessages(library(tibble))
suppressMessages(library(scales))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(garnett))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(reticulate))
suppressMessages(library(sctransform))
suppressMessages(library(RColorBrewer))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))


####################################################################################
# For details and treshold values look at
#    https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html
#    https://github.com/satijalab/seurat/issues/1679
#    https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html
#    https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/cellbrowser.html
#    https://github.com/satijalab/seurat-wrappers
####################################################################################

# To use Garnett classifiers trained with Ensembl gene names
SPECIES_DATA <- list(
    "hs"="org.Hs.eg.db",
    "mm"="org.Mm.eg.db",
    "none"="none"          # Garnett will not check/convert gene IDs, so your CDS and marker file must have the same gene ID type
)

###############################################################################################################
# Modified code from
# https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/R/ExportToCellbrowser-seurat.R
###############################################################################################################


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
exprMatrix="exprMatrix.tsv.gz"
meta="meta.tsv"
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


###############################################################################################################


set_threads <- function (threads) {
    invisible(capture.output(plan("multiprocess", workers=threads)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(threads)))
}


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


load_cell_cycle_data <- function (location) {
    if (!is.null(location)){
        cell_cycle_data <- read.table(
            location,
            sep=get_file_type(location),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        print(paste("Cell cycle data is successfully loaded from ", location))
        return (cell_cycle_data)
    }
    print("Cell cycle data is not provided")
    return (NULL)
}


load_cell_identity_data <- function (location) {
    cell_identity_data <- read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    # prepend with LETTERS, otherwise the order on the plot will be arbitrary sorted
    cell_identity_data <- cell_identity_data %>% mutate("library_id"=paste(LETTERS[1:nrow(cell_identity_data)], .$library_id))
    return (cell_identity_data)
}


load_condition_data <- function(location, cell_identity_data) {
    default_condition_data <- data.frame(
        library_id=cell_identity_data$library_id,
        condition=cell_identity_data$library_id,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    if (!is.null(location)){
        condition_data <- read.table(
            location,
            sep=get_file_type(location),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        # prepend with LETTERS to correspond to the library_id from the cell_identity_data
        condition_data <- condition_data %>% mutate("library_id"=paste(LETTERS[1:nrow(condition_data)], .$library_id))
        if ( (nrow(condition_data) == nrow(cell_identity_data)) && all(is.element(cell_identity_data$library_id, condition_data$library_id)) ){
            print(paste("Condition data is successfully loaded from ", location))
            return (condition_data)
        } else {
            print(paste("Applying defaults - condition data loaded from", location, "is malformed"))
            return (default_condition_data)
        }
    }
    print("Condition data is not provided. Applying defaults")
    return (default_condition_data)
}


load_barcodes_data <- function(location, seurat_data) {
    default_barcodes_data <- Cells(seurat_data)                               # to include all available cells
    if (!is.null(location)){
        barcodes_data <- read.table(
            location,
            sep=get_file_type(location),
            header=FALSE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )[,1]                                                                 # to get it as vector
        print(paste("Barcodes data is successfully loaded from ", location))
        return (barcodes_data)
    }
    print("Barcodes data is not provided. Using all cells")
    return (default_barcodes_data)
}


load_seurat_data <- function(locations, mincells, cell_identity_data, condition_data) {
    if (length(locations) == 1){  # Cell Ranger Aggr or single Cell Ranger Count experiments
        current_location <- locations
        if (length(as.vector(cell_identity_data$library_id)) > 1) {
            print(paste("Read aggregated 10x data from", current_location, "using the original barcode suffixes"))
        } else {
            print(paste("Read not aggregated 10x data from", current_location, "using the original barcode suffixes"))
        }
        print(paste("  applying --mincells filtering parameter", mincells))
        seurat_data <- CreateSeuratObject(
            counts=Read10X(data.dir=current_location),
            min.cells=mincells,
            names.delim="-",  # to get cell identity index from Cell Ranger Aggr or Count output
            names.field=2
        )
        idents <- as.numeric(as.character(Idents(seurat_data)))              # need to properly convert factor to numeric vector
        new_ident <- cell_identity_data$library_id[idents]
        if (sum(is.na(new_ident)) > 0){
            print("Identity file includes less than expected number of rows. Exiting.")
            quit(save="no", status=1, runLast=FALSE)
        }
        seurat_data[["new.ident"]] <- new_ident
        seurat_data[["condition"]] <- condition_data$condition[match(seurat_data$new.ident, condition_data$library_id)]
        Idents(seurat_data) <- "new.ident"
        if ( nrow(cell_identity_data) > length(unique(as.vector(as.character(Idents(seurat_data))))) ){
            print("Identity file includes more than expected number of rows. Exiting.")
            quit(save="no", status=1, runLast=FALSE)
        }
        return (seurat_data)
    } else {                      # multiple Cell Ranger Count experiments
        merged_seurat_data <- NULL
        for (i in 1:length(locations)){
            current_location <- locations[i]
            print(paste("Read not aggregated 10x data from", current_location, "replacing the original barcode suffixes with", i))
            print(paste(
                "  ignoring --mincells filtering parameter as it can't be applied for individual datasets",
                "due to its influence on the genes number"
            ))
            seurat_data <- CreateSeuratObject(
                counts=Read10X(data.dir=current_location, strip.suffix=TRUE),  # removes suffix from barcode
            )
            idents <- i
            new_ident <- cell_identity_data$library_id[idents]
            if (sum(is.na(new_ident)) > 0){
                print("Identity file includes less than expected number of rows. Exiting.")
                quit(save="no", status=1, runLast=FALSE)
            }
            seurat_data[["new.ident"]] <- new_ident
            seurat_data[["condition"]] <- condition_data$condition[match(seurat_data$new.ident, condition_data$library_id)]
            Idents(seurat_data) <- "new.ident"
            seurat_data <- RenameCells(seurat_data, new.names=paste0(Cells(seurat_data), "-", idents))  # to add new barcode suffix
            if (is.null(merged_seurat_data)){
                merged_seurat_data <- seurat_data
            } else {
                merged_seurat_data <- merge(merged_seurat_data, y=seurat_data)
            }
        }
        if ( nrow(cell_identity_data) > length(unique(as.vector(as.character(Idents(merged_seurat_data))))) ){
            print("Identity file includes more than expected number of rows. Exiting.")
            quit(save="no", status=1, runLast=FALSE)
        }
        return (merged_seurat_data)
    }
}


load_classifier <- function(location) {
    if (!is.null(location)){
        classifier <- readRDS(location)
        print(paste("Classifier is successfully loaded from ", location))
        return (classifier)
    }
    print("Classifier is not provided. Skip cell type prediction")
    return (NULL)
}


apply_cell_filters <- function(seurat_data, barcodes_data) {
    filtered_seurat_data <- subset(seurat_data, cells=barcodes_data)
    return (filtered_seurat_data)
}


add_qc_metrics <- function(seurat_data, mitopattern) {
    # Number of genes detected per UMI: this metric with give us an idea of the
    # complexity of our dataset (more genes detected per UMI, more complex our data)
    # Mitochondrial percentage: this metric will give us a percentage of cell reads
    # originating from the mitochondrial genes
    seurat_data$log10_gene_per_log10_umi <- log10(seurat_data$nFeature_RNA) / log10(seurat_data$nCount_RNA)
    seurat_data$mito_percentage <- PercentageFeatureSet(seurat_data, pattern=mitopattern) 
    return (seurat_data)
}


apply_qc_filters <- function(seurat_data, cell_identity_data, args) {
    merged_seurat_data <- NULL
    for (i in 1:length(args$minfeatures)){
        identity <- cell_identity_data$library_id[i]
        minfeatures <- args$minfeatures[i]
        maxfeatures <- args$maxfeatures[i]
        minumi <- args$minumi[i]
        minnovelty <- args$minnovelty[i]
        print(paste("Filtering", identity))
        print(paste(" ", minfeatures, "<= Genes per cell <=", maxfeatures))
        print(paste(" ", "UMIs per cell >=", minumi))
        print(paste(" ", "Novelty score >=", minnovelty))
        print(paste(" ", "Mitochondrial gene expression <=", args$maxmt))
        filtered_seurat_data <- NULL
        tryCatch(
            expr = {
                filtered_seurat_data <- subset(
                    seurat_data,
                    idents=identity,
                    subset=(nFeature_RNA >= minfeatures) &
                        (nFeature_RNA <= maxfeatures) &
                        (nCount_RNA >= minumi) &
                        (log10_gene_per_log10_umi >= minnovelty) &
                        (mito_percentage <= args$maxmt)
                )
            },
            error = function(e){
                print(paste(" ", "Failed to apply QC filters for", identity, "due to", e))
            }
        )
        if (is.null(filtered_seurat_data)){ next }
        if (is.null(merged_seurat_data)){
            merged_seurat_data <- filtered_seurat_data
        } else {
            merged_seurat_data <- merge(merged_seurat_data, y=filtered_seurat_data)
        }
    }
    return (merged_seurat_data)
}


explore_unwanted_variation <- function(seurat_data, cell_cycle_data, args) {
    temp_seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
    tryCatch(
        expr = {
            temp_seurat_data <- CellCycleScoring(
                temp_seurat_data,
                s.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="s", "gene_id"]),
                g2m.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="g2/m", "gene_id"])
            )
        },
        error = function(e){
            print("Failed to run cell cycle scoring when exploring cell cycle phase as a source of unwanted variation")
        }
    )
    tryCatch(
        expr = {
            mito_quartiles <- quantile(temp_seurat_data@meta.data$mito_percentage, c(0.25, 0.5, 0.75))
            temp_seurat_data@meta.data$mito_factor <- cut(
                temp_seurat_data@meta.data$mito_percentage, 
                breaks=c(-Inf, mito_quartiles[1], mito_quartiles[2], mito_quartiles[3], Inf), 
                labels=c("Low", "Medium", "Medium high", "High")
            )
        },
        error = function(e){
            print("Failed to run mito factor scoring when exploring mitochondrial gene expression as a source of unwanted variation")
        }
    )
    temp_seurat_data <- FindVariableFeatures(temp_seurat_data, verbose=FALSE)
    temp_seurat_data <- ScaleData(temp_seurat_data, verbose=FALSE)
    temp_seurat_data <- RunPCA(temp_seurat_data, verbose=FALSE)
    suppressWarnings(
        temp_seurat_data <- RunUMAP(
            temp_seurat_data,
            reduction="pca",
            dims=1:args$ndim,
            verbose=FALSE
        )
    )
    export_dim_plot(
        data=temp_seurat_data,
        reduction="pca",
        plot_title="Split by cell cycle phase PCA of filtered unintegrated/scaled datasets",
        legend_title="Cell cycle phase",
        split_by="Phase",
        group_by="Phase",
        rootname=paste(args$output, "_fltr_pca_spl_by_ph", sep=""),
        palette="Paired",
        pdf=args$pdf
    )
    export_dim_plot(
        data=temp_seurat_data,
        reduction="pca",
        plot_title="Split by level of mitochondrial gene expression PCA of filtered unintegrated/scaled datasets",
        legend_title="Expression level",
        split_by="mito_factor",
        group_by="mito_factor",
        rootname=paste(args$output, "_fltr_pca_spl_by_mito_perc", sep=""),
        palette="Paired",
        pdf=args$pdf
    )
    export_dim_plot(
        data=temp_seurat_data,
        reduction="umap",
        plot_title="Split by identity UMAP projected PCA of filtered unintegrated/scaled datasets",
        legend_title="Identity",
        split_by="new.ident",
        rootname=paste(args$output, "_fltr_umap_spl_by_idnt", sep=""),
        palette="Paired",
        pdf=args$pdf
    )
}


integrate_seurat_data <- function(seurat_data, cell_cycle_data, args) {
    splitted_seurat_data <- SplitObject(seurat_data, split.by="new.ident")
    if (length(splitted_seurat_data) == 1){
        print(paste(
            "Skipping datasets integration as only one identity is present.",
            "Running log-normalization and scalling instead."
        ))
        scaled_norm_seurat_data <- NormalizeData(splitted_seurat_data[[1]], verbose=FALSE)
        tryCatch(
            expr = {
                scaled_norm_seurat_data <- CellCycleScoring(
                    scaled_norm_seurat_data,
                    s.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="s", "gene_id"]),
                    g2m.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="g2/m", "gene_id"]),
                    verbose=FALSE
                )
            },
            error = function(e){
                print("Failed to run cell cycle scoring on log-normalized data")
            }
        )
        scaled_norm_seurat_data <- FindVariableFeatures(
            scaled_norm_seurat_data,
            nfeatures=args$highvarcount,
            verbose=FALSE
        )
        vars_to_regress <- NULL
        if (args$regressmt){ vars_to_regress <- c("mito_percentage") }
        if (args$regresscellcycle){ vars_to_regress <- c("S.Score", "G2M.Score") }
        if (args$regressmt && args$regresscellcycle){ vars_to_regress <- c("mito_percentage", "S.Score", "G2M.Score") }
        scaled_norm_seurat_data <- ScaleData(
            scaled_norm_seurat_data,
            vars.to.regress=vars_to_regress,
            verbose=FALSE
        )
        return (scaled_norm_seurat_data)
    } else {
        for (i in 1:length(splitted_seurat_data)) {
            vars_to_regress <- NULL
            if (args$regressmt){
                vars_to_regress <- c("mito_percentage")
            }
            splitted_seurat_data[[i]] <- SCTransform(
                splitted_seurat_data[[i]],
                assay="RNA",
                new.assay.name="SCT",
                variable.features.n=args$highvarcount,
                vars.to.regress=vars_to_regress,
                verbose=FALSE
            )
            tryCatch(
                expr = {
                    splitted_seurat_data[[i]] <- CellCycleScoring(
                        splitted_seurat_data[[i]],
                        s.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="s", "gene_id"]),
                        g2m.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="g2/m", "gene_id"]),
                        assay="SCT",
                        verbose=FALSE
                    )
                },
                error = function(e){
                    print(paste("Failed to run cell cycle scoring on SCT normalized data for", Idents(splitted_seurat_data[[i]])[1]))
                }
            )
            if (args$regresscellcycle){
                if (!is.null(vars_to_regress)) {
                    vars_to_regress <- append(vars_to_regress, c("S.Score", "G2M.Score"))
                } else {
                    vars_to_regress <- c("S.Score", "G2M.Score")
                }
                tryCatch(
                    expr = {
                        splitted_seurat_data[[i]] <- SCTransform(
                            splitted_seurat_data[[i]],
                            assay="RNA",
                            new.assay.name="SCT",
                            variable.features.n=args$highvarcount,
                            vars.to.regress=vars_to_regress,
                            verbose=FALSE
                        )
                    },
                    error = function(e){
                        print(paste("Failed to regress cell cycle phase as a confounding source of variation for", Idents(splitted_seurat_data[[i]])[1]))
                    }
                )
            }
        }
        set_threads(1)  # need to use one thread, otherwise gets stuck
        integration_features <- SelectIntegrationFeatures(splitted_seurat_data, nfeatures=args$highvarcount)
        splitted_seurat_data <- PrepSCTIntegration(
            splitted_seurat_data, 
            anchor.features=integration_features,
            verbose=FALSE
        )
        integration_anchors <- FindIntegrationAnchors(
            splitted_seurat_data,
            normalization.method="SCT",
            anchor.features=integration_features,
            verbose=FALSE
        )
        integrated_seurat_data <- IntegrateData(
            integration_anchors, 
            new.assay.name="integrated",
            normalization.method="SCT",
            k.weight=min(min(table(Idents(seurat_data))), 100),  # 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
            verbose=FALSE
        )
        DefaultAssay(integrated_seurat_data) <- "integrated"
        set_threads(args$threads)
        return (integrated_seurat_data)
    }
}


get_all_putative_markers <- function(seurat_data, args, assay="RNA", min_diff_pct=-Inf){
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    cluster_prefix <- "RNA"
    if ("integrated" %in% names(seurat_data@assays)) {
        cluster_prefix <- "integrated"
    }
    all_putative_markers <- NULL
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        Idents(seurat_data) <- paste(paste(cluster_prefix, "snn_res", sep="_"), resolution, sep=".")
        markers <- FindAllMarkers(
            seurat_data,
            logfc.threshold=args$logfc,
            min.pct=args$minpct,
            only.pos=args$onlypos,
            test.use=args$testuse,
            min.diff.pct=min_diff_pct,
            verbose=FALSE
        ) %>% relocate(cluster, gene, .before=1)
        if (nrow(markers) > 0) {
            markers <- markers %>% cbind(resolution=resolution, .)
        } else {
            markers <- markers %>% add_column(resolution=numeric(), .before=1)  # safety measure in case markers was empty
        }
        if (!is.null(all_putative_markers)) {
            all_putative_markers <- rbind(all_putative_markers, markers)
        } else {
            all_putative_markers <- markers
        }
        Idents(seurat_data) <- "new.ident"
    }
    DefaultAssay(seurat_data) <- backup_assay
    return (all_putative_markers)
}


get_all_conserved_markers <- function(seurat_data, args, assay="RNA", min_diff_pct=-Inf){
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "RNA"
    cluster_prefix <- "RNA"
    if ("integrated" %in% names(seurat_data@assays)) {
        cluster_prefix <- "integrated"
    }
    all_conserved_markers <- NULL
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        clustering_by <- paste(paste(cluster_prefix, "snn_res", sep="_"), resolution, sep=".")
        Idents(seurat_data) <- clustering_by
        conserved_markers <- map_dfr(
            sort(unique(seurat_data@meta.data[, clustering_by])),
            get_conserved_markers,
            seurat_data,
            "condition",
            resolution,
            args$onlypos,
            args$logfc,
            args$minpct,
            args$testuse,
            min_diff_pct
        )
        if (!is.null(all_conserved_markers)) {
            all_conserved_markers <- rbind(all_conserved_markers, conserved_markers)
        } else {
            all_conserved_markers <- conserved_markers
        }
        Idents(seurat_data) <- "new.ident"
    }
    DefaultAssay(seurat_data) <- backup_assay
    return (all_conserved_markers)
}


get_conserved_markers <- function(cluster, seurat_data, grouping_var, resolution, only_pos=FALSE, logfc_threshold=0.25, min_pct=0.1, test_use="wilcox", min_diff_pct=-Inf){
    conserved_markers <- FindConservedMarkers(
        seurat_data,
        ident.1=cluster,
        grouping.var=grouping_var,
        min.cells.group=0,                 # to catch situation when after splitting by grouping.var we got to few cells in a group
        only.pos=only_pos,
        logfc.threshold=logfc_threshold,
        min.pct=min_pct,
        test.use=test_use,
        min.diff.pct=min_diff_pct,
        verbose=FALSE
    ) %>% rownames_to_column(var="gene")
    if (nrow(conserved_markers) > 0) {
        conserved_markers <- conserved_markers %>% cbind(resolution=resolution, cluster=cluster, .) 
    } else {
        conserved_markers <- conserved_markers %>% add_column(resolution=numeric(), cluster=factor(), .before=1)
    }
    return (conserved_markers)
}


assign_cell_types <- function (seurat_data, classifier, assay, matrix_slot, args){
    cluster_prefix <- "RNA"
    if ("integrated" %in% names(seurat_data@assays)) {
        cluster_prefix <- "integrated"
    }
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        tryCatch(
            expr = {
                monocle_data <- get_monocle_data(
                    seurat_data,
                    features=rownames(seurat_data),                                   # get all features, may not work well for other than "RNA" assays
                    cluster_field=paste(paste(cluster_prefix, "snn_res", sep="_"), resolution, sep="."),
                    assay=assay,
                    matrix_slot=matrix_slot
                )
                monocle_data <- classify_cells(
                    monocle_data,
                    classifier,
                    db=get(as.character(SPECIES_DATA[args$species])),
                    cluster_extend=TRUE,
                    cds_gene_id_type="SYMBOL"
                )
                seurat_data@meta.data[paste("cluster_ext_type_res", resolution, sep=".")] <- monocle_data$cluster_ext_type
            },
            error = function(e){
                print(paste("Failed to make cell type prediction for clusters with resolution ", resolution, sep=""))
            }
        )
    }
    return (seurat_data)
}


get_monocle_data <- function (seurat_data, features, cluster_field, assay="RNA", matrix_slot="counts") {
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    gene_metadata <- data.frame(
        "gene_short_name"=features,
        row.names=features,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    cell_metadata <- data.frame(
        seurat_data@meta.data,
        row.names=colnames(seurat_data),
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    cell_metadata["garnett_cluster"] <- cell_metadata[cluster_field]
    expression_matrix <- GetAssayData(seurat_data, slot=matrix_slot)[features, ]
    monocle_data <- new_cell_data_set(
        expression_matrix,
        cell_metadata=cell_metadata,
        gene_metadata=gene_metadata
    )
    DefaultAssay(seurat_data) <- backup_assay
    return (monocle_data)
}


export_geom_bar_plot <- function(data, rootname, x_axis, color_by, x_label, y_label, legend_title, plot_title, palette="Paired", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                geom_bar(colour="black") +
                geom_text(stat="count", aes(label=..count..), vjust=-1) +
                xlab(x_label) +
                ylab(y_label) +
                guides(fill=guide_legend(legend_title), x=guide_axis(angle=45)) +
                ggtitle(plot_title) +
                scale_fill_brewer(palette=palette)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export geom bar plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export geom bar plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_geom_density_plot <- function(data, rootname, x_axis, color_by, facet_by, x_left_intercept, x_label, y_label, legend_title, plot_title, x_right_intercept=NULL,  scale_x_log10=FALSE, scale_y_log10=FALSE, zoom_on_intercept=FALSE, alpha=0.9, show_ranked=FALSE, ranked_x_label="Ranked cells", palette="Paired", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            palette_colors <- RColorBrewer::brewer.pal(99, palette)  # use 99 to get all available colors. As we added LETTERS, the order will be correct
            intercept_data <- data %>%
                              dplyr::select(all_of(color_by), all_of(facet_by)) %>%
                              distinct() %>%
                              arrange(all_of(color_by)) %>%
                              add_column(color=palette_colors[1:nrow(.)],x_left=x_left_intercept)

            plot <- ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                    geom_density(alpha=alpha) +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(fill=guide_legend(legend_title)) +
                    ggtitle(plot_title) +
                    facet_wrap(as.formula(paste("~", facet_by))) +
                    scale_fill_brewer(palette=palette) +
                    geom_vline(intercept_data, mapping=aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=x_left, y=Inf, label=x_left),
                        color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                        show.legend=FALSE
                    )

            if (!is.null(x_right_intercept)){
                intercept_data <- intercept_data %>% add_column(x_right=x_right_intercept)
                plot <- plot +
                        geom_vline(intercept_data, mapping=aes(xintercept=x_right), color=intercept_data$color, alpha=0.7) +
                        geom_label_repel(
                            intercept_data, mapping=aes(x=x_right, y=Inf, label=x_right),
                            color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                            show.legend=FALSE
                        )
            }

            if (scale_x_log10){ plot <- plot + scale_x_log10() }
            if (scale_y_log10){ plot <- plot + scale_y_log10() }

            if (zoom_on_intercept) {
                zoomed_plot <- ggplot(data, aes_string(x=x_axis, color=color_by)) +
                               geom_density(show.legend=FALSE, size=2) +
                               xlab(x_label) +
                               ylab(y_label) +
                               scale_color_brewer(palette=palette) +
                               geom_vline(intercept_data, mapping=aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                               geom_label_repel(
                                   intercept_data, mapping=aes(x=x_left, y=Inf, label=x_left),
                                   color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                   show.legend=FALSE
                               )
                if (!is.null(x_right_intercept)){
                    zoomed_plot <- zoomed_plot +
                                   coord_cartesian(xlim=c(min(intercept_data$x_left), max(intercept_data$x_right))) +
                                   geom_vline(intercept_data, mapping=aes(xintercept=x_right), color=intercept_data$color, alpha=0.7) +
                                   geom_label_repel(
                                       intercept_data, mapping=aes(x=x_right, y=Inf, label=x_right),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                       show.legend=FALSE
                                   )
                } else {
                    zoomed_plot <- zoomed_plot + coord_cartesian(xlim=c(NA, max(intercept_data$x_left)))
                }
                if (scale_x_log10){ zoomed_plot <- zoomed_plot + scale_x_log10() }
                if (scale_y_log10){ zoomed_plot <- zoomed_plot + scale_y_log10() }
                plot <- plot / zoomed_plot
            }

            if (show_ranked) {
                ranked_plot <- data %>%
                               arrange(get(color_by), get(x_axis)) %>%
                               ggplot(aes_string(x=seq_along(data[[x_axis]]), y=x_axis, color=color_by)) +
                               geom_point(show.legend=FALSE, size=0.5) +
                               xlab(ranked_x_label) +
                               ylab(x_label) +
                               scale_y_log10() +
                               scale_color_brewer(palette=palette) +
                               geom_hline(intercept_data, mapping=aes(yintercept=x_left), color=intercept_data$color, alpha=0.7) +
                               geom_label_repel(
                                   intercept_data, mapping=aes(x=Inf, y=x_left, label=x_left),
                                   color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="x", size=3,
                                   show.legend=FALSE
                               )
                if (!is.null(x_right_intercept)){
                    ranked_plot <- ranked_plot +
                                   geom_hline(intercept_data, mapping=aes(yintercept=x_right), color=intercept_data$color, alpha=0.7) +
                                   geom_label_repel(
                                       intercept_data, mapping=aes(x=-Inf, y=x_right, label=x_right),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="x", size=3,
                                       show.legend=FALSE
                                   )
                }
                plot <- plot / ranked_plot
            }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export geom density plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export geom density plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_geom_point_plot <- function(data, rootname, x_axis, y_axis, facet_by, x_left_intercept, y_low_intercept, y_high_intercept, color_by, colors, color_limits, color_break, x_label, y_label, legend_title, plot_title, scale_x_log10=FALSE, scale_y_log10=FALSE, alpha=0.2, palette="Paired", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            palette_colors <- RColorBrewer::brewer.pal(99, palette)  # use 99 to get all available colors. As we added LETTERS, the order will be correct
            intercept_data <- data %>%
                              dplyr::select(all_of(facet_by)) %>%
                              distinct() %>%
                              arrange(all_of(facet_by)) %>%
                              add_column(color=palette_colors[1:nrow(.)], x_left=x_left_intercept, y_low=y_low_intercept, y_high=y_high_intercept)

            plot <- ggplot(data, aes_string(x=x_axis, y=y_axis, color=color_by)) +
                    geom_point(alpha=alpha) +
                    scale_colour_gradientn(
                        colours=c(colors[1], colors),
                        values=rescale(c(color_limits[1], color_break-0.01*color_break, color_break, color_limits[2])),
                        breaks=c(color_break),
                        limits=color_limits
                    ) +
                    stat_smooth(method=lm) +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(color=guide_colourbar(legend_title)) +
                    ggtitle(plot_title) +
                    facet_wrap(as.formula(paste("~", facet_by))) +
                    geom_vline(intercept_data, mapping=aes(xintercept=x_left), color=intercept_data$color, alpha=0.5) +
                    geom_hline(intercept_data, mapping=aes(yintercept=y_low), color=intercept_data$color, alpha=0.5) +
                    geom_hline(intercept_data, mapping=aes(yintercept=y_high), color=intercept_data$color, alpha=0.5) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=x_left, y=Inf, label=x_left),
                        color="black", fill=intercept_data$color, alpha=0.5, direction="y", size=3,
                        show.legend=FALSE
                    ) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=Inf, y=y_low, label=y_low),
                        color="black", fill=intercept_data$color, alpha=0.5, direction="x", size=3,
                        show.legend=FALSE
                    ) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=Inf, y=y_high, label=y_high),
                        color="black", fill=intercept_data$color, alpha=0.5, direction="x", size=3,
                        show.legend=FALSE
                    )

            if (scale_x_log10){ plot <- plot + scale_x_log10() }
            if (scale_y_log10){ plot <- plot + scale_y_log10() }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export geom point plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export geom point plot to ", rootname, ".(png/pdf)", sep=""))
        }
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


export_expr_heatmap <- function(data, features, rootname, plot_title, matrix_slot="data", palette=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- DoHeatmap(
                        data,
                        features=features,
                        slot=matrix_slot
                    ) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    theme(
                        axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks=element_blank()
                    ) +
                    NoLegend()
          
            if (!is.null(palette)){ plot <- plot + scale_fill_gradientn(colors=palette) }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export expression heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export expression heatmap to ", rootname, ".(png/pdf)", sep=""))
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


export_variable_feature_plot <- function(data, rootname, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(VariableFeaturePlot(data)))
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(VariableFeaturePlot(data)))
                dev.off()
            }
            print(paste("Export variable feature plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export variable feature plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_dim_plot <- function(data, rootname, reduction, plot_title, legend_title, split_by=NULL, group_by=NULL, label=FALSE, palette=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
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


export_pca_heatmap <- function(data, rootname, plot_title, x_label, y_label, dims=NULL, cells=500, nfeatures=30, ncol=NULL, combine_guides=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots <- DimHeatmap(
                        data,
                        dims=dims,
                        nfeatures=nfeatures,
                        reduction="pca",
                        cells=cells,
                        balanced=TRUE,
                        fast=FALSE,
                        combine=FALSE
                    )
            plots <- lapply(seq_along(plots), function(i){
                plots[[i]] +
                ggtitle(paste("PC", i, sep=" ")) +
                theme_gray() +
                xlab(x_label) +
                ylab(y_label) +
                theme(axis.text.x=element_blank(), axis.ticks = element_blank())
            })
            combined_plots <- wrap_plots(plots, guides=combine_guides, ncol=ncol) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export PCA heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export PCA heatmap to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_pca_loadings_plot <- function(data, rootname, plot_title, x_label, y_label, dims=NULL, nfeatures=30, ncol=NULL, combine_guides=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots <- VizDimLoadings(
                        data,
                        dims=dims,
                        nfeatures=nfeatures,
                        reduction="pca",
                        combine=FALSE
                    )
            plots <- lapply(seq_along(plots), function(i){
                plots[[i]] +
                ggtitle(paste("PC", i, sep=" ")) +
                theme_gray() +
                xlab(x_label) +
                ylab(y_label)
            })
            combined_plots <- wrap_plots(plots, guides=combine_guides, ncol=ncol) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export PCA loadings plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export PCA loadings plot to ", rootname, ".(png/pdf)", sep=""))
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


export_all_clustering_plots <- function(seurat_data, suffix, args) {
    cluster_prefix <- "RNA"
    if ("integrated" %in% names(seurat_data@assays)) {
        cluster_prefix <- "integrated"
    }
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        export_dim_plot(
            data=seurat_data,
            reduction="umap",
            plot_title="Clustered UMAP projected PCA of filtered integrated/scaled datasets",
            legend_title="Cluster",
            group_by=paste(paste(cluster_prefix, "snn_res", sep="_"), current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="umap",
            plot_title="Split by condition clustered UMAP projected PCA of filtered integrated/scaled datasets",
            legend_title="Cluster",
            split_by="condition",
            group_by=paste(paste(cluster_prefix, "snn_res", sep="_"), current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "umap_spl_by_cond_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="umap",
            plot_title="Grouped by predicted cell types UMAP projected PCA of filtered integrated/scaled datasets",
            legend_title="Cell type prediction",
            group_by=paste("cluster_ext_type_res", current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "umap_ctype_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="umap",
            plot_title="Split by cell cycle phase clustered UMAP projected PCA of filtered integrated/scaled datasets",
            legend_title="Cluster",
            split_by="Phase",
            group_by=paste(paste(cluster_prefix, "snn_res", sep="_"), current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "umap_spl_by_ph_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        Idents(seurat_data) <- paste(paste(cluster_prefix, "snn_res", sep="_"), current_resolution, sep=".")
        export_feature_plot(
            data=seurat_data,
            features=c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "mito_percentage", "log10_gene_per_log10_umi"),
            labels=c("UMIs", "Genes", "S score", "G to M score", "Mitochondrial %", "Novelty score"),
            reduction="umap",
            plot_title="QC metrics for clustered UMAP projected PCA of filtered integrated/scaled datasets",
            label=TRUE,
            alpha=0.4,
            rootname=paste(args$output, suffix, "qc_mtrcs_res", current_resolution, sep="_"),
            combine_guides="keep",
            pdf=args$pdf
        )
        Idents(seurat_data) <- "new.ident"
    }
}


export_all_expression_plots <- function(seurat_data, suffix, args, assay="RNA") {
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    cluster_prefix <- "RNA"
    if ("integrated" %in% names(seurat_data@assays)) {
        cluster_prefix <- "integrated"
    }
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        Idents(seurat_data) <- paste(paste(cluster_prefix, "snn_res", sep="_"), current_resolution, sep=".")
        export_dot_plot(
            data=seurat_data,
            features=args$features,
            plot_title="Scaled average log normalized gene expression per cluster of filtered integrated/scaled datasets",
            x_label="Genes",
            y_label="Clusters",
            cluster_idents=FALSE,
            rootname=paste(args$output, suffix, "avg_per_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_feature_plot(
            data=seurat_data,
            features=args$features,
            labels=args$features,
            reduction="umap",
            plot_title="Log normalized gene expression per cell of clustered filtered integrated/scaled datasets",
            label=TRUE,
            order=TRUE,
            max_cutoff="q99",  # to prevent cells with overexpressed gene from distoring the color bar
            rootname=paste(args$output, suffix, "per_clst_cell_res", current_resolution, sep="_"),
            combine_guides="keep",
            pdf=args$pdf
        )
        export_expr_heatmap(
            data=seurat_data,
            features=args$features,
            plot_title="Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets",
            matrix_slot="data",  # LogNormalized version of the raw counts
            palette=c("black", "orange"),
            rootname=paste(args$output, suffix, "clst_heatmap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_vln_plot(
            data=seurat_data,
            features=args$features,
            labels=args$features,
            plot_title="Log normalized gene expression densities per cluster of filtered integrated/scaled datasets",
            legend_title="Cluster",
            log=TRUE,
            pt_size=0,
            combine_guides="collect",
            rootname=paste(args$output, suffix, "dnst_per_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (paste("cluster_ext_type_res", current_resolution, sep=".") %in% colnames(seurat_data@meta.data)){
            Idents(seurat_data) <- paste("cluster_ext_type_res", current_resolution, sep=".")
            export_dot_plot(
                data=seurat_data,
                features=args$features,
                plot_title="Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets",
                x_label="Genes",
                y_label="Cell types",
                cluster_idents=FALSE,  # no need to cluster cell types together
                rootname=paste(args$output, suffix, "avg_per_ctype_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            export_feature_plot(
                data=seurat_data,
                features=args$features,
                labels=args$features,
                reduction="umap",
                plot_title="Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types",
                label=TRUE,
                order=TRUE,
                max_cutoff="q99",  # to prevent cells with overexpressed gene from distorting the color bar
                rootname=paste(args$output, suffix, "per_ctype_cell_res", current_resolution, sep="_"),
                combine_guides="keep",
                pdf=args$pdf
            )
            export_expr_heatmap(
                data=seurat_data,
                features=args$features,
                plot_title="Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets with predicted cell types",
                matrix_slot="data",  # LogNormalized version of the raw counts
                palette=c("black", "orange"),
                rootname=paste(args$output, suffix, "ctype_heatmap_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            export_vln_plot(
                data=seurat_data,
                features=args$features,
                labels=args$features,
                plot_title="Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets",
                legend_title="Cell type",
                log=TRUE,
                pt_size=0,
                combine_guides="collect",
                rootname=paste(args$output, suffix, "dnst_per_ctype_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
        Idents(seurat_data) <- "new.ident"
    }
    DefaultAssay(seurat_data) <- backup_assay
}


export_all_dimensionality_plots <- function(seurat_data, suffix, args) {
    export_elbow_plot(
        data=seurat_data,
        ndims=50,
        rootname=paste(args$output, suffix, "elbow", sep="_"),
        plot_title="Elbow plot from PCA of filtered integrated/scaled datasets",
        pdf=args$pdf
    )
    export_dim_plot(
        data=seurat_data,
        reduction="pca",
        plot_title="PCA of filtered integrated/scaled datasets",
        legend_title="Identity",
        rootname=paste(args$output, suffix, "pca", sep="_"),
        palette="Paired",
        pdf=args$pdf
    )
    export_pca_heatmap(
        data=seurat_data,
        dims=1:50,
        rootname=paste(args$output, suffix, "pca_heatmap", sep="_"),
        plot_title="Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated/scaled datasets",
        x_label="Cells",
        y_label="Genes",
        cells=500,
        nfeatures=30,
        ncol=3,
        height=6500,
        combine_guides="collect",
        pdf=args$pdf
    )
    export_pca_loadings_plot(
        data=seurat_data,
        dims=1:50,
        rootname=paste(args$output, suffix, "pca_loadings", sep="_"),
        plot_title="PC scores of the most variant genes from PCA of filtered integrated/scaled datasets",
        x_label="PC scores",
        y_label="Genes",
        nfeatures=30,
        ncol=3,
        height=6500,
        combine_guides="collect",
        pdf=args$pdf
    )
    export_dim_plot(
        data=seurat_data,
        reduction="umap",
        plot_title="Split by identity UMAP projected PCA of filtered integrated/scaled datasets",
        legend_title="Identity",
        split_by="new.ident",
        rootname=paste(args$output, suffix, "umap_spl_by_idnt", sep="_"),
        palette="Paired",
        pdf=args$pdf
    )
}


export_all_qc_plots <- function(seurat_data, suffix, args){
    export_geom_bar_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "cell_count", sep="_"),
        x_axis="new.ident",
        color_by="new.ident",
        x_label="Identity",
        y_label="Cells",
        legend_title="Identity",
        plot_title=paste("Number of cells per dataset (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "umi_dnst_spl_by_cond", sep="_"),
        x_axis="nCount_RNA",
        color_by="new.ident",
        facet_by="condition",
        x_left_intercept=args$minumi,
        x_label="UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Split by condition UMI density per cell (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_dnst_spl_by_cond", sep="_"),
        x_axis="nFeature_RNA",
        color_by="new.ident",
        facet_by="condition",
        x_left_intercept=args$minfeatures,
        x_right_intercept=args$maxfeatures,
        x_label="Genes per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Split by condition gene density per cell (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        pdf=args$pdf
    )
    export_geom_point_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_umi_corr_spl_by_ident", sep="_"),
        x_axis="nCount_RNA",
        y_axis="nFeature_RNA",
        facet_by="new.ident",
        x_left_intercept=args$minumi,
        y_low_intercept=args$minfeatures,
        y_high_intercept=args$maxfeatures,
        color_by="mito_percentage",
        colors=c("lightslateblue", "red", "green"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        x_label="UMIs per cell",
        y_label="Genes per cell",
        legend_title="Mitochondrial %",
        plot_title=paste("Split by identity genes vs UMIs per cell correlation (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "mito_perc_dnst_spl_by_cond", sep="_"),
        x_axis="mito_percentage",
        color_by="new.ident",
        facet_by="condition",
        x_left_intercept=args$maxmt,
        x_label="Percentage of transcripts mapped to mitochondrial genes per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Split by condition density of transcripts mapped to mitochondrial genes per cell (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "nvlt_score_dnst_spl_by_cond", sep="_"),
        x_axis="log10_gene_per_log10_umi",
        color_by="new.ident",
        facet_by="condition",
        x_left_intercept=args$minnovelty,
        x_label="log10 Genes / log10 UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Split by condition novelty score density per cell (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_vln_plot(
        data=seurat_data,
        features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi"),
        labels=c("UMIs", "Genes", "Mitochondrial %", "Novelty score"),
        rootname=paste(args$output, suffix, "qc_mtrcs", sep="_"),
        plot_title=paste("QC metrics densities per cell (", suffix, ")", sep=""),
        legend_title="Identity",
        hide_x_text=TRUE,
        pt_size=0,
        palette="Paired",
        combine_guides="collect",
        pdf=args$pdf
    )
    export_vln_plot(
        data=seurat_data,
        features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi"),
        labels=c("UMIs", "Genes", "Mitochondrial %", "Novelty score"),
        rootname=paste(args$output, suffix, "qc_mtrcs_gr_by_cond", sep="_"),
        plot_title=paste("Grouped by condition QC metrics densities per cell (", suffix, ")", sep=""),
        legend_title="Condition",
        hide_x_text=TRUE,
        group_by="condition",
        pt_size=0,
        palette="Paired",
        combine_guides="collect",
        pdf=args$pdf
    )
}


export_cellbrowser_data <- function(seurat_data, assay, matrix_slot, resolution, features, rootname){
    tryCatch(
        expr = {
            backup_assay <- DefaultAssay(seurat_data)
            DefaultAssay(seurat_data) <- assay
            meta_fields <- c("nCount_RNA", "nFeature_RNA", "log10_gene_per_log10_umi", "mito_percentage", "Phase", "S.Score", "G2M.Score")
            meta_fields_names <- c("UMIs", "Genes", "Novelty score", "Mitochondrial %", "Cell cycle phase", "S score", "G to M score")
            cluster_field <- paste("Clustering (", resolution, ")", sep="")[1]
            cluster_prefix <- "RNA"
            if ("integrated" %in% names(seurat_data@assays)) {
                cluster_prefix <- "integrated"
                cluster_field <- "Identity"
                meta_fields <- c(c("new.ident", "condition"), meta_fields)
                meta_fields_names <- c(c("Identity", "Condition"), meta_fields_names)
            }
            meta_fields <- c(
                meta_fields,
                paste(paste(cluster_prefix, "snn_res", sep="_"), resolution, sep="."),
                paste("cluster_ext_type_res", resolution, sep=".")
            )
            meta_fields_names <- c(
                meta_fields_names,
                paste("Clustering (", resolution, ")", sep=""),
                paste("Cell Type Prediction (", resolution, ")", sep="")
            )
            ExportToCellbrowser(
                seurat_data,
                matrix.slot=matrix_slot,
                dir=rootname,
                cluster.field=cluster_field,
                features=features,
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


export_data <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE){
    tryCatch(
        expr = {
            write.table(
                data,
                file=location,
                sep=get_file_type(location),
                row.names=row_names,
                col.names=col_names,
                quote=quote
            )
            print(paste("Export data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data to", location, sep=" "))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description='Runs Seurat for comparative scRNA-seq analysis of across experimental conditions')

    # Import aggregated MEX data from Cell Ranger Aggregate or from multiple Cell Ranger Count runs
    parser$add_argument(
        "--mex",
        help=paste(
            "Path to the folder with not normalized aggregated feature-barcode matrix",
            "from Cell Ranger Aggregate in MEX format. If multiple locations provided",
            "data is assumed to be not aggregated (outputs from multiple Cell Ranger",
            "Count runs) and will be merged."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--identity",
        help=paste(
            "Path to the metadata TSV/CSV file to set the datasets identities.",
            "If --mex points to the Cell Ranger Aggregate outputs, the aggregation.csv",
            "file can be used as well. If multiple locations were provided through --mex,",
            "the file should include at least one column - 'library_id', and be sorted",
            "based on the the order of locations provided in --mex."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--condition",
        help=paste(
            "Path to the TSV/CSV file to define datasets grouping. First column -",
            "'library_id' with the values provided in the correspondent column of the",
            "--identity file, second column 'condition'. Default: each dataset is",
            "assigned to a separate group."
        ),
        type="character"
    )
    parser$add_argument(
        "--classifier",
        help=paste(
            "Path to the Garnett classifier RDS file for cell type prediction.",
            "Default: skip cell type prediction."
        ),
        type="character"
    )
    parser$add_argument(
        "--cellcycle",
        help=paste(
            "Path to the TSV/CSV file with cell cycle data. First column - 'phase',",
            "second column 'gene_id'. Default: skip cell cycle score assignment."
        ),
        type="character"
    )
    parser$add_argument(
        "--barcodes",
        help=paste(
            "Path to the headerless TSV/CSV file with the list of barcodes to select",
            "cells of interest (one barcode per line). Prefilters input feature-barcode",
            "matrix to include only selected cells. Default: use all cells."
        ),
        type="character"
    )

    # Apply QC filters
    parser$add_argument(
        "--mincells",
        help=paste(
            "Include only features detected in at least this many cells. Applied to",
            "aggregated feature-barcode matrix from Cell Ranger Aggregate. Ignored",
            "when --mex points to the locations of multiple Cell Ranger Count runs.",
            "Default: 5"
        ),
        type="integer", default=5
    )
    parser$add_argument(                                                                        # array
        "--minfeatures",
        help=paste(
            "Include cells where at least this many features are detected. If multiple",
            "values provided each of them will be applied to the correspondent dataset",
            "from the --mex input. Default: 250 (applied to all datasets)"
        ),
        type="integer", default=250, nargs="*"
    )
    parser$add_argument(                                                                        # array
        "--maxfeatures",
        help=paste(
            "Include cells with the number of features not bigger than this value. If",
            "multiple values provided each of them will be applied to the correspondent",
            "dataset from the --mex input. Default: 5000 (applied to all datasets)"
        ),
        type="integer", default=5000, nargs="*"
    )
    parser$add_argument(                                                                        # array
        "--minumi",
        help=paste(
            "Include cells where at least this many UMIs (transcripts) are detected. If",
            "multiple values provided each of them will be applied to the correspondent",
            "dataset from the --mex input. Default: 500 (applied to all datasets)"
        ),
        type="integer", default=500, nargs="*"
    )
    parser$add_argument(                                                                        # array
        "--minnovelty",
        help=paste(
            "Include cells with the novelty score not lower than this value, calculated as",
            "log10(genes)/log10(UMIs). If multiple values provided each of them will be",
            "applied to the correspondent dataset from the --mex input. Default: 0.8 (applied",
            "to all datasets)"
        ),
        type="double", default=0.8, nargs="*"
    )
    parser$add_argument(
        "--maxmt",
        help=paste(
            "Include cells with the percentage of transcripts mapped to mitochondrial genes",
            "not bigger than this value. Default: 5"
        ),
        type="double", default=5
    )
    parser$add_argument(
        "--mitopattern",
        help="Regex pattern to identify mitochondrial genes. Default: '^Mt-'",
        type="character", default="^Mt-"
    )

    # Integration, clustering, cell types and marker genes identification parameters
    parser$add_argument(
        "--features",
        help="Features of interest to evaluate expression. Default: None",
        type="character", nargs="*"
    )
    parser$add_argument(
        "--regresscellcycle",
        help="Regress cell cycle as a confounding source of variation. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--regressmt",
        help=paste(
            "Regress mitochondrial genes expression as a confounding source of variation.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--highvarcount",
        help=paste(
            "Number of highly variable features to detect. Used for datasets integration,",
            "scaling, and dimensional reduction. Default: 3000"
        ),
        type="integer", default=3000
    )
    parser$add_argument(
        "--ndim",
        help=paste(
            "Number of principal components to use in UMAP projection and clustering",
            "(from 1 to 50). Use Elbow plot to adjust this parameter. Default: 10"
        ),
        type="integer", default=10
    )
    parser$add_argument(
        "--resolution",
        help="Clustering resolution. Can be set as an array. Default: 0.4 0.6 0.8 1.0 1.4",
        type="double", default=c(0.4, 0.6, 0.8, 1.0, 1.4), nargs="*"
    )
    parser$add_argument(
        "--logfc",
        help=paste(
            "Include only those genes that on average have log fold change difference in",
            "expression between every tested pair of clusters not lower than this value.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--minpct",
        help=paste(
            "Include only those features that are detected in not lower than this fraction",
            "of cells in either of the two tested clusters. Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--onlypos",
        help="Return only positive markers when running gene markers identification. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--testuse",
        help="Statistical test to use for gene markers identification. Default: wilcox",
        type="character", default="wilcox",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
    )
    parser$add_argument(
        "--species",
        help=paste(
            "Select species for gene name conversion when running cell type prediction",
            "with Garnett classifier. Default: do not convert gene names"
        ),
        type="character", choices=names(SPECIES_DATA), default="none"
    )

    # Export results
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--rds",
        help="Save Seurat data to RDS file. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./seurat",
        type="character", default="./seurat"
    )

    # Performance parameters
    parser$add_argument(
        "--threads",
        help="Threads. Default: 1",
        type="integer", default=1
    )

    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}


args <- get_args()

cat("Step 0: Runtime configuration\n")
print("Used parameters")
print(args)
print(paste("Setting parallelizations threads to", args$threads))
set_threads(args$threads)

cat("\n\nStep 1: Loading raw datasets\n")

print(paste("Loading datasets identity data from", args$identity))
cell_identity_data <- load_cell_identity_data(args$identity)
print("Trying to load cell cycle data")
cell_cycle_data <- load_cell_cycle_data(args$cellcycle)
print("Trying to load condition data")
condition_data <- load_condition_data(args$condition, cell_identity_data)
print("Loading feature-barcode matrices")
seurat_data <- load_seurat_data(args$mex, args$mincells, cell_identity_data, condition_data)

print("Validating filtering thresholds")
# if we put it in a function, args will be passed by value
idents_count <- length(unique(as.vector(as.character(Idents(seurat_data)))))
for (key in names(args)){
    if (key %in% c("minfeatures", "maxfeatures", "minumi", "minnovelty")){
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
if (!is.null(args$features)){
    print("Check genes of interest to include only those that are present in the datasets")
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "RNA"
    args$features <- args$features[args$features %in% as.vector(as.character(rownames(seurat_data)))]
    print(args$features)
    DefaultAssay(seurat_data) <- backup_assay
}
print("Trying to load barcodes of interest to prefilter feature-barcode matrices by cells")
barcodes_data <- load_barcodes_data(args$barcodes, seurat_data)
print("Prefiltering feature-barcode matrices by cells of interest")
seurat_data <- apply_cell_filters(seurat_data, barcodes_data)
print("Trying to load Garnett classifier")
classifier <- load_classifier(args$classifier)

cat("\n\nStep 2: Filtering raw datasets\n")

print("Adding QC metrics to not filtered seurat data")
seurat_data <- add_qc_metrics(seurat_data, args$mitopattern)
export_all_qc_plots(seurat_data, "raw", args)                                                                  # <--- raw
print("Applying QC filters")
seurat_data <- apply_qc_filters(seurat_data, cell_identity_data, args)
export_all_qc_plots(seurat_data, "fltr", args)                                                                 # <--- fltr
print("Evaluating effects of cell cycle and mitochodrial gene expression")
explore_unwanted_variation(seurat_data, cell_cycle_data, args)

cat("\n\nStep 3: Integrating filtered datasets\n")

print("Running dataset integration/scaling")
seurat_data <- integrate_seurat_data(seurat_data, cell_cycle_data, args)                                       # sets "integrated" as a default assay for integrated data, and "RNA" for scaled data

cat("\n\nStep 4: Reducing dimensionality of intergrated/scaled datasets\n")

print("Performing PCA reduction of integrated/scaled data. Use all 50 principal components")
seurat_data <- RunPCA(seurat_data, npcs=50, verbose=FALSE)                                                     # runs on "integrated" assay for integrated data, and on "RNA" assay for scaled data
print(paste("Performing UMAP reduction of integrated/scaled data using", args$ndim, "principal components"))
seurat_data <- RunUMAP(seurat_data, reduction="pca", dims=1:args$ndim, verbose=FALSE)                          # runs on "integrated" assay for integrated data, and on "RNA" assay for scaled data
export_all_dimensionality_plots(seurat_data, "ntgr", args)                                                     # <--- ntgr

cat("\n\nStep 5: Clustering and cell type assignment of intergrated/scaled datasets with reduced dimensionality\n")

print(paste("Clustering integrated/scaled data using", args$ndim, "principal components"))
seurat_data <- FindNeighbors(seurat_data, reduction="pca", dims=1:args$ndim, verbose=FALSE)                    # runs on "integrated" assay for integrated data, and on "RNA" assay for scaled data
seurat_data <- FindClusters(seurat_data, resolution=args$resolution)
print("Assigning cell types for all clusters and all resolutions using only highly variable genes")
seurat_data <- assign_cell_types(seurat_data, classifier, "RNA", "counts", args)                               # uses all features from "counts" slot of "RNA" assay
export_all_clustering_plots(seurat_data, "clst", args)                                                         # <--- clst

cat("\n\nStep 6: Running gene expression analysis\n")

if ("integrated" %in% names(seurat_data@assays)) {                                                             # if we run integration, our counts and data slots in RNA assay remain the same, so we need to normalize counts first
    print("Normalizing counts in RNA assay")
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
    DefaultAssay(seurat_data) <- backup_assay
}
if (args$rds){ export_rds(seurat_data, paste(args$output, "_clst_data.rds", sep="")) }
if (!is.null(args$features)){ export_all_expression_plots(seurat_data, "expr", args, assay="RNA") }            # <--- expr
print("Identifying putative gene markers for all clusters and all resolutions")
all_putative_markers <- get_all_putative_markers(seurat_data, args)
export_data(all_putative_markers, paste(args$output, "_clst_pttv_gene_markers.tsv", sep=""))
print("Identifying conserved gene markers for all clusters and all resolutions")
all_conserved_markers <- get_all_conserved_markers(seurat_data, args)
export_data(all_conserved_markers, paste(args$output, "_clst_csrvd_gene_markers.tsv", sep=""))
print("Exporting UCSC Cellbrowser data")
export_cellbrowser_data(
    seurat_data, assay="RNA", matrix_slot="data",
    resolution=args$resolution, features=args$features,
    rootname=paste(args$output, "_cellbrowser", sep="")
)
