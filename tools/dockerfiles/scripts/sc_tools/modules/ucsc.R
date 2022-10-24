import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("data.table", attach=FALSE)
import("reticulate", attach=FALSE)

export(
    "get_matrix",
    "write_matrix",
    "cb_build",
    "export_cellbrowser"
)


get_matrix <- function(object, slot){
    if (slot == "counts") {
        counts <- SeuratObject::GetAssayData(object=object, slot="counts")
    } else if (slot == "scale.data") {
        counts <- SeuratObject::GetAssayData(object=object, slot="scale.data")
    } else if (slot == "data") {
        counts <- SeuratObject::GetAssayData(object=object)
    } else {
        base::print(base::paste("Slot", slot, "not found. Quit."))
        base::quit(save="no", status=1, runLast=FALSE)
    }
}

write_matrix <- function(inMat, outFname, sliceSize=1000){
    fnames <- c()
    mat <- inMat
    geneCount <- base::nrow(mat)
    startIdx <- 1
    while (startIdx < geneCount){
        endIdx <- min(startIdx + sliceSize - 1, geneCount)
        matSlice <- mat[startIdx:endIdx,]
        denseSlice <- base::as.matrix(matSlice)
        dt <- data.table::data.table(denseSlice)
        dt <- base::cbind(gene=base::rownames(matSlice), dt)
        writeHeader <- FALSE
        if (startIdx == 1) { 
            writeHeader <- TRUE
        }
        sliceFname <- base::paste0("temp", startIdx, ".txt")
        data.table::fwrite(dt, sep="\t", file=sliceFname, quote=FALSE, col.names=writeHeader)
        fnames <- base::append(fnames, sliceFname)
        startIdx <- startIdx + sliceSize
    }
    base::system(base::paste("cat", base::paste(fnames, collapse=" "), "| gzip >", outFname, sep=" "))
    base::unlink(fnames)
}

cb_build <- function(object, slot, short_label, rootname, cluster_field, features, meta_fields, meta_fields_names, is_nested, dot_radius, dot_alpha) {
    if (!base::dir.exists(rootname)) {
        base::dir.create(rootname, recursive=TRUE)
    }
    idents <- SeuratObject::Idents(object)
    meta <- object@meta.data
    cell_order <- base::colnames(object)
    counts <- get_matrix(object, slot)
    genes <- base::rownames(x=object)
    dr <- object@reductions

    gzPath <- base::file.path(rootname, "exprMatrix.tsv.gz")
    if ((((base::ncol(counts)/1000)*(base::nrow(counts)/1000))>2000) && methods::is(counts, "sparseMatrix")){
        write_matrix(counts, gzPath);
    } else {
        mat <- base::as.matrix(counts)
        df <- base::as.data.frame(mat, check.names=FALSE)
        df <- base::data.frame(gene=genes, df, check.names=FALSE)
        z <- base::gzfile(gzPath, "w")
        utils::write.table(x=df, sep="\t", file=z, quote=FALSE, row.names=FALSE)
        base::close(con=z)
    }
    embeddings = names(dr)
    embeddings_conf <- c()
    for (embedding in embeddings) {
        emb <- dr[[embedding]]
        df <- emb@cell.embeddings
        if (base::ncol(df) > 2){
            df <- df[, 1:2]
        }
        base::colnames(df) <- c("x", "y")
        df <- base::data.frame(cellId=base::rownames(df), df, check.names=FALSE)
        utils::write.table(
            df[cell_order,],
            sep="\t",
            file=base::file.path(rootname, base::sprintf("%s.coords.tsv", embedding)),
            quote=FALSE,
            row.names=FALSE
        )
        embeddings_conf <- c(
            embeddings_conf,
            base::sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}', embedding)
        )
    }
    df <- base::data.frame(row.names=cell_order, check.names=FALSE)
    enum_fields <- c()
    for (i in 1:length(meta_fields)){
        field <- meta_fields[i]
        name <- meta_fields_names[i]
        if (field %in% base::colnames(meta)){
            df[[name]] <- meta[[field]]
            if (!is.numeric(df[[name]])) {
                enum_fields <- c(enum_fields, name)
            }
        }
    }
    df <- base::data.frame(Cell=base::rownames(df), df, check.names=FALSE)
    utils::write.table(
        base::as.matrix(df[cell_order,]),
        sep="\t",
        file=base::file.path(rootname, "meta.tsv"),
        quote=FALSE,
        row.names=FALSE
    )

    local_config <- '
name="%s"
shortLabel="%s"
geneLabel="Feature"
exprMatrix="exprMatrix.tsv.gz"
meta="meta.tsv"
radius=%i
alpha=%f
geneIdType="auto"
clusterField="%s"
labelField="%s"
enumFields=%s
showLabels=False
coords=%s'

    local_config <- base::sprintf(
        local_config,
        short_label,
        short_label,
        dot_radius,
        dot_alpha,
        cluster_field,
        cluster_field,
        base::paste0("[", base::paste(base::paste0('"', enum_fields, '"'), collapse=", "), "]"),
        base::paste0("[", base::paste(embeddings_conf, collapse = ",\n"), "]")
    )

    if (!is.null(features)){
        utils::write.table(
            features,
            sep="\t",
            file=base::file.path(rootname, "quickGenes.tsv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE
        )
        local_config <- base::paste(local_config, 'quickGenesFile="quickGenes.tsv"', sep="\n")
    }

    html_data_dir <- base::file.path(rootname, "html_data")
    if (is_nested){
        local_config <- base::paste(local_config, 'dataRoot="../"', sep="\n")
        desc <- 'title="%s"\nabstract=""\nmethods=""\nbiorxiv_url=""\ncustom={}'
        desc <- base::sprintf(desc, short_label)
        base::cat(
            desc,
            file=base::file.path(rootname, "desc.conf")
        )
        base::cat(
            'shortLabel="Multiple datasets"',
            file=base::file.path(rootname, "../cellbrowser.conf")
        )
        html_data_dir <- base::file.path(rootname, "../html_data")
    }

    base::cat(
        local_config,
        file=base::file.path(rootname, "cellbrowser.conf")
    )
    cb <- reticulate::import(module = "cellbrowser")
    cb$cellbrowser$build(rootname, html_data_dir)
}

export_cellbrowser <- function(seurat_data, assay, slot, rootname, dot_radius=3, dot_alpha=0.5, short_label="cellbrowser", is_nested=FALSE, features=NULL, meta_fields=NULL, meta_fields_names=NULL){
    base::tryCatch(
        expr = {
            backup_assay <- SeuratObject::DefaultAssay(seurat_data)
            SeuratObject::DefaultAssay(seurat_data) <- assay                 # safety measure
            SeuratObject::Idents(seurat_data) <- "new.ident"                 # safety measure

            if (is.null(meta_fields) || is.null(meta_fields_names)){
                meta_fields <-       c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "S.Score", "G2M.Score",    "Phase", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment",       "nucleosome_signal", "frip", "blacklist_fraction")
                meta_fields_names <- c("RNA UMIs",   "Genes",        "Mitochondrial %", "Novelty score",            "S score", "G to M score", "Phase", "ATAC UMIs",   "Peaks",         "TSS enrichment score", "Nucleosome signal", "FRiP", "Bl. regions")
            }

            for (i in 1:length(meta_fields)){
                current_field <- meta_fields[i]
                if (current_field %in% base::colnames(seurat_data@meta.data) && length(base::unique(base::as.vector(as.character(seurat_data@meta.data[, current_field])))) > 1){
                    default_cluster_field <- meta_fields_names[i]                 # first field from the meta.data that is not unique for all cells
                    break
                }
            }

            if (length(base::unique(base::as.vector(as.character(seurat_data@meta.data$new.ident)))) > 1){
                meta_fields <- base::append(meta_fields, "new.ident", 0)
                meta_fields_names <- base::append(meta_fields_names, "Identity", 0)
            }
            if (all(base::as.vector(as.character(seurat_data@meta.data$new.ident)) != base::as.vector(as.character(seurat_data@meta.data$condition)))){
                meta_fields <- base::append(meta_fields, "condition", 1)
                meta_fields_names <- base::append(meta_fields_names, "Condition", 1)
            }

            clustering_fields <- base::grep("_res\\.", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            clustering_fields_names <- base::as.vector(
                base::unlist(                                                                              # when nothing found it returns named list that should be unlisted
                    base::sapply(
                        clustering_fields,
                        function(line) {
                            split_line <- base::unlist(base::strsplit(line, split="_res\\."))
                            base::paste("Clustering (", split_line[1], " ", split_line[2], ")", sep="")
                        }
                    )
                )
            )
            meta_fields <- base::append(meta_fields, clustering_fields)
            meta_fields_names <- base::append(meta_fields_names, clustering_fields_names)

            custom_fields <- base::grep("^custom_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            custom_fields_names <- base::gsub("custom_", "Custom ", custom_fields)
            meta_fields <- base::append(meta_fields, custom_fields)
            meta_fields_names <- base::append(meta_fields_names, custom_fields_names)

            quartile_fields <- base::grep("^quartile_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            quartile_fields_names <- base::gsub("quartile_", "Quartile ", quartile_fields)
            meta_fields <- base::append(meta_fields, quartile_fields)
            meta_fields_names <- base::append(meta_fields_names, quartile_fields_names)

            cb_build(
                seurat_data,
                slot=slot,
                short_label=short_label,
                rootname=rootname,
                cluster_field=default_cluster_field,
                features=features,
                meta_fields=meta_fields,
                meta_fields_names=meta_fields_names,
                is_nested=is_nested,
                dot_radius=dot_radius,
                dot_alpha=dot_alpha
            )
            base::print(base::paste("Exporting UCSC Cellbrowser data to", rootname, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export UCSC Cellbrowser data with error -", e))
        },
        finally = {
            SeuratObject::DefaultAssay(seurat_data) <- backup_assay
        }
    )
}