import("dplyr", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("loupeR", attach=FALSE)
import("data.table", attach=FALSE)
import("reticulate", attach=FALSE)
import("tidyselect", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))

export(
    "export_cellbrowser",
    "export_loupe"
)

DEFAULT_META_FIELDS <- c(
    "nCount_RNA",
    "nFeature_RNA",
    "mito_percentage",
    "log10_gene_per_log10_umi",
    "S.Score",
    "G2M.Score",
    "Phase",
    "rna_doublets",
    "atac_doublets",
    "nCount_ATAC",
    "nFeature_ATAC",
    "TSS.enrichment",
    "nucleosome_signal",
    "frip",
    "blacklist_fraction",

    "donor",

    "clonotype_TRA",
    "clonotype_TRB",
    "clonotype_IGH",
    "clonotype_IGL",
    "clonotype_both",

    "clonalFrequency_TRA",
    "clonalFrequency_TRB",
    "clonalFrequency_IGH",
    "clonalFrequency_IGL",
    "clonalFrequency_both",

    "quartile_nCount_RNA",
    "quartile_nFeature_RNA",
    "quartile_mito_percentage",
    "quartile_nCount_ATAC",
    "quartile_nFeature_ATAC",
    "quartile_TSS.enrichment",
    "quartile_nucleosome_signal",
    "quartile_frip",
    "quartile_blacklist_fraction",

    "prediction_cell_type",
    "prediction_mapping_score",
    "prediction_confidence_score"
)

DEFAULT_META_FIELDS_NAMES <- c(
    "RNA reads",
    "Genes",
    "Mitochondrial %",
    "Novelty score",
    "S score",
    "G2M score",
    "Phase",
    "RNA doublets",
    "ATAC doublets",
    "ATAC fragments in peaks",
    "Peaks",
    "TSS enrichment score",
    "Nucleosome signal",
    "FRiP",
    "Bl. regions",

    "Clonotype donor",

    "Clonotype TRA",
    "Clonotype TRB",
    "Clonotype IGH",
    "Clonotype IGL",
    "Clonotype Both",

    "Clonotype frequency TRA",
    "Clonotype frequency TRB",
    "Clonotype frequency IGH",
    "Clonotype frequency IGL",
    "Clonotype frequency Both",

    "Quartiles of RNA reads",
    "Quartiles of Genes",
    "Quartiles of Mitochondrial %",
    "Quartiles of ATAC fragments in peaks",
    "Quartiles of Peaks",
    "Quartiles of TSS enrichment score",
    "Quartiles of Nucleosome signal",
    "Quartiles of FRiP",
    "Quartiles of Bl. regions",

    "Prediction cell type",
    "Prediction mapping score",
    "Prediction confidence score"
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

cb_build <- function(object, slot, short_label, rootname, label_field, show_labels, features, markers, meta_fields, meta_fields_names, is_nested, dot_radius, dot_alpha, color_data) {
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
defQuantPal="tol-sq-blue"
defCatPal="rainbow"
radius=%i
alpha=%f
geneIdType="auto"
clusterField="%s"
labelField="%s"
enumFields=%s
showLabels=%s
coords=%s'

    local_config <- base::sprintf(
        local_config,
        short_label,
        short_label,
        dot_radius,
        dot_alpha,
        label_field,
        label_field,
        base::paste0("[", base::paste(base::paste0('"', enum_fields, '"'), collapse=", "), "]"),
        base::ifelse(show_labels, "True", "False"),
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
        local_config <- base::paste(
            local_config,
            'quickGenesFile="quickGenes.tsv"',
            sep="\n"
        )
    }

    if (!is.null(markers)){
        utils::write.table(
            markers,
            sep="\t",
            file=base::file.path(rootname, "markers.tsv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=TRUE
        )
        local_config <- base::paste(
            local_config,
            'markers=[{"file":"markers.tsv", "shortLabel":"Cluster-specific markers"}]',
            sep="\n"
        )
    }

    if (!is.null(color_data)){
        utils::write.table(
            color_data,
            sep=",",
            file=base::file.path(rootname, "colors.csv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE
        )
        local_config <- base::paste(
            local_config,
            'colors="colors.csv"',
            sep="\n"
        )
    }

    if (SeuratObject::DefaultAssay(object) == "ATAC"){                                        # for searching nearest peaks by gene name
        annotation_data <- base::as.data.frame(Signac::Annotation(object)) %>%                # if we have ATAC assay, Annotation should be also present
                           dplyr::filter(.$type=="transcript") %>%
                           dplyr::group_by(seqnames, strand, gene_id) %>%                     # to account for a replicates in the gene names
                           dplyr::arrange(desc(width), .by_group=TRUE) %>%
                           dplyr::slice_head(n=1) %>%                                         # takes the coordinates of the longest transcript within each group
                           dplyr::ungroup() %>%
                           dplyr::mutate("score"=0) %>%                                       # add empty score column to corresond to BED format
                           dplyr::select(
                               c("seqnames", "start", "end", "gene_id", "score", "strand")
                           )
        gene_symbols <- annotation_data %>%
                        dplyr::select(c("gene_id")) %>%
                        dplyr::mutate("gene_id_copy"=gene_id)

        cb_dir <- base::path.expand("~/cellbrowserData/genes")                                # should be saved in the user's home directory
        if (!base::dir.exists(cb_dir)) {
            base::dir.create(cb_dir, recursive=TRUE)
        }
        utils::write.table(
            annotation_data,
            base::gzfile(base::file.path(cb_dir, "genome.current.bed.gz")),
            sep="\t",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE
        )
        utils::write.table(
            gene_symbols,
            base::gzfile(base::file.path(cb_dir, "current.symbols.tsv.gz")),
            sep="\t",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE
        )
        exit_code <- sys::exec_wait(
            cmd="cbGenes",
            args=c(
                "json",
                "genome",
                "current",
                base::path.expand(base::file.path(cb_dir, "genome.current.json"))
            )
        )
        if (exit_code != 0){
            base::print(
                base::paste0(
                    "Failed to create genome.current.json ",
                    "with exit code ", exit_code, ". Exiting."
                )
            )
            base::quit(save="no", status=1, runLast=FALSE)
        }
        local_config <- base::paste(
            local_config,
            'atacSearch="genome.current"',
            sep="\n"
        )
        base::rm(annotation_data, gene_symbols)
        base::gc(verbose=FALSE)
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

get_color_data <- function(
    seurat_data,
    metadata_column,
    collected_color_data=NULL,
    palette_colors=graphics$D40_COLORS
){
    if (length(metadata_column) > 0){
        for (i in 1:length(metadata_column)){
            current_column <- metadata_column[i]
            categories <- base::levels(seurat_data@meta.data[[current_column]])       # will depend on the order of levels
            if (is.null(categories)) {                                                 # not a factor
                categories <- base::sort(                                              # alphabetically sorted
                    base::unique(
                        base::as.vector(
                            as.character(seurat_data@meta.data[[current_column]])
                        )
                    )
                )
            }
            current_color_data <- base::data.frame(
                category=categories,
                color=palette_colors[1:length(categories)],
                check.names=FALSE,
                stringsAsFactors=FALSE
            )
            if (!is.null(collected_color_data)){
                collected_color_data <- base::rbind(collected_color_data, current_color_data)
            } else {
                collected_color_data <- current_color_data
            }
        }
        collected_color_data <- collected_color_data %>%
                                dplyr::distinct(category, .keep_all=TRUE)
    }
    return (collected_color_data)
}


export_cellbrowser <- function(
    seurat_data, assay, slot,
    rootname,
    label_field=NULL,                        # the default label field (should be one of the metadata columns)
    features=NULL,
    dot_radius=3, dot_alpha=0.5,
    palette_colors=graphics$D40_COLORS,
    short_label="cellbrowser",
    is_nested=FALSE,
    meta_fields=NULL,
    meta_fields_names=NULL
){
    base::tryCatch(
        expr = {
            backup_assay <- SeuratObject::DefaultAssay(seurat_data)
            SeuratObject::DefaultAssay(seurat_data) <- assay                       # safety measure
            SeuratObject::Idents(seurat_data) <- "new.ident"                       # safety measure

            datasets_count <- length(base::unique(base::as.vector(as.character(seurat_data@meta.data$new.ident))))
            conditions_count <- length(base::unique(base::as.vector(as.character(seurat_data@meta.data$condition))))
            not_default_conditions <- all(
                base::as.vector(as.character(seurat_data@meta.data$new.ident)) != base::as.vector(as.character(seurat_data@meta.data$condition))
            )

            color_data <- base::data.frame(                                        # here we will collect all colors
                category=character(),
                color=character(),
                check.names=FALSE,
                stringsAsFactors=FALSE
            ) %>%
            tibble::add_row(category="NA", color=graphics$NA_COLOR) %>%            # color for NA
            tibble::add_row(category="G1", color=graphics$CC_COLORS[1]) %>%        # color for G1 cell cycle phase
            tibble::add_row(category="S", color=graphics$CC_COLORS[2]) %>%         # color for S cell cycle phase
            tibble::add_row(category="G2M", color=graphics$CC_COLORS[3])           # color for G2M cell cycle phase

            if (is.null(meta_fields) || is.null(meta_fields_names)){
                meta_fields <- DEFAULT_META_FIELDS
                meta_fields_names <- DEFAULT_META_FIELDS_NAMES
            }

            if (datasets_count > 1){
                meta_fields <- base::append(meta_fields, "new.ident", 0)
                meta_fields_names <- base::append(meta_fields_names, "Dataset", 0)
                color_data <- get_color_data(seurat_data, "new.ident", color_data, palette_colors)
            }

            if (conditions_count > 1 && not_default_conditions){
                meta_fields <- base::append(meta_fields, "condition", 1)
                meta_fields_names <- base::append(meta_fields_names, "Condition", 1)
                color_data <- get_color_data(seurat_data, "condition", color_data, palette_colors)
            }

            clustering_fields <- base::grep(
                "^(?!custom_).*_res\\.", base::colnames(seurat_data@meta.data),    # we exclude columns with "^custom" because the cell types may be stored in "custom_res."
                value=TRUE, ignore.case=TRUE, perl=TRUE                            # perl is true to have (?!) feature available
            )
            clustering_fields_names <- base::as.vector(
                base::unlist(                                                      # when nothing found it returns named list that should be unlisted
                    base::sapply(
                        clustering_fields,
                        function(line) {
                            base::paste(
                                "Clustering", base::paste0("(", base::gsub("_res.", " ", line), ")")
                            )
                        }
                    )
                )
            )
            meta_fields <- base::append(meta_fields, clustering_fields)
            meta_fields_names <- base::append(meta_fields_names, clustering_fields_names)
            color_data <- get_color_data(seurat_data, clustering_fields, color_data, palette_colors)

            custom_fields <- base::grep("^custom_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            custom_fields_names <- base::gsub("custom_", "Custom ", custom_fields)
            meta_fields <- base::append(meta_fields, custom_fields)
            meta_fields_names <- base::append(meta_fields_names, custom_fields_names)
            color_data <- get_color_data(seurat_data, custom_fields, color_data, palette_colors)

            pseudotime_fields <- base::grep("^ptime_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            pseudotime_fields_names <- base::gsub("ptime_", "Pseudotime from ", pseudotime_fields)
            meta_fields <- base::append(meta_fields, pseudotime_fields)
            meta_fields_names <- base::append(meta_fields_names, pseudotime_fields_names)

            collected_markers <- NULL                                                                     # we collect markers only for clustering, custom, and "prediction_cell_type" fields
            collected_top_features <- c()                                                                 # of the current assay
            updated_metadata_fields <- c()                                                                # to keep track of what metadata fields have been already updated
            if (("markers" %in% names(seurat_data@misc)) && (length(seurat_data@misc$markers) > 0)){
                for (i in 1:length(seurat_data@misc$markers)){                                            # iterating over assays
                    current_assay <- names(seurat_data@misc$markers)[i]
                    if (length(seurat_data@misc$markers[[current_assay]]) > 0){
                        for (j in 1:length(seurat_data@misc$markers[[current_assay]])){                   # iterating over markers sets within the current assay
                            current_field <- names(seurat_data@misc$markers[[current_assay]])[j]
                            current_markers <- NULL                                                       # to collect markers only from the current assay
                            if (current_field %in% clustering_fields){                                    # found markers for the clustering field
                                current_suffix <- base::paste0(                                           # to have "(rna 0.1)" in the cluster name
                                    "(", base::gsub("_res.", " ", current_field), ")"
                                )
                                if (!(current_field %in% updated_metadata_fields)){
                                    seurat_data@meta.data <- seurat_data@meta.data %>%                    # we want to rename labels even if it's not a current assay
                                                             dplyr::mutate(                               # so we can see, for example, RNA clusters labels when exporting
                                                                 !!current_field:=base::paste(            # ATAC. We keep track of which metadata column have been already
                                                                     .[[current_field]],                  # updated, because for some of them, such as wnn_res.0.1, we
                                                                     current_suffix                       # can have both RNA and ATAC markers, where each of them will
                                                                 )                                        # be shown in the correspondent tab of the UCSC Cellbrowser,
                                                             )                                            # however we don't want to rename them twice when iterating over
                                    updated_metadata_fields <- base::append(                              # RNA and ATAC assays.
                                        updated_metadata_fields,
                                        current_field
                                    )
                                }
                                if (current_assay == assay){                                              # we want to collect markers only from the current assay
                                    current_markers <- seurat_data@misc$markers[[current_assay]][[current_field]] %>%
                                                       dplyr::mutate(
                                                           "cluster"=base::paste(                         # changing them the same was as in the Seurat object
                                                                         .$cluster,
                                                                         current_suffix
                                                                     )
                                                       ) %>%
                                                       base::cbind(tmp_group=current_field, .)            # we need it to search for top features within each markers set
                                }
                            } else if (
                                (current_field %in% custom_fields) ||
                                (current_field=="prediction_cell_type")                                   # this is where we save Azimuth's predicted cell types
                            ){
                                if (current_assay == assay){                                              # we want to collect markers only from the current assay
                                    current_markers <- seurat_data@misc$markers[[current_assay]][[current_field]] %>%
                                                       base::cbind(tmp_group=current_field, .)            # we need it to search top features within each markers set
                                }
                            } else {
                                next                                                                      # to skip all unrecognized markers sets
                            }
                            if (!is.null(current_markers)){                                               # can be NULL because, for example, when exporting ATAC assay
                                if (is.null(collected_markers)) {                                         # we still want to update the RNA labels, but we don't need to
                                    collected_markers <- current_markers                                  # collect the markers from the RNA assay
                                } else {
                                    collected_markers <- base::rbind(collected_markers, current_markers)
                                }
                            }
                        }
                    }
                }
            }

            if (!is.null(features)){
                collected_top_features <- base::unique(features)                                                    # safety measure
            } else if (!is.null(collected_markers)){                                                                # will include markers only from the current assay
                collected_top_features <- collected_markers %>%
                                          dplyr::group_by(tmp_group) %>%                                            # to find top features within each markers set separately
                                          dplyr::filter(p_val_adj <= 0.05) %>%                                      # to have only significant gene markers
                                          dplyr::filter(pct.1 >= 0.1) %>%                                           # to have at least 10% of cells expressing this gene
                                          dplyr::group_by(feature) %>%                                              # groups by feature
                                          dplyr::arrange(desc(pct.1), .by_group=TRUE) %>%                           # sort all duplicated features by desc(pct.1)
                                          dplyr::slice_head(n=1) %>%                                                # choose the feature with the highest pct.1
                                          dplyr::ungroup(feature) %>%                                               # removes only the feature group
                                          dplyr::group_by(cluster) %>%
                                          dplyr::arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>%
                                          dplyr::group_modify(~ .x %>%
                                              dplyr::slice_head(n=analyses$get_fraction(.x, 0.0001, min_size=10))   # fraction is very small to have exactly 10 features per cluster
                                          ) %>%
                                          dplyr::ungroup() %>%                                                      # removes all groups
                                          dplyr::distinct(feature) %>%
                                          dplyr::pull(feature)                                                      # to have it as a vector
                if (assay == "ATAC"){                                                                               # when exporting ATAC assay quick genes file
                    collected_top_features <- stringr::str_replace(collected_top_features, "(chr\\w+)-", "\\1:")    # should include the ranges formatted as chr1:start-end
                }
            } else {
                collected_top_features <- NULL
            }

            if (!is.null(collected_markers)){
                collected_markers$tmp_group <- NULL                                  # after we found top features we don't need the tmp_group column anymore
            }

            show_labels <- TRUE                                                      # by default we try to show labels, but we hide them if
            if (is.null(label_field) || !(label_field %in% (meta_fields_names))){    # label_field was either not provided or not present in the metadata
                show_labels <- FALSE                                                 # hide labels
                for (i in 1:length(meta_fields)){                                    # searching the first field from the meta.data that is not unique for all cells
                    current_field <- meta_fields[i]
                    if (
                        (current_field %in% base::colnames(seurat_data@meta.data)) &&
                        (
                            length(
                                base::unique(base::as.vector(as.character(seurat_data@meta.data[, current_field])))
                            ) > 1
                        )
                    ){
                        label_field <- meta_fields_names[i]
                        break
                    }
                }
            }

            cb_build(
                seurat_data,
                slot=slot,
                short_label=short_label,
                rootname=rootname,
                label_field=label_field,
                show_labels=show_labels,
                features=collected_top_features,
                markers=collected_markers,
                meta_fields=meta_fields,
                meta_fields_names=meta_fields_names,
                is_nested=is_nested,
                dot_radius=dot_radius,
                dot_alpha=dot_alpha,
                color_data=color_data
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

export_loupe <- function(seurat_data, assay, rootname, active_cluster=NULL, meta_fields=NULL, meta_fields_names=NULL){
    base::tryCatch(
        expr = {

            datasets_count <- length(base::unique(base::as.vector(as.character(seurat_data@meta.data$new.ident))))
            conditions_count <- length(base::unique(base::as.vector(as.character(seurat_data@meta.data$condition))))
            not_default_conditions <- all(
                base::as.vector(as.character(seurat_data@meta.data$new.ident)) != base::as.vector(as.character(seurat_data@meta.data$condition))
            )

            if (is.null(meta_fields) || is.null(meta_fields_names)){
                meta_fields <- DEFAULT_META_FIELDS
                meta_fields_names <- DEFAULT_META_FIELDS_NAMES
            }

            if (datasets_count > 1){
                meta_fields <- base::append(meta_fields, "new.ident", 0)
                meta_fields_names <- base::append(meta_fields_names, "Dataset", 0)
            }

            if (conditions_count > 1 && not_default_conditions){
                meta_fields <- base::append(meta_fields, "condition", 1)
                meta_fields_names <- base::append(meta_fields_names, "Condition", 1)
            }

            clustering_fields <- base::grep(
                "^(?!custom_).*_res\\.", base::colnames(seurat_data@meta.data),    # we exclude columns with "^custom" because the cell types may be stored in "custom_res."
                value=TRUE, ignore.case=TRUE, perl=TRUE                            # perl is true to have (?!) feature available
            )
            clustering_fields_names <- base::as.vector(
                base::unlist(                                                      # when nothing found it returns named list that should be unlisted
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

            pseudotime_fields <- base::grep("^ptime_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
            pseudotime_fields_names <- base::gsub("ptime_", "Pseudotime from ", pseudotime_fields)
            meta_fields <- base::append(meta_fields, pseudotime_fields)
            meta_fields_names <- base::append(meta_fields_names, pseudotime_fields_names)

            all_colnames <- base::colnames(seurat_data@meta.data)              # need to have it outside of the loop
            for (i in 1:length(all_colnames)){                                 # removes all unused columns and renames remaining ones
                if (!(all_colnames[i] %in% meta_fields)){
                    base::print(base::paste("Removing", all_colnames[i], "column"))
                    seurat_data@meta.data[[all_colnames[i]]] <- NULL
                } else {
                    current_alias <- meta_fields_names[base::match(all_colnames[i], meta_fields)]
                    base::print(base::paste("Renaming", all_colnames[i], "column to", current_alias))
                    base::colnames(seurat_data@meta.data)[base::colnames(seurat_data@meta.data) == all_colnames[i]] <- current_alias
                }
            }

            SeuratObject::DefaultAssay(seurat_data) <- assay

            if(
                !is.null(active_cluster) &&
                (active_cluster %in% meta_fields) &&
                (meta_fields_names[base::match(active_cluster, meta_fields)] %in% base::colnames(seurat_data@meta.data))
            ){
                active_cluster <- meta_fields_names[base::match(active_cluster, meta_fields)]
                base::print(base::paste("Setting active cluster to", active_cluster))
                SeuratObject::Idents(seurat_data) <- active_cluster
            } else {
                base::print("Setting active cluster to Dataset")
                SeuratObject::Idents(seurat_data) <- "Dataset"     # Dataset should be always present, because it's our new.ident
            }

            loupeR::create_loupe_from_seurat(
                obj=seurat_data,                                   # will always use counts slot
                output_name=rootname,                              # should be without extension
                force=TRUE                                         # to overwrite existing file
            )

            base::print(
                base::paste(
                    "Exporting counts from the Seurat",
                    "object as Loupe file to", rootname
                )
            )
        },
        error = function(e){
            base::print(
                base::paste(
                    "Failed to export counts from the Seurat",
                    "object as Loupe file with error -", e
                )
            )
        }
    )
}