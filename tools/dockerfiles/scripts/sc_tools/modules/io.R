import("cmapR", attach=FALSE)
import("dplyr", attach=FALSE)
import("purrr", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("tibble", attach=FALSE)
import("sceasy", attach=FALSE)
import("forcats", attach=FALSE)
import("SeuratDisk", attach=FALSE)
import("SCopeLoomR", attach=FALSE)
import("rtracklayer", attach=FALSE)
import("scRepertoire", attach=FALSE)
import("GenomeInfoDb", attach=FALSE)
import("GenomicRanges", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)

export(
    "get_file_type",
    "export_data",
    "export_markers",
    "export_clonotypes",
    "export_rds",
    "export_gct",
    "export_cls",
    "load_cell_identity_data",
    "extend_metadata",
    "extend_metadata_by_barcode",
    "apply_metadata_filters",
    "refine_metadata_levels",
    "load_grouping_data",
    "load_blacklist_data",
    "assign_identities",
    "load_10x_multiome_data",
    "load_10x_rna_data",
    "load_10x_vdj_data",
    "load_10x_atac_data",
    "load_seqinfo_data",
    "export_h5seurat",
    "export_h5ad",
    "export_scope_loom",
    "export_fragments_coverage",
    "export_coverage",
    "load_cell_cycle_data",
    "load_annotation_data",
    "replace_fragments"
)


get_file_type <- function (filename){
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}

# load_barcodes_data <- function(location, seurat_data){
#     default_barcodes_data <- SeuratObject::Cells(seurat_data)                # to include all available cells
#     if (!is.null(location)){
#         barcodes_data <- utils::read.table(
#             location,
#             sep=get_file_type(location),
#             header=TRUE,
#             check.names=FALSE,
#             stringsAsFactors=FALSE
#         )
#         base::print(
#             base::paste(
#                 "Barcodes data is successfully loaded and from", location,
#                 "and will be used to prefilter feature-barcode matrices by",
#                 "cells of interest")
#         )
#         return (barcodes_data)
#     }
#     base::print("Barcodes data is not provided. Using all cells")
#     return (default_barcodes_data)
# }

load_cell_cycle_data <- function(seurat_data, location){
    base::tryCatch(
        expr = {
            cell_cycle_data <- utils::read.table(
                location,
                sep=get_file_type(location),
                header=TRUE,
                check.names=FALSE,
                stringsAsFactors=FALSE
            )
            if (
                length(
                    base::intersect(
                        base::as.vector(as.character(base::rownames(seurat_data))),         # all genes
                        cell_cycle_data$gene_id                                             # cc genes
                    )
                ) == 0
            ){
                base::print("Loaded cell cycle data has 0 genes present in the dataset")
                return (NULL)
            }
            base::print(
                base::paste(
                    "Cell cycle data is successfully loaded from", location
                )
            )
            return (cell_cycle_data)
        },
        error = function(e){
            base::print(
                base::paste(
                    "Failed to load cell cycle data from", location, "due to", e
                )
            )
            return (NULL)
        }
    )
}

export_fragments_coverage <- function(fragments_data, location){
    base::tryCatch(
        expr = {
            coverage_data <- methods::as(
                base::lapply(
                    GenomicRanges::coverage(fragments_data),
                    function(x) signif(10^6 * x/length(fragments_data), 3)               # scale to RPM mapped
                ),
                "SimpleRleList"
            )
            rtracklayer::export.bw(coverage_data, location)
            base::print(base::paste("Exporting ATAC fragments coverage data to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export ATAC fragments coverage data to", location, "due to", e))
        }
    )
}

export_coverage <- function(coverage_data, location){
    base::tryCatch(
        expr = {
            rtracklayer::export.bw(coverage_data, location)
            base::print(base::paste("Exporting coverage data to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export coverage data to", location, "due to", e))
        }
    )
}

export_data <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE){
    base::tryCatch(
        expr = {
            utils::write.table(
                data,
                file=location,
                sep=get_file_type(location),
                row.names=row_names,
                col.names=col_names,
                quote=quote
            )
            base::print(base::paste("Exporting data to", location))
        },
        error = function(e){
            base::print(base::paste("Failed to export data to", location, "due to", e))
        }
    )
}

export_markers <- function(data, assay, markers_regex, location){
    base::tryCatch(
        expr = {
            markers_fields <- base::grep(
                markers_regex,
                names(data@misc$markers[[assay]]),
                value=TRUE, ignore.case=TRUE
            )
            all_markers <- NULL
            for (i in 1:length(markers_fields)){
                current_markers <- data@misc$markers[[assay]][[markers_fields[i]]] %>%
                                   base::cbind(group=markers_fields[i], .)
                if (is.null(all_markers)) {
                    all_markers <- current_markers
                } else {
                    all_markers <- base::rbind(all_markers, current_markers)
                }
            }
            base::print(
                base::paste(
                    "Collecting", assay, "markers calculated for",
                    paste(markers_fields, collapse=", "), "metadata columns"
                )
            )
            export_data(all_markers, location)
        },
        error = function(e){
            base::print(base::paste("Failed to collect markers due to", e))
        }
    )
}

export_clonotypes <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE){
    base::tryCatch(
        expr = {

            conditions_count <- length(
                base::unique(
                    base::as.vector(as.character(data@meta.data$condition))
                )
            )
            not_default_conditions <- all(
                base::as.vector(as.character(data@meta.data$new.ident)) != base::as.vector(as.character(data@meta.data$condition))
            )

            clonotype_data <- base::data.frame(
                                  barcode=SeuratObject::Cells(data),
                                  dataset=data@meta.data$new.ident,
                                  donor=data@meta.data$donor
                              )

            if (conditions_count > 1 && not_default_conditions){
                clonotype_data <- base::cbind(clonotype_data, data@meta.data$condition)
            }

            clonotype_data <- base::cbind(
                clonotype_data,
                data@meta.data[, base::paste0("clonotype_", data@misc$vdj$chains)]
            )

            clonotype_data <- base::cbind(
                clonotype_data,
                data@meta.data[, base::paste0("clonalFrequency_", data@misc$vdj$chains)]
            )

            ct_gene <- stringr::str_split(data@meta.data$CTgene, "_", simplify=TRUE)    # may return "" values
            ct_aa <- stringr::str_split(data@meta.data$CTaa, "_", simplify=TRUE)        # may return "" values
            ct_nt <- stringr::str_split(data@meta.data$CTnt, "_", simplify=TRUE)        # may return "" values
            ct_gene[ct_gene == ""] <- NA                                                # for consistency with all clonotypes and frequency columns
            ct_aa[ct_aa == ""] <- NA                                                    # for consistency with all clonotypes and frequency columns
            ct_nt[ct_nt == ""] <- NA                                                    # for consistency with all clonotypes and frequency columns

            for (i in 1:length(data@misc$vdj$chains)){
                current_chain <- data@misc$vdj$chains[i]
                if (current_chain == "both") next                                       # we don't want to export both
                chain_position <- switch(
                    current_chain,
                    "TRA" = 1,
                    "TRB" = 2,
                    "IGL" = 1,
                    "IGH" = 2
                )
                clonotype_data <- clonotype_data %>%
                    dplyr::mutate(
                        !!paste(current_chain, "gene", sep="_") := ct_gene[, chain_position],
                        !!paste(current_chain, "aa", sep="_") := ct_aa[, chain_position],
                        !!paste(current_chain, "nt", sep="_") := ct_nt[, chain_position]
                    )
            }

            utils::write.table(
                clonotype_data,
                file=location,
                sep=get_file_type(location),
                row.names=row_names,
                col.names=col_names,
                quote=quote
            )
            base::print(base::paste("Exporting clonotypes to", location))
        },
        error = function(e){
            base::print(base::paste("Failed to export clonotypes to", location, "due to", e))
        }
    )
}

export_rds <- function(data, location){
    base::tryCatch(
        expr = {
            base::saveRDS(data, location)
            base::print(base::paste("Exporting data as RDS to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export data as RDS to", location, sep=" "))
        }
    )
}

export_h5seurat <- function(data, location, overwrite=TRUE){
    base::tryCatch(
        expr = {
            SeuratDisk::SaveH5Seurat(data, location, overwrite=overwrite)
            base::print(base::paste("Exporting data as h5seurat to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export data as h5seurat to", location, sep=" "))
        }
    )
}

export_h5ad <- function(data, location, assay="RNA", slot="counts"){
    base::tryCatch(
        expr = {
            sceasy::convertFormat(
                data,
                from="seurat",
                to="anndata",
                outFile=location,
                assay=assay,
                main_layer=slot
            )
            base::print(base::paste("Exporting data as h5ad to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export data as h5ad to", location, sep=" "))
        }
    )
}

export_scope_loom <- function(data, location, assay="RNA", slot="counts"){
    base::tryCatch(
        expr = {
            SeuratObject::DefaultAssay(data) <- assay
            SCopeLoomR::build_loom(
                file.name=location,
                dgem=SeuratObject::GetAssayData(object=data, slot=slot)
            )
            loom_handler <- SCopeLoomR::open_loom(location, mode="r+")              # open Loom file for adding data

            for (i in 1:length(data@reductions)){                                   # at least one reduction will be present, because
                reduction_name <- names(data@reductions)[i]                         # because we are not supposed to call this function
                SCopeLoomR::add_embedding(                                          # when there is nothing to show
                    loom=loom_handler,
                    embedding=base::as.data.frame(                                  # need as.data.frame, otherwise all.equal fails (in SCopeLoomR)
                        SeuratObject::Embeddings(data, reduction=reduction_name)
                    ),
                    name=reduction_name,
                    is.default=(i == 1)                                             # first found reduction will be default
                )
            }

            meta_fields <- base::append(
                c("new.ident", "condition"),                                        # these two fields should be always present in out Seurat object
                base::grep(                                                         # will be [] if no custom fields are present
                    "^custom_",
                    base::colnames(data@meta.data),
                    value=TRUE,
                    ignore.case=TRUE
                )
            )
            for (i in 1:length(meta_fields)){
                SCopeLoomR::add_col_attr(
                    loom=loom_handler,
                    key=meta_fields[i],
                    value=base::as.vector(
                        as.character(data@meta.data[ , meta_fields[i] ])
                    ),
                    as.annotation=TRUE                                              # adding as categorical value
                )
            }

            cluster_fields <- base::grep(
                "^rna_res\\.|^atac_res\\.|^wsnn_res\\.",                            # we need only columns with numeric cluster values
                base::colnames(data@meta.data),
                value=TRUE,
                ignore.case=TRUE
            )
            for (i in 1:length(cluster_fields)){
                cluster_data <- base::as.vector(as.character(data@meta.data[ , cluster_fields[i] ]))
                names(cluster_data) <- rownames(data@meta.data)
                SCopeLoomR::add_annotated_clustering(
                    loom=loom_handler,
                    group="Clustering",
                    name=cluster_fields[i],
                    clusters=cluster_data,
                    annotation=cluster_data
                )
            }

            SCopeLoomR::close_loom(loom_handler)

            base::print(base::paste("Exporting SCope compatible loom file to", location, sep=" "))
        },
        error = function(e){
            base::print(
                base::paste(
                    "Failed to export SCope compatible loom file to", location, "due to", e
                )
            )
        }
    )
}

load_cell_identity_data <- function(location) {
    default_cell_identity_data <- base::data.frame(
        library_id="Sample",
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    if (!is.null(location)){
        cell_identity_data <- utils::read.table(
            location,
            sep=get_file_type(location),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        if ("sample_id" %in% base::colnames(cell_identity_data)){             # the latest cellranger changed the name from library_id to sample_id
            cell_identity_data <- cell_identity_data %>%
                                  dplyr::rename("library_id"="sample_id")     # we use "library_id" in our code
        }
        # prepend with LETTERS, otherwise the order on the plot will be arbitrary sorted
        cell_identity_data <- cell_identity_data %>%
                              dplyr::mutate("library_id"=base::paste(LETTERS[1:base::nrow(cell_identity_data)], .$library_id))
        base::print(base::paste("Datasets identities were successfully loaded from ", location))
        return (cell_identity_data)
    }
    base::print("Datasets identities are not provided. Applying defaults")
    return (default_cell_identity_data)
}

extend_metadata <- function(seurat_data, location, seurat_ref_column, meta_ref_column, split_by=NULL, seurat_target_columns=NULL){
    metadata <- utils::read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    base::print(base::paste("Metadata is successfully loaded from", location))
    base::print(metadata)
    seurat_ref_values <- base::as.vector(as.character(seurat_data@meta.data[, seurat_ref_column]))
    meta_ref_values <- base::as.vector(as.character(metadata[, meta_ref_column]))
    if (!all(base::sort(meta_ref_values) == base::sort(base::unique(seurat_ref_values)))){       # need to check if metadata is valid before splitting in into the smaller
        base::print("Extra metadata file is malformed. Exiting.")                                # groups, because we use expand_grid that will add all possible values
        base::quit(save="no", status=1, runLast=FALSE)                                           # to meta_ref_values that might not be present in the updated seurat_ref_values
    }
    if(!is.null(split_by)){
        meta_ext_values <- base::unique(base::as.vector(as.character(seurat_data@meta.data[, split_by])))
        base::print(
            base::paste(
                "Splitting metadata by", split_by, "adding",
                paste(meta_ext_values, collapse=","), "suffixes"
            )
        )
        metadata <- metadata %>%
                    tidyr::expand_grid(sc_temp=meta_ext_values) %>%
                    dplyr::mutate(
                        dplyr::across(
                            -tidyselect::one_of("sc_temp"),
                            ~ base::paste(., base::get("sc_temp"), sep="_")
                        )
                    ) %>%
                    dplyr::select(-c("sc_temp")) %>%
                    base::as.data.frame()                                                               # otherwise won't work
        base::print("Updated metadata")
        base::print(metadata)
        seurat_ref_values <- base::paste(                                                               # need to overwrite the previous ones
            base::as.vector(as.character(seurat_data@meta.data[, seurat_ref_column])),
            base::as.vector(as.character(seurat_data@meta.data[, split_by])),
            sep="_"
        )
        meta_ref_values <- base::as.vector(as.character(metadata[, meta_ref_column]))                   # need to overwrite the previous ones
    }
    meta_source_columns <- base::colnames(metadata)[base::colnames(metadata) != meta_ref_column]        # all except ref column from metadata
    if (is.null(seurat_target_columns)){
        seurat_target_columns <- meta_source_columns                                                    # using the names from the metadata file
    } else if (length(seurat_target_columns) != length(meta_source_columns)){
        base::print("Provided target column names cannot be used with the loaded metadata. Exiting.")
        base::quit(save="no", status=1, runLast=FALSE)
    }
    for (i in 1:length(meta_source_columns)){
        current_meta_source_column <- meta_source_columns[i]
        current_seurat_target_column <- seurat_target_columns[i]
        base::print(base::paste("Adding", current_meta_source_column, "as", current_seurat_target_column, "column"))
        seurat_data[[current_seurat_target_column]] <- forcats::fct_infreq(                             # a factor ordered by frequency
            base::as.factor(
                metadata[[current_meta_source_column]][base::match(seurat_ref_values, meta_ref_values)]
            )
        )
    }
    base::rm(metadata, seurat_ref_values, meta_ref_values)
    return (seurat_data)
}

refine_metadata_levels <- function(seurat_data){
    for (current_column in base::colnames(seurat_data@meta.data)){
        if (base::is.factor(seurat_data@meta.data[[current_column]])){
            base::print(base::paste("Re-evaluating levels for a factor column", current_column))
            base::print(
                base::paste(
                    "before:", base::paste(base::levels(seurat_data@meta.data[[current_column]]), collapse=", ")
                )
            )
            seurat_data@meta.data[[current_column]] <- base::droplevels(seurat_data@meta.data[[current_column]])  # need to drop levels of the removed values
            base::print(
                base::paste(
                    "after:", base::paste(base::levels(seurat_data@meta.data[[current_column]]), collapse=", ")
                )
            )
        }
    }

    clustering_columns <- base::grep(                                            # should select only clustering columns
        "^rna_res\\.|^atac_res\\.|^wsnn_res\\.",
        base::colnames(seurat_data@meta.data),
        value=TRUE,
        ignore.case=TRUE
    )
    for (current_column in clustering_columns){                                  # can be either factor or vector
        base::tryCatch(
            expr = {
                base::print(
                    base::paste(
                        "Forcing", current_column, "clustering column",
                        "to be a factor with numerically sorted levels"
                    )
                )
                current_data <- as.numeric(
                    base::as.vector(                                             # in case we are working with factor, we need to take values, not levels
                        as.character(seurat_data@meta.data[[current_column]])
                    )
                )
                if(any(is.na(current_data))){                                    # safety measure to stop when converted non-numerical values
                    base::stop("Failed to convert to numeric")
                }
                seurat_data@meta.data[[current_column]] <- base::factor(
                    seurat_data@meta.data[[current_column]],
                    levels=as.character(                                         # need to be a character, otherwise it's not shown as categorical value on plots
                        base::sort(
                            base::unique(current_data)
                        )
                    )
                )
            },
            error = function(e){
                base::print(
                    base::paste(
                        "Failed to reorder", current_column,
                        "clustering column with error -", e
                    )
                )
            }
        )
    }
    return (seurat_data)
}

apply_metadata_filters <- function(seurat_data, target_column, target_values){
    base::print(
        base::paste(
            "Include only", base::paste(target_values, collapse=", "),
            "values from the", target_column, "metadata column."
        )
    )
    base::print(base::paste("Cells before filtering", base::nrow(seurat_data@meta.data)))
    SeuratObject::Idents(seurat_data) <- target_column
    seurat_data <- base::subset(seurat_data, idents=target_values)
    seurat_data <- refine_metadata_levels(seurat_data)                                    # to drop empty levels
    SeuratObject::Idents(seurat_data) <- "new.ident"
    base::print(base::paste("Cells after filtering", base::nrow(seurat_data@meta.data)))
    return (seurat_data)
}

extend_metadata_by_barcode <- function(seurat_data, location, filter=FALSE){
    metadata <- utils::read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    ) %>% dplyr::rename("barcode"=1)                                                      # rename the first column to barcode

    base::print(base::paste("Barcodes metadata is successfully loaded from ", location))

    if (filter){
        base::print(base::paste("Filtering Seurat data by loaded barcodes"))
        base::print(base::paste("Cells before filtering", base::nrow(seurat_data@meta.data)))
        seurat_data <- base::subset(seurat_data, cells=metadata$barcode)                  # subset to only selected cells
        seurat_data <- refine_metadata_levels(seurat_data)                                # to drop empty levels
        base::print(base::paste("Cells after filtering", base::nrow(seurat_data@meta.data)))
    }
    if (base::ncol(metadata) > 1){
        base::print(
            base::paste(
                "Extending Seurat object metadata with",
                base::paste(base::colnames(metadata)[2:length(base::colnames(metadata))], collapse=", "),
                "column(s)"
            )
        )
        refactored_metadata <- base::data.frame(SeuratObject::Cells(seurat_data)) %>%     # create a dataframe with only one column
                               dplyr::rename("barcode"=1) %>%                             # rename that column to barcode
                               dplyr::left_join(metadata, by="barcode") %>%               # intersect with loaded extra metadata by "barcode"
                               tibble::remove_rownames() %>%
                               tibble::column_to_rownames("barcode") %>%
                               replace(is.na(.), "Unknown")                               # in case Seurat object had more barcodes than loaded extra metadata
        seurat_data <- SeuratObject::AddMetaData(
            seurat_data,
            refactored_metadata[SeuratObject::Cells(seurat_data), , drop=FALSE]           # to guarantee the proper cells order
        )
    }
    SeuratObject::Idents(seurat_data) <- "new.ident"                                      # in case we updated new.ident columns with the new values
    return (seurat_data)
}

load_grouping_data <- function(location, cell_identity_data) {
    default_grouping_data <- base::data.frame(
        library_id=cell_identity_data$library_id,
        condition=cell_identity_data$library_id,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    if (!is.null(location)){
        grouping_data <- utils::read.table(
            location,
            sep=get_file_type(location),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        # prepend with LETTERS to correspond to the library_id from the cell_identity_data
        grouping_data <- grouping_data %>%
                          dplyr::mutate("library_id"=base::paste(LETTERS[1:base::nrow(grouping_data)], .$library_id))
        if ( (base::nrow(grouping_data) == base::nrow(cell_identity_data)) && all(base::is.element(cell_identity_data$library_id, grouping_data$library_id)) ){
            base::print(base::paste("Grouping data is successfully loaded from ", location))
            return (grouping_data)
        } else {
            base::print(base::paste("Applying defaults - grouping data loaded from", location, "is malformed"))
            return (default_grouping_data)
        }
    }
    base::print("Grouping data is not provided. Applying defaults")
    return (default_grouping_data)
}

load_blacklist_data <- function(location) {
    default_blacklist_data <- NULL
    if (!is.null(location)){
        blacklist_data <- rtracklayer::import(location, format="BED")
        base::print(base::paste("Genomic blacklist regions data is successfully loaded from ", location))
        return (blacklist_data)
    }
    base::print("File with the genomic blacklist regions is not provided. Cells won't be filtered by --maxblacklist")
    return (default_blacklist_data)
}

load_seqinfo_data <- function(location) {
    raw_data <- utils::read.table(
        location,
        sep=get_file_type(location),
        header=FALSE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    seqinfo_data <- GenomeInfoDb::Seqinfo(
        seqnames=raw_data[[1]],
        seqlengths=raw_data[[2]],
        isCircular=rep(FALSE, length(raw_data[[1]]))                                      # to be able to use trim on GenomicRanges object
    )
    base::print(base::paste("Chromosome length data is successfully loaded from ", location))
    return (seqinfo_data)
}

load_annotation_data <- function(location) {
    annotation_data <- rtracklayer::import(location, format="GFF")
    if( !("gene_biotype" %in% base::colnames(GenomicRanges::mcols(annotation_data))) ){
        base::print(
            paste(
                "Loaded genome annotation doesn't have gene_biotype column.",
                "Setting to NA"
            )
        )
        annotation_data$gene_biotype <- NA                                                               # some Signac functions fail without this column
    }
    if( !("tx_id" %in% base::colnames(GenomicRanges::mcols(annotation_data))) ){
        base::print(
            base::paste(
                "Loaded genome annotation doesn't have tx_id column.",
                "Setting it from transcript_id"
            )
        )
        annotation_data$tx_id <- annotation_data$transcript_id                                                # https://github.com/stuart-lab/signac/issues/1159
    }
    base::print(base::paste("Genome annotation data is successfully loaded from ", location))
    return (annotation_data)
}

assign_identities <- function(seurat_data, cell_identity_data, grouping_data){
    SeuratObject::Idents(seurat_data) <- "orig.ident"                                     # safety measure to make sure we get correct Idents
    idents <- as.numeric(as.character(SeuratObject::Idents(seurat_data)))                 # need to properly convert factor to numeric vector
    new_ident <- cell_identity_data$library_id[idents]
    if (sum(is.na(new_ident)) > 0){
        base::print("Identity file includes less than expected number of rows. Exiting.")
        base::quit(save="no", status=1, runLast=FALSE)
    }
    seurat_data[["new.ident"]] <- new_ident
    seurat_data[["condition"]] <- grouping_data$condition[base::match(seurat_data$new.ident, grouping_data$library_id)]
    if ("donor" %in% base::colnames(cell_identity_data)){                                 # if we work with aggregated Cell Ranger Multi runs we have this column
        seurat_data[["donor"]] <- cell_identity_data$donor[idents]
    }
    SeuratObject::Idents(seurat_data) <- "new.ident"
    if (base::nrow(cell_identity_data) > length(base::unique(base::as.vector(as.character(SeuratObject::Idents(seurat_data)))))){
        base::print("Identity file includes more than expected number of rows. Exiting.")
        base::quit(save="no", status=1, runLast=FALSE)
    }
    base::rm(idents, new_ident)
    return (seurat_data)
}

load_10x_atac_data <- function(args, cell_identity_data, grouping_data, seqinfo_data, annotation_data) {
    barcodes_location <- base::file.path(args$mex, "barcodes.tsv")
    peaks_location <- base::file.path(args$mex, "peaks.bed")
    all_cells <- utils::read.table(              # a vector of all cells found by Cell Ranger ATAC pipeline
        barcodes_location,
        sep=get_file_type(barcodes_location),
        header=FALSE,
        check.names=FALSE,
        stringsAsFactors=FALSE                   # will be already a vector, not factor
    ) %>% dplyr::pull(1)                         # barcodes.tsv doesn't have header, so we pull it by index
    names(all_cells) <- all_cells                # makes is named vector

    base::print(
        base::paste(
            "Preparing ATAC fragments for",
            length(all_cells), "cells"
        )
    )
    fragments <- Signac::CreateFragmentObject(
        path=args$fragments,
        cells=all_cells,
        validate.fragments=TRUE,
        verbose=FALSE
    )
    peak_coordinates <- rtracklayer::import(
        peaks_location,
        format="BED",
        genome=seqinfo_data                      # probably doesn't influence on anything
    )
    base::print(
        base::paste(
            "Counting ATAC fragments overlapping",
            length(peak_coordinates), "peaks"
        )
    )
    peak_counts <- Signac::FeatureMatrix(                                                 # need to recalculate feature-barcode matrix to include fragments, not reads
        fragments=fragments,
        sep=c("-", "-"),                                                                  # unexplainable bug - fails if use sep=c(":", "-")
        features=peak_coordinates,
        cells=all_cells,
        verbose=FALSE
    )

    seurat_data <- SeuratObject::CreateSeuratObject(
        counts=Signac::CreateChromatinAssay(
            counts=peak_counts,
            sep=c("-", "-"),
            fragments=fragments,
            min.cells=args$atacmincells,
            annotation=annotation_data,
            genome=seqinfo_data                                                           # we will need it when we want to export genome coverage to bigWig files
        ),
        assay="ATAC",
        names.delim="-",
        names.field=2
    )
    base::print("Assigning new dataset identities")
    seurat_data <- assign_identities(seurat_data, cell_identity_data, grouping_data)
    base::rm(all_cells, fragments, peak_coordinates, peak_counts)               # removing unused data
    base::gc(verbose=FALSE)
    return (seurat_data)
}

load_10x_multiome_data <- function(args, cell_identity_data, grouping_data, seqinfo_data, annotation_data) {
    base::suppressMessages(raw_data <- Seurat::Read10X(data.dir=args$mex))
    seurat_data <- SeuratObject::CreateSeuratObject(
        counts=raw_data$`Gene Expression`,
        min.cells=base::ifelse(                # in case we call it from sc_atac_filter.R script
            is.null(args$rnamincells),         # and we don't have rnamincells parameter.
            1,                                 # This influences only on the number of genes and
            args$rnamincells                   # is not used as RNA assay is not needed in sc_atac_filter.R
        ),
        names.delim="-",
        names.field=2
    )
    all_cells <- SeuratObject::Cells(seurat_data)
    names(all_cells) <- all_cells
    base::print(
        base::paste(
            "Preparing ATAC fragments for",
            length(all_cells), "cells"
        )
    )
    fragments <- Signac::CreateFragmentObject(
        path=args$fragments,
        cells=all_cells,
        validate.fragments=TRUE,
        verbose=FALSE
    )
    peak_coordinates <- Signac::StringToGRanges(
        regions=base::rownames(raw_data$Peaks),
        sep=c(":", "-")                            # need to run it explicitely with correct sep, as in FeatureMatrix 1) they have bug, 2) for consistency we'll use ("-", "-")
    )
    base::print(
        base::paste(
            "Counting ATAC fragments overlapping",
            length(peak_coordinates), "peaks"
        )
    )
    peak_counts <- Signac::FeatureMatrix(                                                 # need to recalculate feature-barcode matrix to include fragments, not reads
        fragments=fragments,
        sep=c("-", "-"),                                                                  # unexplainable bug - fails if use sep=c(":", "-")
        features=peak_coordinates,
        cells=all_cells,
        verbose=FALSE
    )
    seurat_data[["ATAC"]] <- Signac::CreateChromatinAssay(
        counts=peak_counts,
        sep=c("-", "-"),
        fragments=fragments,
        min.cells=args$atacmincells,
        annotation=annotation_data,
        genome=seqinfo_data                                                               # we will need it when we want to export genome coverage to bigWig files
    )
    base::print("Assigning new dataset identities")
    seurat_data <- assign_identities(seurat_data, cell_identity_data, grouping_data)
    base::rm(raw_data, all_cells, fragments, peak_coordinates, peak_counts)               # removing unused data
    base::gc(verbose=FALSE)
    return (seurat_data)
}

load_10x_rna_data <- function(args, cell_identity_data, grouping_data) {
    if (length(args$mex) == 1){
        base::print("Single feature-barcode matrix is provided. Using the original barcode suffixes")
        seurat_data <- SeuratObject::CreateSeuratObject(
            counts=Seurat::Read10X(data.dir=args$mex),
            min.cells=args$rnamincells,
            names.delim="-",
            names.field=2
        )
        base::print("Assigning new dataset identities")
        seurat_data <- assign_identities(seurat_data, cell_identity_data, grouping_data)
        return (seurat_data)
    } else {
        base::print("Multiple feature-barcode matrices are provided. Original barcode suffixes will be updated")
        merged_seurat_data <- NULL
        for (i in 1:length(args$mex)){
            current_location <- args$mex[i]
            base::print(
                base::paste(
                    "Reading 10x data from", current_location,
                    "replacing the original barcode suffixes with", i
                )
            )
            seurat_data <- SeuratObject::CreateSeuratObject(
                counts=Seurat::Read10X(
                    data.dir=current_location,
                    strip.suffix=TRUE             # removes suffix from barcode
                )
            )
            idents <- i
            new_ident <- cell_identity_data$library_id[idents]
            base::print(base::paste("Assigning new identity", new_ident))
            if (sum(is.na(new_ident)) > 0){
                base::print("Identity file includes less than expected number of rows. Exiting.")
                base::quit(save="no", status=1, runLast=FALSE)
            }
            seurat_data[["new.ident"]] <- new_ident
            seurat_data[["condition"]] <- grouping_data$condition[base::match(seurat_data$new.ident, grouping_data$library_id)]
            SeuratObject::Idents(seurat_data) <- "new.ident"
            seurat_data <- SeuratObject::RenameCells(
                seurat_data,
                new.names=base::paste0(SeuratObject::Cells(seurat_data), "-", idents)      # to add new barcode suffix
            )
            if (is.null(merged_seurat_data)){
                merged_seurat_data <- seurat_data
            } else {
                merged_seurat_data <- base::merge(merged_seurat_data, y=seurat_data)
            }
            base::rm(seurat_data, idents, new_ident)                                       # remove unused data
        }
        if (base::nrow(cell_identity_data) > length(base::unique(base::as.vector(as.character(SeuratObject::Idents(merged_seurat_data)))))){
            base::print("Identity file includes more than expected number of rows. Exiting.")
            base::quit(save="no", status=1, runLast=FALSE)
        }
        base::gc(verbose=FALSE)
        return (merged_seurat_data)
    }
}

load_10x_vdj_data <- function(seurat_data, args) {
    base::print(
        base::paste(
            "Loading congtigs annotations from", args$contigs
        )
    )

    if (!("donor" %in% base::colnames(seurat_data@meta.data))){
        datasets_count <- length(base::unique(base::as.vector(as.character(seurat_data@meta.data$new.ident))))
        if (datasets_count > 1){
            base::print(
                base::paste(
                    "Loaded data includes more then 1 dataset,",        # aggregated data from multiple Cell Ranger Multi samples will always include donor column
                    "but donor column is not present. Exiting."
                )
            )
            base::quit(save="no", status=1, runLast=FALSE)
        }
        base::print(
            base::paste(
                "Assigning all cells to the same default",
                "donor, as only one dataset is loaded."
            )
        )
        seurat_data@meta.data[["donor"]] <- seurat_data@meta.data$new.ident
    }

    congtigs_data <- scRepertoire::createHTOContigList(            # should work even if seurat_data includes only one dataset
        contig=utils::read.table(
            args$contigs,
            sep=get_file_type(args$contigs),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        ),
        sc=seurat_data,                                            # the order of barcode suffixes should be the same as in contigs
        group.by="donor"                                           # we group by donor to call clonotypes similar to what Cell Ranger Multi does
    )

    filter_by_cells=!is.null(args$filter) && args$filter == "cells"
    filter_by_chains=!is.null(args$filter) && args$filter == "chains"
    base::print(
        base::paste0(
            "Combining contigs.",
            base::ifelse(filter_by_cells, " Filtering by cells.", ""),
            base::ifelse(filter_by_chains, " Filtering by chains.", "")
        )
    )
    detected_chains <- base::sort(                                             # we need to combine chains from all if the list items
        base::unique(
            base::unlist(base::lapply(congtigs_data, function(d) d$chain))
        )
    )
    base::print(
        base::paste(
            "Detected chains:", base::paste(detected_chains, collapse=", ")
        )
    )
    if (base::identical(c("TRA", "TRB"), detected_chains)) {
        congtigs_data <- scRepertoire::combineTCR(                             # combining contigs into clonotypes within each donor, because congtigs_data is a list
            input.data=congtigs_data,
            removeNA=FALSE,                                                    # should be FALSE, because of https://github.com/ncborcherding/scRepertoire/issues/293
            removeMulti=filter_by_cells,
            filterMulti=filter_by_chains
        )
    } else if (base::identical(c("IGH", "IGL"), detected_chains)){
        congtigs_data <- scRepertoire::combineBCR(
            input.data=congtigs_data,
            removeNA=FALSE,                                                    # not tested for BCR, but better to follow the same logic as for combineTCR
            removeMulti=filter_by_cells,
            filterMulti=filter_by_chains
        )
    } else {
        base::print("Not implemented chains detected. Exiting")
        quit(save="no", status=1, runLast=FALSE)
    }

    if (args$removepartial){
        base::print("Removing cells with only one chain detected")
        for(i in seq_along(congtigs_data)) {
            congtigs_data[[i]] <- congtigs_data[[i]] %>%
                tidyr::drop_na(
                    tidyselect::any_of(
                        c("TCR1", "TCR2", "IGH", "IGLC")        # doesn't fail because of the missing columns, so it can be used for both TCR and BCR
                    )
                )
        }
    }

    seurat_data@misc$vdj <- list(
        chains=c(detected_chains, "both")                       # will always be either c("TRA", "TRB", "both") or c("IGH", "IGL", "both")
    )
    for (i in 1:length(seurat_data@misc$vdj$chains)){           # we want to calculate it per chain, because clonalFrequency that we filter by will be different
        current_chain <- seurat_data@misc$vdj$chains[i]

        seurat_data <- scRepertoire::combineExpression(         # creates clonalProportion (n/nrow) and clonalFrequency (n) columns calculated within donor based on the selected chain
            input.data=congtigs_data,                           # this is a list split by donor
            sc.data=seurat_data,
            cloneCall=args$cloneby,
            chain=current_chain,                                # chains are concatenated with _ when current_chain is both
            proportion=TRUE,                                    # we don't use calculated clonalProportion, so this parameter is not important
            group.by=NULL                                       # congtigs_data is already a list, so we don't need to provide group.by
        )
        base::print(utils::head(seurat_data@meta.data, n=5))

        seurat_data@meta.data <- seurat_data@meta.data %>%
                                 dplyr::mutate(                                                              # to have 0 instead of NA for numeric column
                                     clonalFrequency=dplyr::if_else(                                         # calculated per donor
                                         is.na(clonalFrequency),
                                         0,
                                         clonalFrequency
                                     )
                                 ) %>%
                                 dplyr::mutate(                                                              # pretty formating because we show it in the UCSC Cell Browser
                                     !!base::paste0("clonotype_", current_chain):=dplyr::if_else(            # use temp for easy access
                                         clonalFrequency < args$minfrequency,                                # calculated for current_chain
                                         NA,
                                         base::gsub(
                                             ";|_|\\.+", " ",                                                # ; means multiple chains, _ splits chains in "both"
                                             base::gsub(
                                                "NA|None", "",                                               # to remove all NA and None values
                                                .[[scRepertoire:::.convertClonecall(args$cloneby)]]
                                             )
                                         )
                                     )
                                 ) %>%
                                 dplyr::rename(
                                     !!base::paste0("clonalFrequency_", current_chain):="clonalFrequency"    # this column is not filtered
                                 ) %>%
                                 dplyr::select(-c("clonalProportion", "cloneSize"))                          # we also remove cloneSize even if it's used by scRepertoire
                                                                                                             # in some plots (they remove by NA). We want to filter by
                                                                                                             # our clonotype_[chain] column. We can't use multiple cloneSize
                                                                                                             # columns because cloneSize is harcoded in scRepertoire
    }

    base::rm(congtigs_data)
    base::gc(verbose=FALSE)
    return (seurat_data)
}

replace_fragments <- function(location, seurat_data){
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"                                # safety measure
    Signac::Fragments(seurat_data[["ATAC"]]) <- NULL                                 # remove old fragments
    all_cells <- SeuratObject::Cells(seurat_data)
    names(all_cells) <- all_cells
    base::print(base::paste("Preparing ATAC fragments for", length(all_cells), "cells"))
    fragments <- Signac::CreateFragmentObject(
        path=location,
        cells=all_cells,
        validate.fragments=TRUE,
        verbose=FALSE
    )
    Signac::Fragments(seurat_data[["ATAC"]]) <- fragments
    return(seurat_data)
}

export_gct <- function(counts_mat, location, row_metadata=NULL, col_metadata=NULL){
    base::tryCatch(
        expr = {
            if (!is.null(row_metadata)){
                row_metadata <- row_metadata %>%
                                tibble::rownames_to_column("id") %>%
                                dplyr::mutate_at("id", base::as.vector)
                counts_mat <- counts_mat[row_metadata$id, ]                      # to guarantee the order and number of rows
            }
            if (!is.null(col_metadata)){
                col_metadata <- col_metadata %>%
                                tibble::rownames_to_column("id") %>%
                                dplyr::mutate_at("id", base::as.vector)
                counts_mat <- counts_mat[, col_metadata$id]                      # to guarantee the order and number of columns
            }
            gct_data <- methods::new(
                "GCT",
                mat=counts_mat,
                rdesc=row_metadata,                                              # can be NULL
                cdesc=col_metadata                                               # can be NULL
            )
            cmapR::write_gct(
                ds=gct_data,
                ofile=location,
                appenddim=FALSE
            )
            base::print(base::paste("Exporting GCT data to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export GCT data to ", location, "with error - ", e, sep=""))
        }
    )
}

export_cls <- function(categories, location){
    base::tryCatch(
        expr = {
            output_stream <- base::file(location, "w")
            on.exit(base::close(output_stream), add=TRUE)           # can't put it in 'finally' as there is no access to output_stream variable
            base::cat(
                base::paste(
                    length(categories),                             # number of datasets
                    length(base::levels(categories)),               # number of different categories
                    "1",                                            # should be always 1
                    sep="\t"
                ),
                base::paste(
                    "#",
                    base::paste(
                        base::unique(as.character(categories)),     # preserves the order, but removes duplicates
                        collapse="\t"
                    ),
                    sep="\t"
                ),
                base::paste(
                    base::paste(
                        as.character(categories),
                        collapse="\t"
                    ),
                    sep="\t"
                ),
                file=output_stream,
                sep="\n"
            )
            base::print(base::paste("Exporting CLS data to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export CLS data to ", location, "with error - ", e, sep=""))
        }
    )
}