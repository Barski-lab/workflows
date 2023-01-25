import("cmapR", attach=FALSE)
import("dplyr", attach=FALSE)
import("purrr", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("tibble", attach=FALSE)
import("sceasy", attach=FALSE)
import("SeuratDisk", attach=FALSE)
import("rtracklayer", attach=FALSE)
import("GenomicRanges", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)

export(
    "get_file_type",
    "export_data",
    "export_rds",
    "export_gct",
    "export_cls",
    "load_cell_identity_data",
    "extend_metadata",
    "extend_metadata_by_barcode",
    "refine_metadata_levels",
    "load_grouping_data",
    "load_blacklist_data",
    "assign_identities",
    "load_10x_multiome_data",
    "load_10x_rna_data",
    "export_h5seurat",
    "export_h5ad",
    "load_cell_cycle_data",
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

load_cell_cycle_data <- function(location){
    if (!is.null(location)){
        cell_cycle_data <- utils::read.table(
            location,
            sep=get_file_type(location),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        base::print(
            base::paste(
                "Cell cycle data is successfully loaded from ", location
            )
        )
        return (cell_cycle_data)
    }
    base::print("Cell cycle data is not provided.")
    return (NULL)
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
            base::print(base::paste("Exporting data to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export data to", location, sep=" "))
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

export_h5ad <- function(data, location){
    base::tryCatch(
        expr = {
            sceasy::convertFormat(data, from="seurat", to="anndata", outFile=location)
            base::print(base::paste("Exporting data as h5ad to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export data as h5ad to", location, sep=" "))
        }
    )
}

load_cell_identity_data <- function(location) {
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
    return (cell_identity_data)
}

extend_metadata <- function(seurat_data, location, seurat_ref_column, meta_ref_column, seurat_target_columns=NULL){
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
    if (!all(base::sort(meta_ref_values) == base::sort(base::unique(seurat_ref_values)))){
        base::print("Extra metadata file is malformed. Exiting.")
        base::quit(save="no", status=1, runLast=FALSE)
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
        seurat_data[[current_seurat_target_column]] <- base::as.factor(                                  # need as.factor, otherwise the order on the plots might be mismatched
            metadata[[current_meta_source_column]][base::match(seurat_ref_values, meta_ref_values)]
        )
    }
    base::rm(metadata, seurat_ref_values, meta_ref_values)
    return (seurat_data)
}

refine_metadata_levels <- function(seurat_data){
    for (i in base::colnames(seurat_data@meta.data)){
        if (base::is.factor(seurat_data@meta.data[[i]])){
            base::print(base::paste("Re-evaluating levels for a factor column", i))
            base::print(
                base::paste(
                    "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
                )
            )
            seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  # need to drop levels of the removed values
            base::print(
                base::paste(
                    "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
                )
            )
        }
    }
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
    SeuratObject::Idents(seurat_data) <- "new.ident"
    if (base::nrow(cell_identity_data) > length(base::unique(base::as.vector(as.character(SeuratObject::Idents(seurat_data)))))){
        base::print("Identity file includes more than expected number of rows. Exiting.")
        base::quit(save="no", status=1, runLast=FALSE)
    }
    base::rm(idents, new_ident)
    return (seurat_data)
}

load_10x_multiome_data <- function(args, cell_identity_data, grouping_data) {
    base::suppressMessages(raw_data <- Seurat::Read10X(data.dir=args$mex))
    seurat_data <- SeuratObject::CreateSeuratObject(
        counts=raw_data$`Gene Expression`,
        min.cells=args$rnamincells,
        names.delim="-",
        names.field=2
    )
    annotation <- rtracklayer::import(args$annotations, format="GFF")

    if( !("gene_biotype" %in% base::colnames(GenomicRanges::mcols(annotation))) ){
        base::print(
            paste(
                "Loaded genome annotation doesn't have 'gene_biotype' column.",
                "Setting to NA"
            )
        )
        annotation$gene_biotype <- NA                                                               # some Signac functions fail without this column
    }
    if( !("tx_id" %in% base::colnames(GenomicRanges::mcols(annotation))) ){
        base::print(
            base::paste(
                "Loaded genome annotation doesn't have 'tx_id' column.",
                "Setting it from 'transcript_id'"
            )
        )
        annotation$tx_id <- annotation$transcript_id                                                # https://github.com/stuart-lab/signac/issues/1159
    }

    all_cells <- SeuratObject::Cells(seurat_data)
    names(all_cells) <- all_cells
    base::print(
        base::paste(
            "Preparing fragments for", length(all_cells), "cells"
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
            "Counting fragments overlapping", length(peak_coordinates), "peaks"
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
        annotation=annotation
    )
    base::print("Assigning new dataset identities")
    seurat_data <- assign_identities(seurat_data, cell_identity_data, grouping_data)
    base::rm(raw_data, annotation, all_cells, fragments, peak_coordinates, peak_counts)   # removing unused data
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

replace_fragments <- function(location, seurat_data){
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"                                # safety measure
    Signac::Fragments(seurat_data[["ATAC"]]) <- NULL                                 # remove old fragments
    all_cells <- SeuratObject::Cells(seurat_data)
    names(all_cells) <- all_cells
    base::print(base::paste("Preparing fragments for", length(all_cells), "cells"))
    fragments <- Signac::CreateFragmentObject(
        path=location,
        cells=all_cells,
        validate.fragments=TRUE,
        verbose=FALSE
    )
    Signac::Fragments(seurat_data[["ATAC"]]) <- fragments
    return(seurat_data)
}

export_gct <- function(counts_mat, row_metadata, col_metadata, location){
    base::tryCatch(
        expr = {
            row_metadata <- row_metadata %>%
                            tibble::rownames_to_column("id") %>%
                            dplyr::mutate_at("id", base::as.vector)
            col_metadata <- col_metadata %>%
                            tibble::rownames_to_column("id") %>%
                            dplyr::mutate_at("id", base::as.vector)
            gct_data <- methods::new(
                "GCT",
                mat=counts_mat[row_metadata$id, col_metadata$id],       # to guarantee the order and number of row/columns
                rdesc=row_metadata,
                cdesc=col_metadata
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