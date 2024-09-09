import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)

export(
    "print_info"
)


print_info <- function(seurat_data, args) {
    if (!is.null(args$verbose) && args$verbose){
        base::cat("Datasets:\n")
        base::cat(
            base::paste(
                base::unique(base::as.vector(as.character(seurat_data@meta.data$new.ident))),
                collapse="\n"
            )
        )
        base::cat("\n\nConditions:\n")
        base::cat(
            base::paste(
                base::unique(base::as.vector(as.character(seurat_data@meta.data$condition))),
                collapse="\n"
            )
        )
        base::cat("\n\nAssays:\n")
        base::print(seurat_data@assays)
        base::cat("Default assay:\n")
        base::cat(SeuratObject::DefaultAssay(seurat_data))
        base::cat("\n\nReductions:\n")
        base::print(seurat_data@reductions)
        base::cat("Graphs:\n")
        base::cat(
            base::paste(
                base::unique(base::as.vector(as.character(SeuratObject::Graphs(seurat_data)))),
                collapse="\n"
            )
        )
        base::cat("\n\nNeighbors:\n")
        base::cat(
            base::paste(
                base::unique(base::as.vector(as.character(SeuratObject::Neighbors(seurat_data)))),
                collapse="\n"
            )
        )
        base::cat("\n\nMiscellaneous:\n")
        if (("markers" %in% names(seurat_data@misc)) && (length(seurat_data@misc$markers) > 0)){
            for (i in 1:length(seurat_data@misc$markers)){
                current_assay <- names(seurat_data@misc$markers)[i]
                if (length(seurat_data@misc$markers[[current_assay]]) > 0){
                    for (j in 1:length(seurat_data@misc$markers[[current_assay]])){
                        current_field <- names(seurat_data@misc$markers[[current_assay]])[j]
                        base::print(
                            base::paste(
                                current_assay, "markers for",
                                current_field, "metadata column"
                            )
                        )
                        base::print(
                            utils::head(
                                seurat_data@misc$markers[[current_assay]][[current_field]]
                            )
                        )
                    }
                }
            }
        }
        if (("trajectories" %in% names(seurat_data@misc)) && (length(seurat_data@misc$trajectories) > 0)){
            base::print(
                base::paste(
                    "Trajectories for the following reductions:",
                    base::paste(names(seurat_data@misc$trajectories), collapse=", ")
                )
            )
        }
        if (("vdj" %in% names(seurat_data@misc)) && ("chains" %in% names(seurat_data@misc$vdj))){
            base::print(
                base::paste(
                    "VDJ chains:",
                    base::paste(seurat_data@misc$vdj$chains, collapse=", ")
                )
            )
        }
        if (
                ("atac_reduce" %in% names(seurat_data@misc)) &&
                ("first_lsi_removed" %in% names(seurat_data@misc$atac_reduce))
        ){
                base::print(
                    base::paste(
                        "The first LSI component of the ATAC dim. reduction removed:",
                        seurat_data@misc$atac_reduce$first_lsi_removed
                    )
                )
        }
        if ("ATAC" %in% names(seurat_data@assays)) {
            base::cat("\n\nFragments from ATAC assay:\n")
            fragments <- Signac::Fragments(seurat_data[["ATAC"]])
            if (length(fragments) > 0){
                for (i in 1:length(fragments)){
                    fragment <- fragments[[i]]
                    base::cat(
                        base::paste(
                            "Fragment", i, "includes",
                            length(methods::slot(fragment, "cells")), "cells",
                            "loaded from", methods::slot(fragment, "path"),
                            "with hash", methods::slot(fragment, "hash")[1],       # hash 1 is for Fragment file, hash 2 is for index
                            "\n"
                        )
                    )
                }
            }
        }
        base::cat("\nMetadata:\n")
        base::print(utils::head(seurat_data@meta.data))
        base::cat("\nMemory usage:\n")
        base::print(base::gc())
    } else {
        base::gc(verbose=FALSE)                                                    # silently cleans gc
    }
}