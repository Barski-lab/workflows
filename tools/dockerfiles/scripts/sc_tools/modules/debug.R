import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)

export(
    "print_info"
)


print_info <- function(seurat_data, args) {
    if (!is.null(args$verbose) && args$verbose){
        base::print("Assays")
        base::print(seurat_data@assays)
        base::print("Reductions")
        base::print(seurat_data@reductions)
        base::print("Graphs")
        base::print(SeuratObject::Graphs(seurat_data))
        base::print("Neighbors")
        base::print(SeuratObject::Neighbors(seurat_data))
        base::print("Default assay")
        base::print(SeuratObject::DefaultAssay(seurat_data))
        base::print("Metadata")
        base::print(utils::head(seurat_data@meta.data))
        if ("ATAC" %in% names(seurat_data@assays)) {
            base::print("Fragments from ATAC assay")
            fragments <- Signac::Fragments(seurat_data[["ATAC"]])
            if (length(fragments) > 0){
                for (i in 1:length(fragments)){
                    fragment <- fragments[[i]]
                    base::print(
                        base::paste(
                            "Fragment", i, "includes",
                            length(methods::slot(fragment, "cells")), "cells",
                            "loaded from", methods::slot(fragment, "path"),
                            "with hash", methods::slot(fragment, "hash")[1]             # hash 1 is for Fragment file, hash 2 is for index
                        )
                    )
                }
            }
        }
        base::print(base::gc())
    } else {
        base::gc(verbose=FALSE)    # silently clean gc
    }
}