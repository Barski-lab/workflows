#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(escape))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell Gene Set Enrichment Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include genes",
            "expression information stored in the RNA assay. Presence of the ATAC assay is",
            "optional. Additionally, 'rnaumap', and/or 'atacumap', and/or 'wnnumap' dimensionality",
            "reductions should be present."
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
        "--species",
        help=paste(
            "Species type to select gene sets from the Molecular Signatures Database (MSigDB).",
            "Default: 'Homo sapiens'"
        ),
        type="character", default="Homo sapiens",
        choices=c("Homo sapiens", "Mus musculus")
    )
    parser$add_argument(
        "--geneset",
        help=paste(
            "Annotated gene set from the Molecular Signatures Database (MSigDB).",
            "Default: H"
        ),
        type="character", default="H",
        choices=c(
            "H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"
        )
    )
    parser$add_argument(
        "--method",
        help=paste(
            "Single-cell gene signature scoring method.",
            "Default: ssGSEA"
        ),
        type="character", default="ssGSEA",
        choices=c("ssGSEA", "UCell")
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

print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)

if ( !("RNA" %in% names(seurat_data@assays)) ){
    print(
        paste(
            "Loaded Seurat object doesn't include required RNA assay.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

if (!any(c("rnaumap", "atacumap", "wnnumap") %in% names(seurat_data@reductions))){
    print(
        paste(
            "Loaded Seurat object includes neither of the required reductions:",
            "'rnaumap', and/or 'atacumap', and/or 'wnnumap'.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"
debug$print_info(seurat_data, args)

if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)
    debug$print_info(seurat_data, args)
}

print("Running Gene Set Enrichment Analysis")
gene_sets <- getGeneSets(
    species=args$species,
    library=args$geneset
)
gsea_nes_df <- enrichIt(                          # data frame of normalized enrichmenet scores (NES)
    obj=seurat_data, 
    gene.sets=gene_sets,
    method=args$method,
    groups=1000,
    cores=args$cpus, 
    min.size=15
)

gsea_nes_matrix <- t(as.matrix(gsea_nes_df))
seurat_data[["GSEA"]] <- CreateAssayObject(
    counts=gsea_nes_matrix,                       # data slot will be equal to counts
    min.cells=0,                                  # to include all features
    min.features=0                                # to include all cells
)

debug$print_info(seurat_data, args)

if(args$cbbuild){
    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        print("Exporting RNA, ATAC, and GSEA assays to UCSC Cellbrowser jointly")
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/rna", sep=""),
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/atac", sep=""),
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="GSEA",
            slot="counts",
            short_label="GSEA",
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/gsea", sep=""),
        )
    } else {
        print("Exporting RNA and GSEA assays to UCSC Cellbrowser jointly")
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/rna", sep=""),
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="GSEA",
            slot="counts",
            short_label="GSEA",
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/gsea", sep=""),
        )
    }
}

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