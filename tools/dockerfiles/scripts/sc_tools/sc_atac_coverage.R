#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(argparse))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell ATAC-Seq Genome Coverage")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file",
            "should include chromatin accessibility information stored",
            "in the ATAC assay with a proper seqinfo data."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment used in the",
            "loaded Seurat object. File should be saved in TSV format and to be",
            "tbi-indexed."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to split cells into groups.",
            "May be one of the columns added with --metadata or --barcodes",
            "parameters. Default: split by dataset"
        ),
        type="character", nargs="*", default=c("new.ident")
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to optionally extend Seurat object metadata with",
            "categorical values using samples identities. First column - 'library_id'",
            "should correspond to all unique values from the 'new.ident' column of the",
            "loaded Seurat object. If any of the provided in this file columns are already",
            "present in the Seurat object metadata, they will be overwritten. When combined",
            "with --barcodes parameter, first the metadata will be extended, then barcode",
            "filtering will be applied. Default: no extra metadata is added"
        ),
        type="character"
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
        "--flank",
        help=paste(
            "Distance in bp to flank both start and end of the each fragment in both",
            "direction to generate cut sites coverage. Default: 5"
        ),
        type="integer", default=5
    )
    parser$add_argument(
        "--verbose",
        help="Print debug information. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--tmpdir",
        help=paste(
            "Directory to keep temporary files. Default: either /tmp",
            "or defined by environment variables TMPDIR, TMP, TEMP."
        ),
        type="character", default=tempdir()
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
print("Setting default assay to ATAC")
DefaultAssay(seurat_data) <- "ATAC"
debug$print_info(seurat_data, args)

seqinfo_data <- seurat_data[["ATAC"]]@seqinfo
if (is.null(seqinfo_data)){
    print("Loaded Seurat object doesn't include seqinfo data. Exiting.")
    quit(save="no", status=1, runLast=FALSE)
}

if (!is.null(args$metadata)){
    print("Extending Seurat object with the extra metadata fields")
    seurat_data <- io$extend_metadata(
        seurat_data=seurat_data,
        location=args$metadata,
        seurat_ref_column="new.ident",
        meta_ref_column="library_id"
    )
    debug$print_info(seurat_data, args)
}

if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)
}
debug$print_info(seurat_data, args)

if (!all(args$splitby %in% colnames(seurat_data@meta.data))){
    print(
        paste(
            "Loaded Seurat object can't be split by", paste(args$splitby, collapse=", "),
            "column(s). Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

print(paste("Loading fragments data from", args$fragments))
seurat_data <- io$replace_fragments(args$fragments, seurat_data)
debug$print_info(seurat_data, args)

peaks_location <- paste0(args$output, "_peaks.bigBed")
print(paste("Exporting peaks data to", peaks_location))
peaks_data <- seurat_data[["ATAC"]]@ranges
seqlevels(peaks_data) <- seqlevels(seqinfo_data)
seqinfo(peaks_data) <- seqinfo_data
peaks_data$score <- 0
export.bb(peaks_data, peaks_location)

print(
    paste(
        "Spliting cells into groups by", paste(args$splitby, collapse=", "),
        "metadata column(s)."
    )
)
seurat_data@meta.data <- seurat_data@meta.data %>%
                         tidyr::unite(splitby, all_of(args$splitby), remove=FALSE) %>%
                         mutate(splitby=tolower(gsub("#|%|&| ", "_", .$splitby)))
debug$print_info(seurat_data, args)

SplitFragments(
    seurat_data,
    assay="ATAC",
    group.by="splitby",           # we added this column based on args$splitby value
    outdir=args$tmpdir,
    verbose=args$verbose
)

groups <- unique(as.vector(as.character(seurat_data@meta.data$splitby)))
for (i in 1:length(groups)){
    fragments_cov_location <- paste0(args$output, "_", groups[i], "_frg_cov.bigWig")
    cut_sites_cov_location <- paste0(args$output, "_", groups[i], "_cut_cov.bigWig")
    tryCatch(
        expr = {
            fragments_data <- rtracklayer::import(
                paste0(args$tmpdir, "/", groups[i], ".bed"),
                format="BED",
                genome=seqinfo_data
            )
            io$export_fragments_coverage(
                fragments_data=fragments_data,
                location=fragments_cov_location
            )
            fragments_data <- unlist(as(
                list(
                    flank(fragments_data, args$flank, both=TRUE, start=TRUE, ignore.strand=TRUE),
                    flank(fragments_data, args$flank, both=TRUE, start=FALSE, ignore.strand=TRUE)
                ),
                "GRangesList"
            ))
            io$export_fragments_coverage(
                fragments_data=fragments_data,
                location=cut_sites_cov_location
            )
        },
        error = function(e){
            print(
                paste(
                    "Failed to export coverage data with error - ", e, sep=""
                )
            )
        }
    )
}
