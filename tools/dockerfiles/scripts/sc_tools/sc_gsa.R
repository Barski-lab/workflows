#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(escape))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(GSEABase))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


export_all_plots <- function(seurat_data, gene_sets, args) {
    DefaultAssay(seurat_data) <- "GSA"                            # safety measure
    Idents(seurat_data) <- args$groupby
    labels <- sub(".*?_", "", tolower(gene_sets))
    labels <- gsub("_|-", " ", labels)

    graphics$dot_plot(
        data=seurat_data,                                         # will use data slot
        features=gene_sets,
        plot_title=paste("Scaled gene signature scores per group"),
        x_label="Gene set",
        y_label="Group",
        cluster_idents=FALSE,
        theme=args$theme,
        rootname=paste(args$output, "gsa_avg", sep="_"),
        pdf=args$pdf
    )
    for (i in 1:length(gene_sets)){
        current_gene_set <- gene_sets[i]
        current_label <- labels[i]
        graphics$vln_plot(
            data=seurat_data,
            features=current_gene_set,
            labels=current_label,
            plot_title=paste("Gene signature score density per cell group"),
            legend_title="Group",
            log=FALSE,
            pt_size=0,
            combine_guides="collect",
            width=800,
            height=600,
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "gsa_dnst", gsub(" ", "_", current_label), sep="_"),
            pdf=args$pdf
        )
        for (reduction in c("rnaumap", "atacumap", "wnnumap")){
            if (!(reduction %in% names(seurat_data@reductions))) {next}                                  # skip missing reductions
            graphics$feature_plot(
                data=seurat_data,
                features=current_gene_set,
                labels=current_label,
                reduction=reduction,
                plot_title=paste0("Gene signature score on cells UMAP (", reduction, " dim. reduction)"),
                label=TRUE,
                order=TRUE,
                max_cutoff="q99",  # to prevent cells with overexpressed gene from distoring the color bar
                combine_guides="keep",
                width=800,
                height=800,
                theme=args$theme,
                rootname=paste(args$output, "gsa_per_cell_rd", reduction, gsub(" ", "_", current_label), sep="_"),
                pdf=args$pdf
            )
        }
    }

    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
}


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell Gene Set Analysis")                               # https://academic.oup.com/bib/article/17/3/393/1744776
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
        "--groupby",
        help=paste(
            "Column from the metadata of the loaded Seurat object to group cells by.",
            "May be one of the columns added by --barcodes parameter."
        ),
        type="character", required="True"
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
            "Gene signature scoring method.",
            "Default: ssGSEA"
        ),
        type="character", default="ssGSEA",
        choices=c("ssGSEA", "UCell")
    )
    parser$add_argument(
        "--mingenes",
        help=paste(
            "Minimum number of genes from the loaded Seurat object to be present",
            "in the gene set to perform gene set analysis.",
            "Default: 15"
        ),
        type="integer", default=15
    )
    parser$add_argument(
        "--topn",
        help=paste(
            "Show top N the most variable gene sets",
            "Default: 10"
        ),
        type="integer", default=10
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

print("Running Gene Set Analysis")
gene_sets <- getGeneSets(
    species=args$species,
    library=args$geneset
)
print(
    paste(
        "Loaded", length(gene_sets), "genesets for", args$species,
        "from", args$geneset, "library"
    )
)

gsa_data <- enrichIt(                                                          # data frame of normalized (0,1) gene signatures scores
    obj=seurat_data,                                                           # always uses counts slot from the RNA assay
    gene.sets=gene_sets,
    method=args$method,
    groups=1000,                                                               # if method is UCell this parameter shouldn't influence on the results
    cores=args$cpus, 
    min.size=args$mingenes,
    ssGSEA.norm=TRUE,                                                          # used only when method is ssGSEA, needed to have results in (0,1) range
    maxRank=max(sapply(GSEABase::geneIds(gene_sets), function(x) length(x)))   # UCell fails if it's smaller than the longest geneset
)

print("Saving results as GSA assay")
seurat_data[["GSA"]] <- CreateAssayObject(
    counts=t(as.matrix(gsa_data)),                               # data slot will be equal to counts
    min.cells=1,                                                 # to remove gene sets that have all signatures scores equal to 0
    min.features=0                                               # to include all cells
)
debug$print_info(seurat_data, args)

print("Searching the most variable gene sets")
DefaultAssay(seurat_data) <- "GSA"
seurat_data <- FindVariableFeatures(
    seurat_data,                                                 # for "vst" counts slot is used
    selection.method="vst",
    nfeatures=args$topn,
    verbose=FALSE
)

selected_gene_sets <- VariableFeatures(seurat_data)
print(selected_gene_sets)

export_all_plots(
    seurat_data=seurat_data,
    gene_sets=selected_gene_sets,
    args=args
)

if(args$cbbuild){
    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        print("Exporting RNA, ATAC, and GSA assays to UCSC Cellbrowser jointly")
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
            assay="GSA",
            slot="counts",
            short_label="GSA",
            is_nested=TRUE,
            features=selected_gene_sets,
            rootname=paste(args$output, "_cellbrowser/gsa", sep=""),
        )
    } else {
        print("Exporting RNA and GSA assays to UCSC Cellbrowser jointly")
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
            assay="GSA",
            slot="counts",
            short_label="GSA",
            is_nested=TRUE,
            features=selected_gene_sets,
            rootname=paste(args$output, "_cellbrowser/gsa", sep=""),
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