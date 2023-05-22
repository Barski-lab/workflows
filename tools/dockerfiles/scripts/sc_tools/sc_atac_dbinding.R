#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
# suppressMessages(library(forcats))
suppressMessages(library(modules))
suppressMessages(library(argparse))
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))

set.seed(42)                                                                            # otherwise MAnorm returns different results


prepare_fragments_and_peaks <- function(seurat_data, seqinfo_data, args){               # only MACS2 results will be saved to the --output. All other files are in --tmpdir
    print(paste("Splitting ATAC fragments into two BED files by", args$splitby))        # we should always have only two groups because we prefiltered our Seurat object
    SplitFragments(                                                                     # will place BED files into args$tmpdir
        seurat_data,
        assay="ATAC",
        group.by=args$splitby,                                                          # this column should have only two values (--first, --second) by this moment
        outdir=args$tmpdir,
        verbose=TRUE
    )
    if (args$test != "manorm2"){
        print("Exporting original Seurat peaks")
        peaks_data <- seurat_data[["ATAC"]]@ranges
        seqlevels(peaks_data) <- seqlevels(seqinfo_data)
        seqinfo(peaks_data) <- seqinfo_data
        peaks_data$score <- 0
        export.bb(
            peaks_data,
            paste0(args$output, "_seurat_peaks.bigBed")
        )
    }
    tmp_locations <- list(
        fragments=list(
            first=paste0(                                                               # BED file with ATAC fragments for cells from --first group
                args$tmpdir, "/",
                gsub(.Platform$file.sep, "_", gsub(" ", "_", args$first)),              # similar to what SplitFragments does
                ".bed"),
            second=paste0(                                                              # BED file with ATAC fragments for cells from --second group
                args$tmpdir, "/",
                gsub(.Platform$file.sep, "_", gsub(" ", "_", args$second)),             # similar to what SplitFragments does
                ".bed")
        ),
        fragments_coverage=list(
            first=paste0(basename(args$output), "_", "first.bigWig"),                   # bigWig file with ATAC fragments coverage for cells from --first group
            second=paste0(basename(args$output), "_", "second.bigWig")                  # bigWig file with ATAC fragments coverage for cells from --second group
        ),
        macs2=list(
            cut_sites=list(                                                             # doesn't include strand information, should be treated as single-read
                first=paste0(args$tmpdir, "/", "first_tn5ct.bed"),                      # BED file with Tn5 cut sites for cells from --first group
                second=paste0(args$tmpdir, "/", "second_tn5ct.bed")                     # BED file with Tn5 cut sites for cells from --second group
            ),
            cut_sites_coverage=list(                                                    # made of *_treat_pileup.bdg output from MACS2
                first=paste0(basename(args$output), "_", "first_tn5ct.bigWig"),         # bigWig file with Tn5 cut sites coverage for cells from --first group
                second=paste0(basename(args$output), "_", "second_tn5ct.bigWig")        # bigWig file with Tn5 cut sites coverage for cells from --second group
            ),
            peaks_xls=list(                                                             # output from MACS2, will be also used as output from the tool
                first=paste0(basename(args$output), "_", "first_peaks.xls"),
                second=paste0(basename(args$output), "_", "second_peaks.xls")
            ),
            peaks_narrow=list(                                                          # output from MACS2, will be used by MAnorm2 and as output from the tool
                first=paste0(basename(args$output), "_", "first_peaks.narrowPeak"),
                second=paste0(basename(args$output), "_", "second_peaks.narrowPeak")
            ),
            peaks_summit=list(                                                          # output from MACS2 (no need to keep it as output from the tool
                first=paste0(basename(args$output), "_", "first_summits.bed"),          # as they will be used only by MAnorm2)
                second=paste0(basename(args$output), "_", "second_summits.bed")
            )
        )
    )
    for (i in 1:length(tmp_locations$fragments)) {
        print(paste("Loading fragments data from", tmp_locations$fragments[[i]]))
        fragments_data=rtracklayer::import(
            tmp_locations$fragments[[i]],
            format="BED",
            genome=seqinfo_data
        )
        print("Generating coverage from the loaded fragments data")
        io$export_fragments_coverage(
            fragments_data=fragments_data,
            location=tmp_locations$fragments_coverage[[i]]
        )
        if (args$test == "manorm2"){
            print("Extracting 1bp length Tn5 cut sites from the loaded fragments data")
            cut_sites_data <- unlist(as(list(
                resize(fragments_data, 1, fix="start", ignore.strand=TRUE),                    # will be 1bp region that corresponds to the beginning of the fragment
                resize(fragments_data, 1, fix="end", ignore.strand=TRUE)                       # will be 1bp region that corresponds to the end of the fragment
            ), "GRangesList"))
            rtracklayer::export.bed(
                cut_sites_data,
                tmp_locations$macs2$cut_sites[[i]],
                ignore.strand=TRUE                                                             # we need to ignore strand otherwise MACS2 will fail to parse it
            )
            rm(cut_sites_data)
            print(
                paste(
                    "Calling peaks with MACS2 from the extracted Tn5 cut sites",
                    "with the following parameters: --qvalue", args$qvalue,
                    "--shift -25 --extsize 50 --keep-dup all", "--gsize",
                    args$genome
                )
            )
            Signac::CallPeaks(                                                                 # we don't need to keep output as a variable
                object=tmp_locations$macs2$cut_sites[[i]],
                outdir=dirname(args$output),                                                   # we save output files directly to the output folder
                name=paste0(basename(args$output), "_", c("first", "second")[i]),        # all the names will have prefix with group name
                extsize=50,                                                                    # will be always called with --nomodel, so --shift
                shift=-25,                                                                     # and --extsize make difference
                effective.genome.size=args$genome,                                             # looks like it's ok to use character instead of an integer
                additional.args=paste(
                    "-q", args$qvalue,
                    "--keep-dup all",                                                          # our fragments are already deduplicated
                    "-B --SPMR",                                                               # to save signal per million reads as bedgraph file
                    "--seed 42",                                                               # just in case
                    "--tempdir", args$tmpdir
                ),
                cleanup=FALSE,                                                                 # we want to keep outputs after running the script
                verbose=FALSE
            )
            print("Generating coverage from the reported by MACS2 pileup data")
            io$export_coverage(
                coverage_data=rtracklayer::import(
                    paste0(basename(args$output), "_", c("first", "second")[i], "_treat_pileup.bdg"),
                    format="bedGraph",
                    genome=seqinfo_data
                ),
                location=tmp_locations$macs2$cut_sites_coverage[[i]]
            )
        } else if (!is.null(tmp_locations$macs2)) {
            tmp_locations$macs2 <- NULL                           # no need to have macs2 slot if we are not going to use it
        }
        rm(fragments_data)
        gc(verbose=FALSE)
    }
    return (tmp_locations)
}

export_raw_plots <- function(seurat_data, args){
    DefaultAssay(seurat_data) <- "ATAC"                           # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure
    for (reduction in c("rnaumap", "atacumap", "wnnumap")) {
        if (!(reduction %in% names(seurat_data@reductions))) {next}                                  # skip missing reductions
        graphics$dim_plot(
            data=seurat_data,
            reduction=reduction,
            plot_title=paste0(
                "Cells UMAP split by ", args$splitby, " ",
                ifelse(
                    (!is.null(args$groupby) && !is.null(args$subset)),
                    paste("subsetted to", paste(args$subset, collapse=", "), "values from the", args$groupby, "column "),
                    ""
                ),
                "(", reduction, " dim. reduction)"
            ),
            legend_title=ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                args$groupby,
                args$splitby
            ),
            split_by=args$splitby,
            group_by=ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                args$groupby,
                args$splitby
            ),
            highlight_group=if(!is.null(args$groupby) && !is.null(args$subset)) args$subset else NULL,
            label=FALSE,
            label_color="black",
            palette_colors=if(!is.null(args$groupby) && !is.null(args$subset)) c("#4CB6BA", "#FB1C0D") else c("darkred", "darkorchid4"),
            theme=args$theme,
            rootname=paste(args$output, "umap_rd", reduction, sep="_"),
            pdf=args$pdf
        )
    }
    gc(verbose=FALSE)
}

export_browser_tracks <- function(db_results, seqinfo_data, args){
    print(
        paste(
            "Filtering differentially bound sites to include",
            "only regions with adjusted P-values not bigger",
            "than", args$padj
        )
    )
    db_sites <- db_results$db_sites %>%
                dplyr::filter(.$padj <= args$padj) %>%
                dplyr::mutate(
                    "name"=paste0(
                        "padj=", format(padj, digits=3, trim=TRUE),
                        ";log2FC=", format(log2FoldChange, digits=3, trim=TRUE)
                    )
                ) %>%
                dplyr::mutate("score"=0) %>%                                     # need it because export.bb fails without score
                dplyr::select(
                    c("chr", "start", "end", "name", "score", "log2FoldChange")  # we keep log2FoldChange for filtering
                )
    print(head(db_sites))

    if (nrow(db_sites) > 0){
        db_sites <- makeGRangesFromDataFrame(
            db_sites,
            seqinfo=seqinfo_data,
            keep.extra.columns=TRUE                                              # to keep name and score metadata columns
        )
        first_db_sites <- db_sites[db_sites$log2FoldChange <= -args$logfc, ]
        second_db_sites <- db_sites[db_sites$log2FoldChange >= args$logfc, ]
        if (length(first_db_sites) > 0){                                         # we use length because it GRanges
            export.bb(
                first_db_sites,
                paste0(args$output, "_first_enrch.bigBed")
            )
            export.bed(
                first_db_sites,
                paste0(args$output, "_first_enrch.bed"),
                ignore.strand=TRUE
            )
        }
        if (length(second_db_sites) > 0){                                        # we use length because it GRanges
            export.bb(
                second_db_sites,
                paste0(args$output, "_second_enrch.bigBed")
            )
            export.bed(
                second_db_sites,
                paste0(args$output, "_second_enrch.bed"),
                ignore.strand=TRUE
            )
        }
    } else (
        print("No significant differentially bound sites found")
    )
}

export_processed_plots <- function(seurat_data, db_results, args){
    DefaultAssay(seurat_data) <- "ATAC"                                 # safety measure
    Idents(seurat_data) <- "new.ident"                                  # safety measure

    graphics$volcano_plot(
        data=db_results$db_sites,                                       # this is not filtered differentially bound sites
        x_axis="log2FoldChange",
        y_axis="padj",
        x_cutoff=args$logfc,
        y_cutoff=args$padj,
        x_label="log2FC",
        y_label="-log10Padj",
        label_column="chr",                                             # doesn't matter what to set as label because features=NA
        plot_title="Differentially bound sites",
        plot_subtitle=paste0(
            args$second, " vs ", args$first, " for cells split by ", args$splitby, ". ",
            "Use ", args$test, " test. ",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste("Subsetted to", paste(args$subset, collapse=", "), "values from", args$groupby, "column. "),
                ""
            ),
            "Displayed thresholds: Padj <= ", args$padj,
            ", |log2FC| >= ", args$logfc
        ),
        caption=paste(nrow(db_results$db_sites), "peaks"),
        features=NA,                                                    # we don't want to show peak names as they are too long
        theme=args$theme,
        rootname=paste(args$output, "dbnd_vlcn", sep="_"),
        pdf=args$pdf
    )
}

get_args <- function(){
    parser <- ArgumentParser(
        description="Single-cell ATAC-Seq Differential Binding Analysis"
    )
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "chromatin accessibility information stored in the ATAC assay. Additionally",
            "'rnaumap', and/or 'atacumap', and/or 'wnnumap' dimensionality reductions",
            "should be present."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment used in the loaded",
            "Seurat object. File should be saved in TSV format with tbi-index file."
        ),
        type="character", required="True"
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
            "metadata by selected barcodes. First column should be named as 'barcode'.",
            "If file includes any other columns they will be added to the Seurat object",
            "metadata ovewriting the existing ones if those are present.",
            "Default: all cells used, no extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Column from the Seurat object metadata to group cells for optional subsetting",
            "when combined with --subset parameter. May be one of the extra metadata columns",
            "added with --metadata or --barcodes parameters. Ignored if --subset is not set.",
            "Default: do not subset, include all cells into analysis."
        ),
        type="character"
    )
    parser$add_argument(
        "--subset",
        help=paste(
            "Values from the column set with --groupby parameter to subset cells",
            "before running differential binding analysis. Ignored if --groupby",
            "is not provided. Default: do not subset cells, include all of them."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to split cells into two groups to",
            "run --second vs --first differential binding analysis. May be one of",
            "the extra metadata columns added with --metadata or --barcodes parameters."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--first",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby parameter",
            "to define the first group of cells for differential binding analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--second",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby parameter",
            "to define the second group of cells for differential binding analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--test",
        help=paste(
            "Test type to use in differential binding analysis. For all tests except",
            "manorm2, peaks present in the loaded Seurat object will be used. If manorm2",
            "test selected, peaks will be called per group defined by --splitby parameter.",
            "Default: logistic-regression"
        ),
        type="character", default="logistic-regression",
        choices=c(
            "negative-binomial",          # (negbinom) Negative Binomial Generalized Linear Model (use FindMarkers with peaks from Seurat object)
            "poisson",                    # (poisson) Poisson Generalized Linear Model (use FindMarkers with peaks from Seurat object)
            "logistic-regression",        # (LR) Logistic Regression (use FindMarkers with peaks from Seurat object)
            "mast",                       # (MAST) MAST package (use FindMarkers with peaks from Seurat object)
            "manorm2"                      # call peaks for each group with MACS2, run MAnorm2
        )
    )
    parser$add_argument(
        "--genome",
        help=paste(
            "Genome type of the sequencing data loaded from the Seurat",
            "object. It will be used for effective genome size selection",
            "when calling peaks with MACS2. Ignored if --test is not set",
            "to manorm2. Default: hs (2.7e9)"
        ),
        type="character", default="hs",
        choices=c(
            "hs",        # 2.7e9
            "mm"         # 1.87e9
        )
    )
    parser$add_argument(
        "--qvalue",
        help=paste(
            "Minimum FDR (q-value) cutoff for MACS2 peak detection. Ignored",
            "if --test is not set to manorm2. Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--minpeakgap",
        help=paste(
            "If a distance between peaks is smaller than the provided value",
            "they will be merged before splitting them into reference genomic",
            "bins of size --binsize. Ignored if --test is not set to manorm2.",
            "Default: 150"
        ),
        type="integer", default=150
    )
    parser$add_argument(
        "--binsize",
        help=paste(
            "The size of non-overlapping reference genomic bins used by",
            "MAnorm2 when generating a table of reads counts per peaks.",
            "Ignored if --test is not set to manorm2. Default: 1000"
        ),
        type="integer", default=1000
    )
    parser$add_argument(
        "--maxpeaks",
        help=paste(
            "The maximum number of the most significant (based on qvalue)",
            "peaks to keep from each group of cells when constructing",
            "reference genomic bins. Ignored if --test is not set to",
            "manorm2. Default: keep all peaks"
        ),
        type="integer"
    )   
    parser$add_argument(
        "--blacklist",
        help=paste(
            "Path to the optional BED file with the genomic blacklist regions",
            "to be filtered out before running differential binding analysis.",
            "Any reference genomic bin overlapping a blacklist region will be",
            "removed from the output. Ignored if --test is not set to manorm2."
        ),
        type="character"
    )
    parser$add_argument(
        "--padj",
        help=paste(
            "In the exploratory visualization part of the analysis output only",
            "differentially bound peaks with adjusted P-value not bigger than",
            "this value. Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--logfc",
        help=paste(
            "In the exploratory visualization part of the analysis output only",
            "differentially bound peaks with log2 Fold Change not smaller than",
            "this value. Default: 1.0"
        ),
        type="double", default=1
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

print("Adjusting --test parameter")
args$test <- switch(
    args$test,
    "negative-binomial"   = "negbinom",
    "poisson"             = "poisson",
    "logistic-regression" = "LR",
    "mast"                = "MAST",
    "manorm2"             = "manorm2"
)
print(paste("--test was adjusted to", args$test))

print(
    paste(
        "Setting parallelization to", args$cpus, "cores, and", args$memory,
        "GB of memory allowed to be shared between the processes"
    )
)
prod$parallel(args)

print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
debug$print_info(seurat_data, args)

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

if (!("ATAC" %in% names(seurat_data@assays))){
    print(
        paste(
            "Loaded Seurat object doesn't include required ATAC assay.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

print("Setting default assay to ATAC")
DefaultAssay(seurat_data) <- "ATAC"

seqinfo_data <- seurat_data[["ATAC"]]@seqinfo
if (is.null(seqinfo_data)){
    print("Loaded Seurat object doesn't include seqinfo data. Exiting.")
    quit(save="no", status=1, runLast=FALSE)
}

print(paste("Loading fragments data from", args$fragments))
seurat_data <- io$replace_fragments(args$fragments, seurat_data)
debug$print_info(seurat_data, args)

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
    debug$print_info(seurat_data, args)
}

print("Subsetting Seurat object to include only cells from the tested conditions")
seurat_data <- io$apply_metadata_filters(seurat_data, args$splitby, c(args$first, args$second))
debug$print_info(seurat_data, args)

export_raw_plots(seurat_data, args)                                                        # will include only cells from --first and --second

if(!is.null(args$groupby) && !is.null(args$subset)){
    print("Subsetting Seurat object to include only selected groups of cells")
    seurat_data <- io$apply_metadata_filters(seurat_data, args$groupby, args$subset)
    debug$print_info(seurat_data, args)
}

if (!all(table(seurat_data@meta.data[[args$splitby]]) > 0)){                               # check if we accidentally removed cells we want to compare
    print(
        paste(
            "Not enough cells for comparison. Check --groupby",
            "and --subset parameters. Exiting."
        )
    )
    print(table(seurat_data@meta.data[[args$splitby]]))
    quit(save="no", status=1, runLast=FALSE)
}

if (args$test != "manorm2"){                                                               # needed only for FindMarkers (not for MAnorm2)
    print("Normalizing ATAC counts after all filters applied")
    seurat_data <- Signac::RunTFIDF(                                                       # might be redundant as it may not depend on the number of cells
        seurat_data,
        assay="ATAC",
        method=1,                                                                          # using default for RunTFIDF "log-tfidf" method
        verbose=FALSE
    )
}

args$tmp_locations <- prepare_fragments_and_peaks(                                         # will create ATAC fragments and Tn5 cut sites BED files
    seurat_data=seurat_data,                                                               # for --first and --second groups of cells. If --test was
    seqinfo_data=seqinfo_data,                                                             # set to manorm2, will call peaks from Tn5 cut sites.
    args=args
)

db_results <- analyses$atac_dbinding_analyze(seurat_data, args)
print(head(db_results$db_sites))

print("Exporting not filtered differentially bound sites")
io$export_data(
    db_results$db_sites,                                                                      # not filtered na-removed differentially bound sites
    paste(args$output, "_db_sites.tsv", sep="")
)
export_processed_plots(seurat_data, db_results, args)

export_browser_tracks(                                                                           # exporting fragments coverage into bigWigs for visualization purposes
    db_results=db_results,
    seqinfo_data=seqinfo_data,
    args
)
