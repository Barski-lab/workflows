#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(knitr))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(forcats))
suppressMessages(library(stringr))
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
suppressMessages(logger <- modules::use(file.path(HERE, "modules/logger.R")))


## ----
get_ordered_metadata <- function(seurat_data, args){                                # reorder metadata only for better plots
    ordered_metadata <- seurat_data@meta.data %>%
                        dplyr::select(
                            new.ident,                                              # this column is always present
                            tidyselect::all_of(args$splitby)                        # should never be NULL, safe when set to new.ident
                        ) %>%
                        dplyr::mutate(
                            !!args$splitby:=base::factor(
                                .[[args$splitby]],
                                levels=c(args$first, args$second)
                            )
                        ) %>%
                        dplyr::arrange_at(c(args$splitby, "new.ident")) %>%         # safe when splitby is new.ident
                        dplyr::mutate(                                              # if new.ident is already a factor arrange_at will sort by factor
                            new.ident=forcats::fct_inorder(new.ident)
                        )
    return (ordered_metadata)
}

## ----
prepare_fragments_and_peaks <- function(seurat_data, seqinfo_data, args){
    ordered_metadata <- get_ordered_metadata(seurat_data, args)                     # we need it for a proper order of items shown on the plots
    aggr_criteria <- ifelse(
        args$test %in% c("manorm2-full", "manorm2-half"),
        "new.ident",
        args$splitby
    )
    aggr_names <- levels(ordered_metadata[[aggr_criteria]])                         # both "new.ident" and args$splitby are already factors
    aggr_conditions <- as.vector(                                                   # we don't want to have factors here
        as.character(
            ordered_metadata[[args$splitby]][                                       # safe to use the first matched value even for new.ident because we
                match(aggr_names, ordered_metadata[[aggr_criteria]])                # checked that --splitby doesn't divide cells from the same datasets
            ]
        )
    )
    rm(ordered_metadata)                                                            # no reason to keep ordered_metadata anymore

    mapped_counts <- AverageCounts(                                                 # need it for bigWig scaling
        seurat_data, assay="ATAC", group.by=aggr_criteria, verbose=FALSE
    ) * CellsPerGroup(
        seurat_data, group.by=aggr_criteria
    )

    metadata <- NULL                                                                # a dataframe to collect files locations
    for (i in 1:length(aggr_names)){
        current_name <- aggr_names[i]                                               # this is used only for labels
        current_condition <- aggr_conditions[i]                                     # this is used both for labels and file names
        current_suffix <- paste(                                                    # this is used only for file names
            tolower(
                gsub("'|\"|\\s|\\t|#|%|&|-", "_", current_name)
            ),
            ifelse(current_condition == args$first, "first", "second"),            # here we always have either args$first or args$seconds values
            sep="_"
        )
        current_data <- data.frame(
            name=current_name,
            suffix=current_suffix,
            condition=current_condition,

            fragments=paste0(                                                       # temporary files, we use it to search for files
                args$tmpdir, "/",
                gsub(.Platform$file.sep, "_",
                gsub(" ", "_", current_name)),
                ".bed"
            ),
            tn5ct=paste0(args$tmpdir, "/", current_suffix, "_tn5ct.bed"),           # temporary files, we use it to create files
            coverage=paste0(args$output, "_", current_suffix, ".bigWig"),           # will be saved as output, we use it to create files
            peaks=paste0(args$output, "_", current_suffix, "_peaks.narrowPeak"),    # will be saved as output, we use it to search for files
            summits=paste0(args$output, "_", current_suffix, "_summits.bed"),       # will be saved as output, we use it to search for files
            xls=paste0(args$output, "_", current_suffix, "_peaks.xls"),             # will be saved as output, we use it to search for files

            read_cnt=paste(current_suffix, "read_cnt", sep="."),                    # will be used by MAnorm for normalization
            occupancy=paste(current_suffix, "occupancy", sep="."),                  # will be used by MAnorm for normalization
            scaling=10^6/mapped_counts[current_name],                               # for scaling coverage
            stringsAsFactors=FALSE,
            check.names=FALSE,
            row.names=i                                                             # otherwise it puts the current_name as the rowname
        )

        if (is.null(metadata)) {
            metadata <- current_data
        } else {
            metadata <- base::rbind(metadata, current_data)
        }
    }

    print(
        paste(
            "Splitting ATAC fragments data into the",
            "separate BED files by", aggr_criteria,
            "metadata column"
        )
    )
    SplitFragments(                                                                 # will put files to locations defined in fragments
        seurat_data,
        assay="ATAC",
        group.by=aggr_criteria,
        outdir=args$tmpdir,
        verbose=TRUE
    )

    for (i in 1:nrow(metadata)){                                                    # iterating over the collected locations
        current_row <- as.list(metadata[i, ])                                       # to have it as list instead of a data.frame with the single row
        print(paste("Processing", current_row$name, "dataset"))
        print(
            paste(
                "Loading ATAC fragments data",
                "from", current_row$fragments
            )
        )
        fragments_data=rtracklayer::import(
            current_row$fragments,
            format="BED",
            genome=seqinfo_data
        )

        io$export_fragments_coverage(
            fragments_data=fragments_data,
            location=current_row$coverage,
            scaling_coef=current_row$scaling
        )

        if (args$test %in% c("manorm2-full", "manorm2-half")){                      # for MAnorm2 we will need to have Tn5 cut sites for later calling with MACS2
            print(
                paste(
                    "Extracting Tn5 cut sites (1bp length)",
                    "from the loaded ATAC fragments data"
                )
            )
            tn5ct_data <- unlist(as(list(
                resize(fragments_data, 1, fix="start", ignore.strand=TRUE),         # will be 1bp region that corresponds to the beginning of the fragment
                resize(fragments_data, 1, fix="end", ignore.strand=TRUE)            # will be 1bp region that corresponds to the end of the fragment
            ), "GRangesList"))

            rtracklayer::export.bed(
                tn5ct_data,
                current_row$tn5ct,
                ignore.strand=TRUE                                                  # we need to ignore strand otherwise MACS2 will fail to parse it
            )
            rm(tn5ct_data)                                                          # no reason to keep it
        }
        rm(fragments_data)
        gc(verbose=FALSE)
    }

    # we call peaks in a separate for loop because we need to have all
    # tn5ct files already created in case the --test was set to manorm2-half

    if (args$test %in% c("manorm2-full", "manorm2-half")){                          # for MAnorm2 we need to call peaks with MACS2
        for (i in 1:nrow(metadata)){                                                # iterating over the collected locations
            current_row <- as.list(metadata[i, ])                                   # to have it as list instead of a data.frame with the single row
            tn5ct_location <- current_row$tn5ct                                     # default tn5ct location from the current dataset
            if (args$test == "manorm2-half"){
                tn5ct_location <- as.vector(                                        # all tn5ct locations that belong to the current comparison group
                    metadata[                                                       # they are all shoudl have been already created in the previous for loop
                        metadata$condition == current_row$condition,
                        "tn5ct"
                    ]
                )
            }
            tryCatch(
                expr = {
                    print(
                        paste0(
                            "Calling peaks with MACS2 for ", current_row$name,
                            " dataset from the extracted Tn5 cut sites (",
                            paste(basename(tn5ct_location), collapse=", "), ") ",
                            "with the following parameters: --qvalue ", args$qvalue,
                            " --shift -25 --extsize 50 --keep-dup all --gsize ",
                            args$genome
                        )
                    )
                    Signac::CallPeaks(                                                      # we don't need to keep output as a variable
                        object=tn5ct_location,
                        outdir=dirname(args$output),                                        # we save output files directly to the output folder
                        name=paste0(basename(args$output), "_", current_row$suffix),        # set the name so we can use our peaks, summits, and xls locations
                        extsize=50,                                                         # will be always called with --nomodel, so --shift
                        shift=-25,                                                          # and --extsize make difference
                        effective.genome.size=args$genome,
                        additional.args=paste(
                            "-q", args$qvalue,
                            "--keep-dup all",                                               # our fragments are already deduplicated
                            "--seed", args$seed,                                            # just in case
                            "--tempdir", args$tmpdir
                        ),
                        cleanup=FALSE,                                                      # we want to keep outputs after running the script
                        verbose=FALSE
                    )
                },
                error = function(e){
                    logger$info(
                        paste0(
                            "Failed to call peaks with MACS2 for ", current_row$name,
                            " dataset from the extracted Tn5 cut sites (",
                            paste(basename(tn5ct_location), collapse=", "), "). ",
                            "Exiting. ", e
                        )
                    )
                    quit(save="no", status=1, runLast=FALSE)
                }
            )
        }
    }

    print("Exporting original peaks from the loaded Seurat object")
    peaks_data <- seurat_data[["ATAC"]]@ranges
    seqlevels(peaks_data) <- seqlevels(seqinfo_data)
    seqinfo(peaks_data) <- seqinfo_data
    peaks_data$score <- 0                                                           # need it to export in bigBed format, otherwise fails
    export.bb(
        peaks_data,
        paste0(args$output, "_dflt_peaks.bigBed")                                   # only for visualization purposes
    )

    return (metadata)
}

## ----
export_raw_plots <- function(seurat_data, args){
    DefaultAssay(seurat_data) <- "ATAC"                                           # safety measure
    Idents(seurat_data) <- "new.ident"                                            # safety measure

    graphics$dim_plot(
        data=seurat_data,
        reduction=args$reduction,
        plot_title="UMAP colored by selected for analysis cells",
        plot_subtitle=paste0(
            "Split by ", args$splitby, "; ",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste(
                    "only cells with",
                    paste(args$subset, collapse=", "),
                    "values in the", args$groupby,
                    "metadata column are selected."
                ),
                "all cells selected."
            )
        ),
        legend_title="Cells",
        split_by=args$splitby,
        group_by=ifelse(                                                          # highlight_group won't work without it
            (!is.null(args$groupby) && !is.null(args$subset)),
            args$groupby,
            args$splitby
        ),
        highlight_group=if(!is.null(args$groupby) && !is.null(args$subset))
                            args$subset
                        else
                            c(args$first, args$second),                           # to highlight all cells
        palette_colors=if(!is.null(args$groupby) && !is.null(args$subset))
                           c(graphics$NA_COLOR, graphics$HIGHLIGHT_COLOR)
                       else
                           c(graphics$HIGHLIGHT_COLOR, graphics$NA_COLOR),        # need to reverse colors
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_tst", sep="_"),
        pdf=args$pdf
    )
}

## ----
export_processed_plots <- function(db_results, seqinfo_data, args){
    DefaultAssay(seurat_data) <- "ATAC"                                           # safety measure
    Idents(seurat_data) <- "new.ident"                                            # safety measure

    ordered_metadata <- get_ordered_metadata(seurat_data, args)
    graphics$geom_bar_plot(
        data=ordered_metadata,
        x_axis=ifelse(
            args$test %in% c("manorm2-full", "manorm2-half"),
            "new.ident",
            args$splitby
        ),
        color_by=args$splitby,
        x_label=ifelse(
            args$test %in% c("manorm2-full", "manorm2-half"),
            "Dataset",
            args$splitby
        ),
        y_label="Cell counts",
        legend_title="Tested\ncondition",
        plot_title=ifelse(
            args$test %in% c("manorm2-full", "manorm2-half"),
            "Number of cells per dataset",
            paste(
                "Number of cells per tested",
                "condition defined by", args$splitby,
                "metadata column"
            )
        ),
        plot_subtitle=paste0(
            "Colored by ", args$splitby, "; ",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste(
                    "only cells with",
                    paste(args$subset, collapse=", "),
                    "values in the", args$groupby,
                    "metadata column are selected."
                ),
                "all cells selected."
            )
        ),
        palette_colors=c(graphics$TRUE_COLOR, graphics$FALSE_COLOR),
        theme=args$theme,
        rootname=paste(args$output, "cell_cnts", sep="_"),
        pdf=args$pdf
    )

    graphics$volcano_plot(
        data=db_results$db_sites,                                                 # this is not filtered differentially bound sites
        x_axis="log2FoldChange",
        y_axis="padj",
        x_cutoff=args$logfc,
        y_cutoff=args$padj,
        x_label="log2 FC",
        y_label="-log10 Padj",
        label_column="chr",                                                       # doesn't matter what to set as label because features=NA
        plot_title="Differentially accessible regions",
        plot_subtitle=paste0(
            args$second, " vs ", args$first,
            " for cells split by ", args$splitby, ". ",
            "Use ", args$test, " test. ",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste(
                    "Subsetted to", paste(args$subset, collapse=", "),
                    "values from", args$groupby, "column. "),
                ""
            ),
            "Displayed thresholds: Padj <= ", args$padj,
            ", |log2FC| >= ", args$logfc
        ),
        caption=paste(nrow(db_results$db_sites), "peaks"),
        features=NA,                                                              # we don't want to show peak names as they are too long
        theme=args$theme,
        rootname=paste(args$output, "vlcn", sep="_"),
        pdf=args$pdf
    )

    all_db_sites <- db_results$db_sites %>%                                           # all not filtered diff. acc. regions
                    dplyr::mutate(
                        "name"=paste0(
                            "padj=", format(padj, digits=3, trim=TRUE),
                            ";log2FC=", format(log2FoldChange, digits=3, trim=TRUE)
                        ),
                        "itemRgb"=dplyr::case_when(
                                    log2FoldChange >= 0 ~ graphics$FALSE_COLOR,
                                    .default = graphics$TRUE_COLOR
                                )
                    ) %>%
                    dplyr::mutate("score"=-log10(padj)*10) %>%                                             # similar to what MACS2 does
                    dplyr::select(
                        c("chr", "start", "end", "name", "score", "itemRgb", "log2FoldChange", "padj")     # we keep log2FoldChange and padj for row_metadata in morpheus heatmap
                    )
    all_db_ranges <- makeGRangesFromDataFrame(
        all_db_sites,
        seqinfo=seqinfo_data,
        keep.extra.columns=TRUE
    )
    all_db_location <- paste0(args$output, "_all_db_sites.bed")                                                # should be kept as output, the same as all_db_sites.tsv
    export.bed(all_db_ranges, all_db_location, ignore.strand=TRUE)
    writeLines(
        c(
            paste(
                "track name=\"All diff. acc. regions\"",
                "description=\"All diff. acc. regions\"",
                "itemRgb=\"On\""),
            readLines(all_db_location)
        ),
        all_db_location
    )

    print(
        paste(
            "Filtering differentially accessible regions to",
            "include only regions with adjusted p-values",
            "not bigger than", args$padj, "and |log2FoldChange|",
            "bigger or equal to", args$logfc
        )
    )
    filtered_db_sites <- all_db_sites %>%
                         dplyr::filter(.$padj <= args$padj) %>%                                      # to include only significant diff. accessible regions
                         dplyr::filter(abs(.$log2FoldChange) >= args$logfc)                          # to include only |log2FoldChange| >= args$logfc
    print(head(filtered_db_sites))

    if (nrow(filtered_db_sites) > 0){
        filtered_db_location <- paste0(args$output, "_fltr_db_sites.bed")                            # should be kept as output
        first_db_location <- paste0(args$tmpdir, "/", "first_enrch.bed")                             # temporary file, needed only for deeptools
        second_db_location <- paste0(args$tmpdir, "/", "second_enrch.bed")                           # temporary file, needed only for deeptools
        score_matrix_location <- paste0(args$tmpdir, "/", "score_matrix.gz")                         # no reason to keep it

        filtered_db_ranges <- makeGRangesFromDataFrame(                                              # fails when filtered_db_sites is empty
            filtered_db_sites,
            seqinfo=seqinfo_data,
            keep.extra.columns=TRUE                                                                  # to keep the log2FoldChange column
        )
        export.bed(filtered_db_ranges, filtered_db_location, ignore.strand=TRUE)
        writeLines(
            c(
                paste(
                    "track name=\"Filtered diff. acc. regions\"",
                    "description=\"Filtered diff. acc. regions\"",
                    "itemRgb=\"On\""),
                readLines(filtered_db_location)
            ),
            filtered_db_location
        )

        first_db_ranges <- filtered_db_ranges[filtered_db_ranges$log2FoldChange <= -args$logfc, ]    # at least one of these groups won't be empty
        second_db_ranges <- filtered_db_ranges[filtered_db_ranges$log2FoldChange >= args$logfc, ]    # because we checked nrow(filtered_db_sites) > 0)

        regions_locations <- c()
        regions_labels <- c()
        if (length(first_db_ranges) > 0){
            regions_locations <- c(first_db_location, regions_locations)
            regions_labels <- c(args$first, regions_labels)
            export.bed(
                first_db_ranges[order(mcols(first_db_ranges)$log2FoldChange, decreasing=FALSE)],
                first_db_location,
                ignore.strand=TRUE
            )
        }
        if (length(second_db_ranges) > 0){
            regions_locations <- c(second_db_location, regions_locations)
            regions_labels <- c(args$second, regions_labels)
            export.bed(
                second_db_ranges[order(mcols(second_db_ranges)$log2FoldChange, decreasing=TRUE)],
                second_db_location,
                ignore.strand=TRUE
            )
        }

        compute_matrix_args <- c(
            "reference-point",
            "--scoreFileName", args$metadata$coverage,
            "--regionsFileName", regions_locations,
            "--referencePoint", "center",
            "--beforeRegionStartLength", 5000,
            "--afterRegionStartLength", 5000,
            "--binSize", 100,
            "--sortRegions", "descend",                                         # sort by mean of each row in a descending order
            "--sortUsing", "mean",
            "--samplesLabel", args$metadata$name,
            "--outFileName", score_matrix_location,
            "--missingDataAsZero",
            "--numberOfProcessors", args$cpus
        )
        exit_code <- sys::exec_wait(
            cmd="computeMatrix",                                                # if it's not found in PATH, R will fail with error
            args=compute_matrix_args
        )
        if (exit_code != 0){                                                    # we were able to run profile_bins, but something went wrong
            print(                                                              # no reason to print through logger as the workflow shoudn't fail
                paste(
                    "Failed to run computeMatrix command",
                    "with the exit code", exit_code
                )
            )
        } else {
            plot_heatmap_args <- c(
                "--plotTitle", "Tag density around the centers of differentially accessible regions, sorted in descending order by the mean value of each region",
                "--matrixFile", score_matrix_location,
                "--outFileName", paste0(args$output, "_", "tag_dnst_htmp.png"),
                "--plotType", "lines",
                "--sortRegions", "keep",                                        # the regions should be already sorted in the computeMatrix step
                "--averageTypeSummaryPlot", "mean",
                "--whatToShow", "plot, heatmap and colorbar",
                "--refPointLabel", "Center",
                "--regionsLabel", regions_labels,
                "--xAxisLabel", "(bp)",
                "--yAxisLabel", "Signal mean",
                "--plotFileFormat", "png",
                "--legendLocation", "upper-left",
                "--boxAroundHeatmaps", "no",
                "--zMin", "auto",                                               # 1st percentile
                "--zMax", "auto",                                               # 98th percentile
                "--colorList", "white,darkblue"
            )
            exit_code <- sys::exec_wait(
                cmd="plotHeatmap",                                              # if it's not found in PATH, R will fail with error
                args=plot_heatmap_args
            )
            if (exit_code != 0){                                                # we were able to run profile_bins, but something went wrong
                print(                                                          # no reason to print through logger as the workflow shoudn't fail
                    paste(
                        "Failed to run plotHeatmap command",
                        "with the exit code", exit_code
                    )
                )
            }

            print("Exporting Morpheus heatmap")
            score_matrix_data <- read.table(
                score_matrix_location,
                sep="\t",
                skip=1,                                                         # we need to skip the first row because it includes metadata
                header=FALSE,
                check.names=FALSE,
                stringsAsFactors=FALSE,
                quote=""                                                        # safety measure
            ) %>%
            dplyr::mutate(
                peak=base::paste0(V1, ":", V2+1, "-", V3)
            ) %>%
            dplyr::select(
                -c("V1", "V2", "V3", "V4", "V5", "V6")
            ) %>%
            tibble::remove_rownames() %>%
            tibble::column_to_rownames("peak")

            row_metadata <- filtered_db_sites %>%                               # always not empty dataframe
                            dplyr::mutate(
                                "peak"=paste0(chr, ":", start, "-", end)        # to make it correspond to score_matrix_data
                            ) %>%
                            tibble::remove_rownames() %>%
                            tibble::column_to_rownames("peak") %>%
                            dplyr::select("log2FoldChange", "padj") %>%         # we need only these two columns from filtered_db_sites
                            dplyr::mutate(                                      # we need this column only for grouping in Morpheus
                                "regions_group"=dplyr::case_when(
                                                    log2FoldChange >= args$logfc ~ "open",
                                                    log2FoldChange <= -args$logfc ~ "closed",
                                                    .default = "other"         # this group should never exist in the filtered_db_sites, but we add it just in case to avoid NA
                                                )
                            )
            row_metadata <- row_metadata[rownames(score_matrix_data), ]         # we want to have the rows order from deeptools results, should never fail

            col_metadata <- NULL
            n_times <- ncol(score_matrix_data) / length(args$metadata$name)     # shouldn't fail, as it's alsways dividable by the number of coverage files
            for (i in 1:length(args$metadata$name)){
                current_name <- args$metadata$name[i]
                current_condition <- args$metadata$condition[i]
                if (args$test %in% c("manorm2-full", "manorm2-half")){
                    current_col_metadata <- setNames(
                        data.frame(
                            list(
                                rep(current_condition, n_times),
                                rep(current_name, n_times)
                            )
                        ),
                        c(args$splitby, "new.ident")
                    )
                } else {
                    current_col_metadata <- setNames(
                        data.frame(
                            list(
                                rep(current_condition, n_times)
                            )
                        ),
                        c(args$splitby)
                    )
                }
                if (is.null(col_metadata)){
                    col_metadata <- current_col_metadata
                } else {
                    col_metadata <- rbind(col_metadata, current_col_metadata)
                }
            }
            rownames(col_metadata) <- colnames(score_matrix_data)               # score_matrix_data is already sorted by the values from provided as --samplesLabel

            score_matrix_mat <- as.matrix(score_matrix_data)
            score_limits <- stats::quantile(                                    # to exclude outliers
                score_matrix_mat,                                               # no reason to take the absolute value as everything in the score_matrix_mat >= 0
                c(0.01, 0.98),                                                  # will return unnamed vector with 2 values for 1st and 98th percentile
                na.rm=TRUE, names=FALSE
            )

            io$export_gct(
                counts_mat=score_matrix_mat,
                row_metadata=row_metadata,
                col_metadata=col_metadata,
                location=paste0(args$output, "_", "tag_dnst_htmp.gct")
            )
            graphics$morpheus_html_heatmap(
                gct_location=paste0(args$output, "_", "tag_dnst_htmp.gct"),
                rootname=paste0(args$output, "_tag_dnst_htmp"),
                color_scheme=list(
                    scalingMode="fixed",
                    stepped=FALSE,
                    values=as.list(score_limits),
                    colors=c("white", "darkblue")
                )
            )
            io$export_data(
                row_metadata %>% tibble::rownames_to_column(var="feature"),
                paste(args$output, "tag_dnst_htmp.tsv", sep="_")
            )
        }
    }
}

## ----
get_args <- function(){
    parser <- ArgumentParser(
        description="Single-Cell ATAC-Seq Differential Accessibility Analysis"
    )
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should",
            "include chromatin accessibility information stored in the ATAC assay.",
            "The dimensionality reductions selected in the --reduction parameter",
            "should be present in the loaded Seurat object."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--reduction",
        help=paste(
            "Dimensionality reduction to be used for generating UMAP plots."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment used in",
            "the loaded Seurat object. File should be saved in TSV format",
            "with tbi-index file."
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
            "before running differential accessibility analysis. Ignored if --groupby",
            "is not provided. Default: do not subset cells, include all of them."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to split cells into two groups to",
            "run --second vs --first differential accessibility analysis. If --test",
            "parameter is set to manorm2-full or manorm2-half, the --splitby shouldn't",
            "put cells from the same dataset into the different comparison groups.",
            "May be one of the extra metadata columns added with --metadata or",
            "--barcodes parameters."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--first",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby parameter",
            "to define the first group of cells for differential accessibility analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--second",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby parameter",
            "to define the second group of cells for differential accessibility analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--test",
        help=paste(
            "Test type to use in the differential accessibility analysis. For all tests",
            "except manorm2-full and manorm2-half, peaks already present in the loaded",
            "Seurat object will be used. If manorm2-full or manorm2-half test is selected,",
            "reads will be aggregated to pseudo bulk form either by dataset or comparison",
            "group and then peaks will be called with MACS2 per dataset.",
            "Default: logistic-regression"
        ),
        type="character", default="logistic-regression",
        choices=c(
            "negative-binomial",          # (negbinom) Negative Binomial Generalized Linear Model (use FindMarkers with peaks from Seurat object)
            "poisson",                    # (poisson) Poisson Generalized Linear Model (use FindMarkers with peaks from Seurat object)
            "logistic-regression",        # (LR) Logistic Regression (use FindMarkers with peaks from Seurat object)
            "mast",                       # (MAST) MAST package (use FindMarkers with peaks from Seurat object)
            "manorm2-full",               # call peaks for each dataset with MACS2, then run MAnorm2 with datasets
            "manorm2-half"                # call peaks for each comparison group with MACS2, then run MAnorm2 with datasets
        )
    )
    parser$add_argument(
        "--genome",
        help=paste(
            "Genome type of the sequencing data loaded from the Seurat",
            "object. It will be used for effective genome size selection",
            "when calling peaks with MACS2. Ignored if --test is not set",
            "to either manorm2-full or manorm2-half. Default: hs (2.7e9)"
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
            "if --test is not set to either manorm2-full or manorm2-half.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--minpeakgap",
        help=paste(
            "If a distance between peaks is smaller than the provided value",
            "they will be merged before splitting them into reference genomic",
            "bins of size --binsize. Ignored if --test is not set to either",
            "manorm2-full or manorm2-half. Default: 150"
        ),
        type="integer", default=150
    )
    parser$add_argument(
        "--binsize",
        help=paste(
            "The size of non-overlapping reference genomic bins used by",
            "MAnorm2 when generating a table of reads counts per peaks.",
            "Ignored if --test is not set to either manorm2-full or",
            "manorm2-half. Default: 1000"
        ),
        type="integer", default=1000
    )
    parser$add_argument(
        "--minoverlap",
        help=paste(
            "Keep only those reference genomic bins that are present",
            "in at least this fraction of datasets within each of the",
            "comparison groups. Used only when --test is set to",
            "manorm2-full. For manorm2-half this parameter will be",
            "automatically set to 1. Default: 0.5"
        ),
        type="double", default=0.5
    )
    parser$add_argument(
        "--maxpeaks",
        help=paste(
            "The maximum number of the most significant (based on qvalue)",
            "peaks to keep from each group of cells when constructing",
            "reference genomic bins. Ignored if --test is not set to",
            "either manorm2-full or manorm2-half. Default: keep all peaks"
        ),
        type="integer"
    )   
    parser$add_argument(
        "--blacklist",
        help=paste(
            "Path to the optional BED file with the genomic blacklist regions",
            "to be filtered out before running differential accessibility analysis.",
            "Any reference genomic bin overlapping a blacklist region will be",
            "removed from the output. Ignored if --test is not set to either",
            "manorm2-full or manorm2-half."
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
    parser$add_argument(
        "--seed",
        help="Seed number for random values. Default: 42",
        type="integer", default=42
    )
    args <- parser$parse_args(str_subset(commandArgs(trailingOnly=TRUE), "\\.R$", negate=TRUE))  # to exclude itself when executed from the sc_report_wrapper.R
    args$test <- switch(                                                                         # need to adjust --test parameter
        args$test,
        "negative-binomial"   = "negbinom",
        "poisson"             = "poisson",
        "logistic-regression" = "LR",
        "mast"                = "MAST",
        "manorm2-full"        = "manorm2-full",
        "manorm2-half"        = "manorm2-half"
    )
    if (args$test == "manorm2-half"){                                                            # when --test is set to "manorm2-half" all datasets within
        args$minoverlap <- 1                                                                     # the same comparison group have identical peaks so better
    }                                                                                            # to set --minoverlap to 1 to avoid confusion.
    logger$setup(
        file.path(dirname(ifelse(args$output == "", "./", args$output)), "error_report.txt"),
        header="Single-Cell ATAC-Seq Differential Accessibility Analysis (sc_atac_dbinding.R)"
    )
    print(args)
    return (args)
}

## ----
args <- get_args()
prod$parallel(args)

## ----
print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
debug$print_info(seurat_data, args)

## ----
if (!("ATAC" %in% names(seurat_data@assays))){
    logger$info(
        paste(
            "Loaded Seurat object doesn't include",
            "the required ATAC assay. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (!(args$reduction %in% names(seurat_data@reductions))){
    logger$info(
        paste0(
            "Loaded Seurat object doesn't include selected ",
            "reduction ", args$reduction, ". Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
print("Setting default assay to ATAC")
DefaultAssay(seurat_data) <- "ATAC"

## ----
seqinfo_data <- seurat_data[["ATAC"]]@seqinfo
if (is.null(seqinfo_data)){
    logger$info(
        "Loaded Seurat object doesn't include seqinfo data. Exiting."
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
print(paste("Loading ATAC fragments data from", args$fragments))
seurat_data <- io$replace_fragments(args$fragments, seurat_data)
debug$print_info(seurat_data, args)

## ----
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

## ----
if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)         # sets identities to new.ident
    debug$print_info(seurat_data, args)
}

## ----
if (args$test %in% c("manorm2-full", "manorm2-half")){
    # need to make sure that --splitby doesn't put
    # cells from the same dataset into the different
    # comparison groups, because when we use manorm2-full
    # or manorm2-half we aggregate fragments counts to
    # pseudobulk form per dataset.
    cells_counts <- table(
        seurat_data@meta.data$new.ident,
        seurat_data@meta.data[[args$splitby]]
    )
    if (any(rowSums(cells_counts > 0) > 1)){
        logger$info(
            paste(
                "Dividing cells by", args$splitby, "puts cells",
                "from the same dataset into the different",
                "comparison groups, which is not supported when",
                "--test parameter is set to either manorm2-full",
                "or manorm2-half. Exiting."
            )
        )
        logger$info(cells_counts)
        quit(save="no", status=1, runLast=FALSE)
    }
}

## ----
tryCatch(
    expr = {
        print("Subsetting Seurat object to include only cells from the tested conditions")
        seurat_data <- io$apply_metadata_filters(
            seurat_data,
            args$splitby,
            c(args$first, args$second)
        )
    },
    error = function(e){
        logger$info(
            paste(
                "Failed to filter Seurat object by",
                paste(c(args$first, args$second), collapse=", "),
                "values from the", args$splitby, "column.",
                "Exiting.", e
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
)
debug$print_info(seurat_data, args)

## ----
export_raw_plots(seurat_data, args)                                                        # need to run it befire we subset by --groupby

## ----
if(!is.null(args$groupby) && !is.null(args$subset)){
    tryCatch(
        expr = {
            print("Subsetting Seurat object to include only selected groups of cells")
            seurat_data <- io$apply_metadata_filters(
                seurat_data,
                args$groupby,
                args$subset
            )
        },
        error = function(e){
            logger$info(
                paste(
                    "Failed to filter Seurat object by",
                    paste(args$subset, collapse=", "),
                    "values from the", args$subset, "column.",
                    "Exiting.", e
                )
            )
            quit(save="no", status=1, runLast=FALSE)
        }
    )
    debug$print_info(seurat_data, args)
}

## ----
if (!all(table(seurat_data@meta.data[[args$splitby]]) > 0)){                               # check if we accidentally removed cells we want to compare
    logger$info(
        paste(
            "Not enough cells for comparison. Check --groupby",
            "and --subset parameters. Exiting."
        )
    )
    logger$info(table(seurat_data@meta.data[[args$splitby]]))
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (!(args$test %in% c("manorm2-full", "manorm2-half"))){                                  # needed only for FindMarkers (not for MAnorm2)
    print("Normalizing ATAC counts after all filters applied")
    seurat_data <- Signac::RunTFIDF(                                                       # might be redundant as it may not depend on the number of cells
        seurat_data,
        assay="ATAC",
        method=1,                                                                          # using default for RunTFIDF "log-tfidf" method
        verbose=FALSE
    )
}

## ----
args$metadata <- prepare_fragments_and_peaks(
    seurat_data=seurat_data,
    seqinfo_data=seqinfo_data,
    args=args
)

## ----
db_results <- analyses$atac_dbinding_analyze(seurat_data, args)
print(head(db_results$db_sites))

## ----
io$export_data(
    db_results$db_sites,                                                                   # not filtered na-removed differentially bound sites
    paste0(args$output, "_all_db_sites.tsv")
)

## ----
export_processed_plots(
    db_results=db_results,
    seqinfo_data=seqinfo_data,
    args
)

if (args$test %in% c("manorm2-full", "manorm2-half")){                                     # these are the only two cases when we create peaks, summits, and xls files
    for (i in 1:nrow(args$metadata)){                                                      # we update files after we run MAnorm because it fails with track information
        current_row <- as.list(args$metadata[i, ])
        track_label <- ifelse(
            args$test == "manorm2-full",
            paste(current_row$name, current_row$condition, sep=", "),                      # we include the dataset name only for manorm2-full analysis
            current_row$condition
        )
        track_color <- ifelse(
            current_row$condition == args$first,
            "0,145,106",                                                                   # TRUE_COLOR
            "235,99,49"                                                                    # FALSE_COLOR
        )
        track_info <- paste0(
            "track name=\"", track_label, "\" ",
            "description=\"", track_label, "\" ",
            "color=\"", track_color, "\""
        )

        for (current_ext in c("_peaks.narrowPeak", "_summits.bed", "_peaks.xls")){
            current_source <- paste0(args$output, "_", current_row$suffix, current_ext)    # the same way we define it in the prepare_fragments_and_peaks
            if(file.exists(current_source)){                                               # safety measure
                if (args$test == "manorm2-full"){                                          # keep all files, add track information to BED files
                    if(current_ext != "_peaks.xls"){                                       # we add track information only to peaks and summits files
                        print(paste("Adding track information to", current_source))
                        writeLines(
                            c(track_info, readLines(current_source)),
                            current_source
                        )
                    }
                } else {                                                                   # rename files, remove dublicates, add track information to BED files
                    updated_suffix <- paste(                                               # we want to exclude datasets from the filenames
                        tolower(
                            gsub("'|\"|\\s|\\t|#|%|&|-", "_", current_row$condition)       # the same as we did it in the prepare_fragments_and_peaks but using condition
                        ),
                        ifelse(current_row$condition == args$first, "first", "second"),
                        sep="_"
                    )
                    target_location <- paste0(
                        args$output, "_",
                        updated_suffix,                                                    # this now include only condition information
                        current_ext
                    )
                    if(!file.exists(target_location)){                                     # we write the file only if it has not been already created
                        print(paste("Copying", current_source, "to", target_location))
                        current_data <- readLines(current_source)
                        if(current_ext != "_peaks.xls"){                                   # we add track information only to peaks and summits files
                            print(paste("Adding", track_info))
                            current_data <- c(track_info, current_data)
                        }
                        writeLines(current_data, target_location)
                    }
                    print(paste("Removing", current_source))
                    file.remove(current_source)                                            # we always remove the source file
                }
            }
        }
    }
}