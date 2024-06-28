#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(cmapR))
suppressMessages(library(dplyr))
suppressMessages(library(limma))
suppressMessages(library(DESeq2))
suppressMessages(library(hopach))
suppressMessages(library(sva))
suppressMessages(library(Glimma))
suppressMessages(library(argparse))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(htmlwidgets))
suppressMessages(library(BiocParallel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(SummarizedExperiment))


# Useful links
# 1. https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
# 2. https://www.statlect.com/glossary/design-matrix#:~:text=A%20design%20matrix%20is%20a,each%20column%20to%20a%20characteristic.
# 3. https://paasp.net/accurate-design-of-in-vitro-experiments-why-does-it-matter/
# 4. https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
# 5. https://github.com/tavareshugo/tutorial_DESeq2_contrasts


COUNTS_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")
D40_COLORS <- c("#FF6E6A", "#71E869", "#6574FF", "#F3A6B5", "#FF5AD6", "#6DDCFE", "#FFBB70", "#43A14E", "#D71C7C", "#E1E333", "#8139A8", "#00D8B6", "#B55C00", "#7FA4B6", "#FFA4E3", "#B300FF", "#9BC4FD", "#FF7E6A", "#9DE98D", "#BFA178", "#E7C2FD", "#8B437D", "#ADCDC0", "#FE9FA4", "#FF53D1", "#D993F9", "#FF47A1", "#FFC171", "#625C51", "#4288C9", "#9767D4", "#F2D61D", "#8EE6FD", "#B940B1", "#B2D5F8", "#9AB317", "#C70000", "#AC8BAC", "#D7D1E4", "#9D8D87")

set_cpus <- function (cpus) {
    register(MulticoreParam(cpus))
}


export_volcano_plot <- function(data, rootname, x_axis, y_axis, x_cutoff, y_cutoff, x_label, y_label, plot_title, plot_subtitle, caption, features=NULL, label_column="feature", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- EnhancedVolcano(
                        data,
                        lab=data[,label_column],
                        x=x_axis,
                        y=y_axis,
                        FCcutoff=x_cutoff,
                        pCutoff=y_cutoff,
                        xlab=x_label,
                        ylab=y_label,
                        selectLab=features,
                        title=plot_title,
                        subtitle=plot_subtitle,
                        caption=caption,
                        labSize=4,
                        labFace="bold",
                        labCol="red4",
                        colAlpha=0.6,
                        col=c("grey30", "forestgreen", "royalblue", "red"),
                        drawConnectors=TRUE,
                        widthConnectors=0.2
                    ) +
                    scale_y_log10() +
                    theme_classic() +
                    theme(legend.position="none", plot.subtitle=element_text(size=8, face="italic", color="gray30"))

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export volcano plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(err){print(paste("Called  dev.off() with error -", err))})
            print(paste0("Failed to export volcano plot to ", rootname, ".(png/pdf) with error - ", e))
        }
    )
}


export_pca_plot <- function(norm_counts_data, rootname, intgroup, plot_title, plot_subtitle, ntop=500, palette_colors=D40_COLORS, alpha=0.5, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            pca_data <- plotPCA(
                norm_counts_data,
                intgroup=intgroup,
                ntop=ntop,               # ntop the most variable features
                returnData=TRUE
            )
            percentVar <- round(100 * attr(pca_data, "percentVar"))
            plot <- ggplot(pca_data, aes(PC1, PC2, color=group)) +
                    geom_point(size=4, shape=19, alpha=alpha) +
                    xlab(paste0("PC1: ",percentVar[1], "% variance")) +
                    ylab(paste0("PC2: ",percentVar[2], "% variance")) + 
                    ggtitle(plot_title, subtitle=plot_subtitle) +
                    geom_text_repel(
                        aes(label=name),
                        point.padding=0.5,
                        box.padding=0.5,
                        check_overlap=TRUE,
                        show.legend=FALSE
                    ) +
                    theme_classic() +
                    theme(plot.subtitle=element_text(size=8, face="italic", color="gray30")) +
                    ggplot2::scale_color_manual(values=palette_colors)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export PCA plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(err){print(paste("Called  dev.off() with error -", err))})
            print(paste0("Failed to export PCA plot to ", rootname, ".(png/pdf) with error - ", e))
        }
    )
}


export_mds_html_plot <- function(norm_counts_data, location){
    tryCatch(
        expr = {
            htmlwidgets::saveWidget(
                glimmaMDS(
                    x=assay(norm_counts_data),
                    groups=as.data.frame(SummarizedExperiment::colData(norm_counts_data)),
                    labels=rownames(SummarizedExperiment::colData(norm_counts_data))
                ),
                file=location
            )
        },
        error = function(e){
            print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
        }
    )
}


get_highlight_features <- function(diff_expr_features, args){
    print(
        paste(
            "Filtering DESeq results to include only",
            "features with padj <= ", args$padj,
            "and |log2FoldChange| >= ", args$logfc
        )
    )
    filt_diff_expr_features <- diff_expr_features %>%
                               filter(.$padj<=args$padj, abs(.$log2FoldChange)>=args$logfc)
    filt_diff_expr_features <- filt_diff_expr_features %>% arrange(desc(log2FoldChange))
    print(paste("Number of significantly differentially expressed features:", nrow(filt_diff_expr_features)))

    topn_diff_expr_features <- filt_diff_expr_features %>% filter(row_number() > max(row_number()) - 10 | row_number() <= 10)
    highlight_features <- as.vector(as.character(topn_diff_expr_features[, "feature"]))                      # default features to highlight
    if (!is.null(args$label)){
        print("Check features of interest to include only those that are differentially expressed regardless of significance tresholds")
        args$label <- unique(args$label)
        args$label <- args$label[args$label %in% as.vector(as.character(diff_expr_features[, "feature"]))]
        highlight_features <- args$label
    }
    print(paste("Features to highlight", paste(highlight_features, collapse=", ")))
    return (highlight_features)
}


export_plots <- function(diff_expr_data, norm_counts_data, metadata, args){
    print("Identifying features to highlight")
    highlight_features <- get_highlight_features(diff_expr_data$res, args)
    export_volcano_plot(
        data=diff_expr_data$res,                                   # this is not filtered differentially expressed features
        rootname=paste(args$output, "volcano_plot", sep="_"),
        x_axis="log2FoldChange",
        y_axis="padj",
        x_cutoff=args$logfc,
        y_cutoff=args$padj,
        x_label="log2FoldChange",
        y_label="-log10 Padj",
        plot_title="Differentially expressed features",
        plot_subtitle=paste0(
            "Differentially expressed features with padj", "<=", args$padj,
            " and |log2FoldChange| >= ", args$logfc,
            ifelse(
                is.null(args$contrast),
                " from multiple contrasts",
                paste0(" from ", args$contrast, " contrast")
            )
        ),
        caption=paste(nrow(diff_expr_data$res), "features"),
        features=highlight_features,
        pdf=args$pdf
    )
    export_pca_plot(
        norm_counts_data=norm_counts_data,
        rootname=paste(args$output, "pca_plot", sep="_"),
        intgroup=colnames(metadata),
        plot_title=paste0("PCA plot of ", args$norm, "-normalized read counts"),
        plot_subtitle=paste(
            "Based on the top 500 features.",
            ifelse(!is.null(args$remove), paste("Remove the effect of", args$remove), "")
        ),
        ntop=500,
        pdf=args$pdf
    )
    export_mds_html_plot(
        norm_counts_data=norm_counts_data,
        location=paste(args$output, "mds_plot.html", sep="_")
    )
}


export_gct <- function(counts_mat, row_metadata, col_metadata, location){
    tryCatch(
        expr = {
            row_metadata <- row_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            col_metadata <- col_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            gct_data <- new(
                "GCT",
                mat=counts_mat[row_metadata$id, col_metadata$id],       # to guarantee the order and number of row/columns
                rdesc=row_metadata,
                cdesc=col_metadata
            )
            write_gct(
                ds=gct_data,
                ofile=location,
                appenddim=FALSE
            )
            print(paste("Exporting GCT data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export GCT data to", location, sep=" "))
        }
    )
}


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


load_metadata <- function(args){
    metadata <- read.table(
        args$metadata,
        sep=get_file_type(args$metadata),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )  %>% remove_rownames() %>% column_to_rownames("sample") %>% mutate_at(colnames(.), factor)
    if (!is.null(args$base)){
        print(
            paste(
                "Attempting to relevel metadata table based on",
                paste(args$base, collapse=", "), "values"
            )
        )
        for (i in 1:length(args$base)) {
            current_base_level <- args$base[i]
            tryCatch(
                expr = {
                    current_column <- colnames(metadata)[i]
                    metadata[[current_column]] <- relevel(metadata[[current_column]], current_base_level)
                    print(paste("Setting", current_base_level, "as a base level for", current_column))
                },
                error = function(e){
                    print(paste("Failed to set", current_base_level, "as a base level"))
                }
            )
        }
    }
    print("Loaded metadata")
    print(metadata)
    return (metadata)
}


load_expression_data <- function(args, counts_colname=COUNTS_COL, rpkm_colname=RPKM_COL, intersect_by=INTERSECT_BY) {
    collected_expression_data <- NULL
    for (i in 1:length(args$expression)) {
        location <- args$expression[i]
        alias <- args$aliases[i]
        expression_data <- read.table(location, sep=get_file_type(location), header=TRUE, stringsAsFactors=FALSE)
        print(paste("Loading", nrow(expression_data), "rows from", location, "as", alias))
        colnames(expression_data)[colnames(expression_data) == counts_colname] <- paste(alias, counts_colname, sep=" ")
        colnames(expression_data)[colnames(expression_data) == rpkm_colname] <- paste(alias, rpkm_colname, sep=" ")
        if (is.null(collected_expression_data)){
            collected_expression_data <- expression_data
        } else {
            collected_expression_data <- merge(collected_expression_data, expression_data, by=intersect_by, sort = FALSE)
        }
    }
    print(paste("Number of rows common for all loaded files ", nrow(collected_expression_data), sep=""))
    collected_expression_data <- collected_expression_data %>% remove_rownames()
    if (args$type == "gene"){
        collected_expression_data <- collected_expression_data %>% column_to_rownames("GeneId")
    } else {
        collected_expression_data$RefseqId <- make.unique(names=collected_expression_data$RefseqId, sep="_")     # in case we have duplicate RefseqId
        collected_expression_data <- collected_expression_data %>% column_to_rownames("RefseqId")
    }
    all_features <- as.vector(as.character(rownames(collected_expression_data)))
    selected_features <- all_features[!all_features %in% args$exclude]
    excluded_features <- all_features[all_features %in% args$exclude]   # not used elsewhere, only to print to the console
    if (length(excluded_features) > 0){
        print(paste("Following features will be excluded from the differential expression analysis:", paste(excluded_features, collapse=", ")))
    }
    collected_expression_data <- collected_expression_data[selected_features, ]

    print("Keeping only those features which total counts for all samples are bigger then 0")
    counts_columns <- grep(
        paste(COUNTS_COL, sep=""),
        colnames(collected_expression_data),
        value=TRUE,
        ignore.case=TRUE
    )
    keep <- rowSums(collected_expression_data[counts_columns]) > 0
    collected_expression_data <- collected_expression_data[keep,]
    print(paste("Number of the remaining features: ", nrow(collected_expression_data), sep=""))
    return (collected_expression_data)
}


assert_args <- function(args){
    if (length(args$expression) != length(args$aliases)){
            print("Exiting: --expression and --aliases have different number of values")
            quit(save = "no", status = 1, runLast = FALSE)
    }
    return (args)
}


get_contrast <- function(design_formula, deseq_data, metadata, contrast_formula, selected_samples=NULL){
    model <- model.matrix(design_formula, metadata)
    print("Model matrix")
    print(model)
    all_design_vars <- unique(all.vars(as.formula(design_formula)))     # need to use design variable instead of metadata columns
    for (i in 1:length(all_design_vars)){                               # because with --remove design formula get updates but
        current_design_var <- all_design_vars[i]                        # metadata remains the same
        unique_keys <- unique(metadata[[current_design_var]])
        for (j in 1:length(unique_keys)){
            key <- as.character(unique_keys[j])
            print(
                paste0(
                    "Filtering model matrix by ",
                    current_design_var, ": ", key
                )
            )
            if (!is.null(selected_samples)) {                           # need to subset by sample and key in one step
                print(
                    paste0(
                        "Filtering model matrix by sample: ",
                        paste(selected_samples, collapse=", ")
                    )
                )
                subset_model <- model[deseq_data[[current_design_var]] == key & rownames(model) %in% selected_samples, ]
            } else {
                subset_model <- model[deseq_data[[current_design_var]] == key, ]
            }
            print(subset_model)
            assign(key, colMeans(subset_model))
        }
    }
    print(paste("Evaluating the contrast", contrast_formula))
    contrast <- eval(parse(text=contrast_formula))
    print(contrast)
    return (contrast)
}


get_diff_expr_data <- function(expression_data, metadata, args){

    print("Selecting all columns with raw read counts data.")
    counts_columns <- grep(
        paste(COUNTS_COL, sep=""),
        colnames(expression_data),
        value=TRUE,
        ignore.case=TRUE
    )
    counts_data = expression_data[counts_columns]
    colnames(counts_data) <- lapply(
        colnames(counts_data),
        function(s){
            paste(head(unlist(strsplit(s, " ", fixed=TRUE)), -1), collapse=" ")
        }
    )
    print(
        paste(
            "Reordering read counts data columns based on the",
            "row names from the provided metadata file."
        )
    )
    counts_data <- counts_data[, as.vector(rownames(metadata))]

    if ( all(colnames(counts_data) != rownames(metadata)) ){              # safety measure
        print(
            paste(
                "Assert failed: columns order of the counts data",
                "should be the same as rows order in metadata."
            )
        )
        quit(save = "no", status = 1, runLast = FALSE)
    }

    if(!is.null(args$remove)){
        print(paste("Removing the effect of", args$remove, "using ComBat-Seq"))
        counts_data <- ComBat_seq(
            as.matrix(counts_data),
            batch=metadata[[args$remove]],       # columns from counts data are alread ordered by rownames from metadata
            group=NULL
        )
        print(
            paste(
                "Removing", args$remove, "item from design and",
                "reduced formulas if it was used there"
            )
        )
        args$design <- paste0(
            "~",
            paste(
                grep(
                    args$remove,
                    unlist(
                        strsplit(
                            gsub("~| ", "", args$design),
                            "\\+"
                        )
                    ),
                    value=TRUE,
                    ignore.case=TRUE,
                    invert=TRUE
                ),
                collapse="+"
            )
        )
        print(paste("Updated design formula", args$design))
        if (!is.null(args$reduced)){
            args$reduced <- paste0(
                "~",
                paste(
                    grep(
                        args$remove,
                        unlist(
                            strsplit(
                                gsub("~| ", "", args$reduced),
                                "\\+"
                            )
                        ),
                        value=TRUE,
                        ignore.case=TRUE,
                        invert=TRUE
                    ),
                    collapse="+"
                )
            )
            print(paste("Updated reduced formula", args$reduced))
        }
    }

    print("Loadind data to DESeq2")
    deseq_data <- DESeqDataSetFromMatrix(
        countData=counts_data,
        colData=metadata,
        design=as.formula(args$design)
    )

    if (!is.null(args$reduced)){
        print("Using LRT test to calculate p-values")
        deseq_data <- DESeq(
            deseq_data,
            test="LRT",
            reduced=as.formula(args$reduced),
            quiet=TRUE,
            parallel=TRUE,
            BPPARAM=MulticoreParam(args$cpus)  # add it here as well just in case
        )
    } else {
        print("Using Wald test to calculate p-values")
        deseq_data <- DESeq(
            deseq_data,
            quiet=TRUE,
            parallel=TRUE,
            BPPARAM=MulticoreParam(args$cpus)  # add it here as well just in case
        )
    }
    print("Estimated effects")
    print(resultsNames(deseq_data))

    if (is.null(args$contrast)){                                               # will have to create multiple contrasts
        print(
            paste(
                "Contrast is not set. Union of contrasts for all",
                "pairwise combinations of metadata values",
                "from the columns used in the design formula",
                "will be used instead."
            )
        )
        all_contrasts <- c()                                                   # to keep all collected contrasts
        all_design_vars <- unique(all.vars(as.formula(args$design)))
        for (i in 1:length(all_design_vars)) {
            current_design_var <- all_design_vars[i]
            remaining_design_vars <- setdiff(                                  # to define filtering groups
                all_design_vars,                                               # but it might be emply vector as well
                current_design_var
            )
            all_contrasts <- append(
                all_contrasts,
                unlist(
                    lapply(
                        combn(                                         # will return all combinations of size 2
                            x=unique(
                                metadata[[current_design_var]]         # all metadata columns are already factors
                            ),
                            m=2,
                            simplify=FALSE
                        ),
                        function(v){
                            current_pair <- as.vector(sort(v))                             # sorting by levels and converting to vector
                            sample_groups <- metadata %>%
                                filter(!!sym(current_design_var) %in% current_pair) %>%
                                select(all_of(remaining_design_vars)) %>%                  # may select 0 columns when remaining_design_vars is empty
                                rownames_to_column(var="sample") %>%
                                group_by(across(all_of(remaining_design_vars))) %>%        # may create one group when remaining_design_vars is empty
                                mutate(
                                    sample=paste(sample, collapse="@")
                                ) %>%
                                ungroup() %>%
                                select(sample) %>%
                                distinct() %>%
                                pull(sample)
                            lapply(
                                sample_groups,
                                function(l){
                                    selected_samples <- unlist(strsplit(l, "@"))
                                    alias <- metadata %>%
                                        rownames_to_column(var="sample") %>%
                                        filter(sample %in% selected_samples) %>%
                                        select(all_of(remaining_design_vars)) %>%          # these columns should be identical for all selected samples
                                        mutate(
                                            alias = if (ncol(.) == 0)
                                                        ""
                                                    else
                                                        paste0("_for_", do.call(paste, c(., sep="_")))
                                        ) %>%
                                        select(alias) %>%
                                        distinct() %>%
                                        pull(alias)                                        # should be always a vector with length 1
                                    if (length(alias) > 1){
                                        print("Exiting: alias can't have length more than 1")
                                        quit(save = "no", status = 1, runLast = FALSE)
                                    }
                                    list(
                                        contrast=paste0(
                                            current_pair[2], "-", current_pair[1]      # the order of comparison is always higher level to lower level
                                        ),
                                        samples=selected_samples,                      # may include all samples filtered only by current_design_var
                                        alias=gsub(
                                            "'|\"|\\s|\\t|#|%|&|-",
                                            "_",
                                            paste0(
                                                current_pair[2],
                                                "_vs_",
                                                current_pair[1],
                                                alias                                  # can be empty string
                                            )
                                        )
                                    )
                                }
                            )
                        }
                    ),
                recursive=FALSE
                )
            )
        }
        print("Collected contrasts")
        print(all_contrasts)
        if (length(all_contrasts) == 0){
            print("Exiting: no contrasts collected")
            quit(save = "no", status = 1, runLast = FALSE)
        }
        deseq_results <- list()
        for (i in 1:length(all_contrasts)) {
            current_contrast <- all_contrasts[[i]]
            print(
                paste0(
                    "Processing contrast ", current_contrast$alias,
                    " (", paste(current_contrast$samples, collapse=", "), ")"
                )
            )
            current_deseq_results <- results(
                deseq_data,
                contrast=get_contrast(
                    design_formula=as.formula(args$design),
                    deseq_data=deseq_data,
                    metadata=metadata,
                    contrast_formula=current_contrast$contrast,
                    selected_samples=current_contrast$samples
                ),
                parallel=TRUE,
                BPPARAM=MulticoreParam(args$cpus)         # add it here as well just in case
            )
            print("Current results description")
            print(mcols(current_deseq_results))
            print(mcols(current_deseq_results)$description)
            deseq_results[[current_contrast$alias]] <- as.data.frame(current_deseq_results) %>%
                rownames_to_column(var="feature")
            print(head(deseq_results[[current_contrast$alias]]))
        }
        print("Merging results from multiple contrasts")
        deseq_results <- bind_rows(deseq_results, .id="contrast") %>%
            dplyr::mutate(
                "padj_th"=ifelse(                      # to be able to sort by logFC within significant padj
                    .$padj <= args$padj,
                    args$padj,
                    .$padj
                )
            ) %>%
            dplyr::group_by(feature) %>%
            dplyr::arrange(
                padj_th,                               # sort primarily by |log2FoldChange| within significant padj
                desc(abs(log2FoldChange)),             # then, if padj is not signif. anymore, sort by padj and |log2FoldChange|
                .by_group=TRUE
            ) %>%
            dplyr::slice_head(n=1) %>%                 # all NA should be at the bottom
            dplyr::ungroup() %>%
            dplyr::arrange(
                contrast, desc(log2FoldChange)
            ) %>%
            dplyr::select(-c("padj_th")) %>%           # drop temporary column
            relocate(contrast, .after=last_col())
        print(deseq_results)
    } else {
        print(
            paste(
                "Using user-provided contrast:",
                args$contrast
            )
        )
        deseq_results <- results(
            deseq_data,
            contrast=get_contrast(
                design_formula=as.formula(args$design),
                deseq_data=deseq_data,
                metadata=metadata,
                contrast_formula=args$contrast
            ),
            parallel=TRUE,
            BPPARAM=MulticoreParam(args$cpus)         # add it here as well just in case
        )
        print("Results description")
        print(mcols(deseq_results))
        print(mcols(deseq_results)$description)
        deseq_results <- as.data.frame(deseq_results) %>%
            rownames_to_column(var="feature") %>%
            dplyr::arrange(desc(log2FoldChange))
    }

    print("Adding extra columns to the DESeq output")
    rpkm_columns <- grep(                                                    # for column ordering only
        paste(RPKM_COL, sep=""),
        colnames(expression_data),
        value=TRUE,
        ignore.case=TRUE
    )

    deseq_results <- expression_data[as.vector(as.character(deseq_results$feature)), ] %>%    # need to reorder rows in expression_data to correspond to the deseq_results
        bind_cols(deseq_results) %>%
        na.omit() %>%
        relocate(any_of(rpkm_columns), .after=last_col()) %>%                    # move all rpkm columns to the end
        relocate(any_of(counts_columns), .after=last_col())                      # move all read counts columns to the end
    print(
        paste(
            "Number of differentially expressed features",
            "after excluding NA:", nrow(deseq_results)
        )
    )
    print("DESeq2 results")
    print(head(deseq_results))

    if (!base::identical(rownames(deseq_results), deseq_results$feature)){                    # safety measure
        print("Exiting: genes order in the DESeq results is not correct")
        quit(save = "no", status = 1, runLast = FALSE)
    }

    return (list(
        res=deseq_results,
        raw=deseq_data
    ))
}


get_norm_counts_data <- function(deseq_data, metadata, args){
    if (args$norm == "vst"){
        print("Applying vst transformation (not blind to the experimental design)")
        norm_counts_data <- DESeq2::vst(deseq_data, blind=FALSE)
    } else {
        print("Applying rlog transformation (not blind to the experimental design)")
        norm_counts_data <- DESeq2::rlog(deseq_data, blind=FALSE)
    }

    # if(!is.null(args$remove)){
    #     print(
    #         paste("Removing the effect of", args$remove, "from the normalized counts")
    #     )
    #     keep_formula <- paste(
    #         grep(
    #             args$remove,
    #             unlist(strsplit(args$design, "\\+")),
    #             value=TRUE,
    #             ignore.case=TRUE,
    #             invert=TRUE
    #         ),
    #         collapse="+"
    #     )
    #     print(paste("Formula to include conditions to be preserved", keep_formula))
    #     assay(norm_counts_data) <- limma::removeBatchEffect(
    #         assay(norm_counts_data),
    #         batch=norm_counts_data[[args$remove]],
    #         design=stats::model.matrix(stats::as.formula(keep_formula), metadata)
    #     )
    # }

    print("Normalized read counts")
    print(head(assay(norm_counts_data)))
    print(dim(assay(norm_counts_data)))
    print(SummarizedExperiment::colData(norm_counts_data))
    return (norm_counts_data)
}


get_clustered_data <- function(expression_data, center, dist, transpose) {

    if (transpose){
        print("Transposing expression data")
        expression_data = t(expression_data)
    }
    if (!is.null(center)) {
        print(paste("Centering expression data by ", center, sep=""))
        if (center == "mean"){
            expression_data = expression_data - rowMeans(expression_data)    
        } else {
            expression_data = expression_data - rowMedians(data.matrix(expression_data))    
        }
    }
    print("Creating distance matrix")
    distance_matrix <- distancematrix(expression_data, dist)
    print("Running HOPACH")
    hopach_results <- hopach(expression_data, dmat=distance_matrix)

    if (transpose){
        print("Transposing expression data")
        expression_data = t(expression_data)
    }

    print("Parsing cluster labels")
    options(scipen=999)                                           # need to temporary disable scientific notation, because nchar gives wrong answer
    clusters = as.data.frame(hopach_results$clustering$labels)
    colnames(clusters) = "label"
    clusters = cbind(
        clusters,
        "HCL"=outer(
            clusters$label,
            10^c((nchar(trunc(clusters$label))[1]-1):0),
            function(a, b) {
                paste0("c", a %/% b)
            }
        )
    )
    clusters = clusters[, c(-1), drop=FALSE]
    options(scipen=0)                                            # setting back to the default value

    return (
        list(
            order=as.vector(hopach_results$clustering$order),
            expression=expression_data,
            clusters=clusters
        )
    )
}


export_data <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE, digits=NULL){
    tryCatch(
        expr = {
            if (!is.null(digits)){
                data <- format(data, digits=digits)
            }
            write.table(
                data,
                file=location,
                sep=get_file_type(location),
                row.names=row_names,
                col.names=col_names,
                quote=quote
            )
            print(paste("Export data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data to", location, sep=" "))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description="DESeq2 Multi-factor Analysis")
    parser$add_argument(
        "--expression",
        help=paste(
            "Path to the TSV/CSV files with expression data.",
            "All files should have the following header:",
            "RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm"
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--aliases",
        help=paste(
            "Unique names for files provided in --expression,",
            "no special characters or spaces are allowed.",
            "Number and order of the names should corresponds",
            "to values from --expression."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to provide metadata for the",
            "samples from --expression. First column should have",
            "the name 'sample', other columns may have arbitrary names.",
            "The values from the 'sample' column should correspond to",
            "the values provided in --aliases. For a proper --contrast",
            "intepretation, values defined in each column should not be",
            "used in other columns. All metadata columns are treated as",
            "factors (no covariates are supported)."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--design",
        help=paste(
            "Design formula. Should start with ~ and include terms from",
            "the --metadata table."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--reduced",
        help=paste(
            "Reduced formula with the term(s) of interest removed.",
            "Should start with ~. If provided, force DESeq2 to run",
            "LRT test instead of the Wald."
        ),
        type="character"
    )
    parser$add_argument(
        "--contrast",
        help=paste(
            "Contrast to be be applied for the output, formatted as",
            "a mathematical formula of values from the --metadata table.",
            "If not provided, all possible combinations of values from",
            "the metadata columns present in the --design will be used",
            "(results will be merged giving the priority to significantly",
            "differentially expressed genes with higher absolute",
            "log2FoldChange values)."
        ),
        type="character"
    )
    parser$add_argument(
        "--base",
        help=paste(
            "Value(s) from each metadata file column(s) to be set as",
            "the base level(s). Number and order of provided values should",
            "correspond the order of columns in --metadata file. Default:",
            "define base levels alphabetically for each metadata column."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--type",
        help=paste(
            "Feature type to use for differential expression.",
            "If set to 'gene', use 'GeneId' column from the provided in --expression files.",
            "If set to 'transcript', use 'RefseqId' from the provided in --expression files.",
            "Default: gene"
        ),
        type="character", default="gene",
        choices=c("gene", "transcript")
    )
    parser$add_argument(
        "--exclude",
        help=paste(
            "Features to be excluded from the differential expression analysis.",
            "Default: include all features"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--norm",
        help=paste(
            "Read counts normalization for the exploratory visualization analysis.",
            "Use 'vst' for medium-to-large datasets (n > 30) and 'rlog' for",
            "small datasets (n < 30), when there is a wide range of sequencing",
            "depth across samples.",
            "Default: vst"
        ),
        type="character", default="vst",
        choices=c("vst", "rlog")
    )
    parser$add_argument(
        "--remove",
        help=paste(
            "Column from the metadata file to remove batch effect",
            "before running differential expression analysis. If",
            "present, all components that include this term will be",
            "removed from the design and reduced formulas.",
            "Default: do not remove batch effect"
        ),
        type="character"
    )
    parser$add_argument(
        "--cluster",
        help=paste(
            "Hopach clustering method to be run on normalized read counts for the",
            "exploratory visualization analysis. Default: do not run clustering"
        ),
        type="character",
        choices=c("row", "column", "both")
    )
    parser$add_argument(
        "--rowdist",
        help=paste(
            "Distance metric for HOPACH row clustering. Ignored if --cluster is not",
            "provided. Default: cosangle"
        ),
        type="character", default="cosangle",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--columndist",
        help=paste(
            "Distance metric for HOPACH column clustering. Ignored if --cluster is not",
            "provided. Default: euclid"
        ),
        type="character", default="euclid",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--center",
        help=paste(
            "Apply mean centering for feature expression prior to running",
            "clustering by row. Ignored when --cluster is not row or both.",
            "Default: do not centered"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--label",
        help=paste(
            "Features of interest to label on the generated volcanot plot. Default:",
            "top 10 features with the highest and the lowest log2 fold change",
            "expression values."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--padj",
        help=paste(
            "In the exploratory visualization analysis output only features with",
            "adjusted P-value not bigger than this value. Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--logfc",
        help=paste(
            "In the exploratory visualization analysis output only features with",
            "absolute log2FoldChange bigger or equal to this value. Default: 0"
        ),
        type="double", default=0
    )
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help=paste(
            "Output prefix for generated files"
        ),
        type="character", default="./deseq"
    )
    parser$add_argument(
        "--cpus",
        help="Number of cores/cpus to use. Default: 1",
        type="integer", default=1
    )
    args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
    return (args)
}


# Parse arguments
args <- get_args()

print("Used parameters")
print(args)

print(paste("Setting parallelizations to", args$cpus, "cores"))
set_cpus(args$cpus)

print("Loading expression data")
expression_data <- load_expression_data(args)

print(paste("Loading metadata from", args$metadata))
metadata <- load_metadata(args)

print("Identifying differentially expressed features")
diff_expr_data <- get_diff_expr_data(expression_data, metadata, args)

print("Normalizing read count data")
norm_counts_data <- get_norm_counts_data(diff_expr_data$raw, metadata, args)     # includes all genes even those that have NA in diff_expr_data$res

export_plots(diff_expr_data, norm_counts_data, metadata, args)

print(
    paste(
        "Filtering normalized read counts matrix to include",
        "only differentially expressed features with padj <=", args$padj,
        "and |log2FoldChange| >=", args$logfc
    )
)

row_metadata <- diff_expr_data$res %>%
                remove_rownames() %>%
                column_to_rownames("feature") %>%
                dplyr::select(
                    any_of(
                        c("log2FoldChange", "pvalue", "padj", "contrast")        # need any_of because "contrast" might not exist
                    )
                ) %>%
                filter(.$padj<=args$padj, abs(.$log2FoldChange)>=args$logfc)

col_metadata <- metadata %>%
                mutate_at(colnames(.), as.vector)                                # need to convert to vector, because in our metadata everything was a factor

norm_counts_mat <- assay(norm_counts_data)[as.vector(rownames(row_metadata)), ]
print("Size of the normalized read counts matrix after filtering")
print(dim(norm_counts_mat))

if (!is.null(args$cluster)){
    if (args$cluster == "column" || args$cluster == "both") {
        print("Clustering filtered read counts by columns")
        clustered_data = get_clustered_data(
            expression_data=norm_counts_mat,
            center=NULL,                                              # centering doesn't influence on the samples order
            dist=args$columndist,
            transpose=TRUE
        )
        col_metadata <- cbind(col_metadata, clustered_data$clusters)  # adding cluster labels
        col_metadata <- col_metadata[clustered_data$order, ]          # reordering samples order based on the HOPACH clustering resutls
        print("Reordered samples")
        print(col_metadata)
    }
    if (args$cluster == "row" || args$cluster == "both") {
        print("Clustering filtered normalized read counts by rows")
        clustered_data = get_clustered_data(
            expression_data=norm_counts_mat,
            center=if(args$center) "mean" else NULL,                  # about centering normalized data https://www.biostars.org/p/387863/
            dist=args$rowdist,
            transpose=FALSE
        )
        norm_counts_mat <- clustered_data$expression                  # can be different because of centering by rows mean
        row_metadata <- cbind(row_metadata, clustered_data$clusters)  # adding cluster labels
        row_metadata <- row_metadata[clustered_data$order, ]          # reordering features order based on the HOPACH clustering results
        print("Reordered features")
        print(head(row_metadata))
        cluster_columns <- grep(
            "HCL",
            colnames(row_metadata),
            value=TRUE,
            ignore.case=TRUE
        )
        if (length(cluster_columns) > 0){                             # check the length just in case
            diff_expr_data$res <- merge(
                diff_expr_data$res,
                row_metadata[, cluster_columns, drop=FALSE] %>% rownames_to_column(var="feature"),    # need drop=FALSE in case only one HCL column present
                by="feature",
                all.x=TRUE,
                sort=FALSE
            )
        }
    }
}

print("Exporting differentially expressed features")
export_data(                                                          # may include HCL columns from clustering by row
    diff_expr_data$res,                                                          # this is not filtered differentially expressed features
    location=paste(args$output, "diff_expr_features.tsv", sep="_"),
    digits=5
)

# we do not reorder norm_counts_mat based on the clustering order
# because when exportin to GCT we use row_metadata and col_metadata
# to force the proper order of rows and columns

print("Exporting normalized read counts to GCT format")
export_gct(
    counts_mat=norm_counts_mat,
    row_metadata=row_metadata,                                        # includes features as row names
    col_metadata=col_metadata,                                        # includes samples as row names
    location=paste(args$output, "_norm_read_counts.gct", sep="")
)