import("dplyr", attach=FALSE)
import("purrr", attach=FALSE)
import("tidyr", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("tibble", attach=FALSE)
import("tidyselect", attach=FALSE)
import("data.table", attach=FALSE)
import("scDblFinder", attach=FALSE)
import("bestNormalize", attach=FALSE)
import("GenomicRanges", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)
import("SingleCellExperiment", attach=FALSE)

export(
    "qc_metrics_pca",
    "counts_pca",
    "estimate_doublets",
    "add_rna_qc_metrics",
    "add_atac_qc_metrics",
    "add_peak_qc_metrics",
    "quartile_qc_metrics",
    "add_gene_expr_percentage"
)


qc_metrics_pca <- function(seurat_data, qc_columns, qc_labels, orq_transform=FALSE){
    base::tryCatch(
        expr = {
            base::print("Computing PCA for the following QC metrics")
            base::print(base::paste(qc_labels, collapse=", "))
            qc_columns_corrected <- c()
            qc_labels_corrected <- c()
            for (i in 1:length(qc_columns)){
                if (qc_columns[i] %in% base::colnames(seurat_data@meta.data)){
                    qc_columns_corrected <- c(qc_columns_corrected, qc_columns[i])
                    qc_labels_corrected <- c(qc_labels_corrected, qc_labels[i])
                } else {
                    base::print(
                        base::paste(
                            "Column", qc_columns[i], "was not found,",
                            "skipping", qc_labels[i]
                        )
                    )
                }
            }
            target_data <- base::as.data.frame(seurat_data[[qc_columns_corrected]]) %>%
                           tidyr::drop_na() %>%
                           dplyr::filter_all(dplyr::all_vars(!is.infinite(.)))
            base::print(
                base::paste(
                    "Cells removed due to having infinite values",
                    "in any of the selected column -",
                    base::nrow(seurat_data@meta.data) - base::nrow(target_data)
                )
            )
            if (!is.null(orq_transform) && orq_transform){
                base::print("Running Ordered Quantile (ORQ) normalization transformation")
                target_data <- target_data %>%
                               dplyr::mutate_all(function(x){return (bestNormalize::orderNorm(x)$x.t)})
            }
            pca_raw <- stats::prcomp(
                base::t(target_data),
                center=!orq_transform,         # no need to center or scale when data is already ORQ-transformed
                scale.=!orq_transform
            )
            pca_scores <- base::as.data.frame(pca_raw$x)
            pca_scores$labels <- qc_labels_corrected
            pca_variance <- round(pca_raw$sdev / sum(pca_raw$sdev) * 100, 2)
            return (list(scores=pca_scores, variance=pca_variance))
        },
        error = function(e){
            base::print(base::paste("Failed to compute PCA for QC metrics due to", e))
        }
    )
}

counts_pca <- function(counts_data){
    base::tryCatch(
        expr = {
            base::print("Computing PCA for counts data")
            target_data <- base::as.data.frame(counts_data) %>%
                           tidyr::drop_na() %>%
                           dplyr::filter_all(dplyr::all_vars(!is.infinite(.))) %>%
                           dplyr::filter_all(dplyr::any_vars(. != 0))                       # remove rows with only zeros, otherwise prcomp fails
            pca_raw <- stats::prcomp(
                base::t(target_data),
                center=TRUE,
                scale.=TRUE
            )
            pca_scores <- base::as.data.frame(pca_raw$x) %>%
                          tibble::rownames_to_column(var="group")
            pca_variance <- round(pca_raw$sdev / sum(pca_raw$sdev) * 100, 2)
            return (list(scores=pca_scores, variance=pca_variance))
        },
        error = function(e){
            base::print(base::paste("Failed to compute PCA for counts data due to", e))
        }
    )
}

estimate_doublets <- function(seurat_data, assay, target_column, dbl_rate=NULL, dbl_rate_sd=NULL){    # the logic depends on the value of assay input ("RNA"/"ATAC")
    backup_assay <- SeuratObject::DefaultAssay(seurat_data)
    SeuratObject::DefaultAssay(seurat_data) <- assay                                        # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                        # safety measure

    doublets_data <- NULL                                                                   # to collect doublets
    identities <- base::unique(
        base::as.vector(as.character(SeuratObject::Idents(seurat_data)))
    )
    for (i in 1:length(identities)){                                                        # to run scDblFinder for each sample completely independently
        base::print(
            base::paste(
                "Searching for doublets in", assay, "assay of",
                identities[i], "dataset"
            )
        )
        if (assay == "RNA"){
            base::print(base::paste("Subsetting to nCount_RNA >= 200"))
            subsetted_data <- base::subset(
                seurat_data,
                idents=identities[i],
                subset=(nCount_RNA >= 200)                                                  # scDblFinder says that we need to have this filter for RNA
            )
        } else {
            subsetted_data <- base::subset(
                seurat_data,
                idents=identities[i]
            )
        }
        sce_data <- scDblFinder::scDblFinder(
            Seurat::as.SingleCellExperiment(subsetted_data, assay=assay),
            aggregateFeatures=ifelse(assay == "ATAC", TRUE, FALSE),                         # perform feature aggregation (recommended for ATAC)
            processing=ifelse(assay == "ATAC", "normFeatures", "default"),                  # using normFeatures for ATAC
            dbr=dbl_rate,                                                                   # if NULL scDblFinder uses default parameters
            dbr.sd=dbl_rate_sd,                                                             # if NULL scDblFinder uses default parameters
            verbose=TRUE
        )
        sce_data <- base::as.data.frame(SingleCellExperiment::colData(sce_data)) %>%        # somehow fails if not explicitely converted to data.frame
                    dplyr::select(scDblFinder.class) %>%
                    dplyr::rename(!!target_column:="scDblFinder.class") %>%
                    tibble::rownames_to_column(var="barcode")                               # need to have this column for left_join
        if (is.null(doublets_data)) {
            doublets_data <- sce_data
        } else {
            doublets_data <- base::rbind(doublets_data, sce_data)
        }
    }
    doublets_data <- base::data.frame(SeuratObject::Cells(seurat_data)) %>%                 # create a dataframe with only one column and all cells
                     dplyr::rename("barcode"=1) %>%                                         # rename that column to "barcode"
                     dplyr::left_join(doublets_data, by="barcode") %>%                      # intersect by "barcode"
                     tibble::remove_rownames() %>%
                     tibble::column_to_rownames("barcode") %>%
                     replace(is.na(.), "singlet") %>%                                       # all not doublets (including cells with nCount_RNA < 200 for RNA) will be "singlet"
                     dplyr::mutate(                                                         # need to use factor for a proper order of levels on the plots
                         !!target_column:=base::factor(
                             .[[target_column]],                                            # to take only one column by the name from the variable
                             levels=base::sort(                                             # we can have only one category, that's why we need sort
                                 base::unique(.[[target_column]]),                          # to take only one column by the name from the variable
                                 decreasing=TRUE                                            # will result in "doublet", "singlet" order if both are present
                             )
                         )
                     )
    seurat_data <- SeuratObject::AddMetaData(
        seurat_data,
        doublets_data[SeuratObject::Cells(seurat_data), , drop=FALSE]                       # to guarantee the proper cells order
    )
    base::rm(doublets_data, identities)
    base::gc(verbose=FALSE)
    SeuratObject::DefaultAssay(seurat_data) <- backup_assay
    return (seurat_data)
}

add_rna_qc_metrics <- function(seurat_data, args){
    backup_assay <- SeuratObject::DefaultAssay(seurat_data)
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                        # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                        # safety measure
    seurat_data$log10_gene_per_log10_umi <- log10(seurat_data$nFeature_RNA) / log10(seurat_data$nCount_RNA)
    seurat_data$mito_percentage <- Seurat::PercentageFeatureSet(seurat_data, pattern=args$mitopattern)
    SeuratObject::DefaultAssay(seurat_data) <- backup_assay
    return (seurat_data)
}

add_atac_qc_metrics <- function(seurat_data, args){
    backup_assay <- SeuratObject::DefaultAssay(seurat_data)
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"
    seurat_data <- Signac::NucleosomeSignal(seurat_data, verbose=FALSE)
    # if 'gene_biotype' are all NAs, then annotation doesn't have real gene_biotype data and we need to use NULL
    tss_positions <- Signac::GetTSSPositions(
        ranges=Signac::Annotation(seurat_data[["ATAC"]]),
        biotypes=if(all(is.na(Signac::Annotation(seurat_data[["ATAC"]])$gene_biotype))) NULL else "protein_coding"
    )
    seurat_data <- Signac::TSSEnrichment(
        seurat_data,
        tss.positions=tss_positions,
        fast=FALSE,                                                                         # set fast=FALSE, because we want to build TSS Enrichment plot later
        verbose=FALSE
    )
    base::gc(verbose=FALSE)
    SeuratObject::DefaultAssay(seurat_data) <- backup_assay
    return (seurat_data)
}

add_peak_qc_metrics <- function(seurat_data, blacklist_data, args){
    backup_assay <- SeuratObject::DefaultAssay(seurat_data)
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"
    fragments_data <- Signac::CountFragments(
        fragments=args$fragments,
        cells=base::colnames(seurat_data),                                                     # limit it to only those cells that are present in seurat_data
        verbose=FALSE
    ) %>% tibble::column_to_rownames("CB")                                                     # for easy access to cell barcodes
    seurat_data$fragments <- fragments_data[base::colnames(seurat_data), "frequency_count"]    # select by rownames to make sure the cells order wasn't accidentally changed
    base::rm(fragments_data)                                                                   # remove unused data
    seurat_data <- Signac::FRiP(
        seurat_data,
        assay="ATAC",                                                                          # FRiP can't take the default assay, so we set it explicitly
        total.fragments="fragments",
        col.name="frip",
        verbose=FALSE
    )
    if (!is.null(blacklist_data)){
        seurat_data$blacklist_fraction <- Signac::FractionCountsInRegion(
            seurat_data,
            assay="ATAC",
            regions=blacklist_data
        )
    } else {
        seurat_data$blacklist_fraction <- 0                                                  # blacklist regions file wasn't provided, so we set everything to 0
    }
    SeuratObject::DefaultAssay(seurat_data) <- backup_assay
    base::gc(verbose=FALSE)
    return (seurat_data)
}

add_gene_expr_percentage <- function(seurat_data, target_genes){
    backup_assay <- SeuratObject::DefaultAssay(seurat_data)
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"
    for (i in 1:length(target_genes)){
        current_gene <- target_genes[i]
        seurat_data[[base::paste("perc", current_gene, sep="_")]] <- Seurat::PercentageFeatureSet(
            seurat_data,
            pattern=base::paste0("^", current_gene, "$")
        )
    }
    SeuratObject::DefaultAssay(seurat_data) <- backup_assay
    base::gc(verbose=FALSE)
    return (seurat_data)
}

# DEPRECATED as we now use standard Signac::GetTSSPositions function with biotypes=NULL when all gene_biotype are NA
get_tss_positions <- function(annotation_ranges){
    # Based on GetTSSPositions function from signac/R/utilities.R
    # adapted to work with refgene GTF annotations file
    annotation_df <- data.table::as.data.table(x=annotation_ranges)
    annotation_df$strand <- as.character(x=annotation_df$strand)
    annotation_df$strand <- base::ifelse(
        test = annotation_df$strand == "*",
        yes = "+",
        no = annotation_df$strand
    )
    collapsed_annotation_df <- annotation_df[
        , .(base::unique(seqnames),
            min(start),
            max(end),
            strand[[1]],
            gene_name[[1]]),
        "gene_id"
    ]
    base::colnames(x=collapsed_annotation_df) <- c(
        "gene_id", "seqnames", "start", "end", "strand", "gene_name"
    )
    collapsed_annotation_df$gene_name <- base::make.unique(names=collapsed_annotation_df$gene_name)
    collapsed_ranges <- GenomicRanges::makeGRangesFromDataFrame(
        df=collapsed_annotation_df,
        keep.extra.columns=TRUE
    )
    tss_positions <- IRanges::resize(collapsed_ranges, width=1, fix="start")
    return (tss_positions)
}

quartile_qc_metrics <- function(seurat_data, features, prefix="quartile"){
    for (i in 1:length(features)){
        current_feature <- features[i]
        base::tryCatch(
            expr = {
                quartiles <- stats::quantile(seurat_data@meta.data[, current_feature], c(0.25, 0.5, 0.75))
                seurat_data <- SeuratObject::AddMetaData(
                    object=seurat_data,
                    metadata=base::cut(
                        seurat_data@meta.data[, current_feature],
                        breaks=c(-Inf, quartiles[1], quartiles[2], quartiles[3], Inf), 
                        labels=c("Low", "Medium", "Medium high", "High")
                    ),
                    col.name=base::paste(prefix, current_feature, sep="_")
                )
            },
            error = function(e){
                base::print(base::paste("Failed to quantile ", current_feature, " with error - ", e, sep=""))
            }
        )
    }
    return (seurat_data)
}