import("dplyr", attach=FALSE)
import("limma", attach=FALSE)
import("DAseq", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("sceasy", attach=FALSE)
import("hopach", attach=FALSE)
import("DESeq2", attach=FALSE)
import("MAnorm2", attach=FALSE)
import("harmony", attach=FALSE)
import("tibble", attach=FALSE)
import("glmGamPoi", attach=FALSE)  # safety measure. we don't use it directly, but SCTransform with method="glmGamPoi" needs it
import("S4Vectors", attach=FALSE)
import("reticulate", attach=FALSE)
import("tidyselect", attach=FALSE)
import("rtracklayer", attach=FALSE)
import("BiocParallel", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)
import("SummarizedExperiment", attach=FALSE)

export(
    "rna_analyze",
    "add_clusters",
    "integrate_labels",
    "rna_preprocess",
    "rna_log_single",
    "rna_sct_single",
    "rna_log_integrated",
    "rna_sct_integrated",
    "get_vars_to_regress",
    "get_cell_cycle_scores",
    "get_min_ident_size",
    "get_markers_by_res",
    "get_markers",
    "atac_preprocess",
    "atac_analyze",
    "add_wnn_clusters",
    "rna_de_analyze",
    "atac_dbinding_analyze",
    "da_analyze",
    "get_de_sample_data",
    "get_de_cell_data",
    "get_bulk_counts_data",
    "get_clustered_data",
    "get_aggregated_expession"
)

get_tf_idf_method <- function(method_name){
    return (
        switch(
            method_name,
            "log-tfidf"    = 1,
            "tf-logidf"    = 2,
            "logtf-logidf" = 3,
            "idf"          = 4
        )
    )
}

get_cluster_algorithm <- function(algorithm_name){
    return (
        switch(
            algorithm_name,
            "louvain"      = 1,
            "mult-louvain" = 2,
            "slm"          = 3,
            "leiden"       = 4
        )
    )
}

get_vars_to_regress <- function(seurat_data, args, exclude_columns=NULL) {
    vars_to_regress <- NULL
    arguments <- c(args$regressmt, args$regressccfull, args$regressccdiff)                  # any of then can be also NULL
    metadata_columns <- c("mito_percentage", "S.Score&G2M.Score", "CC.Difference")
    for (i in 1:length(arguments)) {
        current_argument <- arguments[i]
        current_column <- metadata_columns[i]
        if ( is.null(current_argument) || (current_column %in% exclude_columns) ) {
            next
        }
        current_column <- base::unlist(base::strsplit(metadata_columns[i], "&"))
        if ( !all(current_column %in% base::colnames(seurat_data@meta.data)) ){             # the column doesn't exists in metadata
            next
        }
        if (current_argument) {
            if (is.null(vars_to_regress)) {
                vars_to_regress <- current_column
            } else {
                vars_to_regress <- base::append(vars_to_regress, current_column)
            }
        }
    }
    if (!is.null(args$regressgenes) && length(args$regressgenes) > 0){                      # easier to process regressgenes separately
        for (i in 1:length(args$regressgenes)){
            current_column <- base::paste("perc", args$regressgenes[i], sep="_")
            if (is.null(vars_to_regress)) {
                vars_to_regress <- current_column
            } else {
                vars_to_regress <- base::append(vars_to_regress, current_column)
            }
        }
    }
    return (vars_to_regress)
}

clean_cell_cycle_scores <- function(seurat_data){
    for (target_column in c("S.Score", "G2M.Score", "CC.Difference", "Phase")){
        if (target_column %in% base::colnames(seurat_data@meta.data)){
            base::print(
                base::paste(
                    "Removing", target_column, "from",
                    "the Seurat object metadata"
                )
            )
            seurat_data[[target_column]] <- NULL
        }
    }
    return (seurat_data)
}

get_cell_cycle_scores <- function(seurat_data, assay, cell_cycle_data){       # we need this function to fail if something went wrong, so no tryCatch inside
    SeuratObject::DefaultAssay(seurat_data) <- assay                          # safety measure
    seurat_data <- Seurat::CellCycleScoring(                                  # calls AddModuleScore that uses data slot
        seurat_data,
        s.features=base::as.vector(cell_cycle_data[base::tolower(cell_cycle_data$phase)=="s", "gene_id"]),
        g2m.features=base::as.vector(cell_cycle_data[base::tolower(cell_cycle_data$phase)=="g2/m", "gene_id"]),
        assay=assay,
        verbose=FALSE
    )
    seurat_data[["CC.Difference"]] <- seurat_data[["S.Score"]] - seurat_data[["G2M.Score"]]   # for softer cell cycle removal (https://satijalab.org/seurat/articles/cell_cycle_vignette.html)
    return (seurat_data)
}

rna_log_single <- function(seurat_data, args, cell_cycle_data=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                                # safety measure
    base::print("Applying LogNormalize")
    seurat_data <- Seurat::NormalizeData(seurat_data, verbose=FALSE)
    if(!is.null(cell_cycle_data)){
        base::tryCatch(
            expr = {
                base::print("Trying to assign cell cycle scores for RNA assay")
                seurat_data <- clean_cell_cycle_scores(seurat_data)                # removes S.Score, G2M.Score, and CC.Difference columns
                seurat_data <- get_cell_cycle_scores(                              # if succeded adds S.Score and G2M.Score, and CC.Difference columns
                    seurat_data,
                    "RNA",
                    cell_cycle_data
                )
            },
            error = function(e){
                base::print(base::paste("Failed to run cell cycle scoring for RNA assay with error - ", e))
            }
        )
    }
    seurat_data <- Seurat::FindVariableFeatures(
        seurat_data,
        nfeatures=args$highvargenes,
        verbose=FALSE
    )
    vars_to_regress <- get_vars_to_regress(seurat_data, args)                                   # may or may not include S.Score, G2M.Score, and CC.Difference columns
    base::print(base::paste0("Regressing out [", paste(vars_to_regress, collapse=", "), "]"))
    seurat_data <- Seurat::ScaleData(
        seurat_data,
        vars.to.regress=vars_to_regress,
        verbose=FALSE
    )
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"
    return (seurat_data)
}

sc_transform_helper <- function(seurat_data, args, vars_to_regress, method){
    base::print(
        base::paste0(
            "Applying SCTransform using ", method,
            " method for initial parameter estimation. ",
            "Regressing out [", paste(vars_to_regress, collapse=", "), "]"
        )
    )
    seurat_data <- Seurat::SCTransform(
        seurat_data,
        assay="RNA",
        new.assay.name="SCT",
        variable.features.n=args$highvargenes,
        method=method,
        vars.to.regress=vars_to_regress,
        conserve.memory=args$lowmem,
        verbose=FALSE
    )
    return (seurat_data)
}

rna_sct_single <- function(seurat_data, args, cell_cycle_data=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                # safety measure
    method <- base::ifelse(args$norm=="sctglm", "glmGamPoi", "poisson")
    if(!is.null(cell_cycle_data)){                                                  # if we are planing to assign/re-assign cell cycle scores, the first
        seurat_data <- clean_cell_cycle_scores(seurat_data)                         # SCTransform should not attempt to regress it. That's why we remove
    }                                                                               # S.Score, G2M.Score, and CC.Difference columns from the very beginning
    vars_to_regress <- get_vars_to_regress(seurat_data, args)                       # may or may not include S.Score, G2M.Score, and CC.Difference columns
    seurat_data <- sc_transform_helper(
        seurat_data=seurat_data,
        args=args,
        vars_to_regress=vars_to_regress,
        method=method
    )
    if(!is.null(cell_cycle_data)){
        base::tryCatch(
            expr = {
                base::print("Trying to assign cell cycle scores for SCT assay")
                seurat_data <- get_cell_cycle_scores(                               # if succeded adds S.Score and G2M.Score, and CC.Difference columns
                    seurat_data,                                                    # this object has cell cycle scores columns already removed above
                    "SCT",
                    cell_cycle_data
                )
                vars_to_regress <- get_vars_to_regress(seurat_data, args)           # may or may not include S.Score, G2M.Score, and CC.Difference
                if (                                                                # need to rerun SCTransform to regress all variables at once
                    all( c("S.Score", "G2M.Score") %in% vars_to_regress ) ||        # complete cell cycle genes removal
                    "CC.Difference" %in% vars_to_regress                            # partial cell cycle genes removal
                ){
                    seurat_data <- sc_transform_helper(
                        seurat_data=seurat_data,
                        args=args,
                        vars_to_regress=vars_to_regress,
                        method=method
                    )
                }
            },
            error = function(e){
                base::print(base::paste("Failed to run cell cycle scoring/regressing for SCT assay with error - ", e))
            }
        )
    }
    SeuratObject::DefaultAssay(seurat_data) <- "SCT"
    return (seurat_data)
}

get_min_ident_size <- function(splitted_seurat_data){
    min_ident_size <- min(
        table(
            unlist(
                lapply(
                    lapply(
                        lapply(
                            splitted_seurat_data,
                            SeuratObject::Idents
                        ),
                        as.character
                    ),
                    base::as.vector
                )
            )
        )
    )
    return (min_ident_size)
}

rna_log_integrated <- function(splitted_seurat_data, args, cell_cycle_data=NULL){
    failed_cell_cycle_scoring <- FALSE

    for (i in 1:length(splitted_seurat_data)){
        base::print(base::paste("Processing", SeuratObject::Idents(splitted_seurat_data[[i]])[1], "dataset"))

        SeuratObject::DefaultAssay(splitted_seurat_data[[i]]) <- "RNA"                            # safety measure
        base::print(base::paste("Applying LogNormalize"))
        splitted_seurat_data[[i]] <- Seurat::NormalizeData(
            splitted_seurat_data[[i]],
            verbose=FALSE
        )
        if(!is.null(cell_cycle_data)){
            if (!failed_cell_cycle_scoring){
                base::tryCatch(
                    expr = {
                        base::print(
                            base::paste(
                                "Trying to assign cell cycle scores for RNA assay of",
                                SeuratObject::Idents(splitted_seurat_data[[i]])[1], "dataset"
                            )
                        )
                        splitted_seurat_data[[i]] <- clean_cell_cycle_scores(splitted_seurat_data[[i]])
                        splitted_seurat_data[[i]] <- get_cell_cycle_scores(                       # if succeded adds S.Score and G2M.Score, and CC.Difference columns
                            splitted_seurat_data[[i]],
                            "RNA",
                            cell_cycle_data
                        )
                    },
                    error = function(e){
                        base::print(base::paste("Failed to run cell cycle scoring for RNA assay with error - ", e))
                        failed_cell_cycle_scoring <- TRUE
                    }
                )
            } else {
                base::print(
                    base::paste(
                        "Skip cell cycle scores assignment for RNA assay of",
                        SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                        "dataset due to failure in one of the other datasets."
                    )
                )
            }
        }
    }

    if (failed_cell_cycle_scoring){
        base::print(
            base::paste(
                "At least one of the datasets failed when",
                "running cell cycle score assignment.",
                "Cleaning up cell cycle scores."
            )
        )
        for (i in 1:length(splitted_seurat_data)){
            splitted_seurat_data[[i]] <- clean_cell_cycle_scores(splitted_seurat_data[[i]])
        }
    }

    for (i in 1:length(splitted_seurat_data)){                          # don't know if cell cycle scores may somehow influence on the Variable features, that's
        splitted_seurat_data[[i]] <- Seurat::FindVariableFeatures(      # why we run in only after either all datasets have cell cycle scores or neither of them
            splitted_seurat_data[[i]],
            nfeatures=args$highvargenes,
            verbose=FALSE
        )
    }

    integration_features <- Seurat::SelectIntegrationFeatures(
        splitted_seurat_data,
        nfeatures=args$highvargenes,
        verbose=FALSE
    )
    integration_anchors <- Seurat::FindIntegrationAnchors(
        splitted_seurat_data,
        normalization.method="LogNormalize",
        anchor.features=integration_features,
        verbose=FALSE
    )
    integrated_seurat_data <- Seurat::IntegrateData(
        integration_anchors, 
        new.assay.name="rna_integrated",
        normalization.method="LogNormalize",
        k.weight=min(get_min_ident_size(splitted_seurat_data), 100),                            # k.weight 100 by default, but shouldn't be bigger than 
        verbose=FALSE                                                                           # the min number of cells among all identities after filtering
    )
    vars_to_regress <- get_vars_to_regress(integrated_seurat_data, args)                        # may or may not include S.Score, G2M.Score, and CC.Difference columns
    base::print(base::paste0("Regressing out [", paste(vars_to_regress, collapse=", "), "]"))
    integrated_seurat_data <- Seurat::ScaleData(
        integrated_seurat_data,
        vars.to.regress=vars_to_regress,
        verbose=FALSE
    )
    SeuratObject::DefaultAssay(integrated_seurat_data) <- "rna_integrated"
    base::rm(integration_features, integration_anchors)                                 # remove unused data
    base::gc(verbose=FALSE)
    return (integrated_seurat_data)
}

rna_sct_integrated <- function(splitted_seurat_data, args, cell_cycle_data=NULL){
    method <- base::ifelse(args$norm=="sctglm", "glmGamPoi", "poisson")
    failed_cell_cycle_scoring <- FALSE
    for (i in 1:length(splitted_seurat_data)) {
        base::print(base::paste("Processing", SeuratObject::Idents(splitted_seurat_data[[i]])[1], "dataset"))

        SeuratObject::DefaultAssay(splitted_seurat_data[[i]]) <- "RNA"                            # safety measure
        if(!is.null(cell_cycle_data)){                                                            # if we are planing to assign/re-assign cell cycle scores, the first
            splitted_seurat_data[[i]] <- clean_cell_cycle_scores(splitted_seurat_data[[i]])       # SCTransform should not attempt to regress it. That's why we remove
        }                                                                                         # S.Score, G2M.Score, and CC.Difference columns from the very beginning
        vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args)                   # may or may not include S.Score, G2M.Score, and CC.Difference columns
        splitted_seurat_data[[i]] <- sc_transform_helper(
            seurat_data=splitted_seurat_data[[i]],
            args=args,
            vars_to_regress=vars_to_regress,
            method=method
        )
        if(!is.null(cell_cycle_data)){
            if (!failed_cell_cycle_scoring){
                base::tryCatch(
                    expr = {
                        base::print(
                            base::paste(
                                "Trying to assign cell cycle scores for SCT assay of",
                                SeuratObject::Idents(splitted_seurat_data[[i]])[1], "dataset"
                            )
                        )
                        splitted_seurat_data[[i]] <- get_cell_cycle_scores(                              # if succeded adds S.Score and G2M.Score, and CC.Difference columns
                            splitted_seurat_data[[i]],                                                   # this object has cell cycle scores columns already removed above
                            "SCT",
                            cell_cycle_data
                        )
                        vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args)          # may or may not include S.Score, G2M.Score, and CC.Difference
                        if (                                                                # need to rerun SCTransform to regress all variables at once
                            all( c("S.Score", "G2M.Score") %in% vars_to_regress ) ||        # complete cell cycle genes removal
                            "CC.Difference" %in% vars_to_regress                            # partial cell cycle genes removal
                        ){
                            splitted_seurat_data[[i]] <- sc_transform_helper(
                                seurat_data=splitted_seurat_data[[i]],
                                args=args,
                                vars_to_regress=vars_to_regress,
                                method=method
                            )
                        }
                    },
                    error = function(e){
                        base::print(base::paste("Failed to run cell cycle scoring for SCT assay with error - ", e))
                        failed_cell_cycle_scoring <- TRUE
                    }
                )
            } else {
                base::print(
                    base::paste(
                        "Skip cell cycle scores assignment for SCT assay of",
                        SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                        "dataset due to failure in one of the other datasets."
                    )
                )
            }
        }
    }

    if (failed_cell_cycle_scoring){
        base::print(
            base::paste(
                "At least one of the datasets failed when",
                "running cell cycle score assignment.",
                "Cleaning up cell cycle scores."
            )
        )
        for (i in 1:length(splitted_seurat_data)){
            vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args)               # may or may not include S.Score, G2M.Score, and CC.Difference
            splitted_seurat_data[[i]] <- clean_cell_cycle_scores(splitted_seurat_data[[i]])       # remove cell cycle score columns after we know what we attempted to regress
            if (                                                                                  # it means we attempted to regress cell cycle scores for this dataset
                all( c("S.Score", "G2M.Score") %in% vars_to_regress ) ||                          # complete cell cycle genes removal
                "CC.Difference" %in% vars_to_regress                                              # partial cell cycle genes removal
            ){
                vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args)           # here it won't include S.Score, G2M.Score, and CC.Difference columns anymore
                splitted_seurat_data[[i]] <- sc_transform_helper(
                    seurat_data=splitted_seurat_data[[i]],
                    args=args,
                    vars_to_regress=vars_to_regress,
                    method=method
                )
            }
        }
    }

    integration_features <- Seurat::SelectIntegrationFeatures(
        splitted_seurat_data,
        nfeatures=args$highvargenes,
        verbose=FALSE
    )
    splitted_seurat_data <- Seurat::PrepSCTIntegration(
        splitted_seurat_data, 
        anchor.features=integration_features,
        verbose=FALSE
    )
    integration_anchors <- Seurat::FindIntegrationAnchors(
        splitted_seurat_data,
        normalization.method="SCT",
        anchor.features=integration_features,
        verbose=FALSE
    )
    integrated_seurat_data <- Seurat::IntegrateData(
        integration_anchors, 
        new.assay.name="rna_integrated",
        normalization.method="SCT",
        k.weight=min(get_min_ident_size(splitted_seurat_data), 100),        # k.weight 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
        verbose=FALSE
    )
    SeuratObject::DefaultAssay(integrated_seurat_data) <- "rna_integrated"
    base::rm(integration_features, integration_anchors)                         # remove unused data
    base::gc(verbose=FALSE)
    return (integrated_seurat_data)
}

rna_preprocess <- function(seurat_data, args, cell_cycle_data=NULL) {
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                                        # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                                        # safety measure
    splitted_seurat_data <- Seurat::SplitObject(seurat_data, split.by="new.ident")                          # to check if we have aggregated datasets
    if (args$ntgr == "none" | args$ntgr == "harmony" | length(splitted_seurat_data) == 1){
        base::print(
            base::paste(
                "Skipping datasets integration (either forced to skip, or will be attempted",
                "to run later with harmony, or only one identity is present). Using the",
                "original not splitted seurat data."
            )
        )
        if (args$norm == "log"){
            processed_seurat_data <- rna_log_single(seurat_data, args, cell_cycle_data)                     # sets default assay to RNA
        } else {
            processed_seurat_data <- rna_sct_single(seurat_data, args, cell_cycle_data)                     # sets default assay to SCT
        }
    } else {
        base::print("Running datasets integration using Seurat on splitted data.")
        if (args$norm == "log"){
            processed_seurat_data <- rna_log_integrated(splitted_seurat_data, args, cell_cycle_data)        # sets default assay to rna_integrated
        } else {
            processed_seurat_data <- rna_sct_integrated(splitted_seurat_data, args, cell_cycle_data)        # sets default assay to rna_integrated
        }
    }
    base::rm(splitted_seurat_data)
    base::gc(verbose=FALSE)
    return (processed_seurat_data)
}

rna_analyze <- function(seurat_data, args, cell_cycle_data=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    backup_reductions <- c()                                                    # RNA integration main remove atac related reductions so we need to back them up
    reduction_names <- c(
        "atac_lsi",
        "atacumap",
        "wnnumap",
        "gene_rnaumi",
        "rnaumi_atacfrgm",
        "tss_atacfrgm"
    )
    if (is.null(cell_cycle_data)){                                              # we are not planning to run cell cycle score assignment, so we want to keep "ccpca"
        reduction_names <- append(reduction_names, "ccpca")
    }
    for (reduction_name in reduction_names){
        if (reduction_name %in% names(seurat_data@reductions)){
            base::print(base::paste("Backing up reduction", reduction_name))
            backup_reductions[[reduction_name]] <- seurat_data[[reduction_name]]
        }
    }
    seurat_data <- rna_preprocess(seurat_data, args, cell_cycle_data)           # sets "rna_integrated" as a default assay for integrated data, and either "RNA" or "SCT" for not integrated data
    if (length(backup_reductions) > 0){                                         # restoring backed up reductions
        for (reduction_name in names(backup_reductions)){
            base::print(base::paste("Restoring reduction", reduction_name, "from backup"))
            seurat_data[[reduction_name]] <- backup_reductions[[reduction_name]]
        }
    }
    if (!is.null(cell_cycle_data) && all(c("S.Score", "G2M.Score", "CC.Difference", "Phase") %in% base::colnames(seurat_data@meta.data))){
        used_features <- base::unique(
            cell_cycle_data$gene_id[
                cell_cycle_data$gene_id %in% base::rownames(SeuratObject::GetAssayData(seurat_data, slot="scale.data"))
            ]
        )
        base::print(
            base::paste(
                "Performing PCA reduction on",
                SeuratObject::DefaultAssay(seurat_data), "assay using 50 principal",
                "components and", length(used_features), "out of", length(cell_cycle_data$gene_id),
                "cell cycle genes"
            )
        )
        seurat_data <- Seurat::RunPCA(            # PCA will be run on the default assay, which can be one of "rna_integrated", "RNA", or "SCT"
            seurat_data,
            reduction.name="ccpca",               # add "ccpca" reduction that can be used for cell cycle evaluation
            npcs=50,
            features=cell_cycle_data$gene_id,     # these features will be overlaped with the features from the scale.data slot
            verbose=FALSE,                        # which means they can be either highly variable genes or integration features if we use rna_integrated assay
            approx=FALSE                          # to avoid warning message
        )
    }
    base::print(
        base::paste(
            "Performing PCA reduction on", SeuratObject::DefaultAssay(seurat_data),
            "assay using 50 principal components for Elbow and QC correlation plots, and",
            max(args$dimensions), "principal components for UMAP projection"
        )
    )
    seurat_data <- Seurat::RunPCA(                # add "qcpca" reduction to be used in elbow and QC correlation plots
        seurat_data,
        npcs=50,
        reduction.name="qcpca",
        verbose=FALSE
    )
    seurat_data <- Seurat::RunPCA(                # add "pca" reduction to be used in UMAP and Harmony integration
        seurat_data,
        npcs=max(args$dimensions),                # need to take the max because dimensions was already adjusted to array of values
        reduction.name="pca",
        verbose=FALSE
    )
    if (args$ntgr == "harmony"){
        if (is.null(args$ntgrby) || length(base::unique(base::as.vector(as.character(SeuratObject::Idents(seurat_data))))) == 1){
            base::print(
                base::paste(
                    "Skipping datasets integration with Harmony. Either --ntgrby",
                    "wasn't provided or data included only a single dataset."
                )
            )
        } else {
            base::print(
                base::paste(
                    "Running datasets integration with harmony using",
                    SeuratObject::DefaultAssay(seurat_data), "assay.",
                    "Integrating over", base::paste(args$ntgrby, collapse=", "),
                    "covariates. Dimensions used:", base::paste(args$dimensions, collapse=", ")
                )
            )
            seurat_data <- harmony::RunHarmony(
                object=seurat_data,
                group.by.vars=args$ntgrby,
                reduction="pca",
                reduction.save="pca",                                  # overwriting old pca reduction
                dims.use=args$dimensions,                              # keep it, but it looks like it doesn't influence anything https://github.com/immunogenomics/harmony/issues/82
                assay.use=SeuratObject::DefaultAssay(seurat_data),     # can be both RNA or SCT depending on --norm parameter
                verbose=FALSE
            )
        }
    }
    base::print(
        base::paste(
            "Performing UMAP projection on 'pca' reduction.",
            "Dimensions used:", base::paste(args$dimensions, collapse=", ")
        )
    )
    base::suppressWarnings(                                            # to not show warning about python version of UMAP
        seurat_data <- Seurat::RunUMAP(
            seurat_data,
            reduction="pca",
            dims=args$dimensions,
            reduction.name="rnaumap",
            reduction.key="RNAUMAP_",
            spread=base::ifelse(is.null(args$uspread), 1, args$uspread),
            min.dist=base::ifelse(is.null(args$umindist), 0.3, args$umindist),
            n.neighbors=base::ifelse(is.null(args$uneighbors), 30, args$uneighbors),
            metric=base::ifelse(is.null(args$umetric), "cosine", args$umetric),
            umap.method=base::ifelse(is.null(args$umethod), "uwot", args$umethod),
            verbose=FALSE
        )
    )
    return (seurat_data)
}

add_clusters <- function(seurat_data, assay, graph_name, reduction, args){
    SeuratObject::DefaultAssay(seurat_data) <- assay                                  # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                  # safety measure
    seurat_data <- Seurat::FindNeighbors(
        seurat_data,
        annoy.metric=base::ifelse(is.null(args$ametric), "euclidean", args$ametric),
        reduction=reduction,
        dims=args$dimensions,
        graph.name=base::paste(graph_name, c("_nn", ""), sep=""),
        verbose=FALSE
    )
    seurat_data <- Seurat::FindClusters(
        seurat_data,
        resolution=args$resolution,
        graph.name=graph_name,
        algorithm=get_cluster_algorithm(args$algorithm),
        verbose=FALSE
    )
    return (seurat_data)
}

integrate_labels <- function(seurat_data, source_columns, args){
    base::print(
        base::paste(
            "Running scTriangulate for", base::paste(source_columns, collapse=", "), "columns.",
            base::ifelse(
                !is.null(args$target),
                base::paste("The results will be saved into the columns with the suffix", args$target),
                ""
            )
        )
    )
    temporary_file <- base::tempfile(                                      # will be automatically removed when R exits
        pattern="sctri",
        tmpdir=base::tempdir(),
        fileext=".h5ad"
    )
    base::print(base::paste("Saving temporary h5ad file to", temporary_file))
    sceasy::convertFormat(seurat_data, from="seurat", to="anndata", outFile=temporary_file)
    script_file <- base::tempfile(                                         # will be automatically removed when R exits
        pattern="sctri",
        tmpdir=base::tempdir(),
        fileext=".py"
    )
    output_stream <- base::file(script_file)
    base::writeLines(
        c(
            "import os",
            "import scanpy",
            "import resource",
            "import sctriangulate",
            "try:",
            "    R_MAX_VSIZE = int(os.getenv('R_MAX_VSIZE'))",
            "    resource.setrlimit(resource.RLIMIT_AS, (R_MAX_VSIZE, R_MAX_VSIZE))",         # ignored if run not on Linux
            "    print(f'Attempting to set the maximum memory limits to {R_MAX_VSIZE}')",
            "except Exception:",
            "    print('Failed to set maximum memory limits')",
            "def sc_triangulate(location, clustering_columns, tmp_dir, cores=None):",
            "    cores = 1 if cores is None else int(cores)",
            "    sctri_data = sctriangulate.ScTriangulate(",
            "        dir=tmp_dir,",
            "        adata=scanpy.read(location),",
            "        query=clustering_columns,",
            "        predict_doublet=False",
            "    )",
            "    sctri_data.compute_metrics(",
            "        cores=cores,",
            "        scale_sccaf=True",
            "    )",
            "    sctri_data.compute_shapley(cores=cores)",
            "    sctri_data.prune_result()",
            "    return sctri_data.adata.obs"
        ),
        output_stream
    )
    base::close(output_stream)
    reticulate::source_python(script_file)
    pruned_clusters <- sc_triangulate(
                           temporary_file,
                           source_columns,
                           base::tempdir(),
                           args$cpus
                       )[, c("pruned", "confidence", "final_annotation")] %>%
                       dplyr::mutate("pruned"=base::gsub("@", "__", .$pruned)) %>%  # @ is not good if it somehow appears in any of the filenames
                       dplyr::rename(
                           !!tidyselect::all_of(base::paste("custom", args$target, "pruned", sep="_")):="pruned"       # need "custom" prefix for UCSC Browser
                       ) %>%
                       dplyr::rename(
                           !!tidyselect::all_of(base::paste("custom", args$target, "confidence", sep="_")):="confidence"
                       ) %>%
                       dplyr::rename(
                           !!tidyselect::all_of(base::paste("custom", args$target, "final_annotation", sep="_")):="final_annotation"
                       )
    base::print(utils::head(pruned_clusters))
    seurat_data <- SeuratObject::AddMetaData(
        seurat_data,
        pruned_clusters[SeuratObject::Cells(seurat_data), , drop=FALSE]    # to guarantee the proper cells order
    )
    base::rm(pruned_clusters)  # remove unused data
    base::gc(verbose=FALSE)
    return (seurat_data)
}

add_wnn_clusters <- function(seurat_data, graph_name, reductions, dimensions, args){
    SeuratObject::Idents(seurat_data) <- "new.ident"                                   # safety measure
    seurat_data <- Seurat::FindMultiModalNeighbors(
        seurat_data,
        reduction.list=reductions,                                       # list("pca", "atac_lsi"),
        dims.list=dimensions,
        snn.graph.name=graph_name,                                       # "wsnn"
        weighted.nn.name="weighted.nn",
        verbose=FALSE
    )
    seurat_data <- Seurat::RunUMAP(
        seurat_data,
        nn.name="weighted.nn",
        reduction.name="wnnumap",
        reduction.key="WNNUMAP_",
        spread=base::ifelse(is.null(args$uspread), 1, args$uspread),
        min.dist=base::ifelse(is.null(args$umindist), 0.3, args$umindist),
        n.neighbors=base::ifelse(is.null(args$uneighbors), 30, args$uneighbors),
        metric=base::ifelse(is.null(args$umetric), "cosine", args$umetric),
        umap.method=base::ifelse(is.null(args$umethod), "uwot", args$umethod),
        verbose=FALSE
    )
    seurat_data <- Seurat::FindClusters(
        seurat_data,
        graph.name=graph_name,
        algorithm=get_cluster_algorithm(args$algorithm),
        resolution=args$resolution,
        verbose=FALSE
    )
    return (seurat_data)
}

get_markers <- function(seurat_data, assay, group_by, args, latent_vars=NULL, min_diff_pct=-Inf){
    SeuratObject::DefaultAssay(seurat_data) <- assay                            # safety measure
    SeuratObject::Idents(seurat_data) <- group_by
    markers <- NULL
    base::tryCatch(
        expr = {
            markers <- Seurat::FindAllMarkers(
                seurat_data,
                logfc.threshold=base::ifelse(is.null(args$logfc), 0.25, args$logfc),
                min.pct=base::ifelse(is.null(args$minpct), 0.1, args$minpct),
                only.pos=base::ifelse(is.null(args$onlypos), FALSE, args$onlypos),
                test.use=base::ifelse(is.null(args$testuse), "wilcox", args$testuse),
                min.diff.pct=min_diff_pct,
                latent.vars=latent_vars,
                verbose=FALSE
            ) %>% dplyr::relocate(cluster, gene, .before=1) %>% dplyr::rename("feature"="gene")
            if (base::nrow(markers) == 0){
                return (NULL)                                                   # to be able to check later just with is.null
            }
        },
        error = function(e){
            base::print(base::paste("Failed to identify markers for cells grouped by", group_by, "with error -", e))
        },
        finally = {
            SeuratObject::Idents(seurat_data) <- "new.ident"                    # safety measure
        }
    )
    return (markers)
}

get_markers_by_res <- function(seurat_data, assay, resolution_prefix, args, latent_vars=NULL, min_diff_pct=-Inf){
    SeuratObject::DefaultAssay(seurat_data) <- assay                # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                # safety measure
    all_putative_markers <- NULL
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        markers <- get_markers(
            seurat_data=seurat_data,
            assay=assay,
            group_by=base::paste(resolution_prefix, resolution, sep="."),
            args=args,
            latent_vars=latent_vars,
            min_diff_pct=min_diff_pct
        )
        if (!is.null(markers)){
            markers <- markers %>% base::cbind(resolution=resolution, .)
            if (!is.null(all_putative_markers)) {
                all_putative_markers <- base::rbind(all_putative_markers, markers)
            } else {
                all_putative_markers <- markers
            }
            base::rm(markers)                                       # remove unused data
        }
    }
    base::gc(verbose=FALSE)
    return (all_putative_markers)
}

atac_preprocess <- function(seurat_data, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"                                           # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                            # safety measure

    base::print(
        base::paste(
            "Applying TF-IDF normalization using", args$norm, "method.",
            "Searching for top highly variable features using", args$minvarpeaks,
            "as a lower percentile bound. Analyzing all datasets jointly."
        )
    )
    processed_seurat_data <- Signac::RunTFIDF(
        seurat_data,
        assay="ATAC",                                                                           # safety measure
        method=get_tf_idf_method(args$norm),
        verbose=FALSE
    )
    processed_seurat_data <- Signac::FindTopFeatures(
        processed_seurat_data,
        assay="ATAC",                                                                           # safety measure
        min.cutoff=args$minvarpeaks,
        verbose=FALSE
    )
    processed_seurat_data <- Signac::RunSVD(
        processed_seurat_data,
        n=50,
        reduction.name="qclsi",                                                                 # adding "qclsi" for QC correlation plots
        verbose=FALSE
    )

    splitted_seurat_data <- Seurat::SplitObject(seurat_data, split.by="new.ident")              # need to use original seurat_data for possible integration below
    if (args$ntgr == "none" | args$ntgr == "harmony" | length(splitted_seurat_data) == 1){
        base::print(
            base::paste(
                "Skipping datasets integration (either forced to skip, or will be attempted",
                "to run later with harmony, or only one identity is present). Using the",
                "original not splitted seurat data."
            )
        )
        processed_seurat_data <- Signac::RunSVD(
            processed_seurat_data,
            n=max(args$dimensions),                                                             # using the max because args$dimensions is already array
            reduction.name="atac_lsi",                                                          # adding "atac_lsi" for UMAP
            verbose=FALSE
        )
        if (args$ntgr == "harmony"){
            if (is.null(args$ntgrby) || length(splitted_seurat_data) == 1){
                base::print(
                    base::paste(
                        "Skipping datasets integration with Harmony. Either --ntgrby",
                        "wasn't provided or data included only a single dataset."
                    )
                )
            } else {
                base::print(
                    base::paste(
                        "Running datasets integration with harmony using",
                        SeuratObject::DefaultAssay(processed_seurat_data), "assay.",
                        "Integrating over", base::paste(args$ntgrby, collapse=", "),
                        "covariates. Dimensions used:", base::paste(args$dimensions, collapse=", ")
                    )
                )
                processed_seurat_data <- harmony::RunHarmony(
                    object=processed_seurat_data,
                    group.by.vars=args$ntgrby,
                    reduction="atac_lsi",
                    reduction.save="atac_lsi",                                    # overwriting old atac_lsi reduction
                    dims.use=args$dimensions,                                     # keed it, but it doesn't actually used anywhere inside RunHarmony function
                    assay.use=SeuratObject::DefaultAssay(processed_seurat_data),
                    project.dim=FALSE,                                            # for ATAC we don't need to project
                    verbose=FALSE
                )
            }
        }
    } else {
        base::print("Running datasets integration using Signac on splitted data.")
        processed_seurat_data <- Signac::RunSVD(
            processed_seurat_data,
            n=max(args$dimensions),                                                             # using the max because args$dimensions is already array
            reduction.name="lsi",                                                               # adding "lsi" as it will be used in integration
            verbose=FALSE
        )
        for (i in 1:length(splitted_seurat_data)){                                              # it was splitted from not updated seurat_data
            SeuratObject::DefaultAssay(splitted_seurat_data[[i]]) <- "ATAC"                     # safety measure
            base::print(
                base::paste(
                    "Applying TF-IDF normalization using", args$norm, "method.",
                    "Searching for top highly variable features using", args$minvarpeaks,
                    "as a lower percentile bound. Analyzing",
                    SeuratObject::Idents(splitted_seurat_data[[i]])[1], "dataset."
                )
            )
            splitted_seurat_data[[i]] <- Signac::RunTFIDF(
                splitted_seurat_data[[i]],
                assay="ATAC",                                                                   # safety measure
                method=get_tf_idf_method(args$norm),
                verbose=FALSE
            )
            splitted_seurat_data[[i]] <- Signac::FindTopFeatures(
                splitted_seurat_data[[i]],
                assay="ATAC",                                                                   # safety measure
                min.cutoff=args$minvarpeaks,
                verbose=FALSE
            )
            splitted_seurat_data[[i]] <- Signac::RunSVD(
                splitted_seurat_data[[i]],
                n=max(args$dimensions),                                                         # using the max because args$dimensions is already array
                reduction.name="lsi",                                                           # adding "lsi" to be used in FindIntegrationAnchors
                verbose=FALSE
            )
        }
        integration_anchors <- Seurat::FindIntegrationAnchors(
            splitted_seurat_data,
            anchor.features=base::rownames(seurat_data),                                        # peaks from our ATAC assay
            reduction="rlsi",                                                                   # will always search for "lsi" reductions in splitted_seurat_data
            dims=args$dimensions,                                                               # don't need to use more than we are going to use in UMAP
            verbose=FALSE
        )
        integrated_seurat_data <- Seurat::IntegrateEmbeddings(
            anchorset=integration_anchors,
            reductions=processed_seurat_data[["lsi"]],
            new.reduction.name="atac_lsi",                                                      # adding "atac_lsi" for consistency
            k.weight=min(get_min_ident_size(splitted_seurat_data), 100),                        # k.weight 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
            dims.to.integrate=args$dimensions                                                   # even if we include 1 dimension it won't influence on the results because we didn't use it in FindIntegrationAnchors
        )
        qclsi_reduction <- processed_seurat_data[["qclsi"]]                                     # need it for QC correlation plots
        processed_seurat_data <- integrated_seurat_data
        processed_seurat_data[["qclsi"]] <- qclsi_reduction
        base::rm(integration_anchors, integrated_seurat_data, qclsi_reduction)                  # remove unused data
    }
    base::rm(splitted_seurat_data)
    base::gc(verbose=FALSE)
    return (processed_seurat_data)
}

atac_analyze <- function(seurat_data, args){
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"                           # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    backup_reductions <- c()                                                    # ATAC integration main remove RNA related reductions so we need to back them up
    reduction_names <- c(
        "ccpca",
        "pca",
        "rnaumap",
        "wnnumap",
        "gene_rnaumi",
        "rnaumi_atacfrgm",
        "tss_atacfrgm"
    )
    for (reduction_name in reduction_names){
        if (reduction_name %in% names(seurat_data@reductions)){
            base::print(base::paste("Backing up reduction", reduction_name))
            backup_reductions[[reduction_name]] <- seurat_data[[reduction_name]]
        }
    }
    seurat_data <- atac_preprocess(seurat_data, args)                            # adds "atac_lsi" reduction
    if (length(backup_reductions) > 0){                                          # restoring backed up reductions
        for (reduction_name in names(backup_reductions)){
            base::print(base::paste("Restoring reduction", reduction_name, "from backup"))
            seurat_data[[reduction_name]] <- backup_reductions[[reduction_name]]
        }
    }
    seurat_data <- Seurat::RunUMAP(
        seurat_data,
        reduction="atac_lsi",
        dims=args$dimensions,
        reduction.name="atacumap",
        reduction.key="ATACUMAP_",
        spread=base::ifelse(is.null(args$uspread), 1, args$uspread),
        min.dist=base::ifelse(is.null(args$umindist), 0.3, args$umindist),
        n.neighbors=base::ifelse(is.null(args$uneighbors), 30, args$uneighbors),
        metric=base::ifelse(is.null(args$umetric), "cosine", args$umetric),
        umap.method=base::ifelse(is.null(args$umethod), "uwot", args$umethod),
        verbose=FALSE
    )
    return (seurat_data)
}

get_aggregated_expession <- function(seurat_data, group_by, selected_genes=NULL, slot="counts"){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"        # safety measure
    SeuratObject::Idents(seurat_data) <- group_by           # will be used by AggregateExpression because by default it's called with group.by="ident"
    aggregated_seurat_data <- Seurat::AggregateExpression(
        seurat_data,
        assays="RNA",                                       # need only RNA assay
        slot=slot,                                          # for slot="counts" no exponentiation is performed prior to aggregating
        features=selected_genes,                            # if NULL use all genes
        return.seurat=TRUE,                                 # summed values are saved in "counts", log-normalized - in "data", and scaled - in "scale.data"
        verbose=FALSE
    )
    SeuratObject::Idents(seurat_data) <- "new.ident"
    return (aggregated_seurat_data)
}

get_de_cell_data <- function(seurat_data, args, samples_order=NULL){
    if (!is.null(samples_order)){
        cell_annotations <- c("new.ident", args$splitby)
    } else {
        cell_annotations <- c(args$splitby, "new.ident")
    }
    if(!is.null(args$batchby)){
        cell_annotations <- c(cell_annotations, args$batchby)
    }
    if(
        all(base::as.vector(as.character(seurat_data@meta.data$new.ident)) != base::as.vector(as.character(seurat_data@meta.data$condition))) &&
        length(base::unique(base::as.vector(as.character(seurat_data@meta.data$condition)))) > 1
    ){
        cell_annotations <- c(cell_annotations, "condition")                                      # several conditions found
    }
    custom_fields <- base::grep("^custom_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
    if(length(custom_fields) > 0){
        cell_annotations <- c(cell_annotations, custom_fields)                                    # adding all custom fields
    }
    cell_annotations <- base::unique(cell_annotations)                                            # in case any of the found columns were defined in splitby or batchby
    cell_data <- seurat_data@meta.data %>%
                 dplyr::select(tidyselect::all_of(cell_annotations))
    if (!is.null(samples_order)){
        cell_data <- cell_data %>%
                     dplyr::mutate(new.ident=base::factor(new.ident, levels=samples_order))       # setting levels for new.ident from samples_order
    }
    cell_data[[args$splitby]] <- base::factor(                                                    # to have --second on the left side from --first (on the plots)
                                     cell_data[[args$splitby]],
                                     levels=c(args$second, args$first)
                                 )
    cell_data <- cell_data %>%
                 dplyr::arrange_at(cell_annotations) %>%                                          # instead of dplyr::arrange(dplyr::across(tidyselect::all_of(cell_annotations)))
                 dplyr::mutate_at(base::colnames(.), base::as.vector)                             # no reasons to have it as factor
    base::print(utils::head(cell_data))
    return (cell_data)
}

get_de_sample_data <- function(seurat_data, samples_order, args){                                 # splitby and batchby should be concordant with new.ident
    sample_data <- seurat_data@meta.data %>%
                   dplyr::select(
                       new.ident,
                       tidyselect::all_of(args$splitby),
                       tidyselect::any_of(args$batchby)                                           # use any_of(args$batchby) because it can be NULL
                   ) %>%
                   dplyr::distinct() %>%
                   dplyr::mutate(new.ident=base::factor(new.ident, levels=samples_order)) %>%     # setting levels for new.ident from samples_order
                   dplyr::arrange(new.ident) %>%                                                  # sorting by levels defined from samples_order
                   tibble::remove_rownames() %>%
                   tibble::column_to_rownames("new.ident") %>%
                   dplyr::mutate_at(base::colnames(.), base::factor)                              # DEseq prefers factors
    sample_data[[args$splitby]] <- stats::relevel(sample_data[[args$splitby]], args$first)        # relevel to have args$first as a base for DESeq comparison
    base::print(sample_data)
    return (sample_data)
}

get_bulk_counts_data <- function(deseq_data, sample_data){
    base::print("Applying rlog transformation (not blind to the experimental design)")
    bulk_counts_data <- DESeq2::rlog(deseq_data, blind=FALSE)                              # no reason to use VST as we have less than 30 datasets

    # if (args$norm == "vst"){
    #     base::print("Applying vst transformation (not blind to the experimental design)")
    #     bulk_counts_data <- DESeq2::vst(deseq_data, blind=FALSE)
    # }
    # if(!is.null(args$batchby) && !is.null(args$remove) && args$remove){
    #     base::print("Removing batch effect from the normalized counts")
    #     SummarizedExperiment::assay(bulk_counts_data) <- limma::removeBatchEffect(
    #         SummarizedExperiment::assay(bulk_counts_data),
    #         batch=bulk_counts_data[[args$batchby]],
    #         design=stats::model.matrix(stats::as.formula(base::paste("~",args$splitby)), sample_data)  # should include only splitby
    #     )
    # }

    base::print("Normalized read counts")
    base::print(utils::head(SummarizedExperiment::assay(bulk_counts_data)))
    base::print(dim(SummarizedExperiment::assay(bulk_counts_data)))
    base::print(SummarizedExperiment::colData(bulk_counts_data))
    return (bulk_counts_data)
}

get_clustered_data <- function(expression_data, center_row, dist, transpose) {
    if (transpose){
        base::print("Transposing expression data")
        expression_data = t(expression_data)
    }
    if (center_row) {
        base::print("Centering expression data by row means")
        expression_data = expression_data - base::rowMeans(expression_data)    
    }
    base::print("Creating distance matrix")
    distance_matrix <- hopach::distancematrix(expression_data, dist)
    base::print("Running HOPACH")
    hopach_results <- hopach::hopach(expression_data, dmat=distance_matrix)

    if (transpose){
        base::print("Transposing expression data back")
        expression_data = base::t(expression_data)
    }

    base::print("Parsing cluster labels")
    clusters = base::as.data.frame(hopach_results$clustering$labels)
    base::colnames(clusters) = "label"
    clusters = base::cbind(
        clusters,
        "HCL"=outer(
            clusters$label,
            10^c((base::nchar(trunc(clusters$label))[1]-1):0),
            function(a, b) {
                base::paste0("c", a %/% b %% 10)
            }
        )
    )
    clusters = clusters[, c(-1), drop=FALSE]
    return (
        list(
            order=base::as.vector(hopach_results$clustering$order),
            expression=expression_data,
            clusters=clusters
        )
    )
}

atac_dbinding_analyze <- function(seurat_data, args){
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"                                   # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                    # safety measure

    base::print(
        base::paste0(
            "Running ", args$second, " vs ", args$first,
            " differential binding analysis using ", args$test,
            " test for cells split by ", args$splitby,
            base::ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste(
                    " subsetted to", base::paste(args$subset, collapse=", "),
                    "values from", args$groupby, "column."
                ),
                "."
            )
        )
    )

    results <- list(db_sites=NULL)                                          # to collect all outputs

    if(args$test == "manorm2"){
        base::print("Counting reads in reference genomic bins")
        profile_bins_args <- c(
            base::paste0(
                "--peaks=",
                args$tmp_locations$macs2$peaks_narrow$second, ",",
                args$tmp_locations$macs2$peaks_narrow$first
            ),
            base::paste0(
                "--reads=",
                args$tmp_locations$macs2$cut_sites$second, ",",
                args$tmp_locations$macs2$cut_sites$first
            ),
            base::paste0(
                "--summits=",
                args$tmp_locations$macs2$peaks_summit$second, ",",
                args$tmp_locations$macs2$peaks_summit$first
            ),
            base::paste0("--typical-bin-size=", args$binsize),
            base::paste0("--min-peak-gap=", args$minpeakgap),
            "--labs=second,first",
            "--shiftsize=0",
            "--keep-dup=all",
            "-n", "peak"
        )
        if (!is.null(args$blacklist)){
            profile_bins_args <- c(profile_bins_args, base::paste0("--filter=", args$blacklist))
        }
        if (!is.null(args$maxpeaks)){
            profile_bins_args <- c(profile_bins_args, base::paste0("--keep-peaks=", args$maxpeaks))
        }
        exit_code <- sys::exec_wait(
            cmd="profile_bins",                                             # if it's not found in PATH, R will fail with error
            args=profile_bins_args
        )
        if (exit_code != 0){                                                # we were able to run profile_bins, but something went wrong
            base::print(
                base::paste0(
                    "Failed to count reads in reference genomic ",
                    "bins with exit code ", exit_code, ". Exiting."
                )
            )
            base::quit(save="no", status=1, runLast=FALSE)                  # force R to exit with error
        }

        base::print("Loadind read counts data")
        read_counts_data <- utils::read.table(
            "peak_profile_bins.xls",
            sep="\t",
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        base::print(utils::head(read_counts_data))

        base::print(
            base::paste(
                "Normalizing read counts data excluding chrX and chrY",
                "from being searched for common reference genomic bins."
            )
        )
        read_counts_data <- MAnorm2::normalize(                                              # baseline should be selected automaticaly so the order of columns doesn't matter
            read_counts_data,
            count=c("second.read_cnt", "first.read_cnt"),                                    # columns with read counts
            occupancy=c("second.occupancy", "first.occupancy"),                              # shows which reference genomic bins are overlapped by peaks
            common.peak.regions=!(read_counts_data$chrom %in% c("chrX", "chrY", "X", "Y"))   # to exclude chrX and Y from being searched for common regions
        )
        base::print(
            base::paste(
                "Automatically selected as a baseline:",
                attr(read_counts_data, "baseline")
            )
        )
        base::print(utils::head(read_counts_data))

        base::print("Adding biological conditions")
        bio_conditions <- list(
            second = MAnorm2::bioCond(
                norm.signal=read_counts_data$second.read_cnt,
                occupancy=read_counts_data$second.occupancy,
                meta.info=read_counts_data[, c("chrom", "start", "end")],                    # to have reference genomic bins coordinates embedded (not used by MAnorm2)
                name="second"
            ),
            first = MAnorm2::bioCond(
                norm.signal=read_counts_data$first.read_cnt,
                occupancy=read_counts_data$first.occupancy,
                name="first"
            ),
            common = MAnorm2::bioCond(
                norm.signal=read_counts_data[, c("second.read_cnt", "first.read_cnt")],
                occupancy=read_counts_data[, c("second.occupancy", "first.occupancy")],
                occupy.num=2,                                                                # include only reference genomic bins occupied by all conditions
                name="common"
            )
        )
        base::rm(read_counts_data)
        base::print("Fitting mean-variance curve based on the common bins")
        bio_conditions <- MAnorm2::fitMeanVarCurve(                                          # will add fit.info field to the bioCond objects
            bio_conditions,                                                                  # only "common" contains replicates, so it will be used for fitting MVC
            method="parametric fit",
            occupy.only=TRUE,
            init.coef=c(0.1, 10)                                                             # per manual it is expected to suit most practical datasets
        )

        base::print("Running differential test")
        db_sites <- MAnorm2::diffTest(                                                       # data frame without coordinates, Mval is calculated as Y/X
                        x=bio_conditions$first,
                        y=bio_conditions$second,
                    ) %>%
                    tibble::rownames_to_column(var="rowname") %>%                            # we need it to join with coordinates
                    dplyr::left_join(
                        bio_conditions$second$meta.info %>%                                  # we take the coordinates from the meta.info of the second bio condition
                        tibble::rownames_to_column(var="rowname"),
                        by="rowname"
                    ) %>%
                    stats::na.omit() %>%                                                     # we shouldn't have any NAs, but filter just in case 
                    dplyr::rename(
                        "pvalue"="pval",
                        "log2FoldChange"="Mval",
                        "lfcSE"="Mval.se",                                                   # to have the name similar to what DESeq2 reports
                        "padj"="padj",
                        "chr"="chrom"                                                        # for consistency
                    ) %>%
                    dplyr::mutate(
                        "baseMean"=(first.mean + second.mean) / 2                            # to have the same column name as DESeq2 reports
                    ) %>%
                    dplyr::select(-c("rowname", "Mval.t", "first.mean", "second.mean")) %>%  # removing not relevant columns
                    dplyr::select(                                                           # to have a proper coliumns order
                        c("chr", "start", "end", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
                    )
        base::print(
            base::paste(
                "Number of differentially bound sites:",
                base::nrow(db_sites)
            )
        )
        base::rm(bio_conditions)
        results$db_sites <- db_sites                               # not filtered differentialy bound sites
    } else {
        SeuratObject::Idents(seurat_data) <- args$splitby
        db_sites <- Seurat::FindMarkers(
                        seurat_data,
                        assay="ATAC",
                        slot="data",                               # using normalized counts
                        ident.1=args$second,
                        ident.2=args$first,                        # this is the reference as logFC = ident.1 / ident.2
                        logfc.threshold=0.25,                      # only to make it less computationaly heavy
                        min.pct=0.05,                              # for ATAC we need smaller value
                        min.diff.pct=-Inf,                         # using default value
                        only.pos=FALSE,                            # we want to have positive and negative results
                        test.use=args$test,                        # this will always be something that supports latent.vars
                        latent.vars="nCount_ATAC",                 # recommended for ATAC
                        base=2,                                    # to make sure we use log2 scale
                        verbose=FALSE
                    ) %>%
                    stats::na.omit() %>%                           # we shouldn't have any NAs, but filter just in case 
                    tibble::rownames_to_column(var="peak") %>%     # peak, p_val, avg_log2FC, pct.1, pct.2, p_val_adj
                    tidyr::separate(
                        col="peak",
                        into=c("chr", "start", "end"),
                        sep="-",
                        remove=TRUE
                    ) %>%
                    dplyr::rename(                                 # peak, pvalue, log2FoldChange, pct.1, pct.2, padj
                        "pvalue"="p_val",
                        "log2FoldChange"="avg_log2FC",
                        "padj"="p_val_adj"
                    ) %>%
                    dplyr::select(                                 # to have a proper coliumns order
                        c("chr", "start", "end", "log2FoldChange", "pct.1", "pct.2", "pvalue", "padj")
                    )
        SeuratObject::Idents(seurat_data) <- "new.ident"
        base::print(
            base::paste(
                "Number of differentially bound sites:",
                base::nrow(db_sites)
            )
        )
        results$db_sites <- db_sites                                                    # not filtered differentialy bound sites
    }

    base::print("Adding closest genes for differentially bound sites")
    results$db_sites <- results$db_sites %>%
                        tidyr::unite(
                            query_region,                                               # we will left join by this column
                            chr:start:end,
                            remove=FALSE,
                            sep="-"
                        ) %>%
                        dplyr::left_join(
                            Signac::ClosestFeature(                                     # will add query_region and gene_name columns
                                seurat_data,
                                regions=GenomicRanges::makeGRangesFromDataFrame(
                                    results$db_sites,
                                    seqinfo=seurat_data[["ATAC"]]@seqinfo               # need to have seqinfo in Seurat object
                                )
                            ) %>% dplyr::select(c("query_region", "gene_name")),
                            by="query_region"
                        ) %>%
                        dplyr::select(-c("query_region")) %>%
                        dplyr::relocate(gene_name, .after=last_col())                   # move gene_name to the end

    return (results)
}

rna_de_analyze <- function(seurat_data, args, excluded_genes=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                    # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                    # safety measure

    selected_genes <- base::as.vector(as.character(base::rownames(seurat_data)))        # all available genes
    if (!is.null(excluded_genes) && length(excluded_genes) > 0){
        base::print(base::paste("Excluding", length(excluded_genes), "genes"))
        selected_genes <- selected_genes[!(selected_genes %in% excluded_genes)]
    }

    base::print(
        base::paste0(
            "Running ", args$second, " vs ", args$first,
            " differential expression analysis using ",
            args$test, " test for cells",
            base::ifelse(
                (args$test %in% c("deseq", "lrt")),
                " aggregated to pseudobulk form by dataset ",
                " "
            ),
            "split by ", args$splitby,
            base::ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste(
                    " subsetted to", base::paste(args$subset, collapse=", "),
                    "values from", args$groupby, "column. "
                ),
                ". "
            ),
            base::ifelse(
                !is.null(args$batchby),
                base::paste0("Model batch effect by ", args$batchby, " column."),
                ""
            )
        )
    )

    results <- list(de_genes=NULL, bulk=NULL, cell=NULL)                                # to collect all outputs

    if(args$test %in% c("deseq", "lrt")){                                               # aggregating to the pseudobulk form
        aggregated_seurat_data <- get_aggregated_expession(
            seurat_data,
            group_by="new.ident",                                                       # aggregating by sample, because we want to run pseudobulk
            selected_genes=selected_genes,
            slot="counts"                                                               # using not normalized counts, because DESeq needs raw counts
        )
        raw_counts <- SeuratObject::GetAssayData(                                       # dgCMatrix
            aggregated_seurat_data,
            assay="RNA",                                                                # set assay in case GetAssayData won't take the default value
            slot="counts"                                                               # we need not normalized counts, because DESeq needs raw counts
        )
        sample_data <- get_de_sample_data(                                              # rows will be sorted by columns from raw_counts
            seurat_data=seurat_data,
            samples_order=base::unname(raw_counts@Dimnames[[2]]),                       # column names from dgCMatrix
            args=args
        )
        design_formula <- stats::as.formula(
            base::paste0(
                "~",
                base::ifelse(
                    is.null(args$batchby),
                    "",
                    base::paste0(args$batchby, "+")
                ),
                args$splitby                                                            # safer to have the condition of interest on the last position
            )
        )
        deseq_data <- DESeq2::DESeqDataSetFromMatrix(
            countData=raw_counts,
            colData=sample_data,
            design=design_formula
        )

        base::rm(aggregated_seurat_data, raw_counts)                                   # to free up some memory

        if (args$test == "lrt"){
            reduced_formula <- stats::as.formula("~1")
            if(!is.null(args$batchby)){
                reduced_formula <- stats::as.formula(base::paste0("~", args$batchby))
            }
            base::print(                                                                # see this post for details https://support.bioconductor.org/p/95493/#95572
                base::paste(
                    "Using LRT test with the design formula",
                    base::paste(design_formula, collapse=""),
                    "and the reduced formula",
                    base::paste(reduced_formula, collapse=""),
                    "to calculate p-values."
                )
            )
            deseq_data <- DESeq2::DESeq(
                deseq_data,
                test="LRT",
                reduced=reduced_formula,
                quiet=TRUE,
                parallel=TRUE,
                BPPARAM=BiocParallel::MulticoreParam(args$cpus)                         # add it here as well just in case
            )
        } else {
            base::print(
                base::paste(
                    "Using Wald test with the design formula",
                    base::paste(design_formula, collapse=""),
                    "to calculate p-values."
                )
            )
            deseq_data <- DESeq2::DESeq(
                deseq_data,
                quiet=TRUE,
                parallel=TRUE,
                BPPARAM=BiocParallel::MulticoreParam(args$cpus)                         # add it here as well just in case
            )
        }

        print("Estimated effects")
        base::print(DESeq2::resultsNames(deseq_data))

        de_genes <- DESeq2::results(
            deseq_data,
            contrast=c(args$splitby, args$second, args$first),            # we are interested in seconds vs first fold change values
            alpha=base::ifelse(is.null(args$padj), 0.1, args$padj),       # recommended to set to our FDR threshold https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
            parallel=TRUE,
            BPPARAM=BiocParallel::MulticoreParam(args$cpus)               # add it here as well just in case
        )

        base::print("Results description")
        base::print(S4Vectors::mcols(de_genes))
        base::print(utils::head(de_genes))

        de_genes <- base::as.data.frame(de_genes) %>%
                    stats::na.omit() %>%                                  # exclude all rows where NA is found in any column. See http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
                    tibble::rownames_to_column(var="gene")                # gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj

        base::print(
            base::paste(
                "Number of differentially expressed genes after",
                "excluding 'NA':", base::nrow(de_genes)
            )
        )

        base::print("Normalizing read counts data aggregated to pseudobulk form")
        counts_data <- get_bulk_counts_data(deseq_data, sample_data)

        row_metadata <- de_genes %>%
                        tibble::remove_rownames() %>%
                        tibble::column_to_rownames("gene") %>%
                        dplyr::select(log2FoldChange, pvalue, padj)  %>%
                        dplyr::filter(.$padj<=args$padj) %>%
                        dplyr::arrange(desc(log2FoldChange))

        column_metadata <- sample_data %>%          # we keep the original order without changing it by --splitby
                           dplyr::mutate_at(        # because we can mess up when merging it with clusters
                               base::colnames(.),
                               base::as.vector      # need to convert to vector, because in our sample_data everything was a factor
                           )

        if (!is.null(args$cluster)){
            base::print(
                base::paste(
                    "Filtering normalized read counts data to include",
                    "only differentially expressed features with padj",
                    "<=", args$padj, "(for clustering)"
                )
            )
            # need it only for clustering to define rows and columns order
            filtered_counts_mat <- SummarizedExperiment::assay(counts_data)[base::as.vector(base::rownames(row_metadata)),]
            base::print("Size after filtering")
            base::print(dim(filtered_counts_mat))

            if (args$cluster == "column" || args$cluster == "both") {
                base::print("Clustering filtered normalized read counts data by columns")
                clustered_data = get_clustered_data(
                    expression_data=filtered_counts_mat,
                    center_row=FALSE,                                                       # centering doesn't influence on the samples order
                    dist=args$columndist,
                    transpose=TRUE
                )
                column_metadata <- base::cbind(column_metadata, clustered_data$clusters)    # adding cluster labels
                column_metadata <- column_metadata[clustered_data$order, ]                  # reordering samples order based on the HOPACH clustering resutls
                base::print("Reordered samples")
                base::print(column_metadata)
            }
            if (args$cluster == "row" || args$cluster == "both") {
                base::print("Clustering filtered normalized read counts data by rows")
                clustered_data = get_clustered_data(
                    expression_data=filtered_counts_mat,
                    center_row=base::ifelse(is.null(args$center), FALSE, args$center),   # about centering normalized data https://www.biostars.org/p/387863/
                    dist=args$rowdist,
                    transpose=FALSE
                )
                row_metadata <- base::cbind(row_metadata, clustered_data$clusters)       # adding cluster labels
                row_metadata <- row_metadata[clustered_data$order, ]                     # reordering features order based on the HOPACH clustering results
                base::print("Reordered features")
                base::print(utils::head(row_metadata))
            }
            base::rm(clustered_data, filtered_counts_mat)                                # these two object will be always defined if we got here
        }

        results$de_genes <- de_genes                                                     # not filtered differentialy expressed genes
        results$bulk <- list(
            column_metadata=column_metadata,                                             # column metadata per dataset always either default or clustered by new.ident
            counts_data=counts_data                                                      # not filtered normalized counts in a form of SummarizedExperiment
        )
        base::rm(column_metadata)                                                        # remove it just in case
    } else {
        SeuratObject::Idents(seurat_data) <- args$splitby
        de_genes <- Seurat::FindMarkers(
                        seurat_data,
                        assay="RNA",
                        slot="data",                               # using normalized counts
                        ident.1=args$second,
                        ident.2=args$first,                        # this is the reference as logFC = ident.1 / ident.2
                        features=selected_genes,
                        logfc.threshold=0.25,                      # using default value
                        min.pct=0.1,                               # using default value
                        min.diff.pct=-Inf,                         # using default value
                        only.pos=FALSE,                            # we want to have both up and dowregulated genes
                        test.use=args$test,                        # at this moment it should one of the supported types
                        latent.vars=args$batchby,                  # will work only for negbinom, poisson, LR, and MAST, otherwise should be ignored
                        base=2,                                    # to make sure we use log2 scale
                        verbose=FALSE
                    ) %>%
                    stats::na.omit() %>%                           # we shouldn't have any NAs, but filter just in case 
                    tibble::rownames_to_column(var="gene") %>%     # gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj
                    dplyr::rename(                                 # gene, pvalue, log2FoldChange, pct.1, pct.2, padj
                        "pvalue"="p_val",
                        "log2FoldChange"="avg_log2FC",
                        "padj"="p_val_adj"
                    )
        SeuratObject::Idents(seurat_data) <- "new.ident"
        base::print(base::paste("Number of differentially expressed genes:", base::nrow(de_genes)))

        row_metadata <- de_genes %>%
                        tibble::remove_rownames() %>%
                        tibble::column_to_rownames("gene") %>%
                        dplyr::select(log2FoldChange, pct.1, pct.2, pvalue, padj)  %>%
                        dplyr::filter(.$padj<=args$padj) %>%
                        dplyr::arrange(desc(log2FoldChange))

        results$de_genes <- de_genes                                             # not filtered differentialy expressed genes
    }

    counts_mat <- base::as.matrix(                                               # will be either sorted by log2FoldChange or clustered based on pseudobulk expression
        SeuratObject::GetAssayData(                                              # dgCMatrix
            seurat_data,
            assay="RNA",                                                         # set assay in case GetAssayData won't take the default value
            slot="data"                                                          # these are normalized counts
        )                                                                        # to subset only to sorted by log2FoldChange and filtered by args$padj genes
    )[base::as.vector(base::rownames(row_metadata)), ]                           # row_metadata can come either from pseudobulk or by cell parts

    base::print("Size of the normalized cell read counts matrix after filtering")
    base::print(dim(counts_mat))
    base::print(counts_mat[1:5, 1:5])

    column_metadata <- get_de_cell_data(
                           seurat_data=seurat_data,
                           samples_order=if (!is.null(results$bulk) && !is.null(args$cluster) && (args$cluster %in% c("both", "column")))   # we keep the clustered datasets order if available
                                             base::rownames(results$bulk$column_metadata)
                                         else
                                            NULL,
                           args=args
                       )

    if (is.null(results$bulk) && !is.null(args$cluster) && args$cluster == "row"){      # we first check if it was first processed by pseudobulk
        base::print("Clustering filtered normalized cell read counts by rows")          # if yes, then to row order should be be changed (it's either
        clustered_data = get_clustered_data(                                            # already clustered or sorted by log2FoldChange
            expression_data=counts_mat,
            center_row=base::ifelse(is.null(args$center), FALSE, args$center),          # about centering normalized data https://www.biostars.org/p/387863/
            dist=args$rowdist,
            transpose=FALSE
        )
        counts_mat <- clustered_data$expression                                         # can be different because of optional centering by rows mean
        row_metadata <- base::cbind(row_metadata, clustered_data$clusters)              # adding cluster labels
        row_metadata <- row_metadata[clustered_data$order, ]                            # reordering features order based on the HOPACH clustering results
        base::print("Reordered features")
        base::print(utils::head(row_metadata))
        base::rm(clustered_data)
    }

    results$cell <- list(
        row_metadata=row_metadata,                          # either sorted by log2FoldChange or clustered (both from pseudobulk or from cell levels)
        column_metadata=column_metadata,                    # column metadata per cell with optionally clusterd new.ident (comes from pseudobulk)
        counts_mat=counts_mat                               # not ordered cell counts matrix, filtered by padj
    )
    if (!is.null(results$bulk)){
        results$bulk$row_metadata <- row_metadata           # adding for easy access, but it should be always the same as results$cell$row_metadata
    }

    return (results)

}

da_analyze <- function(seurat_data, args){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                # safety measure
    base::print("Defining experimental design")
    idents <- base::as.vector(as.character(SeuratObject::Idents(seurat_data)))
    sample_data <- get_de_sample_data(
        seurat_data=seurat_data,
        samples_order=base::unique(idents),                                               # we don't really care about idents order here
        args=args
    )
    first_group <- base::as.vector(
        as.character(
            rownames(sample_data[sample_data[[args$splitby]] == args$first, , drop=FALSE])
        )
    )
    second_group <- base::as.vector(
        as.character(
            rownames(sample_data[sample_data[[args$splitby]] == args$second, , drop=FALSE])
        )
    )
    base::print(
        base::paste(
            "First group of cells identities:", base::paste(first_group, collapse=", ")
        )
    )
    base::print(
        base::paste(
            "Second group of cells identities:", base::paste(second_group, collapse=", ")
        )
    )
    embeddings <- SeuratObject::Embeddings(
        seurat_data,
        reduction=args$reduction
    )
    embeddings <- embeddings[, args$dimensions]   # subset to specific dimensions
    base::print("Selected embeddings")
    base::print(utils::head(embeddings))

    da_cells <- DAseq::getDAcells(
        X=embeddings,
        cell.labels=idents,
        labels.1=first_group,
        labels.2=second_group,
        pred.thres=args$ranges,                   # if NULL, will be calculated automatically
        k.vector=args$knn,                        # if NULL, will be calculated based on the cells number
        n.runs=5,                                 # the same as default value
        n.rand=5,                                 # instead of default 2
        do.plot=TRUE                              # will save only rand.plot
    )

    # DA-seq computes for each cell a score based on the relative prevalence of cells from both biological
    # states in the cell’s neighborhood. DA score measures of how much a cell’s neighborhood is dominated
    # by cells from one of the biological states
    da_sufix <- base::paste0(args$second, "_vs_", args$first)
    seurat_data[[base::paste("custom", "da_score", da_sufix, sep="_")]] <- da_cells$da.pred

    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        base::print(base::paste("Identifying DA subpopulations using resolution", current_resolution))
        # DA-seq clusters the cells whose DA measure is above or below a certain threshold
        da_regions <- DAseq::getDAregion(
            X=embeddings,
            da.cells=da_cells,
            cell.labels=idents,
            labels.1=first_group,
            labels.2=second_group,
            resolution=current_resolution
        )
        seurat_data[[base::paste0("da_", da_sufix, "_res.", current_resolution)]] <- da_regions$da.region.label
        base::print(da_regions$DA.stat)
        base::rm(da_regions)
    }

    base::rm(idents, sample_data, first_group, second_group, embeddings)         # remove unused data
    base::gc(verbose=FALSE)
    return (
        list(
            seurat_data=seurat_data,
            da_cells=da_cells,
            thresholds=c(max(unlist(da_cells$rand.pred)), min(unlist(da_cells$rand.pred)))
        )
    )
}