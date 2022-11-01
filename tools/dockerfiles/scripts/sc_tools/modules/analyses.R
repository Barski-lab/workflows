import("dplyr", attach=FALSE)
import("limma", attach=FALSE)
import("DAseq", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("sceasy", attach=FALSE)
import("hopach", attach=FALSE)
import("DESeq2", attach=FALSE)
import("harmony", attach=FALSE)
import("tibble", attach=FALSE)
import("glmGamPoi", attach=FALSE)  # safety measure. we don't use it directly, but SCTransform with method="glmGamPoi" needs it
import("S4Vectors", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)
import("reticulate", attach=FALSE)
import("BiocParallel", attach=FALSE)
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
    "get_putative_markers",
    "atac_preprocess",
    "atac_analyze",
    "add_wnn_clusters",
    "rna_de_analyze",
    "da_analyze",
    "get_de_sample_data",
    "get_norm_counts_data",
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
    arguments <- c(args$regressmt, args$regresscellcycle)   # any of then can be also NULL
    metadata_columns <- c("mito_percentage", "S.Score&G2M.Score")
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

get_cell_cycle_scores <- function(seurat_data, assay, cell_cycle_data){   # we need this function to fail if something went wrong, so no tryCatch inside
    SeuratObject::DefaultAssay(seurat_data) <- assay                      # safety measure
    seurat_data <- Seurat::CellCycleScoring(
        seurat_data,
        s.features=base::as.vector(cell_cycle_data[base::tolower(cell_cycle_data$phase)=="s", "gene_id"]),
        g2m.features=base::as.vector(cell_cycle_data[base::tolower(cell_cycle_data$phase)=="g2/m", "gene_id"]),
        assay=assay,
        verbose=FALSE
    )
    # seurat_data[["CC.Diff"]] <- seurat_data[["S.Score"]] - seurat_data[["G2M.Score"]]   # for softer cell cycle removal (https://satijalab.org/seurat/articles/cell_cycle_vignette.html)
    return (seurat_data)
}

rna_log_single <- function(seurat_data, args, cell_cycle_data=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                                # safety measure
    base::print("Applying LogNormalize")
    scaled_norm_seurat_data <- Seurat::NormalizeData(seurat_data, verbose=FALSE)
    if(!is.null(cell_cycle_data)){
        base::tryCatch(
            expr = {
                base::print("Trying to assign cell cycle scores for RNA assay")
                scaled_norm_seurat_data <- get_cell_cycle_scores(                                       # if succeded adds S.Score and G2M.Score columns to medadata
                    scaled_norm_seurat_data,
                    "RNA",
                    cell_cycle_data
                )
            },
            error = function(e){
                base::print(base::paste("Failed to run cell cycle scoring for RNA assay with error - ", e))
            }
        )
    }
    scaled_norm_seurat_data <- Seurat::FindVariableFeatures(
        scaled_norm_seurat_data,
        nfeatures=args$highvargenes,
        verbose=FALSE
    )
    vars_to_regress <- get_vars_to_regress(scaled_norm_seurat_data, args)                           # may or may not include S.Score and G2M.Score columns
    base::print(base::paste0("Regressing out [", paste(vars_to_regress, collapse=", "), "]"))
    scaled_norm_seurat_data <- Seurat::ScaleData(
        scaled_norm_seurat_data,
        vars.to.regress=vars_to_regress,
        verbose=FALSE
    )
    SeuratObject::DefaultAssay(scaled_norm_seurat_data) <- "RNA"
    return (scaled_norm_seurat_data)
}

rna_sct_single <- function(seurat_data, args, cell_cycle_data=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                # safety measure
    method <- base::ifelse(args$norm=="sctglm", "glmGamPoi", "poisson")
    vars_to_regress <- get_vars_to_regress(seurat_data, args)                       # will never include "S.Score", "G2M.Score" columns as cell cycle scores are not assigned
    base::print(
        base::paste0(
            "Applying SCTransform using ", method,
            " method for initial parameter estimation. ",
            "Regressing out [", paste(vars_to_regress, collapse=", "), "]"
        )
    )
    scaled_norm_seurat_data <- Seurat::SCTransform(
        seurat_data,                                                                # use not splitted seurat_data
        assay="RNA",
        new.assay.name="SCT",
        variable.features.n=args$highvargenes,
        method=method,
        vars.to.regress=vars_to_regress,                                            # first portion of variables to regress. will never include "S.Score" and "G2M.Score"
        conserve.memory=args$lowmem,
        verbose=FALSE                                                               # too many stdout
    )
    if(!is.null(cell_cycle_data)){
        base::tryCatch(
            expr = {
                base::print("Trying to assign cell cycle scores for SCT assay")
                scaled_norm_seurat_data <- get_cell_cycle_scores(                       # if succeded adds S.Score and G2M.Score columns to medadata
                    scaled_norm_seurat_data,
                    "SCT",
                    cell_cycle_data
                )
                vars_to_regress <- get_vars_to_regress(scaled_norm_seurat_data, args)   # may include S.Score and G2M.Score if run with --regresscellcycle
                if (all(c("S.Score", "G2M.Score") %in% vars_to_regress)){               # need to rerun SCTransform to regress all variables at once
                    base::print(
                        base::paste0(
                            "Re-applying SCTransform using ", method,
                            " method for initial parameter estimation. ",
                            "Regressing out [", base::paste(vars_to_regress, collapse=", "), "]"
                        )
                    )
                    scaled_norm_seurat_data <- Seurat::SCTransform(
                        scaled_norm_seurat_data,
                        assay="RNA",
                        new.assay.name="SCT",
                        variable.features.n=args$highvargenes,
                        method=method,
                        vars.to.regress=vars_to_regress,
                        conserve.memory=args$lowmem,
                        verbose=FALSE
                    )
                }
            },
            error = function(e){
                base::print(base::paste("Failed to run cell cycle scoring/regressing for SCT assay with error - ", e))
            }
        )
    }
    SeuratObject::DefaultAssay(scaled_norm_seurat_data) <- "SCT"
    return (scaled_norm_seurat_data)
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
        SeuratObject::DefaultAssay(splitted_seurat_data[[i]]) <- "RNA"                            # safety measure
        base::print(
            base::paste(
                "Applying LogNormalize for", SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                "dataset"
            )
        )
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
                        splitted_seurat_data[[i]] <- get_cell_cycle_scores(                       # if succeded adds S.Score and G2M.Score columns to medadata
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
        splitted_seurat_data[[i]] <- Seurat::FindVariableFeatures(
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
        k.weight=min(get_min_ident_size(splitted_seurat_data), 100),        # k.weight 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
        verbose=FALSE
    )
    if (failed_cell_cycle_scoring){
        base::print(
            base::paste(
                "At least one of the datasets failed in cell cycle score assignment.",
                "Removing S.Score and G2M.Score columns from metadata."
            )
        )
        integrated_seurat_data[["S.Score"]] <- NULL
        integrated_seurat_data[["G2M.Score"]] <- NULL
    }
    vars_to_regress <- get_vars_to_regress(integrated_seurat_data, args)                # may or may not include S.Score and G2M.Score columns
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
        SeuratObject::DefaultAssay(splitted_seurat_data[[i]]) <- "RNA"                            # safety measure
        vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args)                   # will never include S.Score and G2M.Score
        base::print(
            base::paste0(
                "Applying SCTransform for ", SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                " dataset using ", method, " method for initial parameter estimation. ",
                "Regressing out [", base::paste(vars_to_regress, collapse=", "), "]"
            )
        )
        splitted_seurat_data[[i]] <- Seurat::SCTransform(
            splitted_seurat_data[[i]],
            assay="RNA",
            new.assay.name="SCT",
            variable.features.n=args$highvargenes,
            method=method,
            vars.to.regress=vars_to_regress,
            conserve.memory=args$lowmem,
            verbose=FALSE
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
                        splitted_seurat_data[[i]] <- get_cell_cycle_scores(                              # if succeded adds S.Score and G2M.Score columns to medadata
                            splitted_seurat_data[[i]],
                            "SCT",
                            cell_cycle_data
                        )
                        vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args)          # may include S.Score and G2M.Score if run with --regresscellcycle
                        if (all(c("S.Score", "G2M.Score") %in% vars_to_regress)){                        # need to rerun SCTransform to regress all variables at once
                            base::print(
                                base::paste0(
                                    "Re-applying SCTransform for ", SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                                    " dataset using ", method, " method for initial parameter estimation. ",
                                    "Regressing out [", base::paste(vars_to_regress, collapse=", "), "]"
                                )
                            )
                            splitted_seurat_data[[i]] <- Seurat::SCTransform(
                                splitted_seurat_data[[i]],
                                assay="RNA",
                                new.assay.name="SCT",
                                variable.features.n=args$highvargenes,
                                method=method,
                                vars.to.regress=vars_to_regress,
                                conserve.memory=args$lowmem,
                                verbose=FALSE
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
        base::print(base::paste("At least one of the datasets failed in cell cycle score assignment."))
        for (i in 1:length(splitted_seurat_data)){
            if(!is.null(args$regresscellcycle) && args$regresscellcycle){
                vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args, "S.Score&G2M.Score")    # force to exclude "S.Score&G2M.Score"
                base::print(
                    base::paste0(
                        "Re-applying SCTransform for ", SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                        " dataset using ", method, " method for initial parameter estimation. ",
                        "Regressing out [", base::paste(vars_to_regress, collapse=", "), "]. ",
                        "Skipping cell cycle score assignment."
                    )
                )
                splitted_seurat_data[[i]] <- Seurat::SCTransform(
                    splitted_seurat_data[[i]],
                    assay="RNA",
                    new.assay.name="SCT",
                    variable.features.n=args$highvargenes,
                    method=method,
                    vars.to.regress=vars_to_regress,
                    conserve.memory=args$lowmem,
                    verbose=FALSE
                )
            } else {
                base::print(
                    base::paste(
                        "Removing S.Score and G2M.Score columns from the metadata of",
                        SeuratObject::Idents(splitted_seurat_data[[i]])[1]
                    )
                )
                splitted_seurat_data[[i]][["S.Score"]] <- NULL
                splitted_seurat_data[[i]][["G2M.Score"]] <- NULL
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
    for (reduction_name in c("atac_lsi", "atacumap", "wnnumap")){
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
    base::print(
        base::paste(
            "Performing PCA reduction on", SeuratObject::DefaultAssay(seurat_data),
            "assay using 50 principal components"
        )
    )
    seurat_data <- Seurat::RunPCA(seurat_data, npcs=50, verbose=FALSE)          # add "pca" reduction that should be used in UMAP
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
                dims.use=args$dimensions,
                assay.use=SeuratObject::DefaultAssay(seurat_data),     # can be both RNA or SCT depending on --norm parameter
                verbose=FALSE
            )
        }
    }
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

get_putative_markers <- function(seurat_data, assay, resolution_prefix, args, latent_vars=NULL, min_diff_pct=-Inf){
    SeuratObject::DefaultAssay(seurat_data) <- assay                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    all_putative_markers <- NULL
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        base::tryCatch(
            expr = {
                SeuratObject::Idents(seurat_data) <- paste(resolution_prefix, resolution, sep=".")
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

                if (base::nrow(markers) > 0) {
                    markers <- markers %>% base::cbind(resolution=resolution, .)
                } else {
                    markers <- markers %>% tibble::add_column(resolution=base::numeric(), .before=1)  # safety measure in case markers was empty
                }
                if (!is.null(all_putative_markers)) {
                    all_putative_markers <- base::rbind(all_putative_markers, markers)
                } else {
                    all_putative_markers <- markers
                }
                base::rm(markers)                                                             # remove unused data
            },
            error = function(e){
                base::print(base::paste("Failed to identify putative markers for resolution", resolution, "with error -", e))
            },
            finally = {
                SeuratObject::Idents(seurat_data) <- "new.ident"
            }
        )
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
            n=50,                                                                               # by default computes 50 singular values
            reduction.name="atac_lsi",                                                          # adding "atac_lsi" for consistency
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
                    dims.use=args$dimensions,
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
            n=50,                                                                               # by default computes 50 singular values
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
                n=50,                                                                           # by default computes 50 singular values
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
            dims.to.integrate=args$dimensions                                                   # don't need to use more than we are going to use in UMAP
        )
        processed_seurat_data <- integrated_seurat_data
        base::rm(integration_anchors, integrated_seurat_data)                                   # remove unused data
    }
    base::rm(splitted_seurat_data)
    base::gc(verbose=FALSE)
    return (processed_seurat_data)
}

atac_analyze <- function(seurat_data, args){
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"                           # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    backup_reductions <- c()                                                    # ATAC integration main remove RNA related reductions so we need to back them up
    for (reduction_name in c("pca", "rnaumap", "wnnumap")){
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

get_de_sample_data <- function(seurat_data, samples_order, args){
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

get_norm_counts_data <- function(deseq_data, sample_data, args){
    if (args$norm == "vst"){
        base::print("Applying vst transformation (not blind to the experimental design)")
        norm_counts_data <- DESeq2::vst(deseq_data, blind=FALSE)
    } else {
        base::print("Applying rlog transformation (not blind to the experimental design)")
        norm_counts_data <- DESeq2::rlog(deseq_data, blind=FALSE)
    }
    if(!is.null(args$batchby) && !is.null(args$remove) && args$remove){
        base::print("Removing batch effect from the normalized counts")
        SummarizedExperiment::assay(norm_counts_data) <- limma::removeBatchEffect(
            SummarizedExperiment::assay(norm_counts_data),
            batch=norm_counts_data[[args$batchby]],
            design=stats::model.matrix(stats::as.formula(base::paste("~",args$splitby)), sample_data)  # should include only splitby
        )
    }
    base::print("Normalized read counts")
    base::print(utils::head(SummarizedExperiment::assay(norm_counts_data)))
    base::print(dim(SummarizedExperiment::assay(norm_counts_data)))
    base::print(SummarizedExperiment::colData(norm_counts_data))
    return (norm_counts_data)
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

rna_de_analyze <- function(seurat_data, args, excluded_genes=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                # safety measure
    base::print(base::paste("Aggregating raw gene expression counts by dataset"))
    selected_genes <- base::as.vector(as.character(base::rownames(seurat_data)))    # all available genes
    if (!is.null(excluded_genes) && length(excluded_genes) > 0){
        base::print(base::paste("Excluding", length(excluded_genes), "genes"))
        selected_genes <- selected_genes[!(selected_genes %in% excluded_genes)]
    }
    aggregated_seurat_data <- get_aggregated_expession(
        seurat_data,
        group_by="new.ident",                                                       # aggregating by sample
        selected_genes=selected_genes,
        slot="counts"                                                               # use not normalized counts
    )
    raw_counts <- SeuratObject::GetAssayData(                                       # dgCMatrix
        aggregated_seurat_data,
        assay="RNA",                                                                # set assay in case GetAssayData won't take the default value
        slot="counts"                                                               # we need not normalized counts
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
    if (args$lrt){
        reduced_formula <- stats::as.formula("~1")
        if(!is.null(args$batchby)){
            reduced_formula <- stats::as.formula(base::paste0("~", args$batchby))
        }
        base::print(                                                               # see this post for details https://support.bioconductor.org/p/95493/#95572
            base::paste(
                "Using LRT test with the design formula", base::paste(design_formula, collapse=""),
                "and the reduced formula", base::paste(reduced_formula, collapse=""),
                "to calculate p-values."
            )
        )
        deseq_data <- DESeq2::DESeq(
            deseq_data,
            test="LRT",
            reduced=reduced_formula,
            quiet=TRUE,
            parallel=TRUE,
            BPPARAM=BiocParallel::MulticoreParam(args$cpus)  # add it here as well just in case
        )
    } else {
        base::print(
            base::paste(
                "Using Wald test with the design formula", base::paste(design_formula, collapse=""),
                "to calculate p-values."
            )
        )
        deseq_data <- DESeq2::DESeq(
            deseq_data,
            quiet=TRUE,
            parallel=TRUE,
            BPPARAM=BiocParallel::MulticoreParam(args$cpus)  # add it here as well just in case
        )
    }
    print("Estimated effects")
    base::print(DESeq2::resultsNames(deseq_data))

    de_genes <- DESeq2::results(
        deseq_data,
        contrast=c(args$splitby, args$second, args$first),            # we are interested in seconds vs first fold change values
        alpha=base::ifelse(is.null(args$padj), 0.1, args$padj),       # recommended to set to our FDR threshold https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
        parallel=TRUE,
        BPPARAM=BiocParallel::MulticoreParam(args$cpus)                             # add it here as well just in case
    )
    base::print("Results description")
    base::print(S4Vectors::mcols(de_genes))
    base::print(utils::head(de_genes))

    de_genes <- base::as.data.frame(de_genes) %>%
                stats::na.omit() %>%                           # exclude all rows where NA is found in any column. See http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
                tibble::rownames_to_column(var="gene")
    base::print(
        base::paste(
            "Number of DE genes after excluding 'NA':", base::nrow(de_genes)
        )
    )

    base::print("Normalizing read count data")
    norm_counts_data <- get_norm_counts_data(deseq_data, sample_data, args)

    base::print(
        base::paste(
            "Creating filtered normalized read counts matrix to include",
            "only differentially expressed features with padj <= ", args$padj
        )
    )
    row_metadata <- de_genes %>%
                    tibble::remove_rownames() %>%
                    tibble::column_to_rownames("gene") %>%
                    dplyr::select(log2FoldChange, pvalue, padj)  %>%                 # we are interested only in these three columns
                    dplyr::filter(.$padj<=args$padj) %>%
                    dplyr::arrange(desc(log2FoldChange))
    col_metadata <- sample_data %>%
                    dplyr::mutate_at(base::colnames(.), base::as.vector)             # need to convert to vector, because in our sample_data everything was a factor
    norm_counts_mat <- SummarizedExperiment::assay(norm_counts_data)[base::as.vector(base::rownames(row_metadata)), ]
    base::print("Size of the normalized read counts matrix after filtering")
    base::print(dim(norm_counts_mat))

    if (!is.null(args$cluster)){
        if (args$cluster == "column" || args$cluster == "both") {
            base::print("Clustering filtered read counts by columns")
            clustered_data = get_clustered_data(
                expression_data=norm_counts_mat,
                center_row=FALSE,                                                    # centering doesn't influence on the samples order
                dist=args$columndist,
                transpose=TRUE
            )
            col_metadata <- base::cbind(col_metadata, clustered_data$clusters)       # adding cluster labels
            col_metadata <- col_metadata[clustered_data$order, ]                     # reordering samples order based on the HOPACH clustering resutls
            base::print("Reordered samples")
            base::print(col_metadata)
        }
        if (args$cluster == "row" || args$cluster == "both") {
            base::print("Clustering filtered normalized read counts by rows")
            clustered_data = get_clustered_data(
                expression_data=norm_counts_mat,
                center_row=base::ifelse(is.null(args$center), FALSE, args$center),   # about centering normalized data https://www.biostars.org/p/387863/
                dist=args$rowdist,
                transpose=FALSE
            )
            norm_counts_mat <- clustered_data$expression                             # can be different because of optional centering by rows mean
            row_metadata <- base::cbind(row_metadata, clustered_data$clusters)       # adding cluster labels
            row_metadata <- row_metadata[clustered_data$order, ]                     # reordering features order based on the HOPACH clustering results
            base::print("Reordered features")
            base::print(utils::head(row_metadata))
        }
        base::rm(clustered_data)
    }

    base::rm(raw_counts, aggregated_seurat_data, selected_genes)                       # remove unused data
    base::gc(verbose=FALSE)

    return (
        list(
            de_genes=de_genes,                     # not filtered differentialy expressed genes
            de_data=deseq_data,                    # raw DESeq output
            sample_data=sample_data,               # we return sample_data even if we can get the same as colData(de_data)
            norm_counts_data=norm_counts_data,     # not filtered trasformed DESeq counts data, includes all genes
            norm_counts_mat=norm_counts_mat,       # filtered to only signif. genes normalized counts matrix (not reordered, for order refer to row/col_metadata)
            row_metadata=row_metadata,             # filtered to only signif. genes ordered based on clusters row metadata for normalized counts matrix
            col_metadata=col_metadata              # filtered to only signif. genes ordered based on clusters column metadata for normalized counts matrix
        )
    )
}

da_analyze <- function(seurat_data, args){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                # safety measure
    base::print("Defining experimental design")
    idents <- base::as.vector(as.character(SeuratObject::Idents(seurat_data)))
    sample_data <- get_de_sample_data(
        seurat_data=seurat_data,
        samples_order=unique(idents),                                               # we don't really care about idents order here
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