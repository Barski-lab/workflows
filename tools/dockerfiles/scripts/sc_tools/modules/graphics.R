import("grid", attach=FALSE)
import("dplyr", attach=FALSE)
import("purrr", attach=FALSE)
import("tidyr", attach=FALSE)
import("knitr", attach=FALSE)
import("scales", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("tibble", attach=FALSE)
import("traviz", attach=FALSE)
import("Glimma", attach=FALSE)
import("ggplot2", attach=FALSE)
import("dynwrap", attach=FALSE)
import("dynplot", attach=FALSE)
import("ggrepel", attach=FALSE)
import("circlize", attach=FALSE)
import("cluster", attach=FALSE)
import("reshape2", attach=FALSE)
import("dittoSeq", attach=FALSE)
import("Nebulosa", attach=FALSE)
import("patchwork", attach=FALSE)
import("ggnewscale", attach=FALSE)
import("tidyselect", attach=FALSE)
import("htmlwidgets", attach=FALSE)
import("scRepertoire", attach=FALSE)
import("RColorBrewer", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)
import("EnhancedVolcano", attach=FALSE)
import("SummarizedExperiment", attach=FALSE)


export(
    "clonotype_quant_plot",
    "clonotype_abundance_plot",
    "clonotype_alluvial_plot",
    "clonotype_feature_plot",
    "clonotype_homeostasis_plot",
    "clonotype_overlap_plot",
#   "clonotype_network_plot",
    "clonotype_chord_plot",
    "clonotype_diversity_plot",
    "filter_by_clonotype",
    "geom_bar_plot",
    "geom_density_plot",
    "geom_point_plot",
    "feature_scatter_plot",
    "vln_plot",
    "dim_plot",
    "elbow_plot",
    "corr_plot",
    "tss_plot",
    "fragments_hist",
    "pca_plot",
    "mds_html_plot",
    "dot_plot",
    "feature_plot",
    "expression_density_plot",
    "dim_heatmap",
    "dim_loadings_plot",
    "silhouette_plot",
    "composition_plot",
    "coverage_plot",
    "volcano_plot",
    "feature_heatmap",
    "daseq_permutations",
    "trajectory_plot",
    "trajectory_graph",
    "trajectory_expression",
    "trajectory_hist",
    "dendro_plot",
    "topology_plot",
    "trajectory_heatmap",
    "expand_qc_suffix",
    "D40_COLORS",
    "D24_COLORS",
    "CC_COLORS",
    "NA_COLOR"
)

# https://sashamaps.net/docs/resources/20-colors/
# https://cran.r-project.org/web/packages/Polychrome/vignettes/testgg.html
D40_COLORS <- c(
  "#FB1C0D", "#0DE400", "#0D00FF", "#E8B4BD", "#FD00EA", "#0DD1FE", "#FF9B0D", "#0D601C",
  "#C50D69", "#CACA16", "#722A91", "#00DEBF", "#863B00", "#5D7C91", "#FD84D8", "#C100FB",
  "#8499FC", "#FD6658", "#83D87A", "#968549", "#DEB6FB", "#832E60", "#A8CAB0", "#FE8F95",
  "#FE1CBB", "#DF7CF8", "#FF0078", "#F9B781", "#4D493B", "#1C5198", "#7C32CE", "#EFBC16",
  "#7CD2DE", "#B30DA7", "#9FC0F6", "#7A940D", "#9B0000", "#946D9B", "#C8C2D9", "#94605A",
  "#34A1EB", "#FFA07A", "#FFD700", "#8A2BE2", "#7FFF00", "#D2691E", "#FF7F50", "#6495ED",
  "#DC143C", "#00FFFF", "#00008B", "#008B8B", "#B8860B", "#A9A9A9", "#006400", "#BDB76B",
  "#8B008B", "#556B2F", "#FF8C00", "#9932CC", "#8B0000", "#E9967A", "#8FBC8F", "#483D8B",
  "#2F4F4F", "#00CED1", "#9400D3", "#FF1493", "#00BFFF", "#696969", "#1E90FF", "#B22222"
)
# D40_COLORS <- c("#00D8B6", "#71E869", "#6574FF", "#F3A6B5", "#FF5AD6", "#6DDCFE", "#FFBB70", "#43A14E", "#D71C7C", "#E1E333", "#8139A8", "#FF6E6A", "#B55C00", "#7FA4B6", "#FFA4E3", "#B300FF", "#9BC4FD", "#FF7E6A", "#9DE98D", "#BFA178", "#E7C2FD", "#8B437D", "#ADCDC0", "#FE9FA4", "#FF53D1", "#D993F9", "#FF47A1", "#FFC171", "#625C51", "#4288C9", "#9767D4", "#F2D61D", "#8EE6FD", "#B940B1", "#B2D5F8", "#9AB317", "#C70000", "#AC8BAC", "#D7D1E4", "#9D8D87")
D24_COLORS <- c(
  "#4E79A7", "#F28E2C", "#4EAE4B", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#1B9E77",
  "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
  "#66A61E", "#E6AB02", "#A6761D", "#1B9E77"
)
NA_COLOR <- "#E1F6FF"
CC_COLORS <- c("#FB1C0D", "#0DE400", "#0D00FF", NA_COLOR)        # we added NA color, because sometimes the cell cycle phase is not being assigned

get_theme <- function(theme){
    return (
        switch(
            theme,
            "gray"     = ggplot2::theme_gray(),
            "bw"       = ggplot2::theme_bw(),
            "linedraw" = ggplot2::theme_linedraw(),
            "light"    = ggplot2::theme_light(),
            "dark"     = ggplot2::theme_dark(),
            "minimal"  = ggplot2::theme_minimal(),
            "classic"  = ggplot2::theme_classic(),
            "void"     = ggplot2::theme_void()
        )
    )
}

expand_qc_suffix <- function(suffix){
    return (
        switch(
            suffix,
            "raw"      = "Unfiltered",
            "mid_fltr" = "Unfiltered, after MACS2 peak calling",
            "fltr"     = "Filtered"
        )
    )
}

filter_by_clonotype <- function(data, clone_by, chain, min_frequency){
    clonotype_column <- scRepertoire:::.convertClonecall(clone_by)
    frequency_column <- base::paste0("clonalFrequency_", chain)
    data@meta.data <- data@meta.data %>%
                      dplyr::mutate(
                          !!clonotype_column:=dplyr::if_else(            # we don't need to update other columns, because we set cloneSize based on it
                              .[[frequency_column]] < min_frequency,
                              NA,
                              .[[clonotype_column]]
                          )
                      ) %>%
                      dplyr::mutate(
                          cloneSize=.[[clonotype_column]]                # scRepertoire uses it for filtering to exclude NA
                      )
    return (data)
}

clonotype_quant_plot <- function(data, rootname, clone_by, chains, group_by, x_label, y_label, legend_title, plot_title, plot_subtitle=NULL, min_frequency=0, palette_colors=D40_COLORS, combine_guides=NULL, ncol=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plots <- list()
            for (i in 1:length(chains)){                                            # should work fine even if chains is not a vector
                current_chain <- chains[i]
                plots[[current_chain]] <- scRepertoire::clonalQuant(
                    input.data=filter_by_clonotype(data, clone_by, current_chain, min_frequency),
                    cloneCall=clone_by,
                    chain=current_chain,
                    group.by=group_by,
                    scale=TRUE                                           # to always show percentage
                ) +
                ggplot2::xlab(x_label) +
                ggplot2::ylab(y_label) +
                ggplot2::ylim(c(0, 100)) +                               # to have the same scale when multiple chains are shown
                ggplot2::guides(
                    fill=ggplot2::guide_legend(legend_title),
                    x=ggplot2::guide_axis(angle=45)
                ) +
                ggplot2::scale_fill_manual(values=palette_colors) +
                ggplot2::ggtitle(current_chain) +
                get_theme(theme)
            }
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol) +
                              patchwork::plot_annotation(title=plot_title, subtitle=plot_subtitle)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting clonotype quant plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export clonotype quant plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

clonotype_abundance_plot <- function(data, rootname, clone_by, chains, group_by, x_label, y_label, legend_title, plot_title, plot_subtitle=NULL, min_frequency=0, scale=TRUE, palette_colors=D40_COLORS, combine_guides=NULL, ncol=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plots <- list()
            for (i in 1:length(chains)){          # should work fine even if chains is not a vector
                current_chain <- chains[i]
                plots[[current_chain]] <- scRepertoire::clonalAbundance(
                    input.data=filter_by_clonotype(data, clone_by, current_chain, min_frequency),
                    cloneCall=clone_by,
                    chain=current_chain,
                    group.by=group_by,
                    scale=scale
                ) +
                ggplot2::xlab(x_label) +
                ggplot2::ylab(y_label) +
                ggplot2::guides(
                    fill=ggplot2::guide_legend(legend_title),
                    color=ggplot2::guide_legend(legend_title),
                    x=ggplot2::guide_axis(angle=45)
                ) +
                ggplot2::scale_fill_manual(values=palette_colors) +
                ggplot2::scale_color_manual(values=palette_colors) +
                ggplot2::ggtitle(current_chain) +
                get_theme(theme)
            }
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol) +
                              patchwork::plot_annotation(title=plot_title, subtitle=plot_subtitle)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting clonotype abundance plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export clonotype abundance plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

clonotype_alluvial_plot <- function(data, rootname, clone_by, chains, group_by, x_label, y_label, legend_title, plot_title, plot_subtitle=NULL, min_frequency=0, top_clones=NULL, legend_position="bottom", palette_colors=D40_COLORS, combine_guides=NULL, ncol=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plots <- list()
            start_from <- 1
            for (i in 1:length(chains)){                                                               # should work fine even if chains is not a vector
                current_chain <- chains[i]
                plots[[current_chain]] <- scRepertoire::clonalCompare(
                    input.data=filter_by_clonotype(data, clone_by, current_chain, min_frequency),
                    cloneCall=clone_by,
                    chain=current_chain,
                    top.clones=top_clones,                                                             # if it's too many, we won't have enough colors and plot will fail
                    group.by=group_by,
                    graph="alluvial"
                ) +
                ggplot2::xlab(x_label) +
                ggplot2::ylab(y_label) +
                ggplot2::ylim(c(0, 1)) +                                                               # to have the same scale when multiple chains are shown
                ggplot2::guides(
                    fill=ggplot2::guide_legend(base::paste(legend_title, current_chain)),
                    x=ggplot2::guide_axis(angle=45)
                ) +
                ggplot2::ggtitle(current_chain) +
                get_theme(theme)

                colors_needed <- length(base::unique(ggplot2::ggplot_build(plots[[current_chain]])$data[[1]]$fill))
                plots[[current_chain]] <- plots[[current_chain]] +
                                          ggplot2::scale_fill_manual(
                                              values=palette_colors[start_from:(start_from+colors_needed-1)]
                                          )
                start_from <- start_from+colors_needed
            }

            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol) +
                              patchwork::plot_annotation(
                                  title=plot_title,
                                  subtitle=plot_subtitle,
                                  theme=ggplot2::theme(legend.position=legend_position)
                              )

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting clonotype alluvial plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export clonotype alluvial plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

clonotype_homeostasis_plot <- function(data, rootname, clone_by, chains, group_by, x_label, y_label, legend_title, plot_title, plot_subtitle=NULL, min_frequency=0, palette_colors=D40_COLORS, combine_guides=NULL, ncol=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plots <- list()
            for (i in 1:length(chains)){          # should work fine even if chains is not a vector
                current_chain <- chains[i]
                plots[[current_chain]] <- scRepertoire::clonalHomeostasis(
                    input.data=filter_by_clonotype(data, clone_by, current_chain, min_frequency),
                    cloneCall=clone_by,
                    chain=current_chain,
                    group.by=group_by
                ) +
                ggplot2::xlab(x_label) +
                ggplot2::ylab(y_label) +
                ggplot2::guides(fill=ggplot2::guide_legend(legend_title), x=ggplot2::guide_axis(angle=45)) +
                ggplot2::scale_fill_manual(values=palette_colors) +
                ggplot2::ggtitle(current_chain) +
                get_theme(theme)
            }
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol) +
                              patchwork::plot_annotation(
                                  title=plot_title,
                                  subtitle=plot_subtitle
                              )

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting clonotype homeostasis plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export clonotype homeostasis plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

clonotype_overlap_plot <- function(data, rootname, clone_by, chains, group_by, x_label, y_label, plot_title, plot_subtitle=NULL, min_frequency=0, methods="morisita", gradient_colors=c("lightgrey", "blue"), na_color="white", combine_guides=NULL, ncol=NULL, nrow=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plots <- list()
            start_from <- 1                                                                      # to have different colors for different methods
            for (i in 1:length(methods)){
                current_method <- methods[i]
                current_range <- NA                                                              # should be independent for each method
                for (j in 1:length(chains)){
                    current_chain <- chains[j]
                    current_plot <- base::paste(current_method, current_chain, sep="_")
                    plots[[current_plot]] <- scRepertoire::clonalOverlap(
                        input.data=filter_by_clonotype(data, clone_by, current_chain, min_frequency),
                        cloneCall=clone_by,
                        method=current_method,
                        chain=current_chain,
                        group.by=group_by
                    ) +
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    ggplot2::guides(
                        x=ggplot2::guide_axis(angle=45),
                        y=ggplot2::guide_axis(angle=90)
                    ) +
                    ggplot2::ggtitle(base::paste(current_chain)) +
                    ggplot2::coord_fixed(ratio=1) +
                    get_theme(theme)

                    if (current_method == "raw"){                                                    # we want don't know the ranges only for raw
                        current_range <- range(
                            as.numeric(                                                              # doesn't fail with NA
                                c(
                                    current_range,                                                   # it's either NA in the first iteration or a vector with the previous min max
                                    ggplot2::ggplot_build(plots[[current_plot]])$data[[3]]$label     # it's ok to take labels only for raw, because labels are rounded and raw are integers
                                )
                            ),
                            na.rm=TRUE
                        )
                    } else if (current_method == "cosine") {
                        current_range <- c(-1, 1)
                    } else {
                        current_range <- c(0, 1)
                    }
                }

                if (length(gradient_colors) >= (2 * length(methods))){                               # check if we have enough colors (2 per method)
                    current_colors <- gradient_colors[start_from:(start_from+1)]
                    start_from <- start_from + 2
                } else {
                    current_colors <- gradient_colors                                                # apply the same colors for all methods
                }

                for (j in 1:length(chains)){                                                         # we need to rescale all plots to the same limits otherwise legend is not combined
                    current_chain <- chains[j]
                    current_plot <- base::paste(current_method, current_chain, sep="_")
                    plots[[current_plot]] <- plots[[current_plot]] +
                                             ggplot2::scale_fill_gradientn(
                                                 colors=current_colors,
                                                 limits=current_range,
                                                 na.value=na_color
                                             ) 
                }
            }

            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol, nrow=nrow) +
                              patchwork::plot_annotation(
                                  title=plot_title,
                                  subtitle=plot_subtitle
                              )

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting clonotype overlap plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export clonotype overlap plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

# clonotype_network_plot <- function(data, reduction, rootname, clone_by, chains, group_by, legend_title, plot_title, plot_subtitle=NULL, min_frequency=0, fixed=TRUE, pt_size=1, alpha=1, legend_position="bottom", palette_colors=D40_COLORS, combine_guides=NULL, ncol=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
#     base::tryCatch(
#         expr = {
#             library("ggraph")                                             # need to include it as library, otherwise clonalNetwork fails
#             groups <- base::sort(
#                 base::unique(
#                     base::as.vector(
#                         as.character(data@meta.data[[group_by]])
#                     )
#                 )
#             )
#             plots <- list()
#             for (i in 1:length(chains)){
#                 current_chain <- chains[i]
#                 current_data <- filter_by_clonotype(                       # may have different number of clonotypes depending on the selected clone_by and chain
#                     data=data,
#                     clone_by=clone_by,
#                     chain=current_chain,
#                     min_frequency=min_frequency
#                 )

#                 for (j in 1:length(groups)){
#                     current_group <- groups[j]
#                     current_plot <- base::paste(current_chain, current_group, sep="_")
#                     plots[[current_plot]] <- scRepertoire::clonalNetwork(
#                         sc.data=current_data,
#                         reduction=reduction,
#                         group.by=group_by,
#                         cloneCall=clone_by,
#                         filter.identity=current_group,
#                         chain=current_chain
#                     ) +
#                     ggplot2::scale_color_manual(
#                         values=palette_colors
#                     ) +
#                     ggplot2::guides(                                          # need to add new guides because they were completely removed in clonalNetwork
#                         colour=ggplot2::guide_legend(
#                             title=legend_title,
#                             title.position="top",
#                             direction="horizontal"
#                         )
#                     ) +
#                     ggplot2::ggtitle(base::paste(current_chain, current_group)) +
#                     get_theme(theme)

#                     plots[[current_plot]]$layers[[1]]$aes_params$size <- pt_size      # we need to add at least some size, because alpha doesn't work without it
#                     plots[[current_plot]]$layers[[1]]$aes_params$alpha <- alpha

#                     if (!is.null(fixed) && fixed){
#                         plots[[current_plot]] <- plots[[current_plot]] + ggplot2::coord_fixed(ratio=1)
#                     }
#                 }
#             }
#             combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol) +
#                               patchwork::plot_annotation(
#                                   title=plot_title,
#                                   subtitle=plot_subtitle,
#                                   theme=ggplot2::theme(legend.position=legend_position)
#                               )

#             grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
#             base::suppressMessages(base::print(combined_plots))
#             grDevices::dev.off()

#             if (!is.null(pdf) && pdf) {
#                 grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
#                 base::suppressMessages(base::print(combined_plots))
#                 grDevices::dev.off()
#             }

#             base::print(base::paste("Exporting clonotype network plot to ", rootname, ".(png/pdf)", sep=""))
#         },
#         error = function(e){
#             base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
#             base::print(base::paste("Failed to export clonotype network plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
#         }
#     )
# }

clonotype_diversity_plot <- function(data, rootname, clone_by, chains, group_by, x_label, y_label, legend_title, plot_title, plot_subtitle=NULL, split_by=NULL, min_frequency=0, palette_colors=D40_COLORS, combine_guides=NULL, ncol=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plots <- list()
            for (i in 1:length(chains)){
                current_chain <- chains[i]
                plots[[current_chain]] <- scRepertoire::clonalDiversity(
                    input.data=filter_by_clonotype(data, clone_by, current_chain, min_frequency),
                    cloneCall=clone_by,
                    chain=current_chain,
                    group.by=group_by,
                    x.axis=split_by
                ) +
                ggplot2::xlab(x_label) +
                ggplot2::ylab(y_label) +
                ggplot2::guides(fill=ggplot2::guide_legend(legend_title), x=ggplot2::guide_axis(angle=45)) +
                ggplot2::scale_fill_manual(values=palette_colors) +
                ggplot2::ggtitle(current_chain) +
                get_theme(theme)

                if (is.null(split_by)){                                       # to have clean x axis if we don't split by anything
                    plots[[current_chain]] <- plots[[current_chain]] +
                        ggplot2::theme(
                            axis.text.x=ggplot2::element_blank(),
                            axis.ticks.x=ggplot2::element_blank()
                        )
                }

            }
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol) +
                              patchwork::plot_annotation(
                                  title=plot_title,
                                  subtitle=plot_subtitle
                              )

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting clonotype diversity plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export clonotype diversity plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

clonotype_feature_plot <- function(data, rootname, clone_by, chains, group_by, x_label, y_label, plot_title, plot_subtitle=NULL, min_frequency=0, order_by="variance", scale=TRUE, combine_guides=NULL, ncol=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plots <- list()
            genes_shown <- c()
            for (i in 1:length(chains)){
                current_chain <- chains[i]
                chain_features <- switch(
                    current_chain,
                    "TRA" = c("TRAV", "TRAJ"),
                    "TRB" = c("TRBV", "TRBD", "TRBJ"),
                    "IGH" = c("IGHV", "IGHJ"),
                    "IGL" = c("IGLV"),
                    "both" = if("TRA" %in% data@misc$vdj$chains)           # we check for "both" in data, because chains input can be just a single value "both"
                             c("TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ")
                             else c("IGHV", "IGHJ", "IGLV")
                )
                current_data <- filter_by_clonotype(                       # may have different number of clonotypes depending on the selected clone_by and chain
                    data=data,
                    clone_by=clone_by,
                    chain=current_chain,
                    min_frequency=min_frequency
                )
                for (j in 1:length(chain_features)){
                    current_feature <- chain_features[j]
                    current_plot <- base::paste(current_chain, current_feature, sep="_")

                    plots[[current_plot]] <- scRepertoire::vizGenes(
                        input.data=current_data,
                        x.axis=current_feature,
                        y.axis=NULL,                                       # will be internally assigned to group_by parameter
                        group.by=group_by,
                        order=order_by,                                    # can be either "variance" or "gene"
                        plot="barplot",
                        scale=scale
                    ) +
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    ggplot2::guides(x=ggplot2::guide_axis(angle=45)) +
                    ggplot2::ggtitle(
                        base::paste0(current_chain, " (", current_feature, ")")
                    ) +
                    get_theme(theme)

                    if (!is.null(scale) && scale){
                        plots[[current_plot]] <- plots[[current_plot]] +    # to have the same scale if we are showing proportion 
                                                 ggplot2::ylim(c(0, 1))
                    }

                    genes_shown <- c(
                        genes_shown,
                        length(
                            ggplot2::ggplot_build(plots[[current_plot]])$layout$panel_params[[1]]$x$get_labels()
                        )
                    )
                }
            }

            combined_plots <- patchwork::wrap_plots(
                                  plots,
                                  guides=combine_guides,
                                  ncol=ncol,
                                  widths=genes_shown[1:ncol]                # safe even when ncol is bigger than length of genes_shown
                              ) +
                              patchwork::plot_annotation(
                                  title=plot_title,
                                  subtitle=plot_subtitle
                              )

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting clonotype feature plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export clonotype feature plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

clonotype_chord_plot <- function(data, rootname, clone_by, group_by, plot_title, proportion=TRUE, min_frequency=0, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=800, height=800, resolution=100){
    base::tryCatch(
        expr = {

            chord_data <- scRepertoire::getCirclize(
                sc=filter_by_clonotype(data, clone_by, "both", min_frequency),   # getCirclize doesn't support chain selection, so we use both here
                cloneCall=clone_by,
                group.by=group_by,
                proportion=proportion
            )
            selected_colors <- palette_colors[1:length(base::unique(base::as.vector(as.character(chord_data$from))))]

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            circlize::chordDiagram(
                chord_data,
                self.link=TRUE,
                directional=TRUE,
                direction.type="arrows",
                link.arr.type="big.arrow",
                annotationTrack=c("name", "grid"),
                grid.col=selected_colors
            )
            graphics::title(main=plot_title, adj=0, line=-0.5)
            circlize::circos.clear()
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                circlize::chordDiagram(
                    chord_data,
                    self.link=TRUE,
                    directional=TRUE,
                    direction.type="arrows",
                    link.arr.type="big.arrow",
                    annotationTrack=c("name", "grid"),
                    grid.col=selected_colors
                )
                graphics::title(main=plot_title, adj=0, line=-0.5)
                circlize::circos.clear()
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                circlize::chordDiagram(
                    chord_data,
                    self.link=TRUE,
                    directional=TRUE,
                    direction.type="arrows",
                    link.arr.type="big.arrow",
                    annotationTrack=c("name", "grid"),
                    grid.col=selected_colors
                )
                graphics::title(main=plot_title, adj=0, line=-0.5)
                circlize::circos.clear()
            }

            base::print(base::paste("Exporting chord diagram to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export chord diagram to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

geom_bar_plot <- function(data, rootname, x_axis, color_by, x_label, y_label, legend_title, plot_title, plot_subtitle=NULL, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plot <- ggplot2::ggplot(data, ggplot2::aes_string(x=x_axis, fill=color_by)) +
                ggplot2::geom_bar(colour="black") +
                ggplot2::geom_text(stat="count", ggplot2::aes(label=..count..), vjust=-1) +
                ggplot2::xlab(x_label) +
                ggplot2::ylab(y_label) +
                ggplot2::guides(fill=ggplot2::guide_legend(legend_title), x=ggplot2::guide_axis(angle=45)) +
                ggplot2::ggtitle(plot_title, subtitle=plot_subtitle) +
                ggplot2::scale_fill_manual(values=palette_colors) +
                get_theme(theme)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting geom bar plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export geom bar plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

geom_density_plot <- function(data, rootname, x_axis, group_by, split_by, x_label, y_label, legend_title, plot_title, plot_subtitle=NULL, x_left_intercept=NULL, x_right_intercept=NULL,  scale_x_log10=FALSE, scale_y_log10=FALSE, alpha=0.7, show_zoomed=FALSE, show_ranked=FALSE, ranked_x_label="Ranked cells", palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            intercept_data <- data %>%
                              dplyr::select(tidyselect::all_of(group_by), tidyselect::all_of(split_by)) %>%
                              dplyr::distinct() %>%
                              dplyr::arrange(dplyr::across(tidyselect::all_of(group_by))) %>%
                              tibble::add_column(color=palette_colors[1:base::nrow(.)])

            plot <- ggplot2::ggplot(data, ggplot2::aes_string(x=x_axis, fill=group_by)) +
                    ggplot2::geom_density(alpha=alpha, trim=TRUE) +                           # set trim to TRUE to make it correspond to what we have in violin plots
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    ggplot2::guides(fill=ggplot2::guide_legend(legend_title)) +
                    ggplot2::ggtitle(plot_title, subtitle=plot_subtitle) +
                    ggplot2::facet_wrap(stats::as.formula(base::paste("~", split_by))) +
                    ggplot2::scale_fill_manual(values=palette_colors) +
                    get_theme(theme)

            if (!is.null(x_left_intercept)){
                intercept_data <- intercept_data %>% tibble::add_column(x_left=x_left_intercept)
                plot <- plot +
                        ggplot2::geom_vline(intercept_data, mapping=ggplot2::aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                        ggrepel::geom_label_repel(
                            intercept_data, mapping=ggplot2::aes(x=x_left, y=Inf, label=x_left),
                            color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                            show.legend=FALSE
                        )
            }

            if (!is.null(x_right_intercept)){
                intercept_data <- intercept_data %>% tibble::add_column(x_right=x_right_intercept)
                plot <- plot +
                        ggplot2::geom_vline(intercept_data, mapping=ggplot2::aes(xintercept=x_right), color=intercept_data$color, alpha=0.7) +
                        ggrepel::geom_label_repel(
                            intercept_data, mapping=ggplot2::aes(x=x_right, y=Inf, label=x_right),
                            color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                            show.legend=FALSE
                        )
            }

            if (scale_x_log10){ plot <- plot + ggplot2::scale_x_log10() + ggplot2::annotation_logticks(sides="b", alpha=0.3) }
            if (scale_y_log10){ plot <- plot + ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides="l", alpha=0.3) }

            if (show_zoomed){
                zoomed_plot <- ggplot2::ggplot(data, ggplot2::aes_string(x=x_axis, color=group_by)) +
                               ggplot2::geom_density(show.legend=FALSE, size=2) +
                               ggplot2::xlab(x_label) +
                               ggplot2::ylab(y_label) +
                               ggplot2::scale_color_manual(values=palette_colors) +
                               get_theme(theme)

                if (!is.null(x_left_intercept)){
                    zoomed_plot <- zoomed_plot +
                                   ggplot2::geom_vline(intercept_data, mapping=ggplot2::aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                                   ggrepel::geom_label_repel(
                                       intercept_data, mapping=ggplot2::aes(x=x_left, y=Inf, label=x_left),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                       show.legend=FALSE
                                   )
                }
                if (!is.null(x_right_intercept)){
                    zoomed_plot <- zoomed_plot +
                                   ggplot2::geom_vline(intercept_data, mapping=ggplot2::aes(xintercept=x_right), color=intercept_data$color, alpha=0.7) +
                                   ggrepel::geom_label_repel(
                                       intercept_data, mapping=ggplot2::aes(x=x_right, y=Inf, label=x_right),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                       show.legend=FALSE
                                   )
                }
                zoomed_plot <- zoomed_plot +
                               ggplot2::coord_cartesian(
                                   xlim=c(
                                       base::ifelse(!is.null(x_left_intercept), min(intercept_data$x_left), NA),
                                       base::ifelse(!is.null(x_right_intercept), max(intercept_data$x_right), NA)
                                   )
                               )

                if (scale_x_log10){ zoomed_plot <- zoomed_plot + ggplot2::scale_x_log10() + ggplot2::annotation_logticks(sides="b", alpha=0.3) }
                if (scale_y_log10){ zoomed_plot <- zoomed_plot + ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides="l", alpha=0.3) }
                plot <- plot / zoomed_plot
            }

            if (show_ranked) {
                ranked_plot <- data %>%
                               dplyr::arrange(dplyr::across(tidyselect::all_of(group_by)), dplyr::across(tidyselect::all_of(x_axis))) %>%
                               ggplot2::ggplot(ggplot2::aes_string(x=seq_along(data[[x_axis]]), y=x_axis, color=group_by)) +
                               ggplot2::geom_point(show.legend=FALSE, size=0.5) +
                               ggplot2::xlab(ranked_x_label) +
                               ggplot2::ylab(x_label) +
                               ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides="l", alpha=0.3) +
                               ggplot2::scale_color_manual(values=palette_colors) +
                               get_theme(theme)

                if (!is.null(x_left_intercept)){
                    ranked_plot <- ranked_plot +
                                   ggplot2::geom_hline(intercept_data, mapping=ggplot2::aes(yintercept=x_left), color=intercept_data$color, alpha=0.7) +
                                       ggrepel::geom_label_repel(
                                       intercept_data, mapping=ggplot2::aes(x=Inf, y=x_left, label=x_left),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="x", size=3,
                                       show.legend=FALSE
                                   )
                }

                if (!is.null(x_right_intercept)){
                    ranked_plot <- ranked_plot +
                                   ggplot2::geom_hline(intercept_data, mapping=ggplot2::aes(yintercept=x_right), color=intercept_data$color, alpha=0.7) +
                                   ggrepel::geom_label_repel(
                                       intercept_data, mapping=ggplot2::aes(x=-Inf, y=x_right, label=x_right),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="x", size=3,
                                       show.legend=FALSE
                                   )
                }
                plot <- plot / ranked_plot
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting geom density plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export geom density plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

geom_point_plot <- function(
    data, rootname,
    x_axis, y_axis,
    split_by, color_by,
    gradient_colors, color_limits, color_break,
    x_label, y_label,
    legend_title, plot_title, plot_subtitle=NULL,
    x_left_intercept=NULL, y_low_intercept=NULL, y_high_intercept=NULL,
    highlight_rows=NULL, highlight_color="black", highlight_shape=4,
    scale_x_log10=FALSE, scale_y_log10=FALSE,
    show_lm=FALSE, show_density=FALSE, density_bins=6,
    alpha=0.2, alpha_intercept=0.5,
    palette_colors=D40_COLORS, theme="classic",
    pdf=FALSE, width=1200, height=800, resolution=100
){
    base::tryCatch(
        expr = {
            intercept_data <- data %>%
                              dplyr::select(tidyselect::all_of(split_by)) %>%
                              dplyr::distinct() %>%
                              dplyr::arrange(dplyr::across(tidyselect::all_of(split_by))) %>%
                              tibble::add_column(color=palette_colors[1:base::nrow(.)])

            if (!is.null(x_left_intercept)){
                intercept_data <- intercept_data %>% tibble::add_column(x_left=x_left_intercept)
            }
            if (!is.null(y_low_intercept)){
                intercept_data <- intercept_data %>% tibble::add_column(y_low=y_low_intercept)
            }
            if (!is.null(y_high_intercept)){
                intercept_data <- intercept_data %>% tibble::add_column(y_high=y_high_intercept)
            }

            plot <- ggplot2::ggplot(data, ggplot2::aes_string(x=x_axis, y=y_axis, color=color_by)) +
                    ggplot2::geom_point(alpha=alpha) +
                    ggplot2::scale_colour_gradientn(
                        colours=c(gradient_colors[1], gradient_colors),
                        values=scales::rescale(c(color_limits[1], color_break-0.01*color_break, color_break, color_limits[2])),
                        breaks=c(color_break),
                        limits=color_limits
                    ) +
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    ggplot2::guides(color=ggplot2::guide_colourbar(legend_title)) +
                    ggplot2::ggtitle(plot_title, subtitle=plot_subtitle) +
                    ggplot2::facet_wrap(stats::as.formula(base::paste("~", split_by))) +
                    get_theme(theme)

            if (show_lm){
                plot <- plot +
                        ggplot2::stat_smooth(
                            method=stats::lm,
                            linetype="dashed",
                            colour="black",
                            size=0.25
                        )
            }
            if (!is.null(highlight_rows) && length(highlight_rows) > 0){
                plot <- plot + ggplot2::geom_point(
                    data=data[highlight_rows, ],
                    shape=highlight_shape,
                    color=highlight_color,
                    alpha=0.5
                )
            }
            if (show_density){                                              # we want to have it shown above the highlighted cells
                plot <- plot +
                        ggplot2::geom_density_2d(                           # this densities can't be used for comparing the number of cell between the facets
                            linetype="solid",
                            bins=density_bins,
                            contour_var="ndensity",                         # scales densities to a maximum of 1 making the number of contours the same in each facet
                            colour="black",                                 # all contours will be black
                            size=0.1
                        )
            }
            if (!is.null(x_left_intercept)){
                plot <- plot +
                        ggplot2::geom_vline(intercept_data, mapping=ggplot2::aes(xintercept=x_left), color=intercept_data$color, alpha=alpha_intercept) +
                        ggrepel::geom_label_repel(
                            intercept_data, mapping=ggplot2::aes(x=x_left, y=Inf, label=x_left),
                            color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="y", size=3,
                            show.legend=FALSE
                        )
            }
            if (!is.null(y_low_intercept)){
                plot <- plot +
                        ggplot2::geom_hline(intercept_data, mapping=ggplot2::aes(yintercept=y_low), color=intercept_data$color, alpha=alpha_intercept) +
                        ggrepel::geom_label_repel(
                            intercept_data, mapping=ggplot2::aes(x=Inf, y=y_low, label=y_low),
                            color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="x", size=3,
                            show.legend=FALSE
                        )
            }
            if (!is.null(y_high_intercept)){
                plot <- plot +
                        ggplot2::geom_hline(intercept_data, mapping=ggplot2::aes(yintercept=y_high), color=intercept_data$color, alpha=alpha_intercept) +
                        ggrepel::geom_label_repel(
                            intercept_data, mapping=ggplot2::aes(x=Inf, y=y_high, label=y_high),
                            color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="x", size=3,
                            show.legend=FALSE
                        )
            }

            if (scale_x_log10){ plot <- plot + ggplot2::scale_x_log10() + ggplot2::annotation_logticks(sides="b", alpha=0.3) }
            if (scale_y_log10){ plot <- plot + ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides="l", alpha=0.3) }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting geom point plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export geom point plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

feature_scatter_plot <- function(data, rootname, x_axis, y_axis, x_label, y_label, split_by, color_by, plot_title, legend_title, combine_guides=NULL, palette_colors=NULL, alpha=NULL, jitter=FALSE, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            SeuratObject::Idents(data) <- split_by
            identities <- base::unique(base::as.vector(as.character(SeuratObject::Idents(data))))
            plots <- list()
            for (i in 1:length(identities)) {
                current_identity <- identities[i]
                filtered_data <- base::subset(data, idents=current_identity)
                plots[[current_identity]] <- Seurat::FeatureScatter(
                    filtered_data,
                    feature1=x_axis,
                    feature2=y_axis,
                    group.by=color_by,
                    plot.cor=FALSE,         # will be overwritten by title anyway
                    jitter=jitter
                )
            }
            SeuratObject::Idents(data) <- "new.ident"
            plots <- base::lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] +
                              ggplot2::ggtitle(identities[i]) +
                              ggplot2::xlab(x_label) +
                              ggplot2::ylab(y_label) +
                              ggplot2::guides(color=ggplot2::guide_legend(legend_title)) +
                              get_theme(theme)
                if (!is.null(palette_colors)) { plots[[i]] <- plots[[i]] + ggplot2::scale_color_manual(values=palette_colors) }
                if (!is.null(alpha)) { plots[[i]]$layers[[1]]$aes_params$alpha <- alpha }
                return (plots[[i]])
            })
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides) + patchwork::plot_annotation(title=plot_title)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting feature scatter plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export feature scatter plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

vln_plot <- function(
    data, features, labels,
    rootname, plot_title, legend_title,
    plot_subtitle=NULL,
    from_meta=FALSE, scale_y_log10=FALSE,
    group_by=NULL, split_by=NULL,
    ncol=NULL, show_box_plots=FALSE, hide_x_text=FALSE, rotate_labels=FALSE,
    legend_position="right", pt_size=NULL, palette_colors=NULL, combine_guides=NULL,
    theme="classic", pdf=FALSE, width=NULL, height=NULL, resolution=100
){
    base::tryCatch(
        expr = {

            features_corrected <- features
            labels_corrected <- labels
            if (length(scale_y_log10) == 1){
                scale_y_log10 <- rep(scale_y_log10, length(features_corrected))
            }
            scale_y_log10_corrected <- scale_y_log10
            if (from_meta){
                features_corrected <- c()
                labels_corrected <- c()
                scale_y_log10_corrected <- c()
                for (i in 1:length(features)){
                    if (features[i] %in% base::colnames(data@meta.data)){
                        features_corrected <- c(features_corrected, features[i])
                        labels_corrected <- c(labels_corrected, labels[i])
                        scale_y_log10_corrected <- c(scale_y_log10_corrected, scale_y_log10[i])
                    } else {
                        base::print(
                            base::paste(
                                "Feature", features[i], "was not found,",
                                "skipping", labels[i]
                            )
                        )
                    }
                }
            }

            plots <- Seurat::VlnPlot(
                         data,
                         features=features_corrected,
                         pt.size=pt_size,
                         group.by=group_by,  # will be used instead of the default identities
                         split.by=split_by,
                         log=FALSE,          # we will scale it later with scale_y_log10 if needed
                         combine=FALSE       # to return a list of gglots
                     )
            plots <- base::lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] +
                              ggplot2::ggtitle(labels_corrected[i]) +
                              get_theme(theme) +
                              ggplot2::theme(
                                  axis.title.x=ggplot2::element_blank(),
                                  legend.position=legend_position
                              ) +
                              ggplot2::guides(fill=ggplot2::guide_legend(legend_title)) +
                              Seurat::RotatedAxis()

                if (show_box_plots){
                    plots[[i]] <- plots[[i]] +
                                  ggplot2::stat_boxplot(width=0.1, geom="errorbar", linewidth=0.25) +
                                  ggplot2::geom_boxplot(width=0.1, outlier.alpha=0, linewidth=0.25) +
                                  ggplot2::stat_summary(                                                 # https://stackoverflow.com/questions/56434187/label-whiskers-on-ggplot-boxplot-when-there-are-outliers
                                      ggplot2::aes(label=scales::comma(round(..y.., digits=2))),
                                      geom="text",
                                      fun.y=function(y) grDevices::boxplot.stats(y)$stats[c(1, 5)],
                                      position=ggplot2::position_nudge(x=0.35), 
                                      size=2.5,
                                      color="darkblue"
                                  )
                }

                if (rotate_labels){                                                                      # removes titles and puts the labels on the Y axis
                    plots[[i]] <- plots[[i]] +
                                  ggplot2::ggtitle(NULL) +
                                  ggplot2::ylab(labels_corrected[i])
                }

                if (scale_y_log10_corrected[i]){
                    plots[[i]] <- plots[[i]] +
                                  ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides="l", alpha=0.3)
                }

                if (!is.null(palette_colors)){
                    plots[[i]] <- plots[[i]] +
                                  ggplot2::scale_fill_manual(values=palette_colors)
                }

                if (hide_x_text){
                    plots[[i]] <- plots[[i]] +
                                  ggplot2::theme(axis.text.x=ggplot2::element_blank())
                }

                return (plots[[i]])
            })

            columns_number <- ifelse(is.null(ncol), ceiling(sqrt(length(plots))), ncol)               # tries to make a square
            combined_plots <- patchwork::wrap_plots(
                                  plots,
                                  guides=combine_guides,
                                  ncol=columns_number
                              ) +
                              patchwork::plot_annotation(
                                  title=plot_title,
                                  subtitle=plot_subtitle,
                                  theme=ggplot2::theme(legend.position=legend_position)
                              )

            if (is.null(width)){
                current_identities <- SeuratObject::Idents(data)
                if (!is.null(group_by)){
                    current_identities <- data@meta.data[[group_by]]                                  # in case we group not by the default identities
                }
                width_scale <- ifelse(show_box_plots, 0.5, 0.3)                                       # need more space when showing stats
                width <- round(
                    length(
                        base::unique(base::as.vector(as.character(current_identities)))
                    ) * columns_number * width_scale * resolution
                )
                width <- base::ifelse(width < 600, 600, width)
            }
            if (is.null(height)){
                height_scale <- ifelse(from_meta, 2, 1.5)                                             # need more space when showing QC metrics from the meta data
                height <- round(
                    ceiling(length(plots)/columns_number) * height_scale * resolution
                )
                height <- base::ifelse(height < 400, 400, height)
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting violin plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export violin plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

dim_plot <- function(data, rootname, reduction, plot_title, legend_title, plot_subtitle=NULL, cells=NULL, split_by=NULL, group_by=NULL, highlight_group=NULL, show_density=FALSE, density_bins=10, ncol=NULL, label=FALSE, label_box=FALSE, label_repel=FALSE, label_color="black", label_size=4, fixed=TRUE, alpha=NULL, palette_colors=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {

            if (
                !is.null(show_density) && show_density &&                             # we want to show density contours
                (!is.null(group_by) || !is.null(split_by))                            # densities will calculated per groups that can be defined by group_by and/or split_by
            ){                                                                        # we need to exclude cells that form groups of only 1 cell
                excluded_groups <- data@meta.data %>%
                                   dplyr::group_by(
                                       dplyr::across(tidyselect::any_of(group_by)),   # at least one of them won't be NULL
                                       dplyr::across(tidyselect::any_of(split_by))
                                   ) %>%
                                   dplyr::summarize(counts=dplyr::n()) %>%            # calculcate cells per group
                                   dplyr::ungroup() %>%
                                   dplyr::filter(.$counts < 2) %>%                    # select groups with only 1 cell
                                   base::as.data.frame()
                excluded_cells <- data@meta.data %>%                                  # vector of barcodes to be excluded, might be character(0)
                                  dplyr::semi_join(
                                      excluded_groups,
                                      by=stats::na.omit(c(split_by, group_by))        # will exclude possible NULL from the vector
                                  ) %>%
                                  tibble::rownames_to_column(var="barcode") %>%
                                  dplyr::pull(barcode)
                data <- base::subset(data, cells=excluded_cells, invert=TRUE)         # safe to run even when excluded_cells is NULL
                cells <- base::subset(cells, !(cells %in% excluded_cells))            # safe to run even when cells and/or excluded_cells are NULL
            }

            highlight_cells <- NULL
            if (!is.null(group_by) && !is.null(highlight_group)){
                SeuratObject::Idents(data) <- group_by
                highlight_cells <- SeuratObject::WhichCells(
                    data,
                    idents=highlight_group                               # highlight_group can be also set as a vector
                )
                SeuratObject::Idents(data) <- "new.ident"                # need to set back to our default identity
            }
            plot <- Seurat::DimPlot(
                        data,
                        reduction=reduction,
                        cells=cells,
                        split.by=split_by,
                        group.by=group_by,
                        label.box=label_box,
                        repel=label_repel,
                        cells.highlight=if(is.null(highlight_cells))         # need to use this trick because ifelse doesn't return NULL, 'Selected' is just a name to display on the plot
                                            NULL
                                        else
                                            list(Selected=highlight_cells),
                        ncol=if(!is.null(split_by) && is.null(ncol))         # attempt to arrage all plots in a square
                                ceiling(
                                    sqrt(
                                        length(base::unique(base::as.vector(as.character(data@meta.data[[split_by]]))))
                                    )
                                )
                             else
                                ncol,
                        label=label,
                        label.color=label_color,
                        label.size=label_size
                    ) +
                    get_theme(theme) +
                    ggplot2::ggtitle(plot_title, subtitle=plot_subtitle) +
                    ggplot2::guides(color=ggplot2::guide_legend(legend_title, override.aes=list(size=3)))

            if (!is.null(fixed) && fixed){
                plot <- plot + ggplot2::coord_fixed(ratio=1)
            }

            if (!is.null(palette_colors)){
                plot <- plot +
                        ggplot2::scale_color_manual(values=palette_colors) +
                        ggplot2::scale_fill_manual(values=palette_colors)     # need it for proper label box background
            }
            if (!is.null(alpha)) { plot$layers[[1]]$aes_params$alpha <- alpha }

            if (!is.null(show_density) && show_density){
                plot <- plot +
                        ggplot2::geom_density_2d(
                            ggplot2::aes_string(
                                x=base::paste0(                             # need to set x and y aesthetics because geom_density_2d can't reach them from the plot
                                    SeuratObject::Key(
                                        data@reductions[[reduction]]
                                    ),
                                    "1"
                                ),
                                y=base::paste0(
                                    SeuratObject::Key(
                                        data@reductions[[reduction]]
                                    ),
                                    "2"
                                ),
                                group=group_by                              # if not NULL, density is calculated per group (a.k.a cluster, cell type)
                            ),
                            bins=density_bins,
                            contour_var="count",                            # contours are made based on the densities * on the cell counts withing each bin
                            linetype="solid",
                            colour="black",                                 # all contours will be black
                            size=0.1
                        )
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            # if (!is.null(htmlwidget) && htmlwidget) {               # in case one day we decide to save html widgets
            #     htmlwidgets::saveWidget(
            #         Seurat::HoverLocator(
            #             plot,
            #             information=SeuratObject::FetchData(
            #                 data,
            #                 vars=c("new.ident", "condition")
            #             ) %>% dplyr::rename("dataset"=new.ident),
            #             axes=FALSE
            #         ),
            #         base::paste(rootname, ".html", sep="")
            #     )
            # }

            base::print(base::paste("Exporting dim plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export dim plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

trajectory_plot <- function(                                # all feature related parameters are skipped
    data, rootname, reduction, plot_title,
    legend_title=NULL,
    color_cells="auto",                                     # one of auto, none, grouping, milestone, pseudotime
    color_density="none",                                   # one of none, grouping    
    alpha=1,
    palette_colors=NULL,
    theme="classic",
    pdf=FALSE,
    width=1200,
    height=800,
    resolution=100
){
    base::tryCatch(
        expr = {

            plot <- dynplot::plot_dimred(
                        trajectory=data@misc$trajectories[[reduction]]$dyno,
                        color_cells=color_cells,
                        color_density=color_density,
                        grouping=if (color_density == "grouping"){
                                    dynwrap::group_onto_nearest_milestones(
                                        data@misc$trajectories[[reduction]]$dyno
                                    )
                                } else { NULL },
                        alpha_cells=alpha
                    ) +
                    ggplot2::ggtitle(plot_title) +
                    get_theme(theme)

            if (!is.null(legend_title)){
                plot <- plot +
                        ggplot2::guides(color=ggplot2::guide_legend(legend_title)) +
                        ggplot2::guides(fill=ggplot2::guide_legend(legend_title))
            }

            if (!is.null(palette_colors) && color_density == "grouping" && color_cells == "grouping"){
                plot <- plot +
                        ggplot2::scale_color_manual(values=palette_colors) +
                        ggplot2::scale_fill_manual(values=palette_colors)
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting trajectory plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export trajectory plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

trajectory_expression <- function(
    data, rootname, reduction, features, plot_title,
    assay="RNA",
    slot="data",
    alpha=0.5,
    combine_guides=NULL,
    theme="classic",
    pdf=FALSE,
    width=1200,
    height=800,
    resolution=100
){
    base::tryCatch(
        expr = {
            plots <- list()
            counts_data <- base::as.matrix(
                SeuratObject::GetAssayData(
                    data,
                    assay=assay,
                    slot=slot
                )
            )
            for (i in 1:length(features)){
                current_feature <- features[i]
                plots[[current_feature]] <- traviz::plotExpression(
                    counts=counts_data,
                    sds=data@misc$trajectories[[reduction]]$slingshot,
                    gene=current_feature,
                    alpha=alpha,
                ) + ggplot2::ggtitle(current_feature)
            }
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides) + patchwork::plot_annotation(title=plot_title)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting expression along pseudotime plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export expression along pseudotime plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

trajectory_hist <- function(
    data, rootname, x_axis, group_by, split_by, x_label, y_label, legend_title, plot_title,
    stack_direction="center", stack_ratio=0.25, method="histodot",
    pt_size=0.5, alpha=0.85, palette_colors=D40_COLORS, theme="classic", pdf=FALSE,
    width=1200, height=800, resolution=100
){
    base::tryCatch(
        expr = {
            plot <- ggplot2::ggplot(data, ggplot2::aes_string(x=x_axis, fill=group_by, color=group_by)) +
                    ggplot2::geom_dotplot(
                        alpha=alpha,
                        method=method,
                        stackdir=stack_direction,
                        stackratio=stack_ratio,
                        dotsize=pt_size,
                        binpositions="all",
                        position=ggplot2::position_jitter(height=0, seed=NULL)
                    ) +
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    ggplot2::guides(
                        fill=ggplot2::guide_legend(legend_title),
                        color=ggplot2::guide_legend(legend_title)
                    ) +
                    ggplot2::ggtitle(plot_title) +
                    ggplot2::facet_wrap(stats::as.formula(base::paste("~", split_by))) +
                    ggplot2::scale_fill_manual(values=palette_colors) +
                    ggplot2::scale_color_manual(values=palette_colors) +
                    get_theme(theme)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting trajectory histogram to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export trajectory histogram to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

trajectory_graph <- function(                               # all feature related parameters are skipped
    data, rootname, reduction, plot_title,
    legend_title=NULL,
    color_cells="auto",                                     # one of auto, none, grouping, milestone, pseudotime
    palette_colors=NULL,
    theme="classic",
    pdf=FALSE,
    width=1200,
    height=800,
    resolution=100
){
    base::tryCatch(
        expr = {

            plot <- dynplot::plot_graph(
                        trajectory=data@misc$trajectories[[reduction]]$dyno,
                        color_cells=color_cells,
                        grouping=if (color_cells == "grouping"){
                                     dynwrap::group_onto_nearest_milestones(
                                         data@misc$trajectories[[reduction]]$dyno
                                     )
                                 } else { NULL }
                    ) +
                    ggplot2::ggtitle(plot_title) +
                    get_theme(theme)

            if (!is.null(legend_title)){
                plot <- plot +
                        ggplot2::guides(color=ggplot2::guide_legend(legend_title)) +
                        ggplot2::guides(fill=ggplot2::guide_legend(legend_title))
            }

            if (!is.null(palette_colors) && color_cells == "grouping"){
                plot <- plot +
                        ggplot2::scale_color_manual(values=palette_colors) +
                        ggplot2::scale_fill_manual(values=palette_colors)
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting trajectory graph to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export trajectory graph to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

dendro_plot <- function(                                    # all feature related parameters are skipped
    data, rootname, reduction, plot_title,
    legend_title=NULL,
    color_cells="auto",                                     # one of auto, none, grouping, milestone, pseudotime
    alpha=1,
    palette_colors=NULL,
    theme="classic",
    pdf=FALSE,
    width=1200,
    height=800,
    resolution=100
){
    base::tryCatch(
        expr = {

            plot <- dynplot::plot_dendro(
                        trajectory=data@misc$trajectories[[reduction]]$dyno,
                        color_cells=color_cells,
                        grouping=if (color_cells == "grouping"){
                                    dynwrap::group_onto_nearest_milestones(
                                        data@misc$trajectories[[reduction]]$dyno
                                    )
                                } else { NULL },
                        alpha_cells=alpha
                    ) +
                    ggplot2::ggtitle(plot_title) +
                    get_theme(theme)

            if (!is.null(legend_title)){
                plot <- plot +
                        ggplot2::guides(color=ggplot2::guide_legend(legend_title)) +
                        ggplot2::guides(fill=ggplot2::guide_legend(legend_title))
            }

            if (!is.null(palette_colors) && color_cells == "grouping"){
                plot <- plot +
                        ggplot2::scale_color_manual(values=palette_colors) +
                        ggplot2::scale_fill_manual(values=palette_colors)
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting dendrogram to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export dendrogram to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

topology_plot <- function(                                    # all feature related parameters are skipped
    data, rootname, reduction, plot_title,
    theme="classic",
    pdf=FALSE,
    width=1200,
    height=800,
    resolution=100
){
    base::tryCatch(
        expr = {

            plot <- dynplot::plot_topology(
                        trajectory=data@misc$trajectories[[reduction]]$dyno
                    ) +
                    ggplot2::ggtitle(plot_title) +
                    get_theme(theme)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting topology plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export topology plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

trajectory_heatmap <- function(                                    # changing the theme makes it look bad
    data, rootname, reduction, plot_title,
    assay="RNA",
    slot="data",
    features=50,                                                   # can be either a number of genes or a list of gene names
    palette_colors=D40_COLORS,                                     # for cluster colors
    pdf=FALSE,
    width=1200,
    height=800,
    resolution=100
){
    base::tryCatch(
        expr = {
            sorted_milestone_ids <- base::sort(data@misc$trajectories[[reduction]]$dyno$milestone_ids)
            plot <- dynplot::plot_heatmap(
                        trajectory=data@misc$trajectories[[reduction]]$dyno,
                        expression_source=base::t(
                            base::as.matrix(
                                SeuratObject::GetAssayData(data, assay=assay, slot=slot)
                            )
                        ),
                        features_oi=features,
                        grouping=dynwrap::group_onto_nearest_milestones(
                            data@misc$trajectories[[reduction]]$dyno
                        ),
                        groups=data.frame(
                            group_id=sorted_milestone_ids,
                            color=palette_colors[0:length(sorted_milestone_ids)],
                            check.names=FALSE
                        )
                    ) +
                    patchwork::plot_annotation(title=plot_title)


            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting trajectory heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export trajectory heatmap to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

elbow_plot <- function(data, rootname, plot_title, reduction="pca", ndims=NULL, x_intercept=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plot <- Seurat::ElbowPlot(
                        data,
                        ndims=ndims,
                        reduction=reduction
                    ) +
                    get_theme(theme) +
                    ggplot2::ggtitle(plot_title)

            if (!is.null(x_intercept)){                                      # won't fail even if x_intercept > ndims
                intercept_data <- base::data.frame(x_coord=x_intercept)      # somehow can't add label without using intercept_data
                plot <- plot +
                        ggplot2::geom_vline(
                            intercept_data,
                            mapping=ggplot2::aes(xintercept=x_coord),
                            color="red",
                            size=1
                        ) +
                        ggrepel::geom_label_repel(
                            intercept_data,
                            mapping=ggplot2::aes(x=x_coord, y=Inf, label=x_coord),
                            color="white",
                            fontface="bold",
                            fill="red",
                            segment.colour=NA,
                            direction="y",
                            size=6,
                            show.legend=FALSE
                        )
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting elbow plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export elbow plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

silhouette_plot <- function(data, rootname, plot_title, legend_title, group_by, dims, downsample=300, reduction="pca", plot_subtitle=NULL, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            SeuratObject::Idents(data) <- group_by
            data <- base::subset(data, downsample=downsample)
            silhouette_data <- cluster::silhouette(
                as.numeric(data@meta.data[, group_by]),
                dist=stats::dist(SeuratObject::Embeddings(data[[reduction]])[, dims])      # we always run PCA with 50 PC, so it's safe to subset based on dims
            )
            data@meta.data$silhouette_score <- silhouette_data[, 3]
            mean_silhouette_score <- base::mean(data@meta.data$silhouette_score)

            plot <- data@meta.data %>%
                    dplyr::mutate(barcode=base::rownames(.)) %>%
                    dplyr::arrange(dplyr::across(tidyselect::all_of(group_by)), -silhouette_score) %>%
                    dplyr::mutate(barcode=base::factor(barcode, levels=barcode)) %>%
                    ggplot2::ggplot() +
                    ggplot2::geom_col(ggplot2::aes_string("barcode", "silhouette_score", fill=group_by)) +
                    ggplot2::geom_hline(yintercept=mean_silhouette_score, color="red", linetype="dashed") +
                    ggplot2::scale_x_discrete(name="Cells") +
                    ggplot2::scale_y_continuous(name="Silhouette score") +
                    ggplot2::scale_fill_manual(values=palette_colors) +
                    get_theme(theme) +
                    ggplot2::theme(
                        axis.text.x=ggplot2::element_blank(),
                        axis.ticks.x=ggplot2::element_blank(),
                        panel.grid.major=ggplot2::element_blank(),
                        panel.grid.minor=ggplot2::element_blank()
                    ) +
                    ggplot2::ggtitle(plot_title, subtitle=plot_subtitle) +
                    ggplot2::guides(fill=ggplot2::guide_legend(legend_title))

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting silhouette plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export silhouette plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

composition_plot <- function(data, rootname, plot_title, legend_title, x_label, y_label, split_by, group_by, bar_position="fill", plot_subtitle=NULL, label=TRUE, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    # bar_position can be one of the following
    #   fill  - stacked, percents are diplayed (default)
    #   dodge - grouped, values are displayed
    base::tryCatch(
        expr = {
            counts_data <- data@meta.data %>%
                           dplyr::group_by(dplyr::across(tidyselect::all_of(split_by)), dplyr::across(tidyselect::all_of(group_by))) %>%
                           dplyr::summarize(counts=dplyr::n()) %>%                                      # uses both groups defined by split_by and group_by
                           tidyr::spread(tidyselect::all_of(group_by), counts, fill=0) %>%              # spreads by the values from group_by
                           dplyr::ungroup() %>%
                           dplyr::mutate(total_counts=base::rowSums(.[c(2:ncol(.))])) %>%               # counts sum from all the columns along the rows
                           dplyr::select(c(split_by, "total_counts", tidyselect::everything())) %>%
                           dplyr::arrange(dplyr::across(tidyselect::all_of(split_by)))                  # sort for consistency
            label_data <- data@meta.data %>%
                          dplyr::group_by(dplyr::across(tidyselect::all_of(split_by))) %>%
                          dplyr::tally() %>%                                                            # calls n()
                          dplyr::ungroup() %>%
                          dplyr::arrange(dplyr::across(tidyselect::all_of(split_by)))                   # sort for consistency

            plot <- counts_data %>%
                    dplyr::select(-c("total_counts")) %>%                                               # removes "total_counts" column
                    reshape2::melt(id.vars=split_by) %>%                                                # creates "variable" and "value" columns
                    ggplot2::ggplot(
                        ggplot2::aes_string(
                            x=split_by,
                            y="value",
                            fill="variable",
                            label="value"
                        )
                    ) +
                    ggplot2::geom_bar(position=bar_position, stat="identity") +
                    ggrepel::geom_label_repel(
                        label_data,
                        mapping=ggplot2::aes_string(x=split_by, y=-Inf, label="n"),
                        color="white",
                        fontface="bold",
                        fill="darkred",
                        segment.colour=NA,
                        direction="y",
                        size=3,
                        show.legend=FALSE
                    ) +
                    ggplot2::scale_fill_manual(values=palette_colors) +
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    get_theme(theme) +
                    ggplot2::ggtitle(plot_title, subtitle=plot_subtitle) +
                    ggplot2::guides(fill=ggplot2::guide_legend(legend_title)) +
                    Seurat::RotatedAxis()

            if (bar_position == "fill"){
                plot <- plot + ggplot2::scale_y_continuous(labels=scales::percent_format(), expand=c(0.05, 0))
            }

            if (!is.null(label) && label){
                if (bar_position == "fill"){
                    plot <- plot +
                            ggplot2::geom_text(
                                position=ggplot2::position_fill(vjust=0.5),
                                angle=0,
                                color="gray95",
                                fontface="bold",
                                size=3
                            )
                } else {
                    plot <- plot +
                            ggplot2::geom_text(
                                position=ggplot2::position_dodge(width=0.9),
                                angle=90,
                                hjust=-0.2,
                                color="gray45",
                                fontface="bold",
                                size=3
                            )
                }
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting composition plot to ", rootname, ".(png/pdf)", sep=""))
            base::print("Composition plot data")
            base::print(base::as.data.frame(counts_data))                                               # counts data used for the plot
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to composition plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

corr_plot <- function(data, reduction, qc_columns, qc_labels, plot_title, rootname, highlight_dims=NULL, combine_guides=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            embeddings <- SeuratObject::Embeddings(data[[reduction]])
            ndims=length(data[[reduction]])
            plots <- list()
            for (i in 1:length(qc_columns)) {
                current_qc_column <- qc_columns[i]
                current_qc_label <- qc_labels[i]
                if ( !(qc_columns[i] %in% base::colnames(data@meta.data)) ){
                    base::print(
                        base::paste(
                            "Column", current_qc_column, "was not found,",
                            "skipping", current_qc_label
                        )
                    )
                    next
                }
                qc_data <- data[[current_qc_column]]
                corr_data <- base::as.data.frame(stats::cor(x=embeddings, y=qc_data))
                corr_data$correlation <- corr_data[, 1]
                corr_data$dimension <- seq_len(length.out=base::nrow(corr_data))
                corr_data$color <- "black"
                if (!is.null(highlight_dims)){
                    corr_data[highlight_dims, "color"] <- "red"
                }
                plots[[current_qc_column]] <- ggplot2::ggplot(corr_data, ggplot2::aes(dimension, correlation)) +
                                              ggplot2::geom_point(color=corr_data$color) +
                                              ggplot2::xlab("Dimension") +
                                              ggplot2::ylab("Correlation") +
                                              ggplot2::xlim(c(0, ndims)) +
                                              ggplot2::ylim(c(-1, 1)) +
                                              get_theme(theme) +
                                              ggplot2::ggtitle(current_qc_label)
            }
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides) + patchwork::plot_annotation(title=plot_title)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting correlation plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export correlation plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

tss_plot <- function(data, rootname, plot_title, split_by, plot_subtitle=NULL, group_by_value=NULL, combine_guides=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            SeuratObject::Idents(data) <- split_by
            identities <- base::unique(base::as.vector(as.character(SeuratObject::Idents(data))))
            plots <- list()
            for (i in 1:length(identities)) {
                current_identity <- identities[i]
                filtered_data <- base::subset(data, idents=current_identity)
                group_by <- NULL
                if (!is.null(group_by_value)){
                    filtered_data$tss_group_by <- base::ifelse(
                        (filtered_data$TSS.enrichment >= group_by_value),
                        base::paste("a) TSS enrichment score >=", group_by_value),
                        base::paste("b) TSS enrichment score <", group_by_value)
                    )
                    group_by <- "tss_group_by"
                }
                plots[[current_identity]] <- Signac::TSSPlot(
                        filtered_data,
                        group.by=group_by
                    ) +
                    ggplot2::ggtitle(current_identity) +
                    get_theme(theme) +
                    Seurat::NoLegend()
            }
            SeuratObject::Idents(data) <- "new.ident"
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides) + patchwork::plot_annotation(title=plot_title, subtitle=plot_subtitle)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting TSS Enrichment plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export TSS Enrichment plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

fragments_hist <- function(data, rootname, plot_title, split_by, plot_subtitle=NULL, group_by_value=NULL, combine_guides=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            SeuratObject::Idents(data) <- split_by
            identities <- base::unique(base::as.vector(as.character(SeuratObject::Idents(data))))
            plots <- list()
            for (i in 1:length(identities)) {
                current_identity <- identities[i]
                filtered_data <- base::subset(data, idents=current_identity)
                group_by <- NULL
                if (!is.null(group_by_value)){
                    filtered_data$ns_group_by <- base::ifelse(
                        (filtered_data$nucleosome_signal <= group_by_value),
                        base::paste("a) Nucl. signal <=", group_by_value),
                        base::paste("b) Nucl. signal >", group_by_value)
                    )
                    group_by <- "ns_group_by"
                }
                plots[[current_identity]] <- Signac::FragmentHistogram(
                        filtered_data,
                        group.by=group_by,
                        region = "chr1-1-100000000"                       # need to set longer region because of https://github.com/stuart-lab/signac/issues/199
                    ) +
                    get_theme(theme) +
                    ggplot2::ggtitle(current_identity) +
                    Seurat::NoLegend()
            }
            SeuratObject::Idents(data) <- "new.ident"
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides) + patchwork::plot_annotation(title=plot_title, subtitle=plot_subtitle)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting ATAC fragments length histogram to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export ATAC fragments length histogram to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

pca_plot <- function(pca_data, pcs, rootname, plot_title, legend_title, plot_subtitle=NULL, color_by="label", label_size=5, pt_size=8, pt_shape=19, alpha=1, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            x_score_column <- base::paste0("PC", pcs[1])
            y_score_column <- base::paste0("PC", pcs[2])
            x_variance <- pca_data$variance[pcs[1]]
            y_variance <- pca_data$variance[pcs[2]]
            plot <- ggplot2::ggplot(
                        pca_data$scores,
                        ggplot2::aes_string(x=x_score_column, y=y_score_column, color=color_by)
                    ) +
                    ggplot2::geom_point(size=pt_size, shape=pt_shape, alpha=alpha) +
                    ggplot2::xlab(base::paste0(x_score_column, ": ", x_variance, "% variance")) +
                    ggplot2::ylab(base::paste0(y_score_column, ": ", y_variance, "% variance")) + 
                    ggrepel::geom_label_repel(
                        ggplot2::aes_string(label=color_by),
                        size=label_size,
                        point.padding=0.5,
                        box.padding=0.5,
                        check_overlap=TRUE,
                        show.legend=FALSE
                    ) +
                    ggplot2::ggtitle(plot_title, subtitle=plot_subtitle) +
                    ggplot2::guides(color=ggplot2::guide_legend(legend_title)) +
                    ggplot2::scale_color_manual(values=palette_colors) +
                    get_theme(theme)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting PCA plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export PCA plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

mds_html_plot <- function(norm_counts_data, rootname){
    tryCatch(
        expr = {
            location <- base::paste0(rootname, ".html")
            htmlwidgets::saveWidget(
                Glimma::glimmaMDS(
                    x=SummarizedExperiment::assay(norm_counts_data),
                    groups=base::as.data.frame(SummarizedExperiment::colData(norm_counts_data)),
                    labels=base::rownames(SummarizedExperiment::colData(norm_counts_data))
                ),
                file=location
            )
            base::print(base::paste("Exporting MDS plot to ", location, sep=""))
        },
        error = function(e){
            print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
        }
    )
}

dot_plot <- function(data, features, rootname, plot_title, x_label, y_label, cluster_idents=FALSE, min_pct=0.01, col_min=-2.5, col_max=2.5, plot_subtitle=NULL, theme="classic", pdf=FALSE, width=1200, height=NULL, resolution=100){
    base::tryCatch(
        expr = {
            plot <- Seurat::DotPlot(
                        data,
                        features=features,
                        cluster.idents=cluster_idents,
                        dot.min=min_pct,
                        col.min=col_min,
                        col.max=col_max,
                        scale=TRUE,
                        scale.by="size"  # for optimal perception
                    ) +
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    get_theme(theme) +
                    ggplot2::ggtitle(plot_title, subtitle=plot_subtitle) +
                    Seurat::RotatedAxis()

            if (is.null(height)){
                height <- round((length(base::unique(base::as.vector(as.character(SeuratObject::Idents(data))))) + 2) * 0.4 * resolution)
                height <- base::ifelse(height < 400, 400, height)
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting dot plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export dot plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}


expression_density_plot <- function(data, features, rootname, reduction, plot_title, joint=FALSE, alpha=NULL, fixed=TRUE, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {

            plots <- Nebulosa::plot_density(
                        data,
                        features=features,
                        reduction=reduction,
                        joint=joint,              # show joint expression density for all features
                        combine=FALSE
                    )
            if (length(features) == 1){
                plots <- list(plots)
                features <- list(features)
            }
            if (joint){
                plots <- list(plots[[length(plots)]])               # get only the joint expression plot
                features <- base::paste(features, collapse="+")
            }
            plots <- base::lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] +
                              ggplot2::ggtitle(features[i]) +
                              get_theme(theme)
                if (!is.null(alpha)) {
                    plots[[i]]$layers[[1]]$aes_params$alpha <- alpha
                }
                if (!is.null(fixed) && fixed){
                    plots[[i]] <- plots[[i]] + ggplot2::coord_fixed(ratio=1)
                }
                return (plots[[i]])
            })

            combined_plots <- patchwork::wrap_plots(plots, guides="keep") + patchwork::plot_annotation(title=plot_title)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting expression density plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export expression density plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}


feature_plot <- function(data, features, labels, rootname, reduction, plot_title, plot_subtitle=NULL, legend_title=NULL, from_meta=FALSE, split_by=NULL, label=FALSE, label_color="black", label_size=4, order=FALSE, color_limits=NULL, color_breaks=NULL, color_scales=NULL, gradient_colors=c("lightgrey", "blue"), min_cutoff=NA, max_cutoff=NA, pt_size=NULL, ncol=NULL, combine_guides=NULL, fixed=TRUE, alpha=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {

            features_corrected <- features
            labels_corrected <- labels
            if (from_meta){
                features_corrected <- c()
                labels_corrected <- c()
                for (i in 1:length(features)){
                    if (features[i] %in% base::colnames(data@meta.data)){
                        features_corrected <- c(features_corrected, features[i])
                        labels_corrected <- c(labels_corrected, labels[i])
                    } else {
                        base::print(
                            base::paste(
                                "Feature", features[i], "was not found,",
                                "skipping", labels[i]
                            )
                        )
                    }
                }
            }
            if (!is.null(split_by)){
                labels_corrected <- rep(      # if splitting, need to repeat labels so we don't display NA
                    labels_corrected,
                    each=length(
                        base::unique(
                            base::as.vector(as.character(data@meta.data[, split_by]))
                        )
                    )
                )
                # avoiding bug of different scales when using split.by https://github.com/satijalab/seurat/issues/5243
                # current fix works only when features_corrected includes only one gene
                if(length(features_corrected) == 1){                 # length of string is always 1
                    feature_expr_data <- SeuratObject::FetchData(
                        object=data,
                        vars=features_corrected,
                        slot="data"                                  # this is the default slot for FeaturePlot
                    )
                    min_feature_value <- min(feature_expr_data)
                    max_feature_value <- max(feature_expr_data)
                }
            }
            # cols are always set here to the default values, otherwise
            # when gradient_colors includes more than 2 colors the plot
            # is rescaled in a way that removes negative values.
            # gradient_colors is used only when user provided both
            # color_scales and color_limits
            plots <- Seurat::FeaturePlot(
                        data,
                        features=features_corrected,
                        pt.size=pt_size,
                        order=order,
                        min.cutoff=min_cutoff,
                        max.cutoff=max_cutoff,
                        reduction=reduction,
                        split.by=split_by,
                        label=label,
                        label.color=label_color,
                        label.size=label_size,
                        combine=FALSE       # to return a list of gglots
                    )
            plots <- base::lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] +
                              ggplot2::ggtitle(labels_corrected[i]) +
                              get_theme(theme)
                if (!is.null(legend_title)){
                    plots[[i]] <- plots[[i]] +
                                  ggplot2::guides(color=ggplot2::guide_colourbar(legend_title))  # the same for all plots, a.k.a "Expression"
                }
                if (!is.null(fixed) && fixed){
                    plots[[i]] <- plots[[i]] + ggplot2::coord_fixed(ratio=1)
                }
                if (!is.null(alpha)) { plots[[i]]$layers[[1]]$aes_params$alpha <- alpha }
                if (!is.null(split_by) && (length(features_corrected) == 1)){                # applying bug fix - redefining gradient limits
                    plots[[i]] <- plots[[i]] +
                                  ggplot2::scale_color_gradientn(
                                      colors=c("lightgrey", "blue"),                         # still using default color values
                                      limits=c(min_feature_value, max_feature_value)
                                  )
                }
                if (!is.null(color_limits) && !is.null(color_scales)){         # overwriting color limits and scales if both are provided
                    plots[[i]] <- plots[[i]] +
                                  ggplot2::scale_colour_gradientn(
                                      colors=gradient_colors,                  # colors can be redefined only when color_limits and color_scales are set
                                      values=scales::rescale(color_scales),
                                      na.value="lightgrey",
                                      limits=color_limits,
                                      breaks=if(!is.null(color_breaks)) color_breaks else ggplot2::waiver()
                                  )
                }
                return (plots[[i]])
            })
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol) +
                              patchwork::plot_annotation(
                                  title=plot_title,
                                  subtitle=plot_subtitle
                              )

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting feature plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export feature plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

dim_heatmap <- function(data, rootname, plot_title, x_label, y_label, reduction="pca", dims=NULL, cells=500, nfeatures=30, ncol=NULL, combine_guides=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plots <- Seurat::DimHeatmap(
                        data,
                        dims=dims,
                        nfeatures=nfeatures,
                        reduction=reduction,
                        cells=cells,
                        balanced=TRUE,
                        fast=FALSE,
                        combine=FALSE
                    )
            plots <- base::lapply(seq_along(plots), function(i){
                plots[[i]] +
                ggplot2::ggtitle(paste("PC", i, sep=" ")) +
                get_theme(theme) +
                ggplot2::xlab(x_label) +
                ggplot2::ylab(y_label) +
                ggplot2::theme(axis.text.x=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
            })
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol) + patchwork::plot_annotation(title=plot_title)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting dimensionality reduction heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export dimensionality reduction heatmap to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

dim_loadings_plot <- function(data, rootname, plot_title, x_label, y_label, reduction="pca", dims=NULL, nfeatures=30, ncol=NULL, combine_guides=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plots <- Seurat::VizDimLoadings(
                        data,
                        dims=dims,
                        nfeatures=nfeatures,
                        reduction=reduction,
                        combine=FALSE
                    )
            plots <- base::lapply(seq_along(plots), function(i){
                plots[[i]] +
                ggplot2::ggtitle(paste("PC", i, sep=" ")) +
                get_theme(theme) +
                ggplot2::xlab(x_label) +
                ggplot2::ylab(y_label)
            })
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides, ncol=ncol) + patchwork::plot_annotation(title=plot_title)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(combined_plots)
            }

            base::print(base::paste("Exporting dimensionality reduction loadings plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export dimensionality reduction loadings plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

coverage_plot <- function(data, assay, region, group_by, plot_title, rootname, plot_subtitle=NULL, idents=NULL, cells=NULL, features=NULL, expression_assay="RNA", expression_slot="data", extend_upstream=0, extend_downstream=0, tile_cells=100, tile_size=100, show_annotation=TRUE, show_peaks=TRUE, show_tile=FALSE, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=NULL, resolution=100){
    base::tryCatch(
        expr = {

            plot <- Signac::CoveragePlot(
                data,
                assay=assay,
                group.by=group_by,
                region=region,
                idents=idents,
                cells=cells,
                features=features,
                expression.assay=expression_assay,
                expression.slot=expression_slot,
                extend.upstream=extend_upstream,
                extend.downstream=extend_downstream,
                annotation=show_annotation,
                peaks=show_peaks,
                links=FALSE,                                       # always FALSE as it requires running aditional function beforehand
                tile=show_tile,
                tile.size=tile_size,
                tile.cells=tile_cells,
                sep=c("-", "-")
            ) + patchwork::plot_annotation(title=plot_title, subtitle=plot_subtitle)

            plot[[1]][[1]] <- plot[[1]][[1]] +                                        # for genome coverage plots
                              get_theme(theme) +
                              ggplot2::scale_fill_manual(values=palette_colors) +
                              Seurat::NoLegend()
            plot[[1]][[2]] <- plot[[1]][[2]] +                                        # for gene expression plots
                              get_theme(theme) +
                              ggplot2::scale_fill_manual(values=palette_colors) +
                              Seurat::NoLegend()

            if (is.null(height)){
                groups = length(base::unique(base::as.vector(as.character(data@meta.data[[group_by]]))))
                height <- round(groups * 0.9 * resolution)           # 0.9 inch per group
                if (show_tile){
                    height <- round(2 * height)                      # to give twice more space for tiles
                    height <- height + round(0.1 * resolution)       # 0.1 inch for space between the plots
                }
                if (show_annotation){
                    height <- height + round(0.5 * resolution)       # 0.5 inch for annotations
                }
                if (show_peaks){
                    height <- height + round(0.3 * resolution)       # 0.3 inch for peaks
                }
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting genome coverage plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export genome coverage plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

volcano_plot <- function(data, rootname, x_axis, y_axis, x_cutoff, y_cutoff, x_label, y_label, plot_title, plot_subtitle, caption, features=NULL, label_column="gene", theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plot <- EnhancedVolcano::EnhancedVolcano(
                        data,
                        x=x_axis,
                        y=y_axis,
                        lab=data[,label_column],
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
                        widthConnectors=0.75
                    ) +
                    ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides="l", alpha=0.3) +
                    get_theme(theme) +
                    ggplot2::theme(legend.position="none", plot.subtitle=ggplot2::element_text(size=8, face="italic", color="gray30"))

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Export volcano plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export volcano plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

feature_heatmap <- function(data, features, rootname, plot_title, assay="RNA", slot="data", cells=NULL, scale_to_max=TRUE, scale="none", color_breaks=NA, highlight_features=NULL, cluster_rows=FALSE, split_rows=NULL, legend_title="Expression", heatmap_colors=c("blue", "black", "yellow"), group_by="new.ident", show_rownames=FALSE, palette_colors=D40_COLORS, pdf=FALSE, width=1200, height=900, resolution=100){
    base::tryCatch(
        expr = {

            # borrowed from the ComplexHeatmap package
            if (scale=="row" && !scale_to_max && is.na(color_breaks)){      # building z-scores but no color breaks provided
                mat <- t(scale(t(base::as.matrix(                           # z-score per row
                           SeuratObject::GetAssayData(
                               object=data,
                               assay=assay,
                               slot=slot                                    # not restricted, but better use counts slot
                           )
                       ))))
                limit <- stats::quantile(abs(mat), 0.99, na.rm=TRUE)        # to exclude outliers
                color_breaks <- base::pretty(c(-limit, limit), n=5)
                base::rm(mat)                                               # remove unused data
                base::gc(verbose=FALSE)
            }

            plot <- dittoSeq::dittoHeatmap(
                data,
                assay=assay,
                slot=slot,
                cells.use=cells,
                genes=features,
                highlight.features=highlight_features,
                scaled.to.max=scale_to_max,
                cluster_rows=cluster_rows,          # will cluster row within each group defined by row_split
                cluster_cols=FALSE,
                show_colnames=FALSE,
                show_rownames=show_rownames,
                main=plot_title,
                heatmap.colors=grDevices::colorRampPalette(heatmap_colors)(50),
                heatmap.colors.max.scaled=grDevices::colorRampPalette(heatmap_colors[1:2])(25),  # only two colors needed
                breaks=color_breaks,
                scale=scale,                        # can be "row"/"column"/"none" but will be forced to "none" if scaled.to.max is TRUE
                annot.by=group_by,
                order.by=group_by,                  # the order of items in group_by will define the levels of ordering
                annot.colors=palette_colors,        # defines colors for the first item set in annot.by
                drop_levels=TRUE,                   # to drop factor levels that are not present in factor values
                use_raster=TRUE,
                silent=TRUE,                        # to prevent saving to file
                complex=TRUE,                       # to use ComplexHeatmap instead of pheatmap
                row_split=split_rows,               # defines rows order and grouping
                cluster_row_slices=FALSE,           # to prevent an additional clustering applied to the mean of row slices
                row_gap = grid::unit(0.1, "mm"),    # instead of the default 1 mm
                name=legend_title                   # to give the heatmap color scale a custom title
            )

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting feature expression heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export feature expression heatmap to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

daseq_permutations <- function(data, rootname, plot_title, x_label, y_label, y_intercepts=NULL, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {

            plot <- data$rand.plot +
                    ggplot2::ggtitle(plot_title) +
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    get_theme(theme)

            if(!is.null(y_intercepts) && length(y_intercepts) > 0){
                intercept_data <- base::data.frame(y_coordinate=y_intercepts) %>%
                                  tibble::add_column(color=palette_colors[1:base::nrow(.)])
                plot <- plot +
                        ggplot2::geom_hline(
                            intercept_data,
                            mapping=ggplot2::aes(yintercept=y_coordinate),
                            color=intercept_data$color
                        ) +
                        ggrepel::geom_label_repel(
                            intercept_data,
                            mapping=ggplot2::aes(x=0, y=y_coordinate, label=y_coordinate),
                            color=intercept_data$color,
                            fill="white",
                            direction="x",
                            size=3,
                            show.legend=FALSE
                        )
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            if (knitr::is_html_output()){
                knitr::opts_chunk$set(
                    fig.width=round(width/resolution),
                    fig.height=round(height/resolution),
                    dpi=resolution
                )
                base::plot(plot)
            }

            base::print(base::paste("Exporting DA permutations plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export DA permutations plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}