import("dplyr", attach=FALSE)
import("purrr", attach=FALSE)
import("tidyr", attach=FALSE)
import("scales", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("tibble", attach=FALSE)
import("Glimma", attach=FALSE)
import("ggplot2", attach=FALSE)
import("ggrepel", attach=FALSE)
import("cluster", attach=FALSE)
import("reshape2", attach=FALSE)
import("dittoSeq", attach=FALSE)
import("Nebulosa", attach=FALSE)
import("patchwork", attach=FALSE)
import("htmlwidgets", attach=FALSE)
import("RColorBrewer", attach=FALSE)
import("magrittr", `%>%`, attach=TRUE)
import("EnhancedVolcano", attach=FALSE)
import("SummarizedExperiment", attach=FALSE)


export(
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
    "D40_COLORS"
)

# https://sashamaps.net/docs/resources/20-colors/
# https://cran.r-project.org/web/packages/Polychrome/vignettes/testgg.html
D40_COLORS <- c("#FB1C0D", "#0DE400", "#0D00FF", "#E8B4BD", "#FD00EA", "#0DD1FE", "#FF9B0D", "#0D601C", "#C50D69", "#CACA16", "#722A91", "#00DEBF", "#863B00", "#5D7C91", "#FD84D8", "#C100FB", "#8499FC", "#FD6658", "#83D87A", "#968549", "#DEB6FB", "#832E60", "#A8CAB0", "#FE8F95", "#FE1CBB", "#DF7CF8", "#FF0078", "#F9B781", "#4D493B", "#1C5198", "#7C32CE", "#EFBC16", "#7CD2DE", "#B30DA7", "#9FC0F6", "#7A940D", "#9B0000", "#946D9B", "#C8C2D9", "#94605A")

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

geom_bar_plot <- function(data, rootname, x_axis, color_by, x_label, y_label, legend_title, plot_title, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plot <- ggplot2::ggplot(data, ggplot2::aes_string(x=x_axis, fill=color_by)) +
                ggplot2::geom_bar(colour="black") +
                ggplot2::geom_text(stat="count", ggplot2::aes(label=..count..), vjust=-1) +
                ggplot2::xlab(x_label) +
                ggplot2::ylab(y_label) +
                ggplot2::guides(fill=ggplot2::guide_legend(legend_title), x=ggplot2::guide_axis(angle=45)) +
                ggplot2::ggtitle(plot_title) +
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

            base::print(base::paste("Exporting geom bar plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export geom bar plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

geom_density_plot <- function(data, rootname, x_axis, color_by, facet_by, x_left_intercept, x_label, y_label, legend_title, plot_title, x_right_intercept=NULL,  scale_x_log10=FALSE, scale_y_log10=FALSE, zoom_on_intercept=FALSE, alpha=0.7, show_ranked=FALSE, ranked_x_label="Ranked cells", palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            intercept_data <- data %>%
                              dplyr::select(tidyselect::all_of(color_by), tidyselect::all_of(facet_by)) %>%
                              dplyr::distinct() %>%
                              dplyr::arrange(dplyr::across(tidyselect::all_of(color_by))) %>%
                              tibble::add_column(color=palette_colors[1:base::nrow(.)],x_left=x_left_intercept)

            plot <- ggplot2::ggplot(data, ggplot2::aes_string(x=x_axis, fill=color_by)) +
                    ggplot2::geom_density(alpha=alpha) +
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    ggplot2::guides(fill=ggplot2::guide_legend(legend_title)) +
                    ggplot2::ggtitle(plot_title) +
                    ggplot2::facet_wrap(stats::as.formula(base::paste("~", facet_by))) +
                    ggplot2::scale_fill_manual(values=palette_colors) +
                    ggplot2::geom_vline(intercept_data, mapping=ggplot2::aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                    ggrepel::geom_label_repel(
                        intercept_data, mapping=ggplot2::aes(x=x_left, y=Inf, label=x_left),
                        color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                        show.legend=FALSE
                    ) +
                    get_theme(theme)

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

            if (scale_x_log10){ plot <- plot + ggplot2::scale_x_log10() }
            if (scale_y_log10){ plot <- plot + ggplot2::scale_y_log10() }

            if (zoom_on_intercept) {
                zoomed_plot <- ggplot2::ggplot(data, ggplot2::aes_string(x=x_axis, color=color_by)) +
                               ggplot2::geom_density(show.legend=FALSE, size=2) +
                               ggplot2::xlab(x_label) +
                               ggplot2::ylab(y_label) +
                               ggplot2::scale_color_manual(values=palette_colors) +
                               ggplot2::geom_vline(intercept_data, mapping=ggplot2::aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                               ggrepel::geom_label_repel(
                                   intercept_data, mapping=ggplot2::aes(x=x_left, y=Inf, label=x_left),
                                   color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                   show.legend=FALSE
                               ) +
                               get_theme(theme)
                if (!is.null(x_right_intercept)){
                    zoomed_plot <- zoomed_plot +
                                   ggplot2::coord_cartesian(xlim=c(min(intercept_data$x_left), max(intercept_data$x_right))) +
                                   ggplot2::geom_vline(intercept_data, mapping=ggplot2::aes(xintercept=x_right), color=intercept_data$color, alpha=0.7) +
                                   ggrepel::geom_label_repel(
                                       intercept_data, mapping=ggplot2::aes(x=x_right, y=Inf, label=x_right),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                       show.legend=FALSE
                                   )
                } else {
                    zoomed_plot <- zoomed_plot + ggplot2::coord_cartesian(xlim=c(NA, max(intercept_data$x_left)))
                }
                if (scale_x_log10){ zoomed_plot <- zoomed_plot + ggplot2::scale_x_log10() }
                if (scale_y_log10){ zoomed_plot <- zoomed_plot + ggplot2::scale_y_log10() }
                plot <- plot / zoomed_plot
            }

            if (show_ranked) {
                ranked_plot <- data %>%
                               dplyr::arrange(dplyr::across(tidyselect::all_of(color_by)), dplyr::across(tidyselect::all_of(x_axis))) %>%
                               ggplot2::ggplot(ggplot2::aes_string(x=seq_along(data[[x_axis]]), y=x_axis, color=color_by)) +
                               ggplot2::geom_point(show.legend=FALSE, size=0.5) +
                               ggplot2::xlab(ranked_x_label) +
                               ggplot2::ylab(x_label) +
                               ggplot2::scale_y_log10() +
                               ggplot2::scale_color_manual(values=palette_colors) +
                               ggplot2::geom_hline(intercept_data, mapping=ggplot2::aes(yintercept=x_left), color=intercept_data$color, alpha=0.7) +
                               ggrepel::geom_label_repel(
                                   intercept_data, mapping=ggplot2::aes(x=Inf, y=x_left, label=x_left),
                                   color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="x", size=3,
                                   show.legend=FALSE
                               ) +
                               get_theme(theme)
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

            base::print(base::paste("Exporting geom density plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export geom density plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

geom_point_plot <- function(data, rootname, x_axis, y_axis, facet_by, x_left_intercept, y_low_intercept, color_by, gradient_colors, color_limits, color_break, x_label, y_label, legend_title, plot_title, y_high_intercept=NULL, scale_x_log10=FALSE, scale_y_log10=FALSE, show_lm=FALSE, show_density=FALSE, alpha=0.2, alpha_intercept=0.5, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            intercept_data <- data %>%
                              dplyr::select(tidyselect::all_of(facet_by)) %>%
                              dplyr::distinct() %>%
                              dplyr::arrange(dplyr::across(tidyselect::all_of(facet_by))) %>%
                              tibble::add_column(color=palette_colors[1:base::nrow(.)], x_left=x_left_intercept, y_low=y_low_intercept)

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
                    ggplot2::ggtitle(plot_title) +
                    ggplot2::facet_wrap(stats::as.formula(base::paste("~", facet_by))) +
                    ggplot2::geom_vline(intercept_data, mapping=ggplot2::aes(xintercept=x_left), color=intercept_data$color, alpha=alpha_intercept) +
                    ggplot2::geom_hline(intercept_data, mapping=ggplot2::aes(yintercept=y_low), color=intercept_data$color, alpha=alpha_intercept) +
                    ggrepel::geom_label_repel(
                        intercept_data, mapping=ggplot2::aes(x=x_left, y=Inf, label=x_left),
                        color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="y", size=3,
                        show.legend=FALSE
                    ) +
                    ggrepel::geom_label_repel(
                        intercept_data, mapping=ggplot2::aes(x=Inf, y=y_low, label=y_low),
                        color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="x", size=3,
                        show.legend=FALSE
                    ) +
                    get_theme(theme)

            if (show_lm){ plot <- plot + ggplot2::stat_smooth(method=stats::lm) }
            if (show_density){ plot <- plot + ggplot2::geom_density_2d() }

            if (!is.null(y_high_intercept)){
                plot <- plot +
                        ggplot2::geom_hline(intercept_data, mapping=ggplot2::aes(yintercept=y_high), color=intercept_data$color, alpha=alpha_intercept) +
                        ggrepel::geom_label_repel(
                            intercept_data, mapping=ggplot2::aes(x=Inf, y=y_high, label=y_high),
                            color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="x", size=3,
                            show.legend=FALSE
                        )
            }

            if (scale_x_log10){ plot <- plot + ggplot2::scale_x_log10() }
            if (scale_y_log10){ plot <- plot + ggplot2::scale_y_log10() }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
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

            base::print(base::paste("Exporting feature scatter plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export feature scatter plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

vln_plot <- function(data, features, labels, rootname, plot_title, legend_title, from_meta=FALSE, log=FALSE, group_by=NULL, split_by=NULL, hide_x_text=FALSE, pt_size=NULL, palette_colors=NULL, combine_guides=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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

            plots <- Seurat::VlnPlot(
                         data,
                         features=features_corrected,
                         pt.size=pt_size,
                         group.by=group_by,
                         split.by=split_by,
                         log=log,
                         combine=FALSE       # to return a list of gglots
                     )
            plots <- base::lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] +
                              ggplot2::ggtitle(labels_corrected[i]) +
                              get_theme(theme) +
                              ggplot2::theme(axis.title.x=ggplot2::element_blank()) +
                              ggplot2::guides(fill=ggplot2::guide_legend(legend_title)) +
                              ggplot2::stat_boxplot(width=0.15, geom="errorbar") +
                              ggplot2::geom_boxplot(width=0.15, outlier.alpha=0) +
                              Seurat::RotatedAxis()
                if (!is.null(palette_colors)){ plots[[i]] <- plots[[i]] + ggplot2::scale_fill_manual(values=palette_colors) }
                if (hide_x_text){ plots[[i]] <- plots[[i]] + ggplot2::theme(axis.text.x=ggplot2::element_blank()) }
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

            base::print(base::paste("Exporting violin plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export violin plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

dim_plot <- function(data, rootname, reduction, plot_title, legend_title, cells=NULL, split_by=NULL, group_by=NULL, highlight_group=NULL, perc_split_by=NULL, perc_group_by=NULL, ncol=NULL, label=FALSE, label_box=FALSE, label_color="black", label_size=4, alpha=NULL, palette_colors=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
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
                    ggplot2::ggtitle(plot_title) +
                    ggplot2::guides(color=ggplot2::guide_legend(legend_title, override.aes=list(size=3)))

            if (!is.null(palette_colors)){
                plot <- plot +
                        ggplot2::scale_color_manual(values=palette_colors) +
                        ggplot2::scale_fill_manual(values=palette_colors)     # need it for proper label box background
            }
            if (!is.null(alpha)) { plot$layers[[1]]$aes_params$alpha <- alpha }

            if(!is.null(perc_split_by) && !is.null(perc_group_by)){
                width <- 2 * width
                perc_data <- data@meta.data %>%
                             dplyr::group_by(dplyr::across(tidyselect::all_of(perc_split_by)), dplyr::across(tidyselect::all_of(perc_group_by))) %>%
                             dplyr::summarise(counts=dplyr::n(), .groups="drop_last") %>%        # drop the perc_group_by grouping level so we can get only groups defined by perc_split_by
                             dplyr::mutate(freq=counts/sum(counts)*100) %>%                      # sum is taken for the group defined by perc_split_by
                             dplyr::arrange(dplyr::across(tidyselect::all_of(perc_split_by)), dplyr::across(tidyselect::all_of(perc_group_by)))        # sort for consistency
                label_data <- data@meta.data %>%
                             dplyr::group_by(dplyr::across(tidyselect::all_of(perc_split_by))) %>%
                             dplyr::summarise(counts=dplyr::n(), .groups="drop_last") %>%                      # drops all grouping as we have only one level
                             dplyr::arrange(dplyr::across(tidyselect::all_of(perc_split_by)))                  # sort for consistency
                perc_plot <- ggplot2::ggplot(perc_data, ggplot2::aes_string(x=perc_split_by, y="freq", fill=perc_group_by)) +
                             ggplot2::geom_col(position="dodge", width=0.9, linetype="solid", color="black", show.legend=FALSE) +
                             ggplot2::xlab("") +
                             ggplot2::ylab("Cells percentage") +
                             get_theme(theme) +
                             ggrepel::geom_label_repel(
                                label_data, mapping=ggplot2::aes(y=-Inf, label=counts),
                                color="black", fill="white", segment.colour=NA,
                                direction="y", size=3, show.legend=FALSE
                             ) +
                             Seurat::RotatedAxis()
                if (!is.null(palette_colors)){ perc_plot <- perc_plot + ggplot2::scale_fill_manual(values=palette_colors) }
                plot <- plot + perc_plot
            }

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
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

elbow_plot <- function(data, rootname, plot_title, reduction="pca", ndims=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {
            plot <- Seurat::ElbowPlot(
                        data,
                        ndims=ndims,
                        reduction=reduction
                    ) +
                    get_theme(theme) +
                    ggplot2::ggtitle(plot_title)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            base::print(base::paste("Exporting elbow plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export elbow plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

silhouette_plot <- function(data, rootname, plot_title, legend_title, group_by, dims, downsample=300, reduction="pca", palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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
                    ggplot2::ggtitle(plot_title) +
                    ggplot2::guides(fill=ggplot2::guide_legend(legend_title))

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            base::print(base::paste("Exporting silhouette plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export silhouette plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

composition_plot <- function(data, rootname, plot_title, legend_title, x_label, y_label, split_by, group_by, bar_position="fill", palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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
                    ggplot2::ggplot(ggplot2::aes_string(x=split_by, y="value")) +
                    ggplot2::geom_bar(ggplot2::aes(fill=variable), position=bar_position, stat="identity") +
                    ggrepel::geom_label_repel(
                        label_data, mapping=ggplot2::aes_string(x=split_by, y="-Inf", label="n"),
                        color="black", fill="white", segment.colour=NA,
                        direction="y", size=3, show.legend=FALSE
                    ) +
                    ggplot2::scale_fill_manual(values=palette_colors) +
                    ggplot2::xlab(x_label) +
                    ggplot2::ylab(y_label) +
                    get_theme(theme) +
                    ggplot2::ggtitle(plot_title) +
                    ggplot2::guides(fill=ggplot2::guide_legend(legend_title)) +
                    Seurat::RotatedAxis()
            if (bar_position == "fill"){
                plot <- plot + ggplot2::scale_y_continuous(labels=scales::percent_format(), expand=c(0.01, 0))
            }
            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            base::print(base::paste("Exporting composition plot to ", rootname, ".(png/pdf)", sep=""))
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

            base::print(base::paste("Exporting correlation plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export correlation plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

tss_plot <- function(data, rootname, plot_title, split_by, group_by_value=NULL, combine_guides=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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
                        filtered_data$TSS.enrichment >= group_by_value ,
                        base::paste("High ", "(bigger or equal to ", group_by_value, ")", sep=""),
                        base::paste("Low ", "(smaller than ", group_by_value, ")", sep="")
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
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides) + patchwork::plot_annotation(title=plot_title)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            base::print(base::paste("Exporting TSS Enrichment plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export TSS Enrichment plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

fragments_hist <- function(data, rootname, plot_title, split_by, group_by_value=NULL, combine_guides=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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
                        filtered_data$nucleosome_signal >= group_by_value ,
                        base::paste("Nucl. signal >= ", group_by_value, sep=""),
                        base::paste("Nucl. signal < ", group_by_value, sep="")
                    )
                    group_by <- "ns_group_by"
                }
                plots[[current_identity]] <- Signac::FragmentHistogram(
                        filtered_data,
                        group.by=group_by
                    ) +
                    get_theme(theme) +
                    ggplot2::ggtitle(current_identity) +
                    Seurat::NoLegend()
            }
            SeuratObject::Idents(data) <- "new.ident"
            combined_plots <- patchwork::wrap_plots(plots, guides=combine_guides) + patchwork::plot_annotation(title=plot_title)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            base::print(base::paste("Exporting fragments length histogram to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export fragments length histogram to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

pca_plot <- function(pca_data, pcs, rootname, plot_title, legend_title, color_by="label", label_size=5, pt_size=8, pt_shape=19, alpha=0.75, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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
                    ggplot2::ggtitle(plot_title) +
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

dot_plot <- function(data, features, rootname, plot_title, x_label, y_label, cluster_idents=FALSE, min_pct=0.01, col_min=-2.5, col_max=2.5, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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
                    ggplot2::ggtitle(plot_title) +
                    Seurat::RotatedAxis()

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            base::print(base::paste("Exporting dot plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export dot plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}


expression_density_plot <- function(data, features, rootname, reduction, plot_title, joint=FALSE, alpha=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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
                plots[[i]] <- plots[[i]] + ggplot2::ggtitle(features[i]) + get_theme(theme)
                if (!is.null(alpha)) { plots[[i]]$layers[[1]]$aes_params$alpha <- alpha }
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

            base::print(base::paste("Exporting expression density plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export expression density plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}


feature_plot <- function(data, features, labels, rootname, reduction, plot_title, from_meta=FALSE, split_by=NULL, label=FALSE, order=FALSE, color_limits=NULL, color_scales=NULL, gradient_colors=c("lightgrey", "blue"), min_cutoff=NA, max_cutoff=NA, pt_size=NULL, combine_guides=NULL, alpha=NULL, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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
                        combine=FALSE       # to return a list of gglots
                    )
            plots <- base::lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] + ggplot2::ggtitle(labels_corrected[i]) + get_theme(theme)
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
                                      limits=color_limits
                                  )
                }
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

            base::print(base::paste("Exporting dimensionality reduction loadings plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export dimensionality reduction loadings plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

coverage_plot <- function(data, assay, region, group_by, plot_title, rootname, idents=NULL, cells=NULL, features=NULL, expression_assay="RNA", expression_slot="data", extend_upstream=0, extend_downstream=0, show_annotation=TRUE, show_peaks=TRUE, palette_colors=D40_COLORS, theme="classic", pdf=FALSE, width=1200, height=800, resolution=100){
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
                tile=FALSE,                                        # always FALSE due to bug in Signac
                sep=c("-", "-")
            ) + patchwork::plot_annotation(title=plot_title)

            plot[[1]][[1]] <- plot[[1]][[1]] + get_theme(theme) + ggplot2::scale_fill_manual(values=palette_colors)    # for genome coverage plots
            plot[[1]][[2]] <- plot[[1]][[2]] + get_theme(theme) + ggplot2::scale_fill_manual(values=palette_colors)    # for gene expression plots

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
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
                    ggplot2::scale_y_log10() +
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

            base::print(base::paste("Export volcano plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export volcano plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

feature_heatmap <- function(data, features, rootname, plot_title, assay="RNA", slot="data", cells=NULL, scale_to_max=TRUE, scale="none", heatmap_colors=c("blue", "black", "yellow"), group_by="new.ident", show_rownames=FALSE, palette_colors=D40_COLORS, pdf=FALSE, width=1200, height=800, resolution=100){
    base::tryCatch(
        expr = {

            plot <- dittoSeq::dittoHeatmap(
                data,
                assay=assay,
                slot=slot,
                cells.use=cells,
                genes=features,
                scaled.to.max=scale_to_max,
                cluster_rows=FALSE,
                cluster_cols=FALSE,
                show_colnames=FALSE,
                show_rownames=show_rownames,
                main=plot_title,
                heatmap.colors=grDevices::colorRampPalette(heatmap_colors)(50),
                heatmap.colors.max.scaled=grDevices::colorRampPalette(heatmap_colors[1:2])(25),  # only two colors needed
                scale=scale,                  # can be "row"/"column"/"none" but will be forced to "none" if scaled.to.max is TRUE
                annot.by=group_by,            # the first item will be used to order the cells
                annot.colors=palette_colors,  # defines colors for the first item set in annot.by
                drop_levels=TRUE,             # to drop factor levels that are not present in factor values
                silent=TRUE                   # to prevent saving to file
            )

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
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

            base::print(base::paste("Exporting DA permutations plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(base::paste("Failed to export DA permutations plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}