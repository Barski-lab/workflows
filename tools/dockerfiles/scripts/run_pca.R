#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(ggrepel))


##########################################################################################
#
# All input CSV/TSV files should have the following header (case-sensitive)
# <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
# <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV
#
# Format of the input files is identified by file extension
# *.csv - CSV
# *.tsv - TSV
# CSV is used by default
#
# v0.0.9
# - Reports loadings, PDF, no 3d plot in PNG anymore
##########################################################################################


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = ","
    if (ext == "tsv"){
        separator = "\t"
    }
    return (separator)
}


load_data_set <- function(filenames, prefixes, target_colname, intersect_by, genelist_data) {
    selected_data <- NULL
    for (i in 1:length(filenames)) {
        raw_data <- read.table(filenames[i], sep=get_file_type(filenames[i]), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
        print(paste("Load ", nrow(raw_data), " rows from ", filenames[i], sep=""))
        colnames(raw_data)[colnames(raw_data) == target_colname] <- prefixes[i]
        if (is.null(selected_data)){
            selected_data <- raw_data
        } else {
            if (is.null(intersect_by)){
                selected_data <- cbind(selected_data, raw_data)
            } else {
                selected_data <- merge(selected_data, raw_data, by=intersect_by, sort = FALSE)
            }
        }
    }
    if (!is.null(genelist_data)){
        print("Apply filter by gene name")
        selected_data <- selected_data[selected_data[,"GeneId"] %in% genelist_data[,1],]
        print(paste("Number of rows after filtering by gene name", nrow(selected_data), sep=" "))
    }
    return (selected_data[, c(intersect_by, prefixes)])
}


export_pca_plot <- function(data, rootname, x, y, x_perc_var, y_perc_var, width=800, height=800, resolution=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                ggplot(data, aes_string(x=x, y=y, color="experiment")) +
                geom_point(size=5, shape=19) +
                xlab(paste0(toString(x), ": ", x_perc_var, "% variance")) +
                ylab(paste0(toString(y), ": ", y_perc_var, "% variance")) + 
                geom_label_repel(aes(label=experiment), point.padding=0.5, box.padding=0.5, check_overlap = TRUE, show.legend = FALSE)
            )
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                ggplot(data, aes_string(x=x, y=y, color="experiment")) +
                geom_point(size=5, shape=19) +
                xlab(paste0(toString(x), ": ", x_perc_var, "% variance")) +
                ylab(paste0(toString(y), ": ", y_perc_var, "% variance")) + 
                geom_label_repel(aes(label=experiment), point.padding=0.5, box.padding=0.5, check_overlap = TRUE, show.legend = FALSE)
            )
            dev.off()

            cat(paste("\nExport PCA-plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export PCA-plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        }
    )
}


export_variance_plot <- function(variance_percentage, rootname, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            
            data = data.frame(variance=variance_percentage)
            data$item = as.numeric(row.names(data))
            
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                ggplot(data, aes(x=item, y=variance)) +
                xlab("") +
                ylab("Variance %") + 
                coord_cartesian(ylim = c(0, 100)) +
                scale_x_continuous(breaks=data$item) +
                geom_line(color="red")+
                geom_point()
            )
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                ggplot(data, aes(x=item, y=variance)) +
                xlab("") +
                ylab("Variance %") + 
                coord_cartesian(ylim = c(0, 100)) +
                scale_x_continuous(breaks=data$item) +
                geom_line(color="red")+
                geom_point()
            )
            dev.off()

            cat(paste("\nExport Variance-plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export Variance-plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        }
    )
}


export_3d_plot_html <- function(data, rootname){
    tryCatch(
        expr = {
            t <- list(
                family = "sans serif",
                size = 14,
                color = toRGB("grey50")
            )
            fig_3d_plot = plot_ly(data, x=~PC1, y=~PC2, z=~PC3, color=~experiment, text=data$experiment)
            fig_3d_plot = fig_3d_plot %>% add_markers()
            fig_3d_plot <- fig_3d_plot %>% layout(
                scene = list(
                    xaxis = list(title = "PC1"),
                    yaxis = list(title = "PC2"),
                    zaxis = list(title = "PC3")
                )
            )
            fig_3d_plot <- fig_3d_plot %>% add_text(textfont = t, textposition = "top right")
            saveWidget(fig_3d_plot, paste(rootname, ".html", sep=""))
            
            cat(paste("\nExport 3D-plot to ", rootname, ".html", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export 3D-plot to ", rootname, ".html", "\n", sep=""))
        }
    )
}


# Parser
parser <- ArgumentParser(description='Run BioWardrobe PCA')
parser$add_argument("-i", "--input",     help='Input CSV/TSV files',                   type="character", required="True", nargs='+')
parser$add_argument("-n", "--name",      help='Input aliases, the order corresponds to --input order. Default: basename of --input files', type="character", nargs='+')
parser$add_argument("-t", "--target",    help='Target column name to be used by PCA',  type="character", choices=c("Rpkm","TotalReads"), default="Rpkm")
parser$add_argument("-c", "--combine",   help='Combine inputs by columns names. Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand', type="character", nargs='+', default=c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand"))
parser$add_argument("-g", "--genelist",  help='Filter genes by the list from the file. Headerless, 1 gene per line', type="character")
parser$add_argument("-o", "--output",    help='Output prefix. Default: ./report_',      type="character", default="./report_")
args <- parser$parse_args(gsub("'", "\"", commandArgs(trailingOnly = TRUE)))

# Set default value for --name if it wasn't provided
if(is.null(args$name)){
    for (i in 1:length(args$input)) {
        args$name = append(args$name, head(unlist(strsplit(basename(args$input[i]), ".", fixed = TRUE)), 1))
    }
}

genelist_data <- NULL
if(!is.null(args$genelist)){
    print(paste("Load gene list from the file", args$genelist, sep=" "))
    genelist_data <- read.table(args$genelist, sep=get_file_type(args$genelist), header=FALSE, stringsAsFactors=FALSE)
    print(paste("PCA will be limited to", nrow(genelist_data), "genes", sep=" "))
}

target_data <- load_data_set(args$input, args$name, args$target, args$combine, genelist_data)
filtered_target_data <- target_data[rowSums(target_data[, args$name]) != 0,]

icolor <- colorRampPalette(c("red", "black", "green", "yellow", "blue", "pink", "brown"))(length(colnames(filtered_target_data[, args$name])))

pca <- prcomp(t(filtered_target_data[, args$name]), cor=TRUE, scale.=T)

loadings <- pca$rotation
loadings_df <- as.data.frame(loadings)
loadings_df <- cbind(filtered_target_data[, !colnames(filtered_target_data) %in% args$name], loadings_df)

write.table(loadings_df,
            file=paste(args$output, "loadings.tsv", sep=""),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

result <- pca$x
result_df <- as.data.frame(result)
result_df <- cbind(experiment=rownames(result_df), result_df)

write.table(result_df,
            file=paste(args$output, "scores.tsv", sep=""),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

variance_percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)

export_pca_plot(result_df, paste(args$output, "pc1_pc2_plot", sep=""), "PC1", "PC2", variance_percentage[1], variance_percentage[2])
export_pca_plot(result_df, paste(args$output, "pc2_pc3_plot", sep=""), "PC2", "PC3", variance_percentage[2], variance_percentage[3])
export_variance_plot(variance_percentage, paste(args$output, "variance_plot", sep=""))
export_3d_plot_html(result_df, paste(args$output, "pca_3d_plot", sep=""))

graphics.off()
