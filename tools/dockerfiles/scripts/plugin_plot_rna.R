#!/usr/bin/env Rscript
options(warn=-1)
options("width"=300)

suppressMessages(library(argparse))
suppressMessages(library(sqldf))
suppressMessages(library(ggplot2))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(BiocParallel))

######################################################################
#
# v0.0.4
# - Support --stranded yes/no/reverse
#   reverse corresponds to the original dUTP mode
#
# v0.0.3
# - Refactored
#
######################################################################


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = ","
    if (ext == "tsv"){
        separator = "\t"
    }
    return (separator)
}


get_coverage <- function(ranges, isoforms, min_rpkm, min_length, bam_file, is_pair, stranded, threads) {
	plt_len <- 99
	is_ranges<-sqldf(
        paste("SELECT ranges.chrom, ranges.strand, ranges.exonCount, ranges.exonStarts, ranges.exonEnds
			   FROM   ranges, isoforms
			   WHERE  ranges.name=isoforms.RefseqId AND
                      ranges.name2=isoforms.GeneId AND
                      ranges.chrom=isoforms.Chrom AND
                      ranges.txStart=isoforms.TxStart AND
					  ranges.txEnd=isoforms.TxEnd AND
					  ranges.strand=isoforms.Strand AND
					  ranges.txEnd-ranges.txStart > ", min_length, " AND
					  isoforms.Rpkm > ", min_rpkm, sep="")
    )
	is_count <- length(is_ranges$exonCount)
	if(is_count < 5) return(rep(0, plt_len + 1))
	namesl <- list()
	irl <- list()
	for(i in seq_len(is_count)){
		st <- as.numeric(strsplit(is_ranges$exonStarts[i], ',')[[1]])
		en <- as.numeric(strsplit(is_ranges$exonEnds[i], ',')[[1]])
		names <- paste0(rep(is_ranges$chrom[i], is_ranges$exonCount[i]), ":", st, "-", en)
		ir <- IRanges(start=st, end=en, names=names)
		metadata(ir) <- list(strand=is_ranges$strand[i])
		namesl <- append(namesl, is_ranges$chrom[i])
		irl <- append(irl, list(ir))
	}
	what <- c("pos","strand")
    which <- IRangesList(irl)
	names(which) <- namesl
	flags <- scanBamFlag(isUnmappedQuery = FALSE,                      # only mapped reads
                         isPaired = if(is_pair) TRUE else NA,          # if it was paired-end data, then we need only reads that have pair, otherwise NA - means any
                         isProperPair = if(is_pair) TRUE else NA,      # if it was paired-end data, then we need only properly paired reads, otherwise NA - means any
                         hasUnmappedMate = if(is_pair) FALSE else NA,  # if it was paired-end data, then the mate should be mapped too, otherwise NA - means any
                         isFirstMateRead = if(is_pair) TRUE else NA    # if it was paired-end data, then we take only the first read in pair, otherwise NA - means any
    )
	param <- ScanBamParam(which = which, what = what, flag = flags)
	bam <- scanBam(bam_file, index = bam_file, param = param)
	cov <- bplapply(irl,
                    function(ir){
                        genewidth <- sum(width(ir))
                        bin <- genewidth / plt_len
                        cov <- c(rep(0, plt_len + 1))
                        for(i in seq_len(length(ir))){
                            abs_pos <- 0
                            idx <- 1
                            cwidth <- 0
                            # Strand specificity
                            if (stranded == "no") {
                                tmp1 <- rle(bam[[names(ir)[i]]]$pos)
                            } else if (stranded == "yes" && metadata(ir)$strand == "+") {
                                tmp1 <- rle(bam[[names(ir)[i]]]$pos[bam[[names(ir)[i]]]$strand == "+"])
                            } else if (stranded == "yes" && metadata(ir)$strand == "-") {
                                tmp1 <- rle(bam[[names(ir)[i]]]$pos[bam[[names(ir)[i]]]$strand == "-"])
                            } else if (stranded == "reverse" && metadata(ir)$strand == "+") {
                                tmp1 <- rle(bam[[names(ir)[i]]]$pos[bam[[names(ir)[i]]]$strand == "-"])
                            } else if (stranded == "reverse" && metadata(ir)$strand == "-") {
                                tmp1 <- rle(bam[[names(ir)[i]]]$pos[bam[[names(ir)[i]]]$strand == "+"])
                            }
                            for(j in seq_len(length(tmp1$lengths))){
                                pos <- tmp1$values[j]
                                if(pos > end(ir)[idx]){
                                    cwidth <- sum(width(ir[end(ir) < pos]))
                                    idx <- sum(start(ir) < pos)
                                    if(idx > length(ir))
                                        break
                                }
                                if(pos - start(ir)[idx] + 1 > width(ir)[idx])
                                    next
                                abs_pos <- (pos - start(ir)[idx]) + cwidth
                                il <- floor(abs_pos / bin)
                                cov[il+1] = cov[il+1] + tmp1$lengths[j]
                            }
                        }
                        cov <- cov / bin
                        if(metadata(ir)$strand == "+") { cov } else { rev(cov) }
                    },
                    BPPARAM = MulticoreParam(workers = threads))
    return(Reduce("+", cov) / is_count)
}


export_tag_density_plot <- function(data, min_rpkm, min_length, rootname, width=800, height=800, resolution=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                ggplot(data, aes(x, y)) +
                ggtitle(paste("Gene body average tag density for RPKM < ", min_rpkm, " and transcript length < ", min_length, sep="")) +
                geom_line() +
                xlab("Gene body percentile (5' -> 3')") +
                ylab("Average Tag Density (per percentile)") +
                scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20))
            )
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                ggplot(data, aes(x, y)) +
                ggtitle("Gene body average tag density") +
                geom_line() +
                xlab("Gene body percentile (5' -> 3')") +
                ylab("Average Tag Density (per percentile)") +
                scale_x_continuous(breaks=seq(0,100,20), labels=seq(0,100,20))
            )
            dev.off()

            cat(paste("\nExport tag density plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export tag density plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        }
    )
}


export_histogram_plot <- function(data, rootname, bins=1000, width=800, height=800, resolution=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                ggplot(data, aes(y)) +
                ggtitle("RPKM distribution") +
                geom_histogram(bins=bins) +
                xlab("RPKM") +
                ylab("Frequency")
            )
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                ggplot(data, aes(y)) +
                ggtitle("Isoforms RPKM distribution") +
                geom_histogram(bins=bins) +
                xlab("RPKM") +
                ylab("Frequency")
            )
            dev.off()

            cat(paste("\nExport histogram to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export histogram to ", rootname, ".(png/pdf)", "\n", sep=""))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description='Gene body average tag density plot and RPKM distribution histogram for isoforms')
    parser$add_argument("--annotation", help="Path to the annotation TSV/CSV file with gene names set in name2 field", type="character", required="True")
    parser$add_argument("--bam",        help="Path to the indexed BAM file", type="character", required="True")
    parser$add_argument("--isoforms",   help="Path to the RPKM isoforms expression TSV/CSV file", type="character", required="True")
    parser$add_argument("--minrpkm",    help="Ignore isoforms with RPKM smaller than --minrpkm. Default: 10", type="double", default=10)
    parser$add_argument("--minlength",  help="Ignore isoforms shorter than --minlength. Default: 1000", type="integer", default=1000)
    parser$add_argument("--mapped",     help="Mapped reads/pairs number", type="integer", required="True")
    parser$add_argument("--pair",       help="Run as paired end. Default: false", action='store_true')
    parser$add_argument("--stranded",   help='Strand specificity. One of "yes", "no" or "reverse". Default: "no"', type="character", choices=c("yes","no","reverse"), default="no")
    #  Whether the data is from a strand-specific assay.
    #  --stranded no      - a read is considered overlapping with a feature regardless of whether
    #                       it is mapped to the same or the opposite strand as the feature.
    #  --stranded yes     - the read has to be mapped to the same strand as the feature.
    #  --stranded reverse - the read has to be mapped to the opposite strand than the feature.
    parser$add_argument("--output",     help="Output prefix. Default: ./coverage", type="character", default="./coverage")
    parser$add_argument("--threads",    help="Threads. Default: 1", type="integer", default=1)
    return (parser$parse_args(commandArgs(trailingOnly = TRUE)))
}


args <- get_args()
if (args$mapped <= 0) {  # safety measure
    args$mapped <- 1
} 

ranges <- read.table(args$annotation, sep=get_file_type(args$annotation), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
isoforms <- read.table(args$isoforms, sep=get_file_type(args$isoforms), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cov_raw <- get_coverage(ranges, isoforms, args$minrpkm, args$minlength, args$bam, args$pair, args$stranded, args$threads)/(args$mapped/1000000)
coverage_df <- data.frame(x=as.numeric(rownames(as.data.frame(cov_raw)))-1, y=as.numeric(cov_raw))

report_filename <- paste(args$output, "_gene_body_report.tsv", sep="")
write.table(coverage_df,
            file = report_filename,
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)
print(paste("Export report to ", report_filename, sep=""))

export_tag_density_plot(coverage_df, args$minrpkm, args$minlength, paste(args$output, "_gene_body_plot", sep=""))

rpkms <- isoforms$Rpkm[isoforms$Rpkm>0 & isoforms$Rpkm<500]
rpkms_df <- data.frame(x=as.numeric(rownames(as.data.frame(rpkms)))-1, y=as.numeric(rpkms))
export_histogram_plot(rpkms_df, paste(args$output, "_rpkm_distribution_plot", sep=""))

graphics.off()