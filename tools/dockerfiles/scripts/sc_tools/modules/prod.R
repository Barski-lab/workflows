import("future", attach=FALSE)
import("data.table", attach=FALSE)
import("BiocParallel", attach=FALSE)

export(
    "parallel"
)


parallel <- function (args) {
    invisible(utils::capture.output(future::plan("multiprocess", workers=args$cpus)))
    invisible(utils::capture.output(future::plan()))
    invisible(utils::capture.output(data.table::setDTthreads(args$cpus)))
    base::options(future.globals.maxSize = args$memory * 1024^3)                      # convert to bytes
    BiocParallel::register(BiocParallel::MulticoreParam(args$cpus, RNGseed=12))       # for DESeq2, RNGseed is hardcoded to make results reproducible
}