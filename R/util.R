#' @importFrom methods setGeneric
setGeneric(".mutateEach", signature=c("grl"),
           function(grl, ...) standardGeneric(".mutateEach"))

#' Efficient Metadata Columns Mutation
#'
#' @param grl a CompressedGRangesList
#' @param ... named list of vectors to insert as metadata columns on each element
#' \code{GRanges}. Each vector length must match the length of the \code{GRangesList}.
#'
#' @return a \code{CompressedGRangesList} with all element \code{GRanges} updated
#' with supplied metadata columns
#'
#' @importFrom methods setMethod slot slot<-
#' @importFrom utils capture.output
#' @importFrom S4Vectors elementNROWS
setMethod(".mutateEach", "CompressedGRangesList",
          function (grl, ...) {
              ## ensure data is a valid length
              inputLengths <- vapply(list(...), length, integer(1), USE.NAMES=TRUE)
              if (!all(inputLengths == length(grl))) {
                stop("Mismatched lengths detected:\n",
                     "\tExpected length ", length(grl),
                     ", but found length(s) ",
                     capture.output(dput(inputLengths)))
              }

              ## expand list to full length
              expandedListData <- lapply(list(...), rep, times=elementNROWS(grl))
              mcols(slot(grl, "unlistData"))[names(list(...))] <- expandedListData
              grl
          }
)

#' @importFrom methods setMethod
setMethod(".mutateEach", "SimpleGRangesList",
          function (grl, ...) .mutateEach(GRangesList(grl, compress=TRUE), ...)
)

#' Compatibility wrapper for `makeTxDbFromGRanges`
#'
#' This internal utility provides a version-safe interface to the
#' `makeTxDbFromGRanges()` function, which has been moved from the
#' \pkg{GenomicFeatures} package to the \pkg{txdbmaker} package in
#' Bioconductor (>= 3.19, GenomicFeatures >= 1.61.1).
#'
#' The function checks at runtime which package provides
#' `makeTxDbFromGRanges()` and dispatches the call accordingly.
#'
#' @param ... Arguments passed directly to
#'   \code{makeTxDbFromGRanges()} from either \pkg{txdbmaker} or
#'   \pkg{GenomicFeatures}.
#'
#' @returns A \link[GenomicFeatures:TxDb-class]{TxDb} object constructed
#'   from the provided \code{GRanges} object(s), consistent with
#'   the behavior of \code{makeTxDbFromGRanges()}.
#'
#' @details
#' If \pkg{txdbmaker} is installed, the function calls
#' \code{txdbmaker::makeTxDbFromGRanges()}.  
#' If not, and if \pkg{GenomicFeatures} provides the function (in versions
#' \eqn{\leq 1.61.1}), it will call
#' \code{GenomicFeatures::makeTxDbFromGRanges()} instead.  
#' If neither package is available, an error is raised.
#'
#' @export
makeTxDbSafe <- function(...) {
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    getFromNamespace("makeTxDbFromGRanges", "txdbmaker")(...)
  } else if (requireNamespace("GenomicFeatures", quietly = TRUE)) {
    getFromNamespace("makeTxDbFromGRanges", "GenomicFeatures")(...)
  } else {
    stop("Neither 'txdbmaker' nor 'GenomicFeatures' packages are available.")
  }
}