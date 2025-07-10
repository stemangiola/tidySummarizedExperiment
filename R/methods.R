#' @importFrom methods getMethod
setClass("tidySummarizedExperiment",
    contains=c("SummarizedExperiment", "RangedSummarizedExperiment"))

#' @name tidy
#' @rdname tidy
#' @title tidy for `Seurat`
#'
#' @param object A `Seurat` object.
#' @return A `tidyseurat` object.
#'
#' @examples
#' data(pasilla)
#' pasilla %>% tidy()
#'
#' @export
tidy <- function(object) {
    UseMethod("tidy", object)
}

#' @importFrom lifecycle deprecate_warn
tidy_ <- function(object) {
    
    # DEPRECATE
    deprecate_warn(
        when = "1.1.1",
        what = "tidy()",
        details = "tidySummarizedExperiment says: tidy() is not needed anymore."
    )
    
    object
}

#' @importFrom methods as
#' @rdname tidy
#' @param object A SummarizedExperiment object
#' @export
tidy.SummarizedExperiment <- tidy_

#' @importFrom methods as
#' @rdname tidy
#' @param object A SummarizedExperiment object
#' @export
tidy.RangedSummarizedExperiment <- tidy_