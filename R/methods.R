#' @importFrom methods getMethod
setMethod(
    f="show",
    signature="SummarizedExperiment",
    definition=function(object) {
        if (isTRUE(x=getOption(x="restore_SummarizedExperiment_show",
            default = FALSE)) |         
            # If the object is a SingleCellExperiment
            # # From BioC 3_14 SingleCellExperiment is SummarizedExperiment and 
            # # we don't want to process with tidySummarizedExperiment
            is(object, "SingleCellExperiment")
        ) {
            f <- getMethod(
                f="show",
                signature="SummarizedExperiment",
                where=asNamespace(ns="SummarizedExperiment")
            )
            f(object=object)
        } else {
            object %>%
                print()
        }
    }
)

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