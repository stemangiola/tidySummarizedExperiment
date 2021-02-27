setMethod(
    f = "show",
    signature = "SummarizedExperiment",
    definition = function(object) {
        if (isTRUE(x = getOption(x = "restore_SummarizedExperiment_show", default = FALSE))) {
            f <- getMethod(
                f = "show",
                signature = "SummarizedExperiment",
                where = asNamespace(ns = "SummarizedExperiment")
            )
            f(object = object)
        } else {
            object %>%
                print()
        }
    }
)


setClass("tidySummarizedExperiment", contains=c("SummarizedExperiment", "RangedSummarizedExperiment"))


#' tidy for SummarizedExperiment
#'
#' @param object A SummarizedExperiment object
#'
#' @return A tidySummarizedExperiment object
#'
#' @name tidy
#'
#' @examples
#'
#' tidySummarizedExperiment::pasilla %>% tidy()
#' @export
tidy <- function(object) {
    UseMethod("tidy", object)
}

#' @importFrom lifecycle deprecate_warn
tidy_ <- function(object) {
    
    # DEPRECATE
    deprecate_warn(
        when = "1.1.0",
        what = "tidy()",
        details = "tidySummarizedExperiment says: tidy() is not needed anymore."
    )
}

#' @importFrom methods as
#'
#' @param object A SummarizedExperiment object
#'
#' @export
tidy.SummarizedExperiment <- tidy_

#' @importFrom methods as
#'
#' @param object A SummarizedExperiment object
#'
#' @export
tidy.RangedSummarizedExperiment <- tidy_