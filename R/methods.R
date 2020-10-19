

setClass("tidySummarizedExperiment", contains=c("SummarizedExperiment", "RangedSummarizedExperiment"))

#' @importFrom methods show
#' @import SummarizedExperiment
#' @importFrom magrittr %>%
setMethod(
    f="show",
    signature="tidySummarizedExperiment",
    definition=function(object) {
        object %>%
            print()
    }
)



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

tidy_ <- function(object) {
    object %>%
        as("RangedSummarizedExperiment") %>%
        as("tidySummarizedExperiment") %>%

        # If there is a column called sample change it's name
        when(
            "sample" %in% colnames(colData(.)) ~ {
                warning("tidySummarizedExperiment says: column sample in your colData have been renamed as sample_, since is a reserved column name.")

                rename(., sample_=sample)
            },
            ~ (.)
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
