

setClass("tidySE", contains=c("SummarizedExperiment", "RangedSummarizedExperiment"))

#' @importFrom methods show
#' @import SummarizedExperiment
#' @importFrom magrittr %>%
setMethod(
    f="show",
    signature="tidySE",
    definition=function(object) {
        x = object %>%
            as_tibble() %>%
            print()


        cli::cat_line(format(x,  n = NULL, width = NULL, n_extra = NULL))

        mat <- trunc_mat(x,  n = NULL, width = NULL, n_extra = NULL)
        format(mat)

        invisible(x)
    }
)



#' tidy for SummarizedExperiment
#'
#' @param object A SummarizedExperiment object
#'
#' @return A tidySE object
#'
#' @examples
#'
#' tidySE::pasilla %>% tidy()
#' @export
tidy <- function(object) {
    UseMethod("tidy", object)
}

tidy_ <- function(object) {
    object %>%
        as("RangedSummarizedExperiment") %>%
        as("tidySE") %>%

        # If there is a column called sample change it's name
        when(
            "sample" %in% colnames(colData(.)) ~ {
                warning("tidySE says: column sample in your colData have been renamed as sample_, since is a reserved column name.")

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
