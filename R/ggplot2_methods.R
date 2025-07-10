#' @name ggplot
#' @rdname ggplot
#' @inherit ggplot2::ggplot
#' @title Create a new \code{ggplot} from a \code{SummarizedExperiment}
#' @return `ggplot`
#'
#' @examples
#' library(ggplot2)
#' data(pasilla)
#' pasilla %>%
#'     ggplot(aes(.sample, counts)) +
#'     geom_boxplot()
#' 
#' @importFrom purrr map
#' @importFrom rlang quo_name
#' @importFrom ggplot2 aes ggplot
#' @export
ggplot.SummarizedExperiment <- function(data=NULL, mapping=aes(),
    ..., environment=parent.frame()) {

    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(data, .cols)) {
        data <- ping_old_special_column_into_metadata(data)
    }
  
    data %>%
        as_tibble() %>%
        ggplot2::ggplot(mapping=mapping)
}

# S7 method for ggplot2 4.0.0+ compatibility
register_s7_ggplot_method <- function() {
    if (requireNamespace("S7", quietly = TRUE)) {
        S7 <- getNamespace("S7")
        # Check if ggplot and SummarizedExperiment are S7 classes
        is_s7_class <- function(x) inherits(x, "S7_class")
        ggplot_s7 <- tryCatch(get("ggplot", envir = asNamespace("ggplot2")), error = function(e) NULL)
        se_s7 <- tryCatch(get("SummarizedExperiment", envir = asNamespace("SummarizedExperiment")), error = function(e) NULL)
        if (!is.null(ggplot_s7) && !is.null(se_s7) && is_s7_class(ggplot_s7) && is_s7_class(se_s7)) {
            S7$method(ggplot_s7, se_s7) <- function(data=NULL, mapping=aes(), ..., environment=parent.frame()) {
                .cols <- enquos(..., .ignore_empty="all") %>% 
                    map(~ quo_name(.x)) %>% unlist()
                if (is_sample_feature_deprecated_used(data, .cols)) {
                    data <- ping_old_special_column_into_metadata(data)
                }
                data %>%
                    as_tibble() %>%
                    ggplot2::ggplot(mapping=mapping)
            }
        }
    }
}

.onLoad <- function(libname, pkgname) {
    register_s7_ggplot_method()
}
