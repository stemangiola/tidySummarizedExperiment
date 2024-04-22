#' @name ggplot
#' @rdname ggplot
#' @inherit ggplot2::ggplot
#' @title Create a new \code{ggplot} from a \code{tidyseurat}
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
