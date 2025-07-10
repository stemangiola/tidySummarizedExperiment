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
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org
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
