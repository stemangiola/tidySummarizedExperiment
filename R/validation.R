#' @importFrom magrittr equals
#' @importFrom dplyr n
is_rectangular <- function(.data) {
    is_rectangular_sample <-
        .data %>%
        count(!!sample_symbol ) %>%
        count(n, name="nn") %>%
        nrow() %>%
        equals(1)

    is_rectangular_transcript <-
        .data %>%
        count(!!feature_symbol) %>%
        count(n, name="nn") %>%
        nrow() %>%
        equals(1)

    is_rectangular_sample & is_rectangular_transcript
}

is_not_duplicated <- function(.data) {
    .data %>%
        count(!!sample_symbol , !!feature_symbol) %>%
        filter(n > 1) %>%
        nrow() %>%
        equals(0)
}
