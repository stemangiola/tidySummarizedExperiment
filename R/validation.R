#' @importFrom magrittr equals
#' @importFrom dplyr n
is_rectangular <- function(.data, se) {
    is_rectangular_sample <-
        .data %>%
        count(!!s_(se)$symbol ) %>%
        count(n, name="nn") %>%
        nrow() %>%
        equals(1)

    is_rectangular_transcript <-
        .data %>%
        count(!!f_(se)$symbol) %>%
        count(n, name="nn") %>%
        nrow() %>%
        equals(1)

    is_rectangular_sample & is_rectangular_transcript
}

is_not_duplicated <- function(.data, se ) {
    .data %>%
        count(!!s_(se)$symbol , !!f_(se)$symbol) %>%
        filter(n > 1) %>%
        nrow() %>%
        equals(0)
}
