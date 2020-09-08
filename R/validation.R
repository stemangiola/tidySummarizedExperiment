#' @importFrom magrittr equals
is_rectangular = function(.data){

  is_rectangular_sample =
    .data %>%
    count(sample) %>%
    count(n, name = "nn") %>%
    nrow %>%
    equals(1)

  is_rectangular_transcript =
    .data %>%
    count(transcript) %>%
    count(n, name = "nn") %>%
    nrow %>%
    equals(1)

  is_rectangular_sample & is_rectangular_transcript

}

is_not_duplicated = function(.data){
  .data %>%
    count( sample, transcript) %>%
    filter(n > 1) %>%
    nrow %>%
    equals(0)
}
