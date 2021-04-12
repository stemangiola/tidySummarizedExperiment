#' Coerce lists, matrices, and more to data frames
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' `as_tibble()` turns an existing object, such as a data frame or
#' matrix, into a so-called tibble, a data frame with class [`tbl_df`]. This is
#' in contrast with [tibble()], which builds a tibble from individual columns.
#' `as_tibble()` is to [`tibble()`] as [base::as.data.frame()] is to
#' [base::data.frame()].
#'
#' `as_tibble()` is an S3 generic, with methods for:
#' * [`data.frame`][base::data.frame()]: Thin wrapper around the `list` method
#'   that implements tibble's treatment of [rownames].
#' * [`matrix`][methods::matrix-class], [`poly`][stats::poly()],
#'   [`ts`][stats::ts()], [`table`][base::table()]
#' * Default: Other inputs are first coerced with [base::as.data.frame()].
#'
#' @importFrom tibble as_tibble
#'
#' @section Row names:
#' The default behaviour is to silently remove row names.
#'
#' New code should explicitly convert row names to a new column using the
#' `rownames` argument.
#'
#' For existing code that relies on the retention of row names, call
#' `pkgconfig::set_config("tibble::rownames"=NA)` in your script or in your
#' package's [.onLoad()]  function.
#'
#' @section Life cycle:
#' Using `as_tibble()` for vectors is superseded as of version 3.0.0,
#' prefer the more expressive maturing `as_tibble_row()` and
#' `as_tibble_col()` variants for new code.
#'
#' @seealso [tibble()] constructs a tibble from individual columns. [enframe()]
#'   converts a named vector to a tibble with a column of names and column of
#'   values. Name repair is implemented using [vctrs::vec_as_names()].
#'
#' @param x A data frame, list, matrix, or other object that could reasonably be
#'   coerced to a tibble.
#' @param ... Unused, for extensibility.
#' @param rownames How to treat existing row names of a data frame or matrix:
#'   * `NULL`: remove row names. This is the default.
#'   * `NA`: keep row names.
#'   * A string: the name of a new column. Existing rownames are transferred
#'     into this column and the `row.names` attribute is deleted.
#'  Read more in [rownames].
#' @param .name_repair see tidyr
#'
#'   For compatibility only, do not use for new code.
#' @return A tibble
#'
#' @rdname tibble-methods
#' @name as_tibble
#'
#' @export
#' @examples
#' tidySummarizedExperiment::pasilla %>%
#'     
#'     as_tibble()
NULL

#' @export
#' @importFrom purrr reduce
#' @importFrom purrr map
#' @importFrom tidyr spread
#' @importFrom tibble enframe
#' @importFrom SummarizedExperiment colData
#'
#'
as_tibble.SummarizedExperiment <- function(x, ...,
    .name_repair=c("check_unique", "unique", "universal", "minimal"),
    rownames=pkgconfig::get_config("tibble::rownames", NULL)) {
  
  .as_tibble_optimised(x = x, ..., .name_repair=.name_repair, rownames=rownames)

}

.as_tibble_optimised = function(x, skip_GRanges = F, .subset = NULL,
                                .name_repair=c("check_unique", "unique", "universal", "minimal"),
                                rownames=pkgconfig::get_config("tibble::rownames", NULL)){
  
  .subset = enquo(.subset)
  
  sample_info <-
    colData(x) %>% 
    
    # If reserved column names are present add .x
    setNames(
      colnames(.) %>% 
        str_replace("^sample$", "sample.x")
    ) %>%
    
    # Convert to tibble
    tibble::as_tibble(rownames="sample")
  
  # range_info =
  #     x@rowRanges %>%
  #     as.data.frame %>%
  #     tibble::as_tibble(rownames="feature")
  range_info <-
    skip_GRanges %>%
    when(
      (.) ~ tibble() %>% list,
      ~  get_special_datasets(x) 
    ) %>%
    reduce(left_join, by="coordinate") %>%
    
    # If reserved column names are present add .x
    setNames(
      colnames(.) %>% 
        str_replace("^coordinate$", "coordinate.x")
    ) 
     
  
  gene_info <-
    rowData(x) %>%
    
    # If reserved column names are present add .x
    setNames(
      colnames(.) %>% 
        str_replace("^feature$", "feature.x")
    ) %>%
    
    # Convert to tibble
    tibble::as_tibble(rownames="feature") 
  
  count_info <- get_count_datasets(x)
  
  # Return 
  .subset %>%
    when(
      quo_is_null(.) ~ 
        count_info %>%
        left_join(sample_info, by="sample") %>%
        left_join(gene_info, by="feature") %>%
        when(nrow(range_info) > 0 ~ (.) %>% left_join(range_info, by="feature"), ~ (.)) ,
      ~ subset_tibble_output(count_info, sample_info, gene_info, range_info, !!.subset)
    )
  
}