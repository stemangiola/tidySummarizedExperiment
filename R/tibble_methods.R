#' Coerce lists, matrices, and more to data frames
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' `as_tibble()` turns a SummarizedExperiment existing object into a so-called tibble, a data frame with class `tbl_df`. 
#'
#'
#' @importFrom tibble as_tibble
#'
#'
#' @param x A `SummarizedExperiment`
#' @param ... This parameter includes `.subset` that can be set to any tidyselect expression. For example .subset = c(sample, type), or .subset = contains("PC").
#'
#' @return A tibble
#'
#' @rdname tibble-methods
#' @name as_tibble
#'
#' @export
#' @examples
#' 
#' tidySummarizedExperiment::pasilla %>%
#'     as_tibble()
#'     
#' tidySummarizedExperiment::pasilla %>%
#'     as_tibble(.subset = -c(condition, type))     
#'     
#'     
#'     
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
                                rownames=pkgconfig::get_config("tibble::rownames", NULL), 
                                use_old_special_names = FALSE){
  
  .subset = enquo(.subset)
  
  if(use_old_special_names){
    sample_name = "sample"
    feature_name = "feature"
  }
  
  sample_info <-
    colData(x) %>% 
    
    # If reserved column names are present add .x
    change_reserved_column_names(feature_name, sample_name) %>%
  
    # Convert to tibble
    tibble::as_tibble(rownames=sample_name)
  
  range_info <-
    skip_GRanges %>%
    when(
      (.) ~ tibble() %>% list,
      ~  get_special_datasets(x) 
    ) %>%
    reduce(left_join, by="coordinate") 
    
  gene_info <-
    rowData(x) %>% 
    
    # If reserved column names are present add .x
    change_reserved_column_names(feature_name, sample_name) %>%
  
    # Convert to tibble
    tibble::as_tibble(rownames=feature_name) 
  
  count_info <- get_count_datasets(x)
  
  # Return 
  .subset %>%
    when(
      quo_is_null(.) ~ 
        count_info %>%
        left_join(sample_info, by=sample_name) %>%
        left_join(gene_info, by=feature_name) %>%
        when(nrow(range_info) > 0 ~ (.) %>% left_join(range_info) %>% suppressMessages(), ~ (.)) ,
      ~ subset_tibble_output(count_info, sample_info, gene_info, range_info, !!.subset)
    )
  
}