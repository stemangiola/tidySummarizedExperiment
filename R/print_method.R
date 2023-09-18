# This file is a replacement of the unexported functions in the tibble
# package, in order to specify "tibble abstraction in the header"

#' @name tbl_format_header
#' @rdname tbl_format_header
#' @inherit pillar::tbl_format_header
#' 
#' @examples
#' # TODO
#' 
#' @importFrom rlang names2
#' @importFrom pillar align
#' @importFrom pillar get_extent
#' @importFrom pillar style_subtle
#' @importFrom pillar tbl_format_header
#' @export
tbl_format_header.tidySummarizedExperiment <- function(x, setup, ...) {
  
    number_of_features <- x |> attr("number_of_features")
    number_of_samples <- x |> attr("number_of_samples")
    named_header <- x |> attr("named_header")
    assay_names <- x |> attr("assay_names")

  
    if (all(names2(named_header) == "")) {
        header <- named_header
    } else {
        header <-
            paste0(
                align(paste0(names2(named_header), ":"), space=NBSP),
                " ",
                named_header
            ) %>%
            # Add further info single-cell
            append(sprintf(
                "\033[90m Features=%s | Samples=%s | Assays=%s\033[39m",
                number_of_features,
                number_of_samples,
                assay_names %>% paste(collapse=", ")
            ), after = 1)
    }
    style_subtle(pillar___format_comment(header, width=setup$width))
}

#' @name formatting
#' @rdname formatting
#' @aliases print
#' @inherit tibble::formatting
#' @return Prints a message to the console describing
#'   the contents of the `tidySummarizedExperiment`.
#' 
#' @param n_extra Number of extra columns to print abbreviated information for,
#'   if the width is too small for the entire tibble. If `NULL`, the default,
#'   will print information about at most `tibble.max_extra_cols` extra columns.
#' 
#' @examples
#' data(pasilla)
#' print(pasilla)
#' 
#' @importFrom vctrs new_data_frame
#' @importFrom SummarizedExperiment assayNames
#' @importFrom stats setNames
#' @export
print.SummarizedExperiment <- function(x, ..., n=NULL,
    width=NULL, n_extra=NULL) {


  # Fix NOTEs
  . <- NULL
  

  # Stop if any column or row names are duplicated
  if (check_if_any_dimnames_duplicated(x, dim = "cols")) {
      stop("tidySummarizedExperiment says: some column names are duplicated")
  }
  if (check_if_any_dimnames_duplicated(x, dim = "rows")) {
      stop("tidySummarizedExperiment says: some row names are duplicated")
  }

  # Stop if column names of assays do not overlap
  if (check_if_assays_are_NOT_overlapped(x, dim = "cols")) { 
      stop( 
          "tidySummarizedExperiment says: the assays in your SummarizedExperiment have column names, 
but they do not completely overlap." 
      )
  }
  if (check_if_assays_are_NOT_overlapped(x, dim = "rows")) { 
      stop( 
          "tidySummarizedExperiment says: the assays in your SummarizedExperiment have row names, 
but they do not completely overlap." 
      )
  }
  
    # reorder assay colnames before printing
    # Rearrange if assays has colnames and rownames
    x <- order_assays_internally_to_be_consistent(x)
    
    my_tibble <-
        x |>
    
    # If I have more than 30 genes select first sample
    when(
      nrow(.) > 30 ~.[1:min(50, nrow(x)), min(1, ncol(x)), drop=FALSE] ,
      ncol(.) == 0 ~ .,
      ~ .[, 1:min(20, ncol(x)), drop=FALSE]
    ) %>%
    
        as_tibble() 
  
    my_tibble |>
        new_data_frame(class=c("tidySummarizedExperiment", "tbl")) %>%
        add_attr(nrow(x),  "number_of_features") %>%
        add_attr(ncol(x),  "number_of_samples") %>%
        add_attr(assays(x) %>% names , "assay_names") %>%
    
    # Set fake dimensions for efficiancy
    add_attr(
        sprintf(
            "%s %s %s", 
            x %>% dim %>% {(.)[1] * (.)[2]} %>%
                format(format="f", big.mark=",", digits=1),
            cli::symbol$times,
            ncol(my_tibble)
        ) %>%
        setNames("A SummarizedExperiment-tibble abstraction"), 
        "named_header"
    ) %>%
    print()
    invisible(x) 
}
