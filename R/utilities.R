#' Get matrix from tibble
#'
#' @keywords internal
#' 
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#' as_matrix(head(dplyr::select(tidybulk::counts_mini, feature, count)), rownames=feature)
#'
#' @noRd
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {
  rownames = enquo(rownames)
  tbl %>%
    
    # Through warning if data frame is not numerical beside the rownames column (if present)
    when(
      do_check &&
        (.) %>%
        # If rownames defined eliminate it from the data frame
        when(!quo_is_null(rownames) ~ (.)[,-1], ~ (.)) %>%
        dplyr::summarise_all(class) %>%
        tidyr::gather(variable, class) %>%
        pull(class) %>%
        unique() %>%
        `%in%`(c("numeric", "integer")) %>% not() %>% any() ~ {
        warning("tidybulk says: there are NON-numerical columns, the matrix will NOT be numerical")
        (.)
      },
      ~ (.)
    ) %>%
    as.data.frame() %>%
    
    # Deal with rownames column if present
    when(
      !quo_is_null(rownames) ~ (.) %>%
        magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
        select(-1),
      ~ (.)
    ) %>%
    
    # Convert to matrix
    as.matrix()
}

#' @importFrom tibble as_tibble
#' @importFrom SummarizedExperiment colData
#'
#' @keywords internal
#'
#' @param .data A tidySummarizedExperiment
#'
#' @noRd
to_tib <- function(.data) {
    colData(.data) %>%
        as.data.frame() %>%
        as_tibble(rownames="cell")
}

# Greater than
gt <- function(a, b) {
    a > b
}

# Smaller than
st <- function(a, b) {
    a < b
}

# Negation
not <- function(is) {
    !is
}

# Raise to the power
pow <- function(a, b) {
    a^b
}

# Equals
eq <- function(a, b) {
    a == b
}

prepend <- function(x, values, before=1) {
    n <- length(x)
    stopifnot(before > 0 && before <= n)
    if (before == 1) {
        c(values, x)
    }
    else {
        c(x[seq_len(before - 1)], values, x[before:n])
    }
}
#' Add class to abject
#'
#'
#' @keywords internal
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
#' 
#' @noRd
add_class <- function(var, name) {
    if (!name %in% class(var)) class(var) <- prepend(class(var), name)

    var
}

#' Remove class to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
#' 
#' @noRd
drop_class <- function(var, name) {
    class(var) <- class(var)[!class(var) %in% name]
    var
}

#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom utils tail
#' @importFrom SummarizedExperiment assays
#'
#' @param .data A tidySummarizedExperiment
#' @param transcripts A character
#' @param all A boolean
#'
#'
#' @return A tidySummarizedExperiment object
#'
#'
#' @noRd
get_abundance_sc_wide <- function(.data, transcripts=NULL, all=FALSE) {

    # Solve CRAN warnings
    . <- NULL

    # For SCE there is not filed for variable features
    variable_feature <- c()

    # Check if output would be too big without forcing
    if (
        length(variable_feature) == 0 &
            is.null(transcripts) &
            all == FALSE
    ) {
        stop("
                Your object does not contain variable feature labels,
                feature argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of feature names
                3. set all=TRUE (this will output a very large object, does your computer have enough RAM?)
                ")
    }

    # Get variable features if existing
    if (
        length(variable_feature) > 0 &
            is.null(transcripts) &
            all == FALSE
    ) {
        variable_genes <- variable_feature
    } # Else
    else {
        variable_genes <- NULL
    }

    # Just grub last assay
    assays(.data) %>%
        as.list() %>%
        tail(1) %>%
        .[[1]] %>%
        when(
            variable_genes %>% is.null() %>% `!`() ~ (.)[variable_genes, , drop=FALSE],
            transcripts %>% is.null() %>% `!`() ~ (.)[transcripts, , drop=FALSE],
            ~ stop("It is not convenient to extract all genes, you should have either variable features or feature list to extract")
        ) %>%
        as.matrix() %>%
        t() %>%
        as_tibble(rownames="cell")
}

#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom purrr when
#' @importFrom purrr map2
#' @importFrom SummarizedExperiment assays
#'
#' @param .data A tidySummarizedExperiment
#' @param transcripts A character
#' @param all A boolean
#' @param exclude_zeros A boolean
#'
#' @return A tidySummarizedExperiment object
#'
#'
#' @noRd
get_abundance_sc_long <- function(.data, transcripts=NULL, all=FALSE, exclude_zeros=FALSE) {

    # Solve CRAN warnings
    . <- NULL

    # For SCE there is not filed for variable features
    variable_feature <- c()

    # Check if output would be too big without forcing
    if (
        length(variable_feature) == 0 &
            is.null(transcripts) &
            all == FALSE
    ) {
        stop("
                Your object does not contain variable feature labels,
                feature argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of feature names
                3. set all=TRUE (this will output a very large object, does your computer have enough RAM?)
                ")
    }


    # Get variable features if existing
    if (
        length(variable_feature) > 0 &
            is.null(transcripts) &
            all == FALSE
    ) {
        variable_genes <- variable_feature
    } # Else
    else {
        variable_genes <- NULL
    }

    assay_names <- assays(.data) %>% names()


    assays(.data) %>%
        as.list() %>%

        # Take active assay
        map2(
            assay_names,

            ~ .x %>%
                when(
                    variable_genes %>% is.null() %>% `!`() ~ .x[variable_genes, , drop=FALSE],
                    transcripts %>% is.null() %>% `!`() ~ .x[toupper(rownames(.x)) %in% toupper(transcripts), , drop=FALSE],
                    all ~ .x,
                    ~ stop("It is not convenient to extract all genes, you should have either variable features or feature list to extract")
                ) %>%

                # Replace 0 with NA
                when(exclude_zeros ~ (.) %>% {
                    x <- (.)
                    x[x == 0] <- NA
                    x
                }, ~ (.)) %>%
                as.matrix() %>%
                data.frame() %>%
                as_tibble(rownames="feature") %>%
                tidyr::pivot_longer(
                    cols=-feature,
                    names_to="cell",
                    values_to="abundance" %>% paste(.y, sep="_"),
                    values_drop_na=TRUE
                )
            # %>%
            # mutate_if(is.character, as.factor) %>%
        ) %>%
        Reduce(function(...) left_join(..., by=c("feature", "cell")), .)
}

#' @importFrom methods .hasSlot
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment rowData<-
#'
#' @keywords internal
#'
#' @param .data A tibble
#' @param SummarizedExperiment_object A tidySummarizedExperiment
#'
#' @noRd
update_SE_from_tibble <- function(.data_mutated, .data, column_belonging = NULL) {

    # Comply to CRAN notes 
    . <- NULL 

    # Get the colnames of samples and feature datasets
    colnames_col <- 
      colnames(colData(.data)) %>% 
      c("sample") %>%
      
      # Forcefully add the column I know the source. This is useful in nesting 
      # where a unique value cannot be linked to sample or feature
      c(names(column_belonging[column_belonging=="sample"]))
    
    colnames_row <- .data %>%
      when(
        .hasSlot(., "rowData") | .hasSlot(., "elementMetadata") ~ colnames(rowData(.)), 
        TRUE ~ c()
      ) %>% 
      c("feature") %>%
      
      # Forcefully add the column I know the source. This is useful in nesting 
      # where a unique value cannot be linked to sample or feature
      c(names(column_belonging[column_belonging=="feature"]))
    
    col_data <-
        .data_mutated %>%
      
        # Eliminate special columns that are read only. Assays
        select_if(!colnames(.) %in% get_special_columns(.data)) %>%

        # Replace for subset
        select(sample, 
            get_subset_columns(., sample) %>%
                 
           # Eliminate feature column
             setdiff(colnames_row)
           
        ) %>%
        distinct() %>%

        # In case unitary SE subset does not work
        select_if(!colnames(.) %in% colnames_row) %>%
        data.frame(row.names=.$sample) %>%
        select(-sample) %>%
        DataFrame()

    row_data <-
        .data_mutated %>%
      
       # Eliminate special columns that are read only 
        select_if(!colnames(.) %in% get_special_columns(.data)) %>%

        # Replace for subset
        select(`feature`, 
               get_subset_columns(., feature) %>%
                 
                 # Eliminate sample column
                 setdiff(colnames_col)
        ) %>%
        distinct() %>%

        # In case unitary SE subset does not work because all same
        select_if(!colnames(.) %in% c(colnames_col, colnames(col_data))) %>%
        data.frame(row.names=.$feature) %>%
        select(-feature) %>%
        DataFrame()

    # Subset if needed. This function is used by many dplyr utilities
    .data <- .data[rownames(row_data), rownames(col_data)]

    # Update
    colData(.data) <- col_data
    rowData(.data) <- row_data
    
    colnames_assay <-
      colnames(.data_mutated) %>% 
      setdiff(c("feature", "sample", get_GRanges_colnames())) %>%
      setdiff(colnames(col_data)) %>% 
      setdiff(colnames(row_data)) %>%
      setdiff(assays(.data) %>% names)
    
    if(length(colnames_assay)>0)
      assays(.data, withDimnames=FALSE) = 
        assays(.data) %>% c(
          .data_mutated %>% 
            
            # Select assays column
            select(feature, sample, colnames_assay) %>% 
            
            # Pivot for generalising to many assays
            pivot_longer(cols = -c(sample, feature)) %>%
            nest(data___ = c(sample, feature, value)) %>%
            
            # Convert to matrix and to named list
            mutate(data___ = map2(
              data___, name,
              ~ .x %>%
                spread(sample, value) %>% 
                as_matrix(rownames = feature)  %>% 
                suppressWarnings() %>%
                list() %>%
                setNames(.y)
            )) %>%
            
            # Create correct list
            pull(data___) %>%
            reduce(c) 
        )
    
    # return
    .data
}


#' @importFrom purrr map_chr
#' @importFrom dplyr select_if
#'
#' @keywords internal
#'
#' @param SummarizedExperiment_object A tidySummarizedExperiment
#'
#' @noRd
#'
get_special_columns <- function(SummarizedExperiment_object) {
    colnames_special <-
        get_special_datasets(SummarizedExperiment_object) %>%

        # In case any of those have feature of sample in column names
        map(
            ~ .x %>%
                select_if(!colnames(.) %in% get_needed_columns()) %>%
                colnames()
        ) %>%
        unlist() %>%
        as.character()

    colnames_counts <-
        get_count_datasets(SummarizedExperiment_object) %>%
        select(-feature, -sample) %>%
        colnames()

    colnames_special %>% c(colnames_counts)
}

#' @importFrom dplyr select
#' @importFrom tidyselect one_of
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom tibble rowid_to_column
#' 
#' @noRd
get_special_datasets <- function(SummarizedExperiment_object) {
  
  SummarizedExperiment_object %>%
    rowRanges() %>%
    when( 
      # If no ranges
      as.data.frame(.) %>%
      nrow() %>%
      equals(0) ~ tibble(),
      
      # If it is a range list (multiple rows per feature)
      class(.) %>% equals("CompressedGRangesList") ~ 
        tibble::as_tibble(.) %>%
        eliminate_GRanges_metadata_columns_also_present_in_Rowdata(SummarizedExperiment_object) %>%
        nest(GenomicRanges = -group_name) %>%
        rename(feature = group_name),
      
      # If standard GRanges (one feature per line)
      ~ {
        transcript_column = 
          rowRanges(SummarizedExperiment_object) %>% 
          as.data.frame() %>% 
          lapply(function(x) rownames(SummarizedExperiment_object)[1] %in% x) %>%
          unlist() %>%
          which() %>%
          names() 
        
        
        # Just rename
        (.) %>%
          
          # If transcript_column exists all good 
          when(
            !is.null(transcript_column) ~  tibble::as_tibble(.) %>%
              eliminate_GRanges_metadata_columns_also_present_in_Rowdata(SummarizedExperiment_object) %>%
              rename(feature := !!transcript_column) ,
            
            # If transcript_column is NULL add numeric column
            ~ tibble::as_tibble(.) %>%
              eliminate_GRanges_metadata_columns_also_present_in_Rowdata(SummarizedExperiment_object) %>%
              rowid_to_column(var = "feature") %>%
              mutate(feature = as.character(feature))
            ) %>%
          
          # Always nest
          nest(GenomicRanges = -feature)
       
      }
    ) %>%
    list()
  
}

#' @importFrom tidyr gather
#' @importFrom dplyr rename
#' @importFrom dplyr left_join
#' @importFrom tibble as_tibble
#' @importFrom purrr reduce
#' @importFrom SummarizedExperiment assays
#' 
#' @noRd
get_count_datasets <- function(SummarizedExperiment_object) {
  map2(
    assays(SummarizedExperiment_object) %>% as.list(),
    names(assays(SummarizedExperiment_object)),
    ~ .x %>%
      tibble::as_tibble(rownames = "feature", .name_repair = "minimal") %>%
      
      # If the matrix does not have sample names, fix column names
      when(colnames(.x) %>% is.null() ~ setNames(., c(
        "feature",  seq_len(ncol(.x)) 
      )),
      ~ (.)
    ) %>%
      
      gather(sample, count,-feature) %>%
      rename(!!.y := count)
  ) %>%
    reduce(left_join, by = c("feature", "sample"))
}

get_needed_columns <- function() {
    c("feature", "sample")
}

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
#' 
#' @noRd
quo_names <- function(v) {
    v <- quo_name(quo_squash(v))
    gsub("^c\\(|`|\\)$", "", v) %>%
        strsplit(", ") %>%
        unlist()
}

#' @importFrom purrr when
#' @importFrom dplyr select
#' @importFrom rlang expr
#' @importFrom tidyselect eval_select
#'
#' @noRd
select_helper <- function(.data, ...) {
    loc <- tidyselect::eval_select(expr(c(...)), .data)

    dplyr::select(.data, loc)
}

outersect <- function(x, y) {
    sort(c(
        setdiff(x, y),
        setdiff(y, x)
    ))
}

#' @importFrom dplyr distinct_at
#' @importFrom dplyr vars
#' @importFrom purrr map
#' @importFrom magrittr equals
#' 
#' @noRd
get_subset_columns <- function(.data, .col) {


    # Comply with CRAN NOTES
    . <- NULL

    # Make col names
    .col <- enquo(.col)

    # x-annotation df
    n_x <- .data %>%
        distinct_at(vars(!!.col)) %>%
        nrow()

    # element wise columns
    .data %>%
        select(-!!.col) %>%
        colnames() %>%
        map(
            ~
            .x %>%
                when(
                    .data %>%
                        distinct_at(vars(!!.col, .x)) %>%
                        nrow() %>%
                        equals(n_x) ~ (.),
                    ~NULL
                )
        ) %>%

        # Drop NULL
        {
            (.)[lengths((.)) != 0]
        } %>%
        unlist()
}

#' @importFrom purrr map_int
#' @importFrom purrr map
is_split_by_transcript = function(.my_data){
    
    
    tot_length = .my_data %>% map(~ pull(.x, feature) ) %>% unlist %>% unique() %>% length
    all_lengths = .my_data %>% map_int(~ pull(.x, feature) %>% unique() %>% length) 
    
    all_lengths %>% unique %>% length() %>% gt(1) |
        (all_lengths != tot_length) %>% any()
    
}

is_split_by_sample = function(.my_data){
    
    
    tot_length = .my_data %>% map(~ pull(.x, sample) ) %>% unlist %>% unique() %>% length
    all_lengths = .my_data %>% map_int(~ pull(.x, sample) %>% unique() %>% length) 
    
    all_lengths %>% unique %>% length() %>% gt(1) |
        (all_lengths != tot_length) %>% any()
    
}


get_GRanges_colnames = function(){
  "GenomicRanges"
}

feature_col_name = "feature"


eliminate_GRanges_metadata_columns_also_present_in_Rowdata = function(.my_data, SummarizedExperiment_object){
  .my_data %>%
    select(-one_of(colnames(rowData(SummarizedExperiment_object)))) %>%
    
    # In case there is not metadata column
    suppressWarnings() 
}

data_frame_returned_message = "tidySummarizedExperiment says: A data frame is returned for independent data analysis."
duplicated_cell_names = "tidySummarizedExperiment says: This operation lead to duplicated feature names. A data frame is returned for independent data analysis."

