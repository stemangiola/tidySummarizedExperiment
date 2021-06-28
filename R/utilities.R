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
                as_tibble(rownames=feature_name) %>%
                tidyr::pivot_longer(
                    cols=-!!feature_symbol,
                    names_to="cell",
                    values_to="abundance" %>% paste(.y, sep="_"),
                    values_drop_na=TRUE
                )
            # %>%
            # mutate_if(is.character, as.factor) %>%
        ) %>%
        Reduce(function(...) left_join(..., by=c(feature_name, "cell")), .)
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
      c(sample_name) %>%
      
      # Forcefully add the column I know the source. This is useful in nesting 
      # where a unique value cannot be linked to sample or feature
      c(names(column_belonging[column_belonging==sample_name]))
    
    colnames_row <- .data %>%
      when(
        .hasSlot(., "rowData") | .hasSlot(., "elementMetadata") ~ colnames(rowData(.)), 
        TRUE ~ c()
      ) %>% 
      c(feature_name) %>%
      
      # Forcefully add the column I know the source. This is useful in nesting 
      # where a unique value cannot be linked to sample or feature
      c(names(column_belonging[column_belonging==feature_name]))
    
    col_data <-
        .data_mutated %>%
      
        # Eliminate special columns that are read only. Assays
        select_if(!colnames(.) %in% get_special_columns(.data)) %>%

        # Replace for subset
        select(!!sample_symbol, 
            get_subset_columns(., !!sample_symbol) %>%
                 
           # Eliminate feature column
             setdiff(colnames_row)
           
        ) %>%
        distinct() %>%

        # In case unitary SE subset does not work
        select_if(!colnames(.) %in% colnames_row) %>%
        data.frame(row.names=pull(., !!sample_symbol)) %>%
        select(-!!sample_symbol) %>%
        DataFrame()

    row_data <-
        .data_mutated %>%
      
       # Eliminate special columns that are read only 
        select_if(!colnames(.) %in% get_special_columns(.data)) %>%

        # Replace for subset
        select(!!feature_symbol, 
               get_subset_columns(., !!feature_symbol) %>%
                 
                 # Eliminate sample column
                 setdiff(colnames_col)
        ) %>%
        distinct() %>%

        # In case unitary SE subset does not work because all same
        select_if(!colnames(.) %in% c(colnames_col, colnames(col_data))) %>%
        data.frame(row.names=pull(.,feature_symbol)) %>%
        select(-!!feature_symbol) %>%
        DataFrame()
    
    # This to avoid the mismatch between missing row names for counts 
    # and numerical row names for rowData
    row_names_row = 
      row_data %>%
      rownames() %>%
      when(
        rownames(.data) %>% is.null ~ as.integer(.),
        ~ (.)
      )
    
    # This to avoid the mismatch between missing column names for counts 
    # and numerical row names for colData
    row_names_col = 
      col_data %>%
      rownames() %>%
      when(
        colnames(.data) %>% is.null ~ as.integer(.),
        ~ (.)
      )
    
    # Subset if needed. This function is used by many dplyr utilities
    .data <- .data[row_names_row, row_names_col]

    # Update
    colData(.data) <- col_data
    rowData(.data) <- row_data
    
    colnames_assay <-
      colnames(.data_mutated) %>% 
      setdiff(c(feature_name, sample_name, get_GRanges_colnames())) %>%
      setdiff(colnames(col_data)) %>% 
      setdiff(colnames(row_data)) %>%
      setdiff(assays(.data) %>% names)
    
    if(length(colnames_assay)>0)
      assays(.data, withDimnames=FALSE) = 
        assays(.data) %>% c(
          .data_mutated %>% 
            
            # Select assays column
            select(!!feature_symbol, !!sample_symbol, colnames_assay) %>% 
            
            # Pivot for generalising to many assays
            pivot_longer(cols = -c(!!sample_symbol, !!feature_symbol)) %>%
            nest(data___ = c(!!sample_symbol, !!feature_symbol, value)) %>%
            
            # Convert to matrix and to named list
            mutate(data___ = map2(
              data___, name,
              ~ .x %>%
                spread(!!sample_symbol, value) %>% 
                as_matrix(rownames = !!feature_symbol)  %>% 
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
        select(-!!feature_symbol, -!!sample_symbol) %>%
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
  
  rr =  SummarizedExperiment_object %>%
    rowRanges() 
  
  rr %>%
    when( 
      
      # If no ranges
      as.data.frame(.) %>%
      nrow() %>%
      equals(0) ~ tibble(),
      
      # If it is a range list (multiple rows per feature)
      class(.) %>% equals("CompressedGRangesList") ~ {
        
        # If GRanges does not have row names
        if(is.null(rr@partitioning@NAMES)) rr@partitioning@NAMES = as.character(1:nrow(SummarizedExperiment_object))
        
        tibble::as_tibble(rr) %>%
          #mutate(!!feature_symbol := rr@partitioning@NAMES)
          eliminate_GRanges_metadata_columns_also_present_in_Rowdata(SummarizedExperiment_object) %>%
          nest(GRangesList = -group_name) %>%
          rename(!!feature_symbol := group_name)
        
      },
      
      # If standard GRanges (one feature per line)
      ~ {
        
        # If GRanges does not have row names
       if(is.null(rr@ranges@NAMES)) rr@ranges@NAMES = as.character(1:nrow(SummarizedExperiment_object))
         
        tibble::as_tibble(rr) %>% 
        mutate(!!feature_symbol := rr@ranges@NAMES) 
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
      tibble::as_tibble(rownames = feature_name, .name_repair = "minimal") %>%
      
      # If the matrix does not have sample names, fix column names
      when(colnames(.x) %>% is.null() ~ setNames(., c(
        feature_name,  seq_len(ncol(.x)) 
      )),
      ~ (.)
    ) %>%
      
      gather(!!sample_symbol, count,-!!feature_symbol) %>%
      rename(!!.y := count)
  ) %>%
    when(
      length(.)>0 ~ reduce(., left_join, by = c(feature_name, sample_name)),
      ~ expand.grid(
        rownames(SummarizedExperiment_object), colnames(SummarizedExperiment_object)
        ) %>% 
        setNames(c(feature_name, sample_name)) %>%
        tibble::as_tibble()
    )
    
}

get_needed_columns <- function() {
    c(feature_name, sample_name)
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
    
    
    tot_length = .my_data %>% map(~ pull(.x, !!feature_symbol) ) %>% unlist %>% unique() %>% length
    all_lengths = .my_data %>% map_int(~ pull(.x, !!feature_symbol) %>% unique() %>% length) 
    
    all_lengths %>% unique %>% length() %>% gt(1) |
        (all_lengths != tot_length) %>% any()
    
}

is_split_by_sample = function(.my_data){
    
    
    tot_length = .my_data %>% map(~ pull(.x, !!sample_symbol) ) %>% unlist %>% unique() %>% length
    all_lengths = .my_data %>% map_int(~ pull(.x, !!sample_symbol) %>% unique() %>% length) 
    
    all_lengths %>% unique %>% length() %>% gt(1) |
        (all_lengths != tot_length) %>% any()
    
}


get_GRanges_colnames = function(){
  "GenomicRanges"
}


eliminate_GRanges_metadata_columns_also_present_in_Rowdata = function(.my_data, SummarizedExperiment_object){
  .my_data %>%
    select(-one_of(colnames(rowData(SummarizedExperiment_object)))) %>%
    
    # In case there is not metadata column
    suppressWarnings() 
}

subset_tibble_output = function(count_info, sample_info, gene_info, range_info, .subset){
  
  .subset = enquo(.subset)
  
  # Build template of the output
  output_colnames = 
    slice(count_info, 0) %>%
    left_join(slice(sample_info, 0), by=sample_name) %>%
    left_join(slice(gene_info, 0), by = feature_name) %>%
    when(nrow(range_info) > 0 ~ (.) %>% left_join(range_info, by=feature_name), ~ (.)) %>%
    select(!!.subset) %>%
    colnames()
  
  
  # Sample table
  sample_info = 
    sample_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(sample_name, output_colnames)) %>%
        suppressWarnings()
    )
  
  # Ranges table
  range_info = 
    range_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(feature_name, output_colnames)) %>%
        suppressWarnings()
    )
  
  # Ranges table
  gene_info = 
    gene_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(feature_name, output_colnames)) %>%
        suppressWarnings()
    )
  
  # Ranges table
  count_info = 
    count_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(feature_name, sample_name, output_colnames)) %>%
        suppressWarnings()
    )
  
  
  if(
    !is.null(sample_info) & !is.null(gene_info) | 
    
    # Make exception for weirs cases (e.g. c(sample, counts))
    (colnames(count_info) %>% outersect(c(feature_name, sample_name)) %>% length() %>% gt(0))
  ) {
    output_df = 
      count_info %>%
      when(!is.null(sample_info) ~ (.) %>% left_join(sample_info, by=sample_name), ~ (.)) %>%
      when(!is.null(gene_info) ~ (.) %>% left_join(gene_info, by=feature_name), ~ (.)) %>%
      when(!is.null(range_info) ~ (.) %>% left_join(range_info, by=feature_name), ~ (.))
  }
  else if(!is.null(sample_info) ){
    output_df = sample_info
  }
  else if(!is.null(gene_info)){
    output_df = gene_info %>%
      
      # If present join GRanges
      when(!is.null(range_info) ~ (.) %>% left_join(range_info, by=feature_name), ~ (.))
  }
  
  output_df %>%
    
    # Cleanup
    select(one_of(output_colnames)) %>%
    suppressWarnings()
  
}

change_reserved_column_names = function(.data){
  
  .data %>%
    
    setNames(
      colnames(.) %>% 
        str_replace(feature_name, sprintf("%s.x", feature_name)) %>% 
        str_replace(sample_name, sprintf("%s.x", sample_name)) %>% 
        str_replace("^coordinate$", "coordinate.x")
    ) 
  
}

#' @importFrom purrr when
join_efficient_for_SE <- function(x, y, by=NULL, copy=FALSE, suffix=c(".x", ".y"), join_function, force_tibble_route = FALSE,
                                  ...) {
  
  
  
  # Comply to CRAN notes 
  . <- NULL 
  
  # Get the colnames of samples and feature datasets
  colnames_col <- get_colnames_col(x)
  colnames_row <- get_rownames_col(x)
  
  # See if join done by sample, feature or both
  columns_query = by %>% when(
    !is.null(.) ~ sapply(., function(.x) ifelse(!is.null(names(.x)), names(.x), .x)), 
    ~ colnames(y)
  )
  
 
columns_query %>%
    when(
      
      # Complex join that it is not efficient yet
      (any(columns_query %in% colnames_row) & any(columns_query %in% colnames_col)) |
        
        # If join is with something else, the inefficient generic solution might work, 
        # or delegate the proper error downstream
        (!any(columns_query %in% colnames_row) & !any(columns_query %in% colnames_col)) |
        
        # Needed for internal recurrency if outcome is not valid
        force_tibble_route ~ {
          
          message("tidySummarizedExperiment says: either your resulting dataset is not a valid SummarizedExperiment or you are joining a dataframe both sample-wise and feature-wise. In the latter case, for efficiency (until further development), it is better to separate your joins and join datasets sample-wise OR feature-wise.")
          x %>%
            as_tibble(skip_GRanges  = T) %>%
            join_function(y, by=by, copy=copy, suffix=suffix, ...) %>%
            when(
              
              # If duplicated sample-feature pair returns tibble
              !is_not_duplicated(.) ~ {
                message(duplicated_cell_names)
                (.)
              },
              
              # Otherwise return updated tidySummarizedExperiment
              ~ update_SE_from_tibble(., x)
            )
          
        },
      
      # Join only feature-wise
      any(columns_query %in% colnames_row) & !any(columns_query %in% colnames_col) ~ {
        
        row_data_tibble = 
          rowData(x) %>% 
          as_tibble(rownames = feature_name) %>%  
          join_function(y, by=by, copy=copy, suffix=suffix, ...) 
        
        # Check if the result is not SE then take the tibble route
        if(
          is.na(pull(row_data_tibble, !!feature_symbol)) %>% any | 
          duplicated(pull(row_data_tibble, !!feature_symbol)) %>% any |
          pull(row_data_tibble, !!feature_symbol) %>% setdiff(rownames(colData(x))) %>% length() %>% gt(0)
        ) return(join_efficient_for_SE(x, y, by=by, copy=copy, suffix=suffix, join_function, force_tibble_route = TRUE, ...))
        
        row_data = 
          row_data_tibble %>% 
          data.frame(row.names=pull(., !!feature_symbol)) %>%
          select(-!!feature_symbol) %>%
          DataFrame()
        
        # Subset in case of an inner join, or a right join
        x = x[rownames(row_data),]  
        
        # Tranfer annotation
        rowData(x) = row_data
        
        # Return
        x
      },
      
      # Join only sample-wise
      any(columns_query %in% colnames_col) & !any(columns_query %in% colnames_row) ~ {
        
        col_data_tibble = 
          colData(x) %>% 
          as_tibble(rownames = sample_name) %>%  
          join_function(y, by=by, copy=copy, suffix=suffix, ...)
        
        # Check if the result is not SE then take the tibble route
        if(
          is.na(pull(col_data_tibble, !!sample_symbol)) %>% any | 
          duplicated(pull(col_data_tibble, !!sample_symbol)) %>% any |
          pull(col_data_tibble, !!sample_symbol) %>% setdiff(rownames(colData(x))) %>% length() %>% gt(0)
        ) return(join_efficient_for_SE(x, y, by=by, copy=copy, suffix=suffix, join_function, force_tibble_route = TRUE, ...))
          
       col_data = 
         col_data_tibble %>% 
          data.frame(row.names=pull(., !!sample_symbol)) %>%
          select(-!!sample_symbol) %>%
          DataFrame()
        
          
        # Subset in case of an inner join, or a right join
        x = x[,rownames(col_data)]  
        
        # Tranfer annotation
        colData(x) = col_data
        
        # Return
        x
      },
      
      ~ stop("tidySummarizedExperiment says: ERROR FOR DEVELOPERS: this option should not exist. In join utility.")
      
    )
  
}

get_ellipse_colnames = function(...){
  
  (enquos(..., .ignore_empty = "all") %>% map(~ quo_name(.x)) %>% unlist)
}

get_colnames_col = function(x){
  colnames(colData(x)) %>% 
    c(sample_name) 
}

get_rownames_col = function(x){
  x %>%
    when(
      .hasSlot(., "rowData") | .hasSlot(., "elementMetadata") ~ colnames(rowData(.)), 
      TRUE ~ c()
    ) %>% 
    c(feature_name) 
}

data_frame_returned_message = "tidySummarizedExperiment says: A data frame is returned for independent data analysis."
duplicated_cell_names = "tidySummarizedExperiment says: This operation lead to duplicated feature names. A data frame is returned for independent data analysis."

# Key column names
feature_name = ".feature"
sample_name = ".sample"
feature_symbol = as.symbol(feature_name)
sample_symbol = as.symbol(sample_name)
