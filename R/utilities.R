#' Drop attribute to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the attribute
#' 
#' @noRd
#'
#' @return A tibble with an additional attribute
drop_attr = function(var, name) {
  attr(var, name) <- NULL
  var
}

#' Drop all attribute to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @noRd
#' 
#' @return A tibble with an additional attribute
drop_all_attr = function(var) {
  
  N = var %>% attributes() %>% names()
  
  for(n in N){
    var = var %>% drop_attr(n)
  }

  var
}

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
                as_tibble(rownames=f_(.data)$name) %>%
                tidyr::pivot_longer(
                    cols=-!!f_(.data)$symbol,
                    names_to="cell",
                    values_to="abundance" %>% paste(.y, sep="_"),
                    values_drop_na=TRUE
                )
            # %>%
            # mutate_if(is.character, as.factor) %>%
        ) %>%
        Reduce(function(...) left_join(..., by=c(f_(.data)$name, "cell")), .)
}

#' @importFrom methods .hasSlot
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom S4Vectors head 
#'
#' @keywords internal
#'
#' @param .data A tibble
#' @param se A tidySummarizedExperiment
#'
#' @noRd
update_SE_from_tibble <- function(.data_mutated, se, column_belonging = NULL) {

    # Comply to CRAN notes 
    . <- NULL 

    # Get the colnames of samples and feature datasets
    colnames_col <- 
      colnames(colData(se)) %>% 
      c(s_(se)$name) %>%
      
      # Forcefully add the column I know the source. This is useful in nesting 
      # where a unique value cannot be linked to sample or feature
      c(names(column_belonging[column_belonging==s_(se)$name]))
    
    colnames_row <- se %>%
      when(
        .hasSlot(., "rowData") | .hasSlot(., "elementMetadata") ~ colnames(rowData(.)), 
        TRUE ~ c()
      ) %>% 
      c(f_(se)$name) %>%
      
      # Forcefully add the column I know the source. This is useful in nesting 
      # where a unique value cannot be linked to sample or feature
      c(names(column_belonging[column_belonging==f_(se)$name]))
    
    secial_columns = get_special_columns(
      
      # Decrease the size of the dataset
      se[1:min(100, nrow(se)), 1:min(20, ncol(se))]
    ) 
    
    new_colnames_col = 
      .data_mutated %>%
      select_if(!colnames(.) %in% setdiff(colnames_col, s_(se)$name)) %>% 
      
      # Eliminate special columns that are read only. Assays
      select_if(!colnames(.) %in% secial_columns) %>%
      select_if(!colnames(.) %in% colnames_row) %>%
      # Replace for subset
      select(!!s_(se)$symbol,     get_subset_columns(., !!s_(se)$symbol)   ) %>% 
      colnames()
     
    col_data <-
        .data_mutated %>%
      
        select(c(colnames_col,new_colnames_col)) %>%
      
        # Make fast distinct()
        filter(pull(., !!s_(se)$symbol) %>% duplicated() %>%  not()) %>% 
      
        # In case unitary SE subset does not work
        data.frame(row.names=pull(., !!s_(se)$symbol), check.names = FALSE) %>%
        select(-!!s_(se)$symbol) %>%
        DataFrame(check.names = FALSE)

    # This to avoid the mismatch between missing column names for counts 
    # and numerical row names for colData
    row_names_col = 
      col_data %>%
      rownames() %>%
      when(
        colnames(se) %>% is.null ~ as.integer(.),
        ~ (.)
      )
    
    row_data <-
        .data_mutated %>%
      
       # Eliminate special columns that are read only 
        select_if(!colnames(.) %in% secial_columns) %>%
      
        #eliminate sample columns directly
        select_if(!colnames(.) %in% c(s_(se)$name, colnames(col_data))) %>%

        # select(one_of(colnames(rowData(se))))
        # Replace for subset
        select(!!f_(se)$symbol,  get_subset_columns(., !!f_(se)$symbol) ) %>%
      
        # Make fast distinct()
        filter(pull(., !!f_(se)$symbol) %>% duplicated() %>%  not()) %>% 

        # In case unitary SE subset does not work because all same
        data.frame(row.names=pull(.,f_(se)$symbol), check.names = FALSE) %>%
        select(-!!f_(se)$symbol) %>%
        DataFrame(check.names = FALSE)
    
    # This to avoid the mismatch between missing row names for counts 
    # and numerical row names for rowData
    row_names_row = 
      row_data %>%
      rownames() %>%
      when(
        rownames(se) %>% is.null ~ as.integer(.),
        ~ (.)
      )
    
    # Subset if needed. This function is used by many dplyr utilities
    se <- se[row_names_row, row_names_col]

    # Update
    colData(se) <- col_data
    rowData(se) <- row_data
    
    # Count-like data that is NOT in the assay slot already 
    colnames_assay <-
      colnames(.data_mutated) %>% 
      setdiff(c(f_(se)$name, s_(se)$name, colnames(as.data.frame(head(rowRanges(se), 1))) )) %>%
      setdiff(colnames(col_data)) %>% 
      setdiff(colnames(row_data)) %>%
      setdiff(assays(se) %>% names)
    
    if(length(colnames_assay)>0)
      assays(se) = #, withDimnames=FALSE) = 
        assays(se, withDimnames = FALSE) %>% c(
          .data_mutated %>% 
            
            # Select assays column
            select(!!f_(se)$symbol, !!s_(se)$symbol, colnames_assay) %>% 
            
            # Pivot for generalising to many assays
            pivot_longer(cols = -c(!!s_(se)$symbol, !!f_(se)$symbol)) %>%
            nest(data___ = c(!!s_(se)$symbol, !!f_(se)$symbol, value)) %>%
            
            # Convert to matrix and to named list
            mutate(data___ = map2(
              data___, name,
              ~ {
                .x = 
                  .x %>%
                  spread(!!s_(se)$symbol, value) %>% 
            
                  
                  as_matrix(rownames = !!f_(se)$symbol)  %>% 
                  suppressWarnings()
                
                # Rearrange if assays has colnames and rownames
                if(!is.null(rownames(se)) & !is.null(rownames(.x))) .x = .x[rownames(se),,drop=FALSE]
                if(!is.null(colnames(se)) & !is.null(colnames(.x))) .x = .x[,colnames(se),drop=FALSE]

                
              .x %>%
                
                list() %>%
                setNames(.y)
              }
            )) %>%
            
            # Create correct list
            pull(data___) %>%
            reduce(c) 
        )
    
    # return
    se
}

#' @importFrom methods is
slice_optimised <- function(.data, ..., .preserve=FALSE) {
  
  # This simulated tibble only gets samples and features so we know those that have been completely omitted already
  # In order to save time for the as_tibble conversion
  simulated_slice = 
    simulate_feature_sample_from_tibble(.data) %>% 
    dplyr::slice(..., .preserve=.preserve)
  
  .data = 
    .data %>%
    
    # Subset the object for samples and features present in the simulated data
    .[rownames(.) %in% simulated_slice[,f_(.)$name], colnames(.) %in% simulated_slice[,s_(.)$name]] %>% 
    inner_join(simulated_slice, by = c(f_(.)$name, s_(.)$name)) 
  
  # If order do not match with the one proposed by slice convert to tibble
  if(.data %>% is("tbl") %>% not()){
    .data_for_match = .data %>%  select(!!f_(.data)$symbol, !!s_(.data)$symbol) %>% as_tibble() 
    
    x = c(pull(.data_for_match, !!f_(.data)$symbol), pull(.data_for_match, !!s_(.data)$symbol))
    y =  c(pull(simulated_slice, !!f_(.data)$symbol), pull(simulated_slice, !!s_(.data)$symbol))
    
    if( identical(x, y) %>% not())
      left_join(simulated_slice, as_tibble(.data), by = c(f_(.data)$name, s_(.data)$name))
    else
      .data
  }
  
 
}


#' @importFrom purrr map_chr
#' @importFrom dplyr select_if
#'
#' @keywords internal
#'
#' @param se A tidySummarizedExperiment
#'
#' @noRd
#'
get_special_columns <- function(.data) {
    colnames_special <-
        get_special_datasets(.data) %>%

        # In case any of those have feature of sample in column names
        map(
            ~ .x %>%
                select_if(!colnames(.) %in% get_needed_columns(.data)) %>%
                colnames()
        ) %>%
        unlist() %>%
        as.character()

    colnames_counts <-
        get_count_datasets(.data) %>%
        select(-!!f_(.data)$symbol, -!!s_(.data)$symbol) %>%
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
get_special_datasets <- function(se) {
  
  rr =  se %>%
    rowRanges() 
  
  rr %>%
    when( 
      
      # If no ranges
      as.data.frame(.) %>%
      nrow() %>%
      equals(0) ~ tibble(),
      
      # If it is a range list (multiple rows per feature)
      is(., "CompressedGRangesList") ~ {
        
        # If GRanges does not have row names
        if(is.null(rr@partitioning@NAMES)) rr@partitioning@NAMES = as.character(1:nrow(se))
        
        tibble::as_tibble(rr) %>%
          eliminate_GRanges_metadata_columns_also_present_in_Rowdata(se) %>%
          nest(GRangesList = -group_name) %>%
          rename(!!f_(se)$symbol := group_name)
        
      },
      
      # If standard GRanges (one feature per line)
      ~ {
        
        # If GRanges does not have row names
       if(is.null(rr@ranges@NAMES)) rr@ranges@NAMES = as.character(1:nrow(se))
         
        tibble::as_tibble(rr) %>% 
          eliminate_GRanges_metadata_columns_also_present_in_Rowdata(se) %>% 
        mutate(!!f_(se)$symbol := rr@ranges@NAMES) 
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
#' @importFrom magrittr equals
#' 
#' @noRd
get_count_datasets <- function(se) { 
  
  # Stop if column names of assays do not overlap
  if( check_if_assays_are_NOT_overlapped(se) ) 
    stop( 
    "tidySummarizedExperiment says: the assays in your SummarizedExperiment have column names, 
but their order is not the same, and they not completely overlap." 
  )

  # Join assays
  map2( 
    assays(se, withDimnames = FALSE) %>% as.list(),
    names(assays(se)),
    ~ {
      
      # If the counts are in a sparse matrix convert to a matrix
      # This might happen because the user loaded tidySummarizedExperiment and is 
      # print a SingleCellExperiment
      if(is(.x, "dgCMatrix") | is(.x, "DelayedArray")) {
        .x = as.matrix(.x) 
      }
      
      # Rearrange if assays has colnames and rownames
      if(!is.null(rownames(se)) & !is.null(rownames(.x))) .x = .x[rownames(se),,drop=FALSE]
      if(!is.null(colnames(se)) & !is.null(colnames(.x))) .x = .x[,colnames(se),drop=FALSE]

      # If I don't have assay colnames and rownames add them
      if(!is.null(rownames(se)) & is.null(rownames(.x))) rownames(.x) = rownames(se) 
      if(!is.null(colnames(se)) & is.null(colnames(.x))) colnames(.x) = colnames(se) 
      
      .x %>%
      # matrix() %>%
      # as.data.frame() %>% 
        
      tibble::as_tibble(rownames = f_(se)$name, .name_repair = "minimal") %>%
      
      # If the matrix does not have sample names, fix column names
      when(colnames(.x) %>% is.null() & is.null(colnames(se)) ~ setNames(., c(
        f_(se)$name,  seq_len(ncol(.x)) 
      )),
      ~ (.)
    ) %>%
      
      gather(!!s_(se)$symbol, !!.y,-!!f_(se)$symbol) 
    
    #%>%
    #  rename(!!.y := count)
  }) %>%
    when(
      length(.)>0 ~ 
        
        # reduce(., left_join, by = c(f_(se)$name, s_(se)$name)),
        bind_cols(.,  .name_repair = c("minimal")) %>% .[!duplicated(colnames(.))], 
      ~ expand.grid(
        rownames(se), colnames(se)
        ) %>% 
        setNames(c(f_(se)$name, s_(se)$name)) %>%
        tibble::as_tibble()
    ) %>% 
    
    # Add dummy sample or feature if we have empty assay. 
    # This is needed for a correct isualisation of the tibble form
    when(
      f_(se)$name %in% colnames(.) %>% not ~ mutate(., !!f_(se)$symbol := as.character(NA)),
      s_(se)$name %in% colnames(.) %>% not ~ mutate(., !!s_(se)$symbol := as.character(NA)),
      ~ (.)
    )
}

get_needed_columns <- function(.data) {
    c(f_(.data)$name, s_(.data)$name)
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
    
  se = .my_data[[1]]
    tot_length = .my_data %>% map(~ pull(.x, !!f_(se)$symbol) ) %>% unlist %>% unique() %>% length
    all_lengths = .my_data %>% map_int(~ pull(.x, !!f_(se)$symbol) %>% unique() %>% length) 
    
    all_lengths %>% unique %>% length() %>% gt(1) |
        (all_lengths != tot_length) %>% any()
    
}

is_split_by_sample = function(.my_data){
    
    se = .my_data[[1]]
    tot_length = .my_data %>% map(~ pull(.x, !!s_(se)$symbol) ) %>% unlist %>% unique() %>% length
    all_lengths = .my_data %>% map_int(~ pull(.x, !!s_(se)$symbol) %>% unique() %>% length) 
    
    all_lengths %>% unique %>% length() %>% gt(1) |
        (all_lengths != tot_length) %>% any()
    
}


get_GRanges_colnames = function(){
  "GenomicRanges"
}


eliminate_GRanges_metadata_columns_also_present_in_Rowdata = function(.my_data, se){
  .my_data %>%
    select(-one_of(colnames(rowData(se)))) %>%
    
    # In case there is not metadata column
    suppressWarnings() 
}

subset_tibble_output = function(.data, count_info, sample_info, gene_info, range_info, .subset){
  # This function outputs a tibble after subsetting the columns
  .subset = enquo(.subset)
  
  # Build template of the output
  output_colnames = 
    slice(count_info, 0) %>%
    left_join(slice(sample_info, 0), by=s_(.data)$name) %>%
    left_join(slice(gene_info, 0), by = f_(.data)$name) %>%
    when(nrow(range_info) > 0 ~ (.) %>% left_join(range_info, by=f_(.data)$name), ~ (.)) %>%
    select(!!.subset) %>%
    colnames()
  
  
  # Sample table
  sample_info = 
    sample_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(s_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )
  
  # Ranges table
  range_info = 
    range_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(f_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )
  
  # Ranges table
  gene_info = 
    gene_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(f_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )
  
  # Ranges table
  count_info = 
    count_info %>%
    when(
      colnames(.) %>% intersect(output_colnames) %>% length() %>% equals(0) ~ NULL,
      select(., one_of(f_(.data)$name, s_(.data)$name, output_colnames)) %>%
        suppressWarnings()
    )
  
  
  if(
    !is.null(count_info) & 
    (
      !is.null(sample_info) & !is.null(gene_info) | 
    
      # Make exception for weirs cases (e.g. c(sample, counts))
      (colnames(count_info) %>% outersect(c(f_(.data)$name, s_(.data)$name)) %>% length() %>% gt(0))
    )
  ) {
    output_df = 
      count_info %>%
      when(!is.null(sample_info) ~ (.) %>% left_join(sample_info, by=s_(.data)$name), ~ (.)) %>%
      when(!is.null(gene_info) ~ (.) %>% left_join(gene_info, by=f_(.data)$name), ~ (.)) %>%
      when(!is.null(range_info) ~ (.) %>% left_join(range_info, by=f_(.data)$name), ~ (.))
  }
  else if(!is.null(sample_info) ){
    output_df = sample_info
  }
  else if(!is.null(gene_info)){
    output_df = gene_info %>%
      
      # If present join GRanges
      when(!is.null(range_info) ~ (.) %>% left_join(range_info, by=f_(.data)$name), ~ (.))
  }
  
  output_df %>%
    
    # Cleanup
    select(one_of(output_colnames)) %>%
    suppressWarnings()
  
}

#' @importFrom stringr str_replace
change_reserved_column_names = function(col_data, .data ){
  
  col_data %>%
    
    setNames(
      colnames(.) %>% 
        sapply(function(x) if(x==f_(.data)$name) sprintf("%s.x", f_(.data)$name) else x) %>% 
        sapply(function(x) if(x==s_(.data)$name) sprintf("%s.x", s_(.data)$name) else x) %>% 
        str_replace("^coordinate$", "coordinate.x")
    ) 
  
}

choose_name_if_present = function(x){
  columns_query = c()
  for(i in 1:length(x)){
    if(is.null(names(x[i]))) columns_query[i] = x[i]
    else columns_query[i] = names(x[i])
  }
  
  columns_query
}

#' @importFrom purrr when
join_efficient_for_SE <- function(x, y, by=NULL, copy=FALSE, suffix=c(".x", ".y"), join_function, force_tibble_route = FALSE,
                                  ...) {
  
  
   
  # Comply to CRAN notes 
  . <- NULL 
  
  # Deprecation of special column names
  if(is_sample_feature_deprecated_used( x, when(by, !is.null(.) ~ by, ~ colnames(y)))){
    x= ping_old_special_column_into_metadata(x)
  }
  
  # Get the colnames of samples and feature datasets
  colnames_col <- get_colnames_col(x)
  colnames_row <- get_rownames_col(x)
  
  # See if join done by sample, feature or both
  columns_query = by %>% when(
    !is.null(.) ~ choose_name_if_present(.), 
    ~ colnames(y) %>% intersect(c(colnames_col, colnames_row))
  )
  
  if(
    
    # Complex join that it is not efficient yet
    (any(columns_query %in% colnames_row) & any(columns_query %in% colnames_col)) |
    
    # If join is with something else, the inefficient generic solution might work, 
    # or delegate the proper error downstream
    (!any(columns_query %in% colnames_row) & !any(columns_query %in% colnames_col)) |
    
    # Needed for internal recurrence if outcome is not valid
    force_tibble_route) {
    
    
    
    # If I have a big dataset
    if(ncol(x)>100) message("tidySummarizedExperiment says: if you are joining a dataframe both sample-wise and feature-wise, for efficiency (until further development), it is better to separate your joins and join datasets sample-wise OR feature-wise.")
    
    x %>%
      as_tibble(skip_GRanges  = TRUE) %>%
      join_function(y, by=by, copy=copy, suffix=suffix, ...) %>%
      when(
        
        # If duplicated sample-feature pair returns tibble
        !is_not_duplicated(., x) | !is_rectangular(., x) ~ {
          message(duplicated_cell_names)
          message(data_frame_returned_message)
          (.)
        },
        
        # Otherwise return updated tidySummarizedExperiment
        ~ update_SE_from_tibble(., x)
      )
    
  }
  
  # Join only feature-wise
  else if(any(columns_query %in% colnames_row) & !any(columns_query %in% colnames_col)) {
    
    row_data_tibble = 
      rowData(x) %>% 
      as_tibble(rownames = f_(x)$name) %>%  
      join_function(y, by=by, copy=copy, suffix=suffix, ...) 
    
    # Check if the result is not SE then take the tibble route
    if(
      is.na(pull(row_data_tibble, !!f_(x)$symbol)) %>% any | 
      duplicated(pull(row_data_tibble, !!f_(x)$symbol)) %>% any |
      pull(row_data_tibble, !!f_(x)$symbol) %>% setdiff(rownames(colData(x))) %>% length() %>% gt(0)
    ) return(join_efficient_for_SE(x, y, by=by, copy=copy, suffix=suffix, join_function, force_tibble_route = TRUE, ...))
    
    row_data = 
      row_data_tibble %>% 
      data.frame(row.names=pull(., !!f_(x)$symbol)) %>%
      select(-!!f_(x)$symbol) %>%
      DataFrame()
    
    # Subset in case of an inner join, or a right join
    x = x[rownames(row_data),]  
    
    # Tranfer annotation
    rowData(x) = row_data
    
    # Return
    x
  }
  
  # Join only sample-wise
  else if(any(columns_query %in% colnames_col) & !any(columns_query %in% colnames_row)) {
    
    col_data_tibble = 
      colData(x) %>% 
      as_tibble(rownames = s_(x)$name) %>%  
      join_function(y, by=by, copy=copy, suffix=suffix, ...)
    
    # Check if the result is not SE then take the tibble route
    if(
      is.na(pull(col_data_tibble, !!s_(x)$symbol)) %>% any | 
      duplicated(pull(col_data_tibble, !!s_(x)$symbol)) %>% any |
      pull(col_data_tibble, !!s_(x)$symbol) %>% setdiff(rownames(colData(x))) %>% length() %>% gt(0)
    ) return(join_efficient_for_SE(x, y, by=by, copy=copy, suffix=suffix, join_function, force_tibble_route = TRUE, ...))
    
    col_data = 
      col_data_tibble %>% 
      data.frame(row.names=pull(., !!s_(x)$symbol)) %>%
      select(-!!s_(x)$symbol) %>%
      DataFrame()
    
    
    # Subset in case of an inner join, or a right join
    x = x[,rownames(col_data)]  
    
    # Tranfer annotation
    colData(x) = col_data
    
    # Return
    x
  }
  
  else stop("tidySummarizedExperiment says: ERROR FOR DEVELOPERS: this option should not exist. In join utility.")
  
  
  
}

get_ellipse_colnames = function(...){
  
  (enquos(..., .ignore_empty = "all") %>% map(~ quo_name(.x)) %>% unlist)
}

get_colnames_col = function(x){
  colnames(colData(x)) %>% 
    c(s_(x)$name) 
}

get_rownames_col = function(x){
  x %>%
    when(
      .hasSlot(., "rowData") | .hasSlot(., "elementMetadata") ~ colnames(rowData(.)), 
      TRUE ~ c()
    ) %>% 
    c(f_(x)$name) 

}

# This function is used for the change of special sample column to .sample
# Check if "sample" is included in the query and is not part of any other existing annotation
#' @importFrom stringr str_detect
#' @importFrom stringr regex
is_sample_feature_deprecated_used = function(.data, user_columns, use_old_special_names = FALSE){
  
  old_standard_is_used_for_sample =  
    (
      ( any(str_detect(user_columns  , regex("\\bsample\\b"))) & !any(str_detect(user_columns  , regex("\\W*(\\.sample)\\W*")))  ) |
        "sample" %in% user_columns 
    ) & 
    !"sample" %in% c(colnames(rowData(.data)), colnames(colData(.data)))
  
  old_standard_is_used_for_feature = 
    (
      ( any(str_detect(user_columns  , regex("\\bfeature\\b"))) & !any(str_detect(user_columns  , regex("\\W*(\\.feature)\\W*")))  ) |
        "feature" %in% user_columns 
    ) & 
    !"feature" %in% c(colnames(rowData(.data)), colnames(colData(.data)))
  
  old_standard_is_used = old_standard_is_used_for_sample | old_standard_is_used_for_feature
  
  if(old_standard_is_used){
    warning("tidySummarizedExperiment says: from version 1.3.1, the special columns including sample/feature id (colnames(se), rownames(se)) has changed to \".sample\" and \".feature\". This dataset is returned with the old-style vocabulary (feature and sample), however we suggest to update your workflow to reflect the new vocabulary (.feature, .sample)")
    
    use_old_special_names = TRUE
  }
  
  use_old_special_names
}
  
data_frame_returned_message = "tidySummarizedExperiment says: A data frame is returned for independent data analysis."
duplicated_cell_names = "tidySummarizedExperiment says: This operation lead to duplicated feature names. A data frame is returned for independent data analysis."

# Key column names
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
ping_old_special_column_into_metadata = function(.data){
  
  metadata(.data)$feature__ = get_special_column_name_symbol("feature")
  metadata(.data)$sample__ = get_special_column_name_symbol("sample")
  
  .data
}

get_special_column_name_symbol = function(name){
  list(name = name, symbol = as.symbol(name))
}

# This function produce artificially the feature ans sample column, 
# to make optimisation before as_tibble is called
# for big datasets
simulate_feature_sample_from_tibble = function(.data){
  
  r = rownames(.data) %>% .[rep(1:length(.), ncol(.data) )]
  c = colnames(.data) %>% .[rep(1:length(.), each=nrow(.data) )]
  
    tibble(!!f_(.data)$symbol := r,  !!s_(.data)$symbol := c)
  
}



feature__ =  get_special_column_name_symbol(".feature")
sample__ = get_special_column_name_symbol(".sample")

#' @importFrom S4Vectors metadata
f_ =  function(x){
  # Check if old deprecated columns are used
  if("feature__" %in% names(metadata(x))) feature__ = metadata(x)$feature__
  return(feature__)
}

#' @importFrom S4Vectors metadata
s_ = function(x){
  if("sample__" %in% names(metadata(x))) sample__ = metadata(x)$sample__
  return(sample__)
}

split_SummarizedExperiment_by_feature_to_list = function(.data){
  if(nrow(.data)>1000)
    message("tidySummarizedExperiment says: grouping a SummarizedExperiment by feature takes 1 minute for ~ 10,000 features.")
  map(1:nrow(.data), ~ .data[.x,])
}

#' Add attribute to abject
#'
#' @keywords internal
#' @noRd
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
  attr(var, name) <- attribute
  var
}

is_filer_columns_in_column_selection = function(.data, ...){
  # columns = enquos(columns)
  tryCatch({
    .data |>
      slice(0) |>
      dplyr::filter(..., .preserve=.preserve)
    TRUE
  },
  error = function(e) FALSE)
}

check_if_assays_are_NOT_consistently_ordered = function(se){
  
  # If I have any assay at all
  assays(se) |> length() |> gt(0) &&
    
    # If I have more than one assay with colnames
    Filter(
      Negate(is.null),
      assays(se, withDimnames = FALSE) |>  
        as.list() |> 
        map(colnames)
    ) |> 
    length() |>
    gt(0) &&
    
  # If I have lack of consistency
  se |> 
    assays(withDimnames = FALSE) |>
    as.list() |>
    purrr::map_dfr(colnames) |>
    apply(1, function(x) x |> unique() |> length()) |>
    equals(1) |>
    all() |>  
    not()
}

check_if_assays_are_NOT_overlapped = function(se){
  
  # If I have any assay at all
  assays(se) |> length() |> gt(0) &&
    
    # If I have more than one assay with colnames
    Filter(
      Negate(is.null),
      assays(se, withDimnames = FALSE) |>  
        as.list() |> 
        map(colnames)
    ) |> 
    length() |>
    gt(0) &&
    
    # If I have lack of consistency
    assays(se, withDimnames = FALSE) |>  
    as.list() |> 
    map(colnames) |> 
    reduce(intersect) |> 
    length() |> 
    equals(ncol(se)) |> 
    not()
}


order_assays_internally_to_be_consistent = function(se){

  se |> 
    assays(withDimnames = FALSE) =
      map2(
        assays(se, withDimnames = FALSE) %>% as.list(),
        names(assays(se)),
        ~ {
          
          if(!is.null(rownames(se)) & !is.null(rownames(.x))) .x = .x[rownames(se),,drop=FALSE]
          if(!is.null(colnames(se)) & !is.null(colnames(.x))) .x = .x[,colnames(se),drop=FALSE]
          
          .x
          
        })
    
    se
}