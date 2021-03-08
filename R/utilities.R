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
#' @keywords internal
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
                Your object does not contain variable transcript labels,
                transcript argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of transcript names
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
            ~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
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
                Your object does not contain variable transcript labels,
                transcript argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of transcript names
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
                    ~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
                ) %>%

                # Replace 0 with NA
                when(exclude_zeros ~ (.) %>% {
                    x <- (.)
                    x[x == 0] <- NA
                    x
                }, ~ (.)) %>%
                as.matrix() %>%
                data.frame() %>%
                as_tibble(rownames="transcript") %>%
                tidyr::pivot_longer(
                    cols=-transcript,
                    names_to="cell",
                    values_to="abundance" %>% paste(.y, sep="_"),
                    values_drop_na=TRUE
                )
            # %>%
            # mutate_if(is.character, as.factor) %>%
        ) %>%
        Reduce(function(...) left_join(..., by=c("transcript", "cell")), .)
}

#' @importFrom methods .hasSlot
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#'
#' @keywords internal
#'
#' @param .data A tibble
#' @param SummarizedExperiment_object A tidySummarizedExperiment
#'
#' @noRd
update_SE_from_tibble <- function(.data_mutated, .data) {

    # Comply to CRAN notes
    . <- NULL

    # Get the colnames of samples and transcript datasets
    colnames_col <- colnames(colData(.data)) %>% c("sample")
    colnames_row <- when(.hasSlot(.data, "rowData") ~ colnames(rowData(.data)), ~ c()) %>% c("transcript")

    col_data <-
        .data_mutated %>%
        select_if(!colnames(.) %in% get_special_columns(.data)) %>%

        # Replace for subset
        select(sample, get_subset_columns(., sample)) %>%
        distinct() %>%

        # In case unitary SE subset does not ork
        select_if(!colnames(.) %in% colnames_row) %>%
        data.frame(row.names=.$sample) %>%
        select(-sample) %>%
        DataFrame()

    row_data <-
        .data_mutated %>%
        select_if(!colnames(.) %in% get_special_columns(.data)) %>%

        # Replace for subset
        select(`transcript`, get_subset_columns(., transcript)) %>%
        distinct() %>%

        # In case unitary SE subset does not work because all same
        select_if(!colnames(.) %in% c(colnames_col, colnames(col_data))) %>%
        data.frame(row.names=.$transcript) %>%
        select(-transcript) %>%
        DataFrame()


    # Subset if needed. This function is used by many dplyr utilities
    .data <- .data[rownames(row_data), rownames(col_data)]

    # Update
    colData(.data) <- col_data
    rowData(.data) <- row_data

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

        # In case any of those have transcript of sample in column names
        map(
            ~ .x %>%
                select_if(!colnames(.) %in% get_needed_columns()) %>%
                colnames()
        ) %>%
        unlist() %>%
        as.character()

    colnames_counts <-
        get_count_datasets(SummarizedExperiment_object) %>%
        select(-transcript, -sample) %>%
        colnames()

    colnames_special %>% c(colnames_counts)
}

#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom SummarizedExperiment rowRanges
get_special_datasets <- function(SummarizedExperiment_object) {
    if (
        "RangedSummarizedExperiment" %in% .class2(SummarizedExperiment_object) &

        rowRanges(SummarizedExperiment_object) %>%
                as.data.frame() %>%
                nrow() %>%
                gt(0)
    ) {
        rowRanges(SummarizedExperiment_object) %>%
            as.data.frame() %>%

            # Take off rowData columns as there is a recursive anomaly within gene ranges
            suppressWarnings(
              select(-one_of(colnames(rowData(SummarizedExperiment_object))))
            ) %>%
            tibble::as_tibble(rownames="transcript") %>%
            list()
    } else {
        tibble() %>% list()
    }
}

#' @importFrom tidyr gather
#' @importFrom dplyr rename
#' @importFrom dplyr left_join
#' @importFrom tibble as_tibble
#' @importFrom purrr reduce
#' @importFrom SummarizedExperiment assays
get_count_datasets <- function(SummarizedExperiment_object) {
    map2(
        assays(SummarizedExperiment_object) %>% as.list(),
        names(assays(SummarizedExperiment_object)),
        ~ .x %>%
            tibble::as_tibble(rownames="transcript") %>%
            gather(sample, count, -transcript) %>%
            rename(!!.y := count)
    ) %>%
        reduce(left_join, by=c("transcript", "sample"))
}

get_needed_columns <- function() {
    c("transcript", "sample")
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

data_frame_returned_message = "tidySummarizedExperiment says: A data frame is returned for independent data analysis."
duplicated_cell_names = "tidySummarizedExperiment says: This operation lead to duplicated transcript names. A data frame is returned for independent data analysis."
