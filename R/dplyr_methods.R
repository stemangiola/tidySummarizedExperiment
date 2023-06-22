#' Efficiently bind multiple data frames by row and column
#'
#' This is an efficient implementation of the common pattern of
#' `do.call(rbind, dfs)` or `do.call(cbind, dfs)` for binding many
#' data frames into one.
#'
#' The output of `bind_rows()` will contain a column if that column
#' appears in any of the inputs.
#'
#' @param ... Data frames to combine.
#'
#'   Each argument can either be a data frame, a list that could be a data
#'   frame, or a list of data frames.
#'
#'   When row-binding, columns are matched by name, and any missing
#'   columns will be filled with NA.
#'
#'   When column-binding, rows are matched by position, so all data
#'   frames must have the same number of rows. To match by value, not
#'   position, see mutate-joins.
#' @param .id Data frame identifier.
#'
#'   When `.id` is supplied, a new column of identifiers is
#'   created to link each row to its original data frame. The labels
#'   are taken from the named arguments to `bind_rows()`. When a
#'   list of data frames is supplied, the labels are taken from the
#'   names of the list. If no names are found a numeric sequence is
#'   used instead.
#' @param add.cell.ids from Seurat 3.0 A character vector of length(x = c(x, y)). Appends the corresponding values to the start of each objects' cell names.
#'
#' @importFrom ttservice bind_rows
#'
#' @return `bind_rows()` and `bind_cols()` return the same type as
#'   the first input, either a data frame, `tbl_df`, or `grouped_df`.
#' @examples
#' 
#' data(se)
#' bind_rows(    se, se  )
#'
#' se_bind = se |> select(dex,  albut)
#' se |> bind_cols(se_bind)
#'
#' 
#' @name bind_rows
NULL

#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @importFrom SummarizedExperiment cbind
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assays<-
#' @importFrom S4Vectors SimpleList
#'
#' @export
#'
bind_rows.SummarizedExperiment <- function(..., .id=NULL, add.cell.ids=NULL) {
    tts <- flatten_if(dots_values(...), is_spliced)

    new_obj <- 
      tts %>%
      when(
        is_split_by_sample(.) & is_split_by_transcript(.) ~ stop("tidySummarizedExperiment says: bind_rows cannot be applied to splits both by sample- and feature-wise information"),
        is_split_by_sample(.) ~ cbind(.[[1]], .[[2]]) ,
        is_split_by_transcript(.) ~ rbind(.[[1]], .[[2]]),
        
        # If there is not split, then bind the samples
        ~ cbind(.[[1]], .[[2]])
      )
    

    # If duplicated sample names
    if (new_obj |> colnames() |> duplicated() |> which() |> length() |> gt(0)) {
        warning("tidySummarizedExperiment says: you have duplicated sample names, they will be made unique.")
    }
    unique_colnames <- make.unique(colnames(new_obj), sep="_")

    colnames(new_obj) <- unique_colnames

    # Change also all assays colnames
    assays(new_obj) <- assays(new_obj)@listData |> map(~ {
        colnames(.x) <- unique_colnames
        .x
    }) |> 
      SimpleList()

    new_obj
}


#'
#' @importFrom ttservice bind_cols
#' 
#' @inheritParams bind_cols
#'
#' @rdname dplyr-methods
NULL

#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
bind_cols_internal = function(..., .id=NULL, column_belonging = NULL) {
    tts <- tts <- flatten_if(dots_values(...), is_spliced)

    tts[[1]] |> 
        as_tibble(skip_GRanges = TRUE) |>
        dplyr::bind_cols(tts[[2]], .id=.id) %>%
        when(

            # If the column added are not sample-wise or feature-wise return tibble
            (colnames(tts[[2]]) %in% c(
                get_subset_columns(., !!s_(tts[[1]])$symbol),
                get_subset_columns(., !!f_(tts[[1]])$symbol)
            )
            ) |> all() ~ update_SE_from_tibble(., tts[[1]], column_belonging = column_belonging),

            # Return tiblle
            ~ {
                warning("tidySummarizedExperiment says: The new columns do not include pure sample-wise or feature-wise. A data frame is returned for independent data analysis.")
                (.)
            }
        )
}

bind_cols_ = function(..., .id=NULL) { bind_cols_internal(..., .id=NULL) }

#' @export
#'
bind_cols.SummarizedExperiment <- bind_cols_

#' @export
#'
bind_cols.RangedSummarizedExperiment <- bind_cols_

#' distinct
#'
#'
#' @param .data A tbl. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#' @param .keep_all If TRUE, keep all variables in .data. If a combination
#'   of ... is not distinct, this keeps the first row of values. (See dplyr)
#'
#' @return A tidySummarizedExperiment object
#'
#' @importFrom dplyr distinct
#'
#' @rdname dplyr-methods
#' @name distinct
#'
#' @examples
#' 
#' data(pasilla)
#' pasilla |> distinct(.sample)
#' 
NULL

#' @inheritParams distinct
#' 
#' 
#' @export
distinct.SummarizedExperiment <- function(.data, ..., .keep_all=FALSE) {
  
  # message(data_frame_returned_message)

  distinct_columns = 
    (enquos(..., .ignore_empty = "all") |> map(~ quo_name(.x)) |> unlist())
  
  # If Ranges column not in query perform fast as_tibble
  skip_GRanges = 
    get_GRanges_colnames() %in% 
    distinct_columns |>
    not()
  
  # Deprecation of special column names
  if(is_sample_feature_deprecated_used(.data, distinct_columns)){
    .data= ping_old_special_column_into_metadata(.data)
    
  }
  
  .data |>
    select(...) |> 
      as_tibble(skip_GRanges = skip_GRanges) |>
      dplyr::distinct(..., .keep_all=.keep_all)
}


#' Subset rows using column values
#'
#' `filter()` retains the rows where the conditions you provide a `TRUE`. Note
#' that, unlike base subsetting with `[`, rows where the condition evaluates
#' to `NA` are dropped.
#'
#' dplyr is not yet smart enough to optimise filtering optimisation
#' on grouped datasets that don't need grouped calculations. For this reason,
#' filtering is often considerably faster on [ungroup()]ed data.
#'
#' @importFrom dplyr filter
#'
#' @section Useful filter functions:
#'
#' * [`==`], [`>`], [`>=`] etc
#' * [`&`], [`|`], [`!`], [xor()]
#' * [is.na()]
#' * [between()], [near()]
#'
#' @section Grouped tibbles:
#'
#' Because filtering expressions are computed within groups, they may
#' yield different results on grouped tibbles. This will be the case
#' as soon as an aggregating, lagging, or ranking function is
#' involved. Compare this ungrouped filtering:
#'
#'
#' The former keeps rows with `mass` greater than the global average
#' whereas the latter keeps rows with `mass` greater than the gender
#'
#' average.
#' @family single table verbs
#' @inheritParams arrange
#' @param .data A tidySummarizedExperiment object or any data frame
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Logical predicates defined in
#'   terms of the variables in `.data`.
#'   Multiple conditions are combined with `&`. Only rows where the
#'   condition evaluates to `TRUE` are kept.
#' @param .preserve when `FALSE` (the default), the grouping structure
#'   is recalculated based on the resulting data, otherwise it is kept as is.
#' @return
#' An object of the same type as `.data`.
#'
#' * Rows are a subset of the input, but appear in the same order.
#' * Columns are not modified.
#' * The number of groups may be reduced (if `.preserve` is not `TRUE`).
#' * Data frame attributes are preserved.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @seealso [filter_all()], [filter_if()] and [filter_at()].
#'
#' @rdname dplyr-methods
#' @name filter
#'
#' @examples
#'
#' data(pasilla)
#' pasilla |>  filter(.sample == "untrt1")
#'
#' # Learn more in ?dplyr_tidy_eval
#' 
NULL

#' @inheritParams filter
#' 
#' @rdname dplyr-methods
#' @name filter
#' 
#' @export
filter.SummarizedExperiment <- function(.data, ..., .preserve=FALSE) {
  
  
  # Deprecation of special column names
  if(is_sample_feature_deprecated_used(
    .data, 
    (enquos(..., .ignore_empty = "all") |> map(~ quo_name(.x)) |> unlist())
  )){
    .data= ping_old_special_column_into_metadata(.data)
  }
  
  # Understand what the filtering is about
  is_filter_from_samples = .data |> colData() |> as_tibble() |>  is_filer_columns_in_column_selection(...)
  is_filter_from_features = .data |> rowData() |> as_tibble() |>  is_filer_columns_in_column_selection(...)
  
  # Get the simpler route if filter is only on samples
  if(is_filter_from_samples & !is_filter_from_features){
    filtered_samples = 
      colData(.data) |> 
      as_tibble(rownames = s_(.data)$name) |> 
      dplyr::filter(..., .preserve=.preserve) |> 
      pull(!!s_(.data)$symbol)
    
    return(.data[,filtered_samples])
  }
  
  # Get the simpler route if filter is only on features
  else if(!is_filter_from_samples & is_filter_from_features){
    filtered_features = 
      rowData(.data) |> 
      as_tibble(rownames = f_(.data)$name) |> 
      dplyr::filter(..., .preserve=.preserve) |> 
      pull(!!f_(.data)$symbol)
    
    return(.data[filtered_features,])
  }  
  
  # If filtering is based on both features and samples
  else{
    
    # Do filtering
    new_meta <- .data |>
      as_tibble(skip_GRanges = TRUE) |>
      dplyr::filter(..., .preserve=.preserve)
    
    # If data cannot be a SummarizedExperiment
    if(!is_rectangular(new_meta, .data)){
      message("tidySummarizedExperiment says: The resulting data frame is not rectangular (all genes for all samples), a tibble is returned for independent data analysis.")
      return(new_meta)
    } else {
      return(.data[  
        unique(pull(new_meta,!!f_(.data)$symbol)), 
        unique(pull(new_meta,!!s_(.data)$symbol)) 
      ])
    }
  }
  
}

#' @param .drop When `.drop=TRUE`, empty groups are dropped. See
#'   [group_by_drop_default()] for what the default value is for this argument.
#' @return A [grouped data frame][grouped_df()], unless the combination of
#'   `...` and `add` yields a non empty set of grouping columns, a
#'   regular (ungrouped) data frame otherwise.
#' @section Methods:
#' These function are **generic**s, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' Methods available in currently loaded packages:
#'
#' @importFrom dplyr group_by_drop_default
#' @importFrom dplyr group_by
#'
#' @rdname dplyr-methods
#' @name group_by
#'
#' @examples
#' 
#' data(pasilla)
#' pasilla  |> group_by(.sample)
#' 
NULL

#' @export
group_by.SummarizedExperiment <- function(.data, ..., .add=FALSE, .drop=group_by_drop_default(.data)) {
    message(data_frame_returned_message)

  # Deprecation of special column names
  if(is_sample_feature_deprecated_used(
    .data, 
    (enquos(..., .ignore_empty = "all") |> map(~ quo_name(.x)) |> unlist())
  )){
    .data= ping_old_special_column_into_metadata(.data)
  }
  
    .data |>
        as_tibble() |>
        dplyr::group_by(..., .add=.add, .drop=.drop)
}


#' Summarise each group to fewer rows
#'
#' @importFrom dplyr summarise
#'
#' @description
#' `summarise()` creates a new data frame. It will have one (or more) rows for
#' each combination of grouping variables; if there are no grouping variables,
#' the output will have a single row summarising all observations in the input.
#' It will contain one column for each grouping variable and one column
#' for each of the summary statistics that you have specified.
#'
#' `summarise()` and `summarize()` are synonyms.
#'
#' @section Useful functions:
#'
#' * Center: [mean()], [median()]
#' * Spread: [sd()], [IQR()], [mad()]
#' * Range: [min()], [max()], [quantile()]
#' * Position: [first()], [last()], [nth()],
#' * Count: [n()], [n_distinct()]
#' * Logical: [any()], [all()]
#'
#' @section Backend variations:
#'
#' The data frame backend supports creating a variable and using it in the
#' same summary. This means that previously created summary variables can be
#' further transformed or combined within the summary, as in [mutate()].
#' However, it also means that summary variables with the same names as previous
#' variables overwrite them, making those variables unavailable to later summary
#' variables.
#'
#' This behaviour may not be supported in other backends. To avoid unexpected
#' results, consider using new names for your summary variables, especially when
#' creating multiple summaries.
#'
#' @param .data A tidySummarizedExperiment object or any data frame
#'
#' @rdname dplyr-methods
#' @name summarise
#'
#' @inheritParams arrange
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Name-value pairs of summary
#'   functions. The name will be the name of the variable in the result.
#'
#'   The value can be:
#'
#'   * A vector of length 1, e.g. `min(x)`, `n()`, or `sum(is.na(y))`.
#'   * A vector of length `n`, e.g. `quantile()`.
#'   * A data frame, to add multiple columns from a single expression.
#' @family single table verbs
#' @return
#' An object _usually_ of the same type as `.data`.
#'
#' * The rows come from the underlying `group_keys()`.
#' * The columns are a combination of the grouping keys and the summary
#'   expressions that you provide.
#' * If `x` is grouped by more than one variable, the output will be another
#'   [grouped_df] with the right-most group removed.
#' * If `x` is grouped by one variable, or is not grouped, the output will
#'   be a [tibble].
#' * Data frame attributes are **not** preserved, because `summarise()`
#'   fundamentally creates a new data frame.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @examples
#' 
#' data(pasilla)
#' pasilla |> summarise(mean(counts))
NULL

#' @export
summarise.SummarizedExperiment <- function(.data, ...) {
    message(data_frame_returned_message)

  # Deprecation of special column names
  if(is_sample_feature_deprecated_used(
    .data, 
    (enquos(..., .ignore_empty = "all") |> map(~ quo_name(.x)) |> unlist())
  )){
    .data= ping_old_special_column_into_metadata(.data)
  }
  
  # If Ranges column not in query perform fast as_tibble
  skip_GRanges = 
    get_GRanges_colnames() %in% 
    (enquos(..., .ignore_empty = "all") |> map(~ quo_name(.x)) |> unlist()) |>
    not()

    .data |>
        as_tibble(skip_GRanges = skip_GRanges) |>
        dplyr::summarise(...)
}


#' Create, modify, and delete columns
#'
#' `mutate()` adds new variables and preserves existing ones;
#' `transmute()` adds new variables and drops existing ones.
#' New variables overwrite existing variables of the same name.
#' Variables can be removed by setting their value to `NULL`.
#'
#' @importFrom dplyr mutate
#'
#' @section Useful mutate functions:
#'
#' * [`+`], [`-`], [log()], etc., for their usual mathematical meanings
#'
#' * [lead()], [lag()]
#'
#' * [dense_rank()], [min_rank()], [percent_rank()], [row_number()],
#'   [cume_dist()], [ntile()]
#'
#' * [cumsum()], [cummean()], [cummin()], [cummax()], [cumany()], [cumall()]
#'
#' * [na_if()], [coalesce()]
#'
#' * [if_else()], [recode()], [case_when()]
#'
#' @section Grouped tibbles:
#'
#' Because mutating expressions are computed within groups, they may
#' yield different results on grouped tibbles. This will be the case
#' as soon as an aggregating, lagging, or ranking function is
#' involved. Compare this ungrouped mutate:
#'
#' @param .data A tidySummarizedExperiment object or any data frame
#'
#' With the grouped equivalent:
#'
#' The former normalises `mass` by the global average whereas the
#' latter normalises by the averages within gender levels.
#'
#' @inheritParams arrange
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Name-value pairs.
#'   The name gives the name of the column in the output.
#'
#'   The value can be:
#'
#'   * A vector of length 1, which will be recycled to the correct length.
#'   * A vector the same length as the current group (or the whole data frame
#'     if ungrouped).
#'   * `NULL`, to remove the column.
#'   * A data frame or tibble, to create multiple columns in the output.
#'
#' @family single table verbs
#' @return
#' An object of the same type as `.data`.
#'
#' For `mutate()`:
#'
#' * Rows are not affected.
#' * Existing columns will be preserved unless explicitly modified.
#' * New columns will be added to the right of existing columns.
#' * Columns given value `NULL` will be removed
#' * Groups will be recomputed if a grouping variable is mutated.
#' * Data frame attributes are preserved.
#'
#' For `transmute()`:
#'
#' * Rows are not affected.
#' * Apart from grouping variables, existing columns will be remove unless
#'   explicitly kept.
#' * Column order matches order of expressions.
#' * Groups will be recomputed if a grouping variable is mutated.
#' * Data frame attributes are preserved.
#' @section Methods:
#' These function are **generic**s, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' Methods available in currently loaded packages:
#'
#' @examples
#' 
#' data(pasilla)
#' pasilla |> mutate(logcounts=log2(counts))
#'
#' @rdname dplyr-methods
#' @name mutate
#' 
NULL

#' @importFrom dplyr mutate
#' @importFrom rlang enquos
#'
#' @export
mutate.SummarizedExperiment <- function(.data, ...) {
    # Check that we are not modifying a key column
    cols <- enquos(...) |> names()
    
    # Deprecation of special column names
    if(is_sample_feature_deprecated_used(
      .data, 
      (enquos(..., .ignore_empty = "all") |> map(~ quo_name(.x)) |> unlist())
    )){
      .data= ping_old_special_column_into_metadata(.data)
    }
    
    secial_columns = get_special_columns(
      
      # Decrease the size of the dataset
      .data[1:min(100, nrow(.data)), 1:min(20, ncol(.data))]
    ) |> 
      c(get_needed_columns(.data))
    
    tst =
        intersect(
            cols,
            secial_columns
        ) |>
        length() |>
        gt(0)


    if (tst) {
        columns =
          secial_columns |>
            paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says: you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData). If you want to mutate a view-only column, make a copy and mutate that one."
        )
    }

    # If Ranges column not in query perform fast as_tibble
    skip_GRanges = 
      get_GRanges_colnames() %in% 
      cols |>
      not()
    
    .data |>
        as_tibble(skip_GRanges= skip_GRanges) |>
        dplyr::mutate(...) |>
        update_SE_from_tibble(.data)
}



#' Rename columns
#'
#' Rename individual variables using `new_name=old_name` syntax.
#'
#' @importFrom dplyr rename
#'
#' @section Scoped selection and renaming:
#'
#' Use the three scoped variants ([rename_all()], [rename_if()], [rename_at()])
#' to renaming a set of variables with a function.
#'
#' @inheritParams arrange
#'
#' @param .data A tidySummarizedExperiment object or any data frame
#' @param ... <[`tidy-select`][dplyr_tidy_select]> Use `new_name=old_name`
#'   to rename selected variables.
#' @return
#' An object of the same type as `.data`.
#' * Rows are not affected.
#' * Column names are changed; column order is preserved
#' * Data frame attributes are preserved.
#' * Groups are updated to reflect new names.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @family single table verbs
#'
#' @rdname dplyr-methods
#' @name rename
#'
#' @examples
#' 
#' data(pasilla)
#' pasilla |> rename(cond=condition)
NULL

#' @importFrom tidyselect eval_select
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment rowData<-
#' 
#' @export
rename.SummarizedExperiment <- function(.data, ...) {

    
    # Cols data frame
    cols_data_frame = 
      bind_cols(
      colData(.data) |> as_tibble() |> slice(0),
      rowData(.data) |> as_tibble() |> slice(0)
    )
  
    # Check that we are not modifying a key column
    cols <- tidyselect::eval_select(  expr(c(...)), cols_data_frame  )

    # Check if column is row-wise of column-wise
    old_names = cols_data_frame[,cols] |> colnames()
    new_names = cols |> names()
    
    # If renaming col and row data at the same time, it is too complicate, so error
    if(
      old_names %in% colnames(colData(.data)) |> any() &
      old_names %in% colnames(rowData(.data)) |> any()
    )
      stop("tidySummarizedExperiment says: renaming columns from both colData and rowData at the same time is an unfeasible abstraction using dplyr. Please run two `rename` commands for sample-wise and feature-wise columns.")
    
    secial_columns = get_special_columns(
      
      # Decrease the size of the dataset
      .data[1:min(100, nrow(.data)), 1:min(20, ncol(.data))]
    ) |> 
      c(get_needed_columns(.data))
    
    tst =
        intersect(
            cols |> names(),
            secial_columns
        ) |>
        length() |>
        gt(0)

    # If column in view-only columns stop
    if (tst) {
        columns =
          secial_columns |>
            paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says: you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData). If you want to mutate a view-only column, make a copy and mutate that one."
        )
    }

    # Rename sample annotation
    if(old_names %in% colnames(colData(.data)) |> any())
      colData(.data) <- dplyr::rename(colData(.data) |> as.data.frame(), ...) |> DataFrame()

    # Rename gene annotation
    if(old_names %in% colnames(rowData(.data)) |> any())
      rowData(.data) <- dplyr::rename(rowData(.data) |> as.data.frame(), ...) |> DataFrame()
    
    .data
}


#' Group input by rows
#'
#'
#' See [this repository](https://github.com/jennybc/row-oriented-workflows)
#' for alternative ways to perform row-wise operations.
#'
#' `rowwise()` is used for the results of [do()] when you
#' create list-variables. It is also useful to support arbitrary
#' complex operations that need to be applied to each row.
#'
#' Currently, rowwise grouping only works with data frames. Its
#' main impact is to allow you to work with list-variables in
#' [summarise()] and [mutate()] without having to
#' use \code{[[1]]}. This makes `summarise()` on a rowwise tbl
#' effectively equivalent to [plyr::ldply()].
#'
#' @importFrom dplyr rowwise
#'
#' @param data Input data frame.
#' @param ... See dplyr::rowwise
#'
#' @return A `tbl`
#'
#'   A `tbl`
#'
#' @rdname dplyr-methods
#' @name rowwise
#'
#' @examples
#' 
#' print("To come...")
#' 
NULL

#' @export
rowwise.SummarizedExperiment <- function(data, ...) {
    message(data_frame_returned_message)

    data |>
        as_tibble() |>
        dplyr::rowwise()
}


#' Left join datasets
#'
#' @importFrom dplyr count
#' @importFrom dplyr left_join
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE,
#'   then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these
#'   suffixes will be added to the output to disambiguate them. Should be a
#'   character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tidySummarizedExperiment object
#'
#' @rdname dplyr-methods
#' @name left_join
#'
#'
#' @examples
#' data(pasilla)
#'
#' tt <- pasilla 
#' tt |> left_join(tt |> distinct(condition) |> mutate(new_column=1:2))
NULL

#' @export
left_join.SummarizedExperiment <- function(x, y, by=NULL, copy=FALSE, suffix=c(".x", ".y"),
    ...) {
  
  join_efficient_for_SE(x, y, by=by, copy=copy, suffix=suffix, dplyr::left_join, ...)
  
}

#' Inner join datasets
#'
#' @importFrom dplyr pull
#' @importFrom dplyr inner_join
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tidySummarizedExperiment object
#'
#' @examples
#' 
#' data(pasilla)
#' 
#' tt <- pasilla 
#' tt |> inner_join(tt |> distinct(condition) |> mutate(new_column=1:2) |> slice(1))
#'
#' @rdname dplyr-methods
#' @name inner_join
#'
NULL

#' @export
inner_join.SummarizedExperiment <- function(x, y, by=NULL, copy=FALSE, suffix=c(".x", ".y"), ...) {
  
  join_efficient_for_SE(x, y, by=by, copy=copy, suffix=suffix, dplyr::inner_join, ...)

}

#' Right join datasets
#'
#' @importFrom dplyr pull
#' @importFrom dplyr right_join
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE,
#'   then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these
#'   suffixes will be added to the output to disambiguate them. Should be a
#'   character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tidySummarizedExperiment object
#'
#' @examples
#' 
#' data(pasilla)
#' 
#' tt <- pasilla
#' tt |> right_join(tt |> distinct(condition) |> mutate(new_column=1:2) |> slice(1))
#'
#' @rdname dplyr-methods
#' @name right_join
#'
NULL

#' @export
right_join.SummarizedExperiment <- function(x, y, by=NULL, copy=FALSE, suffix=c(".x", ".y"),
    ...) {
  
  join_efficient_for_SE(x, y, by=by, copy=copy, suffix=suffix, dplyr::right_join, ...)
}


#' Full join datasets
#'
#' @importFrom dplyr pull
#' @importFrom dplyr full_join
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE,
#'   then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these
#'   suffixes will be added to the output to disambiguate them. Should be a
#'   character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tidySummarizedExperiment object
#'
#' @examples
#' 
#' data(pasilla)
#' 
#' tt <- pasilla
#' tt |> full_join(tibble::tibble(condition="treated", dose=10))
#'
#' @rdname dplyr-methods
#' @name full_join
#'
NULL

#' @export
full_join.SummarizedExperiment <- function(x, y, by=NULL, copy=FALSE, suffix=c(".x", ".y"),
    ...) {
  join_efficient_for_SE(x, y, by=by, copy=copy, suffix=suffix, dplyr::full_join, ...)
}

#' Subset rows using their positions
#'
#' @description
#' `slice()` lets you index rows by their (integer) locations. It allows you
#' to select, remove, and duplicate rows. It is accompanied by a number of
#' helpers for common use cases:
#'
#' * `slice_head()` and `slice_tail()` select the first or last rows.
#' * `slice_sample()` randomly selects rows.
#' * `slice_min()` and `slice_max()` select rows with highest or lowest values
#'   of a variable.
#'
#' If `.data` is a [grouped_df], the operation will be performed on each group,
#' so that (e.g.) `slice_head(df, n=5)` will select the first five rows in
#' each group.
#'
#' @importFrom dplyr slice
#'
#' @details
#' Slice does not work with relational databases because they have no
#' intrinsic notion of row order. If you want to perform the equivalent
#' operation, use [filter()] and [row_number()].
#'
#' @family single table verbs
#' @inheritParams arrange
#' @inheritParams filter
#'
#' @param .data A tidySummarizedExperiment object or any data frame
#' @param ... For `slice()`: <[`data-masking`][dplyr_data_masking]> Integer row
#'   values.
#'
#'   Provide either positive values to keep, or negative values to drop.
#'   The values provided must be either all positive or all negative.
#'   Indices beyond the number of rows in the input are silently ignored.
#'
#'   For `slice_helpers()`, these arguments are passed on to methods.
#'
#'
#'   If `n` is greater than the number of rows in the group (or `prop > 1`),
#'   the result will be silently truncated to the group size. If the
#'   `prop`ortion of a group size is not an integer, it is rounded down.
#' @return
#' An object of the same type as `.data`. The output has the following
#' properties:
#'
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * Data frame attributes are preserved.
#' @section Methods:
#' These function are **generic**s, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' Methods available in currently loaded packages:
#'
#' * `slice()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice")}.
#' * `slice_head()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_head")}.
#' * `slice_tail()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_tail")}.
#' * `slice_min()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_min")}.
#' * `slice_max()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_max")}.
#' * `slice_sample()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_sample")}.
#'
#' @rdname dplyr-methods
#' @name slice
#'
#' 
#' @examples
#'
#' data(pasilla)
#' 
#' pasilla |> slice(1)
NULL

#' @export
slice.SummarizedExperiment <- function(.data, ..., .preserve=FALSE) {
  
  slice_optimised(.data, ..., .preserve=.preserve)
    
    # .data |>
    #     as_tibble(skip_GRanges = T) |>
    #     dplyr::slice(..., .preserve=.preserve) %>%
    #     when(
    # 
    #         # If duplicated sample-feature pair returns tibble
    #         !is_not_duplicated(., .data) | 
    #         !is_rectangular(., .data) ~ {
    #             message(duplicated_cell_names)
    #             (.)
    #         },
    # 
    #         # Otherwise return updated tidySummarizedExperiment
    #         ~ update_SE_from_tibble(., .data)
    #     )
}

#' Subset columns using their names and types
#'
#' @description
#'
#' Select (and optionally rename) variables in a data frame, using a concise
#' mini-language that makes it easy to refer to variables based on their name
#' (e.g. `a:f` selects all columns from `a` on the left to `f` on the
#' right). You can also use predicate functions like is.numeric to select
#' variables based on their properties.
#'
#'
#' @importFrom dplyr select
#'
#' @inheritParams arrange
#'
#' @param .data A tidySummarizedExperiment object or any data frame
#' @param ... <[`tidy-select`][dplyr_tidy_select]> One or more unquoted
#'   expressions separated by commas. Variable names can be used as if they
#'   were positions in the data frame, so expressions like `x:y` can
#'   be used to select a range of variables.
#' @return
#' An object of the same type as `.data`. The output has the following
#' properties:
#'
#' * Rows are not affected.
#' * Output columns are a subset of input columns, potentially with a different
#'   order. Columns will be renamed if `new_name=old_name` form is used.
#' * Data frame attributes are preserved.
#' * Groups are maintained; you can't select off grouping variables.
#'
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("select")}.
#'
#' @examples
#'
#' data(pasilla)
#' 
#' pasilla |> select(.sample, .feature, counts)
#' 
#' @family single table verbs
#'
#' @rdname dplyr-methods
#' @name select
#'
NULL

#' @export
select.SummarizedExperiment <- function(.data, ...) {
   
  
  
  # colnames_col <- get_colnames_col(.data)
  # colnames_row <- get_rownames_col(.data)
  # colnames_assay = .data@assays@data |> names()
  

  
  # Deprecation of special column names
  if(is_sample_feature_deprecated_used(
    .data, 
    (enquos(..., .ignore_empty = "all") |> map(~ quo_name(.x)) |> unlist())
  )){
    .data= ping_old_special_column_into_metadata(.data)
  }
  
  # Warning if column names of assays do not overlap
  if(  check_if_assays_are_NOT_consistently_ordered(.data)  ){
    
    warning(
      "tidySummarizedExperiment says: the assays in your SummarizedExperiment have column names,
  but their order is not the same. Assays were internally reordered to be consistent with each other.
  To avoid unwanted behaviour it is highly reccomended to have assays with the same order of colnames and rownames"
    )
    
    
    # reorder assay colnames before printing
    # Rearrange if assays has colnames and rownames
    .data = order_assays_internally_to_be_consistent(.data)
    
  }

  # See if join done by sample, feature or both
  columns_query = 
    .data %>%
    .[1,1, drop=FALSE] |> 
    as_tibble() |> 
    select_helper(...) |> 
    colnames()
  
  row_data_tibble = 
    rowData(.data) |> 
    as_tibble(rownames = f_(.data)$name) 
  
  row_data_DF =
    row_data_tibble |> 
    select(one_of(columns_query), !!f_(.data)$symbol) |>
    suppressWarnings() %>% 
    data.frame(row.names=pull(., !!f_(.data)$symbol)) |>
    select(-!!f_(.data)$symbol) |>
    DataFrame()
  
  col_data_tibble = 
    colData(.data) |> 
    as_tibble(rownames = s_(.data)$name) 
  
  col_data_DF = 
    col_data_tibble |>  
    select(one_of(columns_query), !!s_(.data)$symbol) |>
    suppressWarnings() %>% 
    data.frame(row.names=pull(., !!s_(.data)$symbol)) |>
    select(-!!s_(.data)$symbol) |>
    DataFrame()

  
  count_data = 
    assays(.data)@listData %>% .[names(assays(.data)@listData) %in% columns_query]
  
  # If it's just col data
  if(ncol(row_data_DF) == 0 & !f_(.data)$name %in% columns_query & length(count_data) == 0 & ( ncol(col_data_DF) > 0 | s_(.data)$name %in% columns_query) ) {
    message("tidySummarizedExperiment says: Key columns are missing. A data frame is returned for independent data analysis.")
    
    col_data_tibble |> 
      select_helper(...) |> 
      slice(rep(1:n(), each=nrow(!!.data) ))
    
  } 
  
  # If it's just row data
  else if(ncol(col_data_DF) == 0 & !s_(.data)$name %in% columns_query & length(count_data) == 0 & ( ncol(row_data_DF) > 0 | f_(.data)$name %in% columns_query)){
    message("tidySummarizedExperiment says: Key columns are missing. A data frame is returned for independent data analysis.")
    
    row_data_tibble |> 
      select_helper(...) |> 
      slice(rep(1:n(), ncol(!!.data)  ))
    
  }
  
  else if(!all(c(get_needed_columns(.data)) %in% columns_query)){
    if(ncol(.data)>100) message("tidySummarizedExperiment says: You are doing a complex selection both sample-wise and feature-wise. In the latter case, for efficiency (until further development), it is better to separate your selects sample-wise OR feature-wise.")
    message("tidySummarizedExperiment says: Key columns are missing. A data frame is returned for independent data analysis.")
    
    .data |>
      as_tibble(skip_GRanges = TRUE) |>
      select_helper(...) 
  }
  
  else {
    rowData(.data) = row_data_DF
    colData(.data) = col_data_DF
    assays(.data) = count_data
    .data
    
  } 
   
}


#' Sample n rows from a table
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("superseded")}
#' `sample_n()` and `sample_frac()` have been superseded in favour of
#' [slice_sample()]. While they will not be deprecated in the near future,
#' retirement means that we will only perform critical bug fixes, so we recommend
#' moving to the newer alternative.
#'
#' These functions were superseded because we realised it was more convenient to
#' have two mutually exclusive arguments to one function, rather than two
#' separate functions. This also made it to clean up a few other smaller
#' design issues with `sample_n()`/`sample_frac`:
#'
#' * The connection to `slice()` was not obvious.
#' * The name of the first argument, `tbl`, is inconsistent with other
#'   single table verbs which use `.data`.
#' * The `size` argument uses tidy evaluation, which is surprising and
#'   undocumented.
#' * It was easier to remove the deprecated `.env` argument.
#' * `...` was in a suboptimal position.
#'
#' @importFrom dplyr sample_n
#'
#' @keywords internal
#' @param tbl A data.frame.
#' @param size <[`tidy-select`][dplyr_tidy_select]>
#'   For `sample_n()`, the number of rows to select.
#'   For `sample_frac()`, the fraction of rows to select.
#'   If `tbl` is grouped, `size` applies to each group.
#' @param replace Sample with or without replacement?
#' @param weight <[`tidy-select`][dplyr_tidy_select]> Sampling weights.
#'   This must evaluate to a vector of non-negative numbers the same length as
#'   the input. Weights are automatically standardised to sum to 1.
#' @param .env DEPRECATED.
#' @param ... ignored
#' @examples
#'
#' data(pasilla)
#' 
#' pasilla |>
#'     
#'     sample_n(50)
#' pasilla |>
#'     
#'     sample_frac(0.1)
#' @return A tidySummarizedExperiment object
#'
#' @rdname dplyr-methods
#' @name sample_n
#'
NULL

#' @export
sample_n.SummarizedExperiment <- function(tbl, size, replace=FALSE,
    weight=NULL, .env=NULL, ...) {
    lifecycle::signal_superseded("1.0.0", "sample_n()", "slice_sample()")

    message(data_frame_returned_message)

    tbl |>
        as_tibble() |>
        dplyr::sample_n(size, replace=replace, weight=weight, .env=.env, ...)
}

#' @importFrom dplyr sample_frac
#'
#' @rdname dplyr-methods
#' @name sample_frac
#'
NULL

#' @export
sample_frac.SummarizedExperiment <- function(tbl, size=1, replace=FALSE,
    weight=NULL, .env=NULL, ...) {
    lifecycle::signal_superseded("1.0.0", "sample_frac()", "slice_sample()")

    message(data_frame_returned_message)

    tbl |>
        as_tibble() |>
        dplyr::sample_frac(size, replace=replace, weight=weight, .env=.env, ...)
}


#' Count observations by group
#'
#' @importFrom dplyr count
#'
#' @description
#' `count()` lets you quickly count the unique values of one or more variables:
#' `df |> count(a, b)` is roughly equivalent to
#' `df |> group_by(a, b) |> summarise(n=n())`.
#' `count()` is paired with `tally()`, a lower-level helper that is equivalent
#' to `df |> summarise(n=n())`. Supply `wt` to perform weighted counts,
#' switching the summary from `n=n()` to `n=sum(wt)`.
#'
#' `add_count()` are `add_tally()` are equivalents to `count()` and `tally()`
#' but use `mutate()` instead of `summarise()` so that they add a new column
#' with group-wise counts.
#'
#' @param x A data frame, data frame extension (e.g. a tibble), or a
#'   lazy data frame (e.g. from dbplyr or dtplyr).
#' @param ... <[`data-masking`][dplyr_data_masking]> Variables to group by.
#' @param wt <[`data-masking`][dplyr_data_masking]> Frequency weights.
#'   Can be `NULL` or a variable:
#'
#'   * If `NULL` (the default), counts the number of rows in each group.
#'   * If a variable, computes `sum(wt)` for each group.
#' @param sort If `TRUE`, will show the largest groups at the top.
#' @param name The name of the new column in the output.
#'
#'   If omitted, it will default to `n`. If there's already a column called `n`,
#'   it will error, and require you to specify the name.
#' @param .drop For `count()`: if `FALSE` will include counts for empty groups
#'   (i.e. for levels of factors that don't exist in the data). Deprecated in
#'   `add_count()` since it didn't actually affect the output.
#' @return
#' An object of the same type as `.data`. `count()` and `add_count()`
#' group transiently, so the output has the same groups as the input.
#' 
#'
#' @rdname dplyr-methods
#' @name count
#'
#' @examples
#'
#' se |> count(dex)
#' 
NULL

#' @export
count.SummarizedExperiment <- function(x, ..., wt=NULL, sort=FALSE, name=NULL, .drop=group_by_drop_default(x)) {
    message(data_frame_returned_message)

  # Deprecation of special column names
  if(is_sample_feature_deprecated_used(
    x, 
    (enquos(..., .ignore_empty = "all") |> map(~ quo_name(.x)) |> unlist())
  )){
    x= ping_old_special_column_into_metadata(x)
  }
  
  
  # If Ranges column not in query perform fast as_tibble
  skip_GRanges = 
    get_GRanges_colnames() %in% 
    (enquos(..., .ignore_empty = "all") |> map(~ quo_name(.x)) |> unlist()) |>
    not()
  
    x |>
        as_tibble(skip_GRanges = skip_GRanges) |>
        dplyr::count(..., wt=!!enquo(wt), sort=sort, name=name, .drop=.drop)
}



#' Extract a single column
#'
#' `pull()` is similar to `$`. It's mostly useful because it looks a little
#' nicer in pipes, it also works with remote data frames, and it can optionally
#' name the output.
#'
#' @importFrom ellipsis check_dots_used
#' @importFrom dplyr pull
#'
#' @inheritParams arrange
#' @inheritParams tidyselect::vars_pull
#'
#' @param .data A tidySummarizedExperiment object or any data frame
#' @param name An optional parameter that specifies the column to be used
#'   as names for a named vector. Specified in a similar manner as \code{var}.
#' @param ... For use by methods.
#' @return A vector the same size as `.data`.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("pull")}.
#'
#' @rdname dplyr-methods
#' @name pull
#'
#' @examples
#'
#' data(pasilla)
#' pasilla |>
#'     
#'     pull(feature)
NULL

#' @export
pull.SummarizedExperiment <- function(.data, var=-1, name=NULL, ...) {
  
  # Fix NOTEs
  . = NULL
  
    var <- enquo(var)
    name <- enquo(name)

    quo_name_name = name %>% when(quo_is_null(.) ~ NULL, quo_name(name))
    
    # Deprecation of special column names
    if(is_sample_feature_deprecated_used(
      .data, 
      quo_name(var)
    )){
      .data= ping_old_special_column_into_metadata(.data)
    }
    
    # If Ranges column not in query perform fast as_tibble
    skip_GRanges = 
      get_GRanges_colnames() %in% 
      quo_name(var) |>
      not()
    
    # Subset column annotation
    if(all(c(quo_names(var), quo_name_name) %in% colnames(colData(.data))))
      return( colData(.data)[,quo_names(var)] %>% .[rep(1:length(.), each=nrow(.data) )])
    
    # Subset row annotation
    if(all(c(quo_names(var), quo_name_name) %in% colnames(rowData(.data))))
      return( colData(.data)[,quo_names(var)] %>% .[rep(1:length(.), ncol(.data) )])
    
    # This returns a vector column wise. With the first sample and all features, 
    # second sample and all features, etc..
    if(all(c(quo_names(var), quo_name_name) %in% names(.data@assays@data))){
      
      # Warning if column names of assays do not overlap
      if( check_if_assays_are_NOT_consistently_ordered(.data) ){
        
        warning( 
          "tidySummarizedExperiment says: the assays in your SummarizedExperiment have column names, 
  but their order is not the same. Pulling assays can return data in a order you don't expect. 
  To avoid unwanted behaviour it is highly reccomended to have assays with the same order of colnames and rownames" 
        )
        
        # reorder assay colnames before printing
        # Rearrange if assays has colnames and rownames
        .data = order_assays_internally_to_be_consistent(.data)
        
      }
       
      
      return(assay(.data, quo_names(var)) |> as.matrix() |> as.vector()) 
      
    }
    
    # Subset rowranges
    if(all(c(quo_names(var), quo_name_name) %in% colnames(as.data.frame(rowRanges(.data)))))
      return( as.data.frame(rowRanges(.data))[,quo_names(var)] %>% .[rep(1:length(.), ncol(.data) )])
    
    # Otherwise (SHOULD NOT HAPPEN) use the long general procedure
    colData(.data) = colData(.data)[,colnames(colData(.data)) %in% c(quo_names(var), quo_name_name), drop=FALSE ]
    rowData(.data) = rowData(.data)[,colnames(rowData(.data)) %in% c(quo_names(var), quo_name_name), drop=FALSE ]
  
    .data |>
        as_tibble(skip_GRanges = skip_GRanges) |>
        dplyr::pull(var=!!var, name=!!name, ...)
}
