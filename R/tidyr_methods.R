#' unnest
#'
#' @importFrom tidyr unnest
#'
#' @param data A tbl. (See tidyr)
#' @param cols <[`tidy-select`][tidyr_tidy_select]> Columns to unnest.
#'   If you `unnest()` multiple columns, parallel entries must be of
#'   compatible sizes, i.e. they're either equal or length 1 (following the
#'   standard tidyverse recycling rules).
#' @param ... <[`tidy-select`][tidyr_tidy_select]> Columns to nest, specified
#'   using name-variable pairs of the form `new_col=c(col1, col2, col3)`.
#'   The right hand side can be any valid tidy select expression.
#'
#'   \Sexpr[results=rd, stage=render]{lifecycle::badge("deprecated")}:
#'   previously you could write `df %>% nest(x, y, z)` and `df %>%
#'   unnest(x, y, z)`. Convert to `df %>% nest(data=c(x, y, z))`.
#'   and `df %>% unnest(c(x, y, z))`.
#'
#'   If you previously created new variable in `unnest()` you'll now need to
#'   do it explicitly with `mutate()`. Convert `df %>% unnest(y=fun(x, y, z))`
#'   to `df %>% mutate(y=fun(x, y, z)) %>% unnest(y)`.
#' @param names_sep If `NULL`, the default, the names will be left
#'   as is. In `nest()`, inner names will come from the former outer names;
#'   in `unnest()`, the new outer names will come from the inner names.
#'
#'   If a string, the inner and outer names will be used together. In `nest()`,
#'   the names of the new outer columns will be formed by pasting together the
#'   outer and the inner column names, separated by `names_sep`. In `unnest()`,
#'   the new inner names will have the outer names (+ `names_sep`) automatically
#'   stripped. This makes `names_sep` roughly symmetric between nesting and unnesting.
#' @param keep_empty See tidyr::unnest
#' @param names_repair See tidyr::unnest
#' @param ptype See tidyr::unnest
#' @param .drop See tidyr::unnest
#' @param .id tidyr::unnest
#' @param .sep tidyr::unnest
#' @param .preserve See tidyr::unnest
#'
#'
#' @return A tidySummarizedExperiment objector a tibble depending on input
#'
#' @examples
#'
#' tidySummarizedExperiment::pasilla %>%
#'
#'     nest(data=-condition) %>%
#'     unnest(data)
#'
#' @rdname tidyr-methods
#' @name unnest
#'
#' @export
NULL


#' @importFrom rlang quo_name
#' @importFrom purrr imap
#'
#' @export
unnest.tidySummarizedExperiment_nested <-
    function(data, cols, ..., keep_empty=FALSE, ptype=NULL, names_sep=NULL, names_repair="check_unique", .drop, .id, .sep, .preserve) {


    # Need this otherwise crashes map
    .data_ <- data

    cols <- enquo(cols)



    .data_ %>%
        when(

            # If my only column to unnest is tidySummarizedExperiment
            pull(., !!cols) %>%
                .[[1]] %>%
                class() %>%
                as.character() %>% 
                eq("SummarizedExperiment") %>%
                any() ~ {

                  se = pull(., !!cols) %>% .[[1]] 
                  
                  # Mark if columns belong to feature or sample
                  my_unnested_tibble =
                    mutate(., !!cols := map(!!cols, ~ as_tibble(.x))) %>%
                    select(-suppressWarnings( one_of(s_(se)$name, f_(se)$name))) %>%
                    unnest(!!cols)

                  # Get which column is relative to feature or sample
                  sample_columns = my_unnested_tibble %>% get_subset_columns(!!s_(se)$symbol)
                  transcript_columns = my_unnested_tibble %>% get_subset_columns(!!f_(se)$symbol)
                  source_column =
                    c(
                      rep(s_(se)$name, length(sample_columns)) %>% setNames(sample_columns),
                      rep(f_(se)$name, length(transcript_columns)) %>% setNames(transcript_columns)
                    )

                  # Do my trick to unnest
                  mutate(., !!cols := imap(
                    !!cols, ~ .x %>%
                      bind_cols_internal(

                        # Attach back the columns used for nesting
                        .data_ %>%
                          select(-!!cols, -suppressWarnings( one_of(s_(se)$name, f_(se)$name))) %>%
                          slice(rep(.y, ncol(.x) * nrow(.x))),

                        # Column sample-wise or feature-wise
                        column_belonging =
                          source_column[
                            .data_ %>%
                              select(-!!cols, -suppressWarnings( one_of(s_(se)$name, f_(se)$name))) %>%
                              colnames()
                          ]
                      )
                  )) %>%
                    pull(!!cols) %>%

                    # See if split by feature or sample
                    when(
                      is_split_by_sample(.) & is_split_by_transcript(.) ~ stop("tidySummarizedExperiment says: for the moment nesting both by sample- and feature-wise information is not possible. Please ask this feature to github/stemangiola/tidySummarizedExperiment"),
                      ~ reduce(., bind_rows)
                    )
                },

            # Else do normal stuff
            ~ (.) %>%
                drop_class("tidySummarizedExperiment_nested") %>%
                tidyr::unnest(!!cols, ..., keep_empty=keep_empty, ptype=ptype, names_sep=names_sep, names_repair=names_repair) %>%
                add_class("tidySummarizedExperiment_nested")
        )
}

#' nest
#'
#' @importFrom tidyr nest
#'
#' @param .data A tbl. (See tidyr)
#' @param ... Name-variable pairs of the form new_col=c(col1, col2, col3) (See tidyr)
#' @param .names_sep See ?tidyr::nest
#'
#' @return A tidySummarizedExperiment objector a tibble depending on input
#'
#' @examples
#'
#' tidySummarizedExperiment::pasilla %>%
#'
#'     nest(data=-condition)
#'
#' @rdname tidyr-methods
#' @name nest
#'
#' @export
NULL

#' @importFrom rlang enquos
#' @importFrom rlang :=
#' @importFrom purrr when
#' @importFrom purrr pmap
#'
#' @export
nest.SummarizedExperiment <- function(.data, ..., .names_sep = NULL) {
    cols <- enquos(...)
    col_name_data <- names(cols)

    # Deprecation of special column names
    if(is_sample_feature_deprecated_used(
      .data, 
      (enquos(..., .ignore_empty = "all") %>% map(~ quo_name(.x)) %>% unlist)
    )){
      .data= ping_old_special_column_into_metadata(.data)
    }
    
    my_data__ <- .data 
    
    
    my_data__nested =
      my_data__ %>%

        # This is needed otherwise nest goes into loop and fails
        as_tibble() %>%
        tidyr::nest(...) %>%

        # Check that sample or feature are in the nesting
        {
            if(c(s_(.data)$name, f_(.data)$name) %>% intersect(colnames(.)) %>% length() %>% `>` (0))
                stop("tidySummarizedExperiment says: You cannot have the columns sample or feature among the nesting")
            (.)
        }

    sample_name = s_(my_data__)$name
    feature_name = f_(my_data__)$name
    sample_symbol = s_(my_data__)$symbol
    feature_symbol = f_(my_data__)$symbol
    
    my_data__nested %>%

        mutate(
            !!as.symbol(col_name_data) := pmap(

              # Add sample feature to map if nesting by those
                list(!!as.symbol(col_name_data)) %>%

                  # Check if nested by sample
                  when(sample_name %in% colnames(my_data__nested) ~ c(., list(!!sample_symbol)), ~ (.)) %>%

                  # Check if nested by feature
                  when(feature_name %in% colnames(my_data__nested) ~ c(., list(!!feature_symbol)), ~ (.)) , ~ {

                    # Check if nested by sample
                    if(sample_name %in% colnames(my_data__nested)) { my_samples=..2 }
                    else {my_samples=pull(..1,!!sample_symbol)}

                    # Check if nested by feature
                    if(sample_name %in% colnames(my_data__nested) & feature_name %in% colnames(my_data__nested)) {my_transcripts=..3}
                    else if(feature_name %in% colnames(my_data__nested)) my_transcripts=..2
                    else my_transcripts=pull(..1,!!feature_symbol)

                  my_data__ %>%

                    # Subset cells
                    filter(!!sample_symbol %in% my_samples & !!feature_symbol %in% my_transcripts) %>%

                    # Subset columns
                    select(colnames(..1) %>% c(sample_name, feature_name) %>% unique)
                }
            )
        ) %>%

        # Coerce to tidySummarizedExperiment_nested for unnesting
        add_class("tidySummarizedExperiment_nested")
}

#' Extract a character column into multiple columns using regular
#' expression groups
#'
#' Given a regular expression with capturing groups, `extract()` turns
#' each group into a new column. If the groups don't match, or the input
#' is NA, the output will be NA.
#'
#' @importFrom tidyr extract
#'
#' @param data A tidySummarizedExperiment object
#' @param col Column name or position. This is passed to
#'   [tidyselect::vars_pull()].
#'
#'   This argument is passed by expression and supports
#'   [quasiquotation][rlang::quasiquotation] (you can unquote column
#'   names or column positions).
#' @param into Names of new variables to create as character vector.
#'    Use `NA` to omit the variable in the output.
#' @param regex a regular expression used to extract the desired values.
#'   There should be one group (defined by `()`) for each element of `into`.
#' @param remove If `TRUE`, remove input column from output data frame.
#' @param convert If `TRUE`, will run [type.convert()] with
#'   `as.is=TRUE` on new columns. This is useful if the component
#'   columns are integer, numeric or logical.
#'
#'   NB: this will cause string `"NA"`s to be converted to `NA`s.
#' @param ... Additional arguments passed on to methods.
#' @seealso [separate()] to split up by a separator.
#'
#' @rdname tidyr-methods
#' @name extract
#'
#' @export
#' @examples
#'
#' tidySummarizedExperiment::pasilla %>%
#'
#'     extract(type, into="sequencing", regex="([a-z]*)_end", convert=TRUE)
#' @return A tidySummarizedExperiment objector a tibble depending on input
#'
#' @importFrom tidyr extract
#'
#' @export
NULL

#' @importFrom rlang enquo
#' @export
extract.SummarizedExperiment <- function(data, col, into, regex="([[:alnum:]]+)", remove=TRUE,
    convert=FALSE, ...) {
    col <- enquo(col)

    # Deprecation of special column names
    if(is_sample_feature_deprecated_used(
      data, 
      c(quo_name(col), into)
    )){
      data= ping_old_special_column_into_metadata(data)
    }
    
    tst =
        intersect(
            into %>% quo_names(),
            get_special_columns(data) %>% c(get_needed_columns(data))
        ) %>%
        length() %>%
        gt(0) &
        remove


    if (tst) {
        columns =
            get_special_columns(data) %>%
            c(get_needed_columns(data)) %>%
            paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says: you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData). If you want to mutate a view-only column, make a copy and mutate that one."
        )
    }

    data %>%
        as_tibble(skip_GRanges = T) %>%
        tidyr::extract(col=!!col, into=into, regex=regex, remove=remove, convert=convert, ...) %>%
        update_SE_from_tibble(data)
}

#' Pivot data from wide to long
#'
#' @description
#'
#' `pivot_longer()` "lengthens" data, increasing the number of rows and
#' decreasing the number of columns. The inverse transformation is
#' `pivot_wider()`
#'
#' Learn more in `vignette("pivot")`.
#'
#' @importFrom ellipsis check_dots_used
#' @importFrom tidyr pivot_longer
#'
#' @details
#' `pivot_longer()` is an updated approach to [gather()], designed to be both
#' simpler to use and to handle more use cases. We recommend you use
#' `pivot_longer()` for new code; `gather()` isn't going away but is no longer
#' under active development.
#'
#' @param data A data frame to pivot.
#' @param cols <[`tidy-select`][tidyr_tidy_select]> Columns to pivot into
#'   longer format.
#' @param names_to A string specifying the name of the column to create
#'   from the data stored in the column names of `data`.
#'
#'   Can be a character vector, creating multiple columns, if `names_sep`
#'   or `names_pattern` is provided. In this case, there are two special
#'   values you can take advantage of:
#'
#'   * `NA` will discard that component of the name.
#'   * `.value` indicates that component of the name defines the name of the
#'     column containing the cell values, overriding `values_to`.
#' @param names_prefix A regular expression used to remove matching text
#'   from the start of each variable name.
#' @param names_sep,names_pattern If `names_to` contains multiple values,
#'   these arguments control how the column name is broken up.
#'
#'   `names_sep` takes the same specification as [separate()], and can either
#'   be a numeric vector (specifying positions to break on), or a single string
#'   (specifying a regular expression to split on).
#'
#'   `names_pattern` takes the same specification as [extract()], a regular
#'   expression containing matching groups (`()`).
#'
#'   If these arguments do not give you enough control, use
#'   `pivot_longer_spec()` to create a spec object and process manually as
#'   needed.
#' @param names_repair What happens if the output has invalid column names?
#'   The default, `"check_unique"` is to error if the columns are duplicated.
#'   Use `"minimal"` to allow duplicates in the output, or `"unique"` to
#'   de-duplicated by adding numeric suffixes. See [vctrs::vec_as_names()]
#'   for more options.
#' @param values_to A string specifying the name of the column to create
#'   from the data stored in cell values. If `names_to` is a character
#'   containing the special `.value` sentinel, this value will be ignored,
#'   and the name of the value column will be derived from part of the
#'   existing column names.
#' @param values_drop_na If `TRUE`, will drop rows that contain only `NA`s
#'   in the `value_to` column. This effectively converts explicit missing values
#'   to implicit missing values, and should generally be used only when missing
#'   values in `data` were created by its structure.
#' @param names_transform,values_transform A list of column name-function pairs.
#'   Use these arguments if you need to change the type of specific columns.
#'   For example, `names_transform=list(week=as.integer)` would convert
#'   a character week variable to an integer.
#' @param names_ptypes,values_ptypes A list of column name-prototype pairs.
#'   A prototype (or ptype for short) is a zero-length vector (like `integer()`
#'   or `numeric()`) that defines the type, class, and attributes of a vector.
#'   Use these arguments to confirm that the created columns are the types that
#'   you expect.
#'
#'   If not specified, the type of the columns generated from `names_to` will
#'   be character, and the type of the variables generated from `values_to`
#'   will be the common type of the input columns used to generate them.
#' @param ... Additional arguments passed on to methods.
#'
#' @return A tidySummarizedExperiment objector a tibble depending on input
#'
#' @rdname tidyr-methods
#' @name pivot_longer
#'
#' @export
#' @examples
#' # See vignette("pivot") for examples and explanation
#'
#' library(dplyr)
#' tidySummarizedExperiment::pasilla %>%
#'
#'     pivot_longer(c(condition, type), names_to="name", values_to="value")
NULL

#' @export
pivot_longer.SummarizedExperiment <- function(data,
    cols,
    names_to="name",
    names_prefix=NULL,
    names_sep=NULL,
    names_pattern=NULL,
    names_ptypes=list(),
    names_transform=list(),
    names_repair="check_unique",
    values_to="value",
    values_drop_na=FALSE,
    values_ptypes=list(),
    values_transform=list(),
    ...) {
    cols <- enquo(cols)

    message(data_frame_returned_message)

    
    # Deprecation of special column names
    if(is_sample_feature_deprecated_used(
      data, 
      c(quo_names(cols))
    )){
      data= ping_old_special_column_into_metadata(data)
    }
    
    data %>%
        as_tibble(skip_GRanges = T) %>%
        tidyr::pivot_longer(!!cols,
            names_to=names_to,
            names_prefix=names_prefix,
            names_sep=names_sep,
            names_pattern=names_pattern,
            names_ptypes=names_ptypes,
            names_transform=names_transform,
            names_repair=names_repair,
            values_to=values_to,
            values_drop_na=values_drop_na,
            values_ptypes=values_ptypes,
            values_transform=values_transform,
            ...
        )
}


#' Pivot data from long to wide
#
#' @description
#' `pivot_wider()` "widens" data, increasing the number of columns and
#' decreasing the number of rows. The inverse transformation is
#' [pivot_longer()].
#'
#' Learn more in `vignette("pivot")`.
#'
#' @details
#' `pivot_wider()` is an updated approach to [spread()], designed to be both
#' simpler to use and to handle more use cases. We recommend you use
#' `pivot_wider()` for new code; `spread()` isn't going away but is no longer
#' under active development.
#'
#' @seealso [pivot_wider_spec()] to pivot "by hand" with a data frame that
#'   defines a pivotting specification.
#' @inheritParams pivot_longer
#' @param id_cols <[`tidy-select`][tidyr_tidy_select]> A set of columns that
#'   uniquely identifies each observation. Defaults to all columns in `data`
#'   except for the columns specified in `names_from` and `values_from`.
#'   Typically used when you have redundant variables, i.e. variables whose
#'   values are perfectly correlated with existing variables.
#' @param names_from,values_from <[`tidy-select`][tidyr_tidy_select]> A pair of
#'   arguments describing which column (or columns) to get the name of the
#'   output column (`names_from`), and which column (or columns) to get the
#'   cell values from (`values_from`).
#'
#'   If `values_from` contains multiple values, the value will be added to the
#'   front of the output column.
#' @param names_sep If `names_from` or `values_from` contains multiple
#'   variables, this will be used to join their values together into a single
#'   string to use as a column name.
#' @param names_prefix String added to the start of every variable name. This is
#'   particularly useful if `names_from` is a numeric vector and you want to
#'   create syntactic variable names.
#' @param names_glue Instead of `names_sep` and `names_prefix`, you can supply
#'   a glue specification that uses the `names_from` columns (and special
#'   `.value`) to create custom column names.
#' @param names_sort Should the column names be sorted? If `FALSE`, the default,
#'   column names are ordered by first appearance.
#' @param values_fill Optionally, a (scalar) value that specifies what each
#'   `value` should be filled in with when missing.
#'
#'   This can be a named list if you want to apply different aggregations
#'   to different value columns.
#' @param values_fn Optionally, a function applied to the `value` in each cell
#'   in the output. You will typically use this when the combination of
#'   `id_cols` and `value` column does not uniquely identify an observation.
#'
#'   This can be a named list if you want to apply different aggregations
#'   to different value columns.
#' @param ... Additional arguments passed on to methods.
#'
#' @importFrom tidyr pivot_wider
#' @rdname tidyr-methods
#' @name pivot_wider
#'
#'
#' @export
#' @examples
#' # See vignette("pivot") for examples and explanation
#'
#' library(dplyr)
#' tidySummarizedExperiment::pasilla %>%
#'
#'     pivot_wider(names_from=feature, values_from=counts)
NULL

#' @export
pivot_wider.SummarizedExperiment <- function(data,
                                   id_cols = NULL,
                                   names_from = name,
                                   names_prefix = "",
                                   names_sep = "_",
                                   names_glue = NULL,
                                   names_sort = FALSE,
                                   names_repair = "check_unique",
                                   values_from = value,
                                   values_fill = NULL,
                                   values_fn = NULL,
                                   ...
) {
  id_cols <- enquo(id_cols)
  name = enquo(names_from)
  value = enquo(values_from)

  message(data_frame_returned_message)

  # Deprecation of special column names
  if(is_sample_feature_deprecated_used(
    data, 
    c(quo_names(id_cols), quo_names(names_from))
  )){
    data= ping_old_special_column_into_metadata(data)
  }
  
  data %>%
    as_tibble(skip_GRanges = T) %>%
    tidyr::pivot_wider( id_cols = !!id_cols,
                        names_from = !!name,
                        names_prefix = names_prefix,
                        names_sep = names_sep,
                        names_glue = names_glue,
                        names_sort = names_sort,
                        names_repair = names_repair,
                        values_from = !!value,
                        values_fill = values_fill,
                        values_fn = values_fn,
                        ...
    )
}



#' Unite multiple columns into one by pasting strings together
#'
#' Convenience function to paste together multiple columns into one.
#'
#' @importFrom ellipsis check_dots_unnamed
#' @importFrom tidyr unite
#'
#' @param data A data frame.
#' @param col The name of the new column, as a string or symbol.
#'
#'   This argument is passed by expression and supports
#'   [quasiquotation][rlang::quasiquotation] (you can unquote strings
#'   and symbols). The name is captured from the expression with
#'   [rlang::ensym()] (note that this kind of interface where
#'   symbols do not represent actual objects is now discouraged in the
#'   tidyverse; we support it here for backward compatibility).
#' @param ... <[`tidy-select`][tidyr_tidy_select]> Columns to unite
#' @param sep Separator to use between values.
#' @param na.rm If `TRUE`, missing values will be remove prior to uniting
#'   each value.
#' @param remove If `TRUE`, remove input columns from output data frame.
#' @seealso [separate()], the complement.
#'
#' @return A tidySummarizedExperiment objector a tibble depending on input
#'
#' @rdname tidyr-methods
#' @name unite
#'
#' @export
#' @examples
#'
#' tidySummarizedExperiment::pasilla %>%
#'
#'     unite("group", c(condition, type))
NULL

#' @export
unite.SummarizedExperiment <- function(data, col, ..., sep="_", remove=TRUE, na.rm=FALSE) {

    # Check that we are not modifying a key column
    cols <- enquo(col)

    # Deprecation of special column names
    # Deprecation of special column names
    if(is_sample_feature_deprecated_used(
      data, 
      (enquos(..., .ignore_empty = "all") %>% map(~ quo_name(.x)) %>% unlist)
    )){
      data= ping_old_special_column_into_metadata(data)
    }
    
    tst =
        intersect(
            cols %>% quo_names(),
            get_special_columns(data) %>% c(get_needed_columns(data))
        ) %>%
        length() %>%
        gt(0) &
        remove


    if (tst) {
        columns =
            get_special_columns(data) %>%
            c(get_needed_columns(data)) %>%
            paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says: you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData). If you want to mutate a view-only column, make a copy and mutate that one."
        )
    }



    data %>%
        as_tibble(skip_GRanges = T) %>%
        tidyr::unite(!!cols, ..., sep=sep, remove=remove, na.rm=na.rm) %>%
        update_SE_from_tibble(data)
}

#' Separate a character column into multiple columns with a regular
#' expression or numeric locations
#'
#' Given either a regular expression or a vector of character positions,
#' `separate()` turns a single character column into multiple columns.
#'
#' @importFrom ellipsis check_dots_used
#' @importFrom tidyr separate
#'
#' @inheritParams extract
#' @param sep Separator between columns.
#'
#'   If character, `sep` is interpreted as a regular expression. The default
#'   value is a regular expression that matches any sequence of
#'   non-alphanumeric values.
#'
#'   If numeric, `sep` is interpreted as character positions to split at. Positive
#'   values start at 1 at the far-left of the string; negative value start at -1 at
#'   the far-right of the string. The length of `sep` should be one less than
#'   `into`.
#' @param extra If `sep` is a character vector, this controls what
#'   happens when there are too many pieces. There are three valid options:
#'
#'   * "warn" (the default): emit a warning and drop extra values.
#'   * "drop": drop any extra values without a warning.
#'   * "merge": only splits at most `length(into)` times
#' @param fill If `sep` is a character vector, this controls what
#'   happens when there are not enough pieces. There are three valid options:
#'
#'   * "warn" (the default): emit a warning and fill from the right
#'   * "right": fill with missing values on the right
#'   * "left": fill with missing values on the left
#' @seealso [unite()], the complement, [extract()] which uses regular
#'   expression capturing groups.
#'
#' @return A tidySummarizedExperiment objector a tibble depending on input
#'
#' @rdname tidyr-methods
#' @name separate
#'
#' @export
#' @examples
#'
#' un <- tidySummarizedExperiment::pasilla %>%
#'
#'     unite("group", c(condition, type))
#' un %>% separate(col=group, into=c("condition", "type"))
NULL

#' @export
separate.SummarizedExperiment <- function(data, col, into, sep="[^[:alnum:]]+", remove=TRUE,
    convert=FALSE, extra="warn", fill="warn", ...) {

    # Check that we are not modifying a key column
    cols <- enquo(col)

    # Deprecation of special column names
    if(is_sample_feature_deprecated_used(
      data, 
      c(quo_names(cols))
    )){
      data= ping_old_special_column_into_metadata(data)
    }
    
    tst =
        intersect(
            cols %>% quo_names(),
            get_special_columns(data) %>% c(get_needed_columns(data))
        ) %>%
        length() %>%
        gt(0) &
        remove


    if (tst) {
        columns =
            get_special_columns(data) %>%
            c(get_needed_columns(data)) %>%
            paste(collapse=", ")
        stop(
            "tidySummarizedExperiment says: you are trying to rename a column that is view only",
            columns,
            "(it is not present in the colData). If you want to mutate a view-only column, make a copy and mutate that one."
        )
    }


    data %>%
        as_tibble(skip_GRanges = T) %>%
        tidyr::separate(!!cols, into=into, sep=sep, remove=remove, convert=convert, extra=extra, fill=fill, ...) %>%
        update_SE_from_tibble(data)
}
