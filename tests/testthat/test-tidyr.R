context("tidyr test")

library(tidySE)

tt <-
    pasilla %>%
    tidy() %>%
    tidySE::mutate(col2 = "other_col")

test_that("nest_unnest", {
    col_names <- colnames(tt@colData) %>% c("sample")
    library(magrittr)

    y <- tibble::tibble(
        sample = c(
            "untrt1", "untrt2", "untrt3", "untrt4",
            "trt1", "trt2", "trt3"
        ),
        counts = c(0L, 0L, 0L, 0L, 0L, 0L, 1L)
    )

    x <- tt %>%
        nest(data = -condition) %>%
        unnest(data) %>%
        head(n = 1) %>%
        select(sample, counts)


    expect_equal(x, y)
})

test_that("unite separate", {
    un <- tt %>% unite("new_col", c(condition, col2), sep = ":")

    expect_equal(un %>% tidySE::select(new_col) %>% slice(1) %>% pull(new_col), "untreated:other_col")

    se <- un %>% separate(col = new_col, into = c("orig.ident", "condition"), sep = ":")

    expect_equal(se %>% tidySE::select(sample) %>% ncol(), 1)
})

test_that("extract", {
    expect_equal(
        tt %>% tidySE::extract(col2, into = "g", regex = "other_([a-z]+)", convert = TRUE) %>% tidySE::pull(g) %>% class(),
        "character"
    )
})

test_that("pivot_longer", {
    expect_equal(
        tt %>% tidySE::pivot_longer(c(sample, condition), names_to = "name", values_to = "value") %>% class() %>% .[1],
        "tbl_df"
    )
})
