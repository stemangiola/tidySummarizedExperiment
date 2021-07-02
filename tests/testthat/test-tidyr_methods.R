context("tidyr test")

library(magrittr)
library(tidySummarizedExperiment)

tt <-
    pasilla %>%
    
    tidySummarizedExperiment::mutate(col2 = "other_col")

test_that("nest_unnest", {

    y <- tibble::tibble(
        .sample = c(
            "untrt1",
            "untrt2",
            "untrt3",
            "untrt4",
            "trt1",
            "trt2",
            "trt3"
        ),
        counts = c(0L, 0L, 0L, 0L, 0L, 0L, 1L)
    )

    x <- tt %>%
        nest(data = -condition) %>%
        unnest(data) %>%
        head(n = 1) %>%
        select(.sample, counts)


    expect_equal(x, y)
})

test_that("nest_unnest_slice_1",{
    
    tt %>%
        nest(data = -condition) %>% 
        slice(1) %>% 
        unnest(data)
    
})

test_that("unite separate", {
    un <- tt %>% unite("new_col", c(condition, col2), sep = ":")

    un %>%
        tidySummarizedExperiment::select(new_col) %>%
        slice(1) %>%
        pull(new_col) %>%
        expect_equal("untreated:other_col")

    se <-
        un %>%
        separate(
            col = new_col,
            into = c("orig.ident", "condition"),
            sep = ":"
        )

    se %>%
        tidySummarizedExperiment::select(.sample) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("extract", {
    tt %>%
        tidySummarizedExperiment::extract(col2,
                        into = "g",
                        regex = "other_([a-z]+)",
                        convert = TRUE) %>%
        tidySummarizedExperiment::pull(g) %>%
        class() %>%
        expect_equal("character")
})

test_that("pivot_longer", {
    tt %>%
        tidySummarizedExperiment::pivot_longer(c(.sample, condition),
                             names_to = "name",
                             values_to = "value") %>%
        class() %>%
        .[1] %>%
        expect_equal("tbl_df")
})
