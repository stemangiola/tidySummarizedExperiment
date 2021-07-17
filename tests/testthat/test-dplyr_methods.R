context("dplyr test")

library(tidySummarizedExperiment)


test_that("bind_rows", {
    pasilla_bind <- tidySummarizedExperiment::bind_rows(pasilla, pasilla)

    pasilla_bind %>%
        count(.sample, .feature) %>%
        dplyr::count(n) %>%
        filter(n > 1) %>%
        nrow() %>%
        expect_equal(0)
})

test_that("distinct", {
    pasilla %>%
        distinct(condition) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("filter", {
    pasilla %>%
        filter(condition == "untreated") %>%
        nrow() %>%
        expect_equal(14599)
})

test_that("group_by", {
    pasilla %>%
        group_by(condition) %>%
        ncol() %>%
        expect_equal(5)
})

test_that("summarise", {
    pasilla %>%
        summarise(mean(counts)) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("mutate", {
    pasilla %>%
        mutate(condition = 1) %>%
        distinct(condition) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("rename", {
    pasilla %>%
        rename(groups = condition) %>%
        select(groups) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("left_join", {
    expect_equal(
        pasilla %>%
            left_join(pasilla %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2)) %>%
            colData() %>%
            ncol(),
        pasilla %>%
            colData() %>%
            ncol() %>%
            sum(1)
    )
})

test_that("inner_join", {
    pasilla %>% inner_join(pasilla %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2) %>%
                          slice(1)) %>%
        ncol() %>%
        expect_equal(4)
})

test_that("right_join", {
    pasilla %>% right_join(pasilla %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2) %>%
                          slice(1)) %>%
        ncol() %>%
        expect_equal(4)
})

test_that("full_join", {
    pasilla %>%
        full_join(tibble::tibble(condition = "A",     other = 1:4)) %>% nrow() %>%
        expect_equal(102197)
})

test_that("slice", {
    pasilla %>%
        slice(1) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("select", {
    pasilla %>%
        select(-condition) %>%
        class() %>%
        as.character() %>%
        expect_equal("SummarizedExperiment")

    pasilla %>%
        select(condition) %>%
        class() %>%
        as.character() %>%
        .[1] %>%
        expect_equal("tbl_df")
})

test_that("sample_n", {
    pasilla %>%
        sample_n(50) %>%
        nrow() %>%
        expect_equal(50)
})

test_that("sample_frac", {
    pasilla %>%
        sample_frac(0.1) %>%
        nrow() %>%
        expect_equal(10219)
})

test_that("count", {
    pasilla %>%
        count(condition) %>%
        nrow() %>%
        expect_equal(2)
})
