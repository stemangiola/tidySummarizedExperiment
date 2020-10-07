context("dplyr test")

library(tidySE)

tt <- pasilla %>% tidy()

test_that("bind_rows", {
    tt_bind <- tidySE::bind_rows(tt, tt)

    tt_bind %>%
        count(sample, transcript) %>%
        dplyr::count(n) %>%
        filter(n > 1) %>%
        nrow() %>%
        expect_equal(0)
})

test_that("distinct", {
    tt %>%
        distinct(condition) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("filter", {
    tt %>%
        filter(condition == "untreated") %>%
        nrow() %>%
        expect_equal(14599)
})

test_that("group_by", {
    tt %>%
        group_by(condition) %>%
        ncol() %>%
        expect_equal(5)
})

test_that("summarise", {
    tt %>%
        summarise(mean(counts)) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("mutate", {
    tt %>%
        mutate(condition = 1) %>%
        distinct(condition) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("rename", {
    tt %>%
        rename(groups = condition) %>%
        select(groups) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("left_join", {
    expect_equal(
        tt %>%
            left_join(tt %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2)) %>%
            `@`(colData) %>%
            ncol(),
        tt %>%
            `@`(colData) %>%
            ncol() %>%
            sum(1)
    )
})

test_that("inner_join", {
    tt %>% inner_join(tt %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2) %>%
                          slice(1)) %>%
        ncol() %>%
        expect_equal(4)
})

test_that("right_join", {
    tt %>% right_join(tt %>%
                          distinct(condition) %>%
                          mutate(new_column = 1:2) %>%
                          slice(1)) %>%
        ncol() %>%
        expect_equal(4)
})

test_that("full_join", {
    tt %>%
        full_join(tibble::tibble(condition = "A",     other = 1:4)) %>% nrow() %>%
        expect_equal(102197)
})

test_that("slice", {
    tt %>%
        slice(1) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("select", {
    tt %>%
        select(-condition) %>%
        class() %>%
        as.character() %>%
        expect_equal("tidySE")

    tt %>%
        select(condition) %>%
        class() %>%
        as.character() %>%
        .[1] %>%
        expect_equal("tbl_df")
})

test_that("sample_n", {
    tt %>%
        sample_n(50) %>%
        nrow() %>%
        expect_equal(50)
})

test_that("sample_frac", {
    tt %>%
        sample_frac(0.1) %>%
        nrow() %>%
        expect_equal(10219)
})

test_that("count", {
    tt %>%
        count(condition) %>%
        nrow() %>%
        expect_equal(2)
})
