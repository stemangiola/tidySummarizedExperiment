context("dplyr test")

library(tidySE)

tt <- pasilla %>% tidy()

test_that("bind_rows", {
    tt_bind <- tidySE::bind_rows(tt, tt )

    expect_equal(
        tt_bind %>% count(sample, transcript) %>% dplyr::count(n) %>% filter(n>1) %>% nrow(),
        0
    )
})

test_that("distinct", {
    expect_equal(tt %>% distinct(condition) %>% ncol(), 1)
})

test_that("filter", {
    expect_equal(tt %>% filter(condition == "untreated") %>% nrow(), 14599)
})

test_that("group_by", {
    expect_equal(tt %>% group_by(condition) %>% ncol(), 5)
})

test_that("summarise", {
    expect_equal(tt %>% summarise(mean(counts)) %>% nrow(), 1)
})

test_that("mutate", {
    expect_equal(tt %>% mutate(condition = 1) %>% distinct(condition) %>% nrow(), 1)
})

test_that("rename", {
    expect_equal(tt %>% rename( groups = condition) %>% select(groups) %>% ncol(), 1)
})

test_that("left_join", {
    expect_equal(
        tt %>% left_join(tt %>% distinct(condition) %>% mutate(new_column = 1:2)) %>% `@`(colData) %>% ncol(),
        tt %>% `@`(colData) %>% ncol() %>% sum(1)
    )
})

test_that("inner_join", {
    expect_equal(tt %>% inner_join(tt %>% distinct(condition) %>% mutate(new_column = 1:2) %>% slice(1)) %>% ncol(), 4)
})

test_that("right_join", {
    expect_equal(tt %>% right_join(tt %>% distinct(condition) %>% mutate(new_column = 1:2) %>% slice(1)) %>% ncol(), 4)
})

test_that("full_join", {
    expect_equal(tt %>% full_join(tibble::tibble(condition = "A", other = 1:4)) %>% nrow(), 102197)
})

test_that("slice", {
    expect_equal(tt %>% slice(1) %>% ncol(), 1)
})

test_that("select", {
    expect_equal(tt %>% select(-condition) %>% class() %>% as.character(), "tidySE")

    expect_equal(tt %>% select(condition) %>% class() %>% as.character() %>% .[1], "tbl_df")
})

test_that("sample_n", {
    expect_equal(tt %>% sample_n(50) %>% nrow(), 50)
})

test_that("sample_frac", {
    expect_equal(tt %>% sample_frac(0.1) %>% nrow(), 10219)
})

test_that("count", {
    expect_equal(tt %>% count(condition) %>% nrow(), 2)
})
