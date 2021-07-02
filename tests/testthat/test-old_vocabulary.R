context("old vocabulary")

library(tidySummarizedExperiment)

warning_message = "the special columns including sample/feature"

test_that("distinct", {
    pasilla %>%
        distinct(sample, condition) %>%
        expect_warning(warning_message)
})

test_that("filter", {
    pasilla %>%
        filter(feature == "FBgn0000003") %>%
        expect_warning(warning_message)
})

test_that("group_by", {
    pasilla %>%
        group_by(sample) %>%
        expect_warning(warning_message)
})

test_that("summarise", {
    pasilla %>%
        summarise(unique(sample )) %>%
        expect_warning(warning_message)
})

test_that("mutate", {
    pasilla %>%
        mutate(condition = sample) %>%
        expect_warning(warning_message)
})

test_that("left_join", {
        pasilla %>%
            left_join(pasilla %>%
                          distinct(sample) %>%
                          mutate(new_column = 1:7)) %>%
            expect_warning(warning_message)
    
    
    pasilla %>%
        left_join(pasilla %>%
                      distinct(feature) %>%
                      mutate(new_column = 1:14599 )) %>%
        expect_warning(warning_message)
    
})

test_that("inner_join", {
    pasilla %>%
        inner_join(pasilla %>%
                      distinct(sample) %>%
                      mutate(new_column = 1:7)) %>%
        expect_warning(warning_message)
})

test_that("right_join", {
    pasilla %>%
        right_join(pasilla %>%
                      distinct(sample) %>%
                      mutate(new_column = 1:7)) %>%
        expect_warning(warning_message)
})

test_that("full_join", {
    pasilla %>%
        full_join(pasilla %>%
                      distinct(sample) %>%
                      mutate(new_column = 1:7)) %>%
        expect_warning(warning_message)
})

test_that("select", {
    pasilla %>%
        select(sample, feature, counts, condition) %>%
        expect_warning(warning_message)

    pasilla %>%
        select(condition) %>%
        class() %>%
        as.character() %>%
        .[1] %>%
        expect_equal("tbl_df")
})

test_that("count", {
    pasilla %>%
        count(sample, condition)  %>%
        expect_warning(warning_message)
})

test_that("pull", {
    pasilla %>%
        pull(sample, condition)  %>%
        expect_warning(warning_message)
})


library(magrittr)
library(tidySummarizedExperiment)

tt <-
    pasilla %>%
    
    tidySummarizedExperiment::mutate(col2 = "other_col")

test_that("nest_unnest", {
    

    
    tt %>%
        nest(data = -sample)  %>% 
        unnest(data) %>% 
        expect_warning(warning_message)

})


test_that("unite separate", {
    un <- 
        tt %>% 
        unite("new_col", c(condition, sample), sep = ":", remove = FALSE) %>% 
        expect_warning(warning_message)
    
        un %>%
        separate(
            col = feature,
            into = c("orig.ident", "condition"),
            sep = ":", remove = FALSE
        ) %>% 
            expect_warning(warning_message)
    

})

test_that("extract", {
    tt %>%
        tidySummarizedExperiment::extract(sample,
                                          into = "g",
                                          regex = "other_([a-z]+)",
                                          convert = TRUE, remove=FALSE) %>% 
        expect_warning(warning_message)
})

test_that("pivot_longer", {
    tt %>%
        tidySummarizedExperiment::pivot_longer(c(sample, condition),
                                               names_to = "name",
                                               values_to = "value") %>%
        class() %>%
        .[1] %>%
        expect_equal("tbl_df")
})
