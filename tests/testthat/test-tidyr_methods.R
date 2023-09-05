context("tidyr test")

library(magrittr)
library(tidySummarizedExperiment)

tt <-
    pasilla %>%
    
    tidySummarizedExperiment::mutate(col2 = "other_col")

# Create SummarizedExperiment object for testing
nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=rowRanges, colData=colData)
rownames(rse) <- sprintf("ID%03d", 1:200)

test_that("RangedSummarizedExperiment_nest_unnest", {
  tryCatch({
    rse_nested <- rse %>%
      nest(data = -.sample)

    rse_unnested <- rse_nested %>%
      unnest(data)
  })
  
  expect_equal(rse@colData, rse_unnested@colData)
  expect_equal(rse@rowRanges, rse_unnested@rowRanges)
})

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

test_that("nest_0_samples",{
  
  rowData(tt)$n = rep(1, nrow(tt))
  
  tt[,0] |> 
    nest(data = -n) 
  
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
            into = c( "condition", "col2"),
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

test_that("nest_unnest_by_feature_chunk", {
  
  chunks =
    tibble::tibble(.feature = rownames(tt)) |>
    mutate(chunk___ = c( 
      rep(1, times = floor(dplyr::n()/2)), 
      rep(2, times = ceiling(dplyr::n()/2))
    ))
  
  statistics_for_features = 
    tibble::tibble(.feature = rownames(tt)) |>
    mutate(pvalue = runif( dplyr::n(), min = 0, max = 1))
  
  tt = tt |> left_join(statistics_for_features)
    
  tt_unnested = 
    tt |> 
    left_join(chunks) |>
    nest(se_chunk = -chunk___) |> 
    unnest(se_chunk) |> 
    select(-chunk___)
  
  identical(tt, tt_unnested)
  
})