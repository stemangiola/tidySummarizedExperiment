context("felix test")

library(magrittr)
library(tidySummarizedExperiment)

# Create dataset
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=rowRanges, colData=colData)



test_that("Example 1 all columns included", {
    
    rse  %>% 
        as_tibble() %>%
        nrow() %>%
        expect_equal(1200)
})


test_that("Example 2", {
    
    colData(rse)$sample <- seq_len(ncol(rse))
    rowData(rse)$transcript <- seq_len(nrow(rse))
    rse  %>% 
        as_tibble() %>%
        nrow() %>%
        expect_equal(1200)
})

test_that("Example 3", {
    
    # rowRanges(rse) <- split(rowRanges(rse),seq_len(nrow(rse)))
    # 
    # rse  %>% 
    #     as_tibble() %>%
    #     nrow() %>%
    #     expect_equal(1200)
    
    rowData(rse)$transcript <- seq_len(nrow(rse))
    
    rse  %>% 
        as_tibble() %>%
        nrow() %>%
        expect_equal(1200)
    
    colnames(rse) <- NULL
    
    rse  %>% 
        as_tibble() %>%
        select(.sample, .feature) %>%
        ncol() %>%
        expect_equal(2) 
})

test_that("Example 4 from tidybulk", {
    
   x = se %>% as_tibble()
})