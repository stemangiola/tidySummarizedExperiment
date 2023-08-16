test_that("get_count_datasets works", {
    # Use testthat edition 3
    local_edition(3)
    
    # Generate test SE
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            mat1 = matrix(seq_len(9), nrow = 3),
            mat2 = matrix(seq(10, 18), nrow = 3),
            mat3 = matrix(seq(19, 27), nrow = 3)
        )
    )
    rownames(se) <- paste0("G", seq_len(3))
    colnames(se) <- paste0("S", seq_len(3))
    SummarizedExperiment::assays(se) <- lapply(
        SummarizedExperiment::assays(se),
        function(x) {
            rownames(x) <- rownames(se)
            colnames(x) <- colnames(se)
            x
        }
    )
    
    # Check that dimnames are as expected
    expect_equal(rownames(se), rownames(assay(se, "mat1", withDimnames = FALSE)))
    expect_equal(rownames(se), rownames(assay(se, "mat2", withDimnames = FALSE)))
    expect_equal(rownames(se), rownames(assay(se, "mat3", withDimnames = FALSE)))
    expect_equal(colnames(se), colnames(assay(se, "mat1", withDimnames = FALSE)))
    expect_equal(colnames(se), colnames(assay(se, "mat2", withDimnames = FALSE)))
    expect_equal(colnames(se), colnames(assay(se, "mat3", withDimnames = FALSE)))
    # Check that matrix values are as expected
    expect_equal(assay(se, "mat1")[, 1], c(1, 2, 3), ignore_attr = TRUE)
    expect_equal(assay(se, "mat2")[, 2], c(13, 14, 15), ignore_attr = TRUE)
    expect_equal(assay(se, "mat3")[, 3], c(25, 26, 27), ignore_attr = TRUE)
    
    # SE and all assays have the same dimnames
    cds <- get_count_datasets(se)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE does not have dimnames, all assays have the same
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    expect_equal(colnames(assay(se1, "mat1", withDimnames = FALSE)), paste0("S", seq_len(3)))
    expect_equal(rownames(assay(se1, "mat1", withDimnames = FALSE)), paste0("G", seq_len(3)))
    expect_null(colnames(se1))
    expect_null(rownames(se1))
    expect_warning(expect_warning(cds <- get_count_datasets(se1), "have column names, but the SummarizedExperiment does not"), "have row names, but the SummarizedExperiment does not")
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE does not have dimnames, two assays have the same, third assay does not have -> error
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- 
        colnames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    expect_error(get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have column names")
    
    # SE does not have dimnames, assays have the same but in different order
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- c("S2", "S3", "S1")
    rownames(assay(se1, "mat3", withDimnames = FALSE)) <- c("G2", "G3", "G1")
    expect_warning(expect_warning(cds <- get_count_datasets(se1), "have column names, but the SummarizedExperiment does not"), "have row names, but the SummarizedExperiment does not")
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, c(16, 17, 18, 10, 11, 12, 13, 14, 15))
    expect_equal(cds$mat3, c(21, 19, 20, 24, 22, 23, 27, 25, 26))
    
    # SE does not have dimnames, assays have nonoverlapping dimnames -> error
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- paste0("A", seq_len(3))
    expect_error(get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have row names")
    
    se1 <- se
    rownames(se1) <- colnames(se1) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- paste0("A", seq_len(3))
    expect_error(get_count_datasets(se1), "at least one of the assays in your SummarizedExperiment have column names")
    
    # Neither SE nor assays have column names
    se1 <- se
    colnames(se1) <- NULL
    colnames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    cds <- get_count_datasets(se1)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(as.character(seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # Neither SE nor assays have row names
    se1 <- se
    rownames(se1) <- NULL
    rownames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    cds <- get_count_datasets(se1)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(as.character(seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # Neither SE nor assays have any dimnames
    se1 <- se
    rownames(se1) <- NULL
    colnames(se1) <- NULL
    rownames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    cds <- get_count_datasets(se1)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(as.character(seq_len(3)), 3))
    expect_equal(cds$.sample, rep(as.character(seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # SE has dimnames, assays have the same dimnames but not overlapping with those of the SE
    se1 <- se
    rownames(se1) <- colnames(se1) <- paste0("A", seq_len(3))
    expect_error(get_count_datasets(se1), "don't agree with the column names of the SummarizedExperiment")
    
    se1 <- se
    rownames(se1) <- paste0("A", seq_len(3))
    expect_error(get_count_datasets(se1), "don't agree with the row names of the SummarizedExperiment")
    
    # SE has dimnames, none of the assays have
    se1 <- se
    rownames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat1", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat2", withDimnames = FALSE)) <- NULL
    rownames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    colnames(assay(se1, "mat3", withDimnames = FALSE)) <- NULL
    cds <- get_count_datasets(se1)
    expect_s3_class(cds, "tbl_df")
    expect_equal(nrow(cds), 9L)
    expect_equal(ncol(cds), 5L)
    expect_named(cds, c(".feature", ".sample", "mat1", "mat2", "mat3"))
    expect_equal(cds$.feature, rep(paste0("G", seq_len(3)), 3))
    expect_equal(cds$.sample, rep(paste0("S", seq_len(3)), each = 3))
    expect_equal(cds$mat1, seq(1, 9))
    expect_equal(cds$mat2, seq(10, 18))
    expect_equal(cds$mat3, seq(19, 27))
    
    # Unnamed assay(s)
    # se1 <- SummarizedExperiment::SummarizedExperiment(
    #     assays = list(
    #         matrix(seq_len(9), nrow = 3),
    #         mat2 = matrix(seq(10, 18), nrow = 3),
    #         matrix(seq(19, 27), nrow = 3)
    #     )
    # )
    # expect_error(cds <- get_count_datasets(se1), "assays must be named")
    # 
    # se1 <- SummarizedExperiment::SummarizedExperiment(
    #     assays = list(
    #         matrix(seq_len(9), nrow = 3),
    #         matrix(seq(10, 18), nrow = 3),
    #         matrix(seq(19, 27), nrow = 3)
    #     )
    # )
    # expect_error(cds <- get_count_datasets(se1), "assays must be named")
})
