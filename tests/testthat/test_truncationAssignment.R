context("truncationAssignment")
library(m6Aboost)

test_that("truncationAssignment works as expected",{
    testpath <- system.file("extdata", package = "m6Aboost")
    test <- readRDS(file.path(testpath, "test.rds"))
    truncationBw_p <- file.path(testpath, "truncation_positive.bw")
    truncationBw_n <- file.path(testpath, "truncation_negative.bw")
    test <- truncationAssignment(test, bw_positive=truncationBw_p,
        bw_negative=truncationBw_n, sampleName = "WT1")

    expect_is(test, "GRanges")
    expect_is(test$WT1, "numeric")
}
)
