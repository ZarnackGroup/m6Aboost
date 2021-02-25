context("truncationAssignment")
library(m6Aboost)

test_that("truncationAssignment works as expected",{
    testpath <- system.file("extdata", package = "m6Aboost")
    test <- readRDS(file.path(testpath, "test.rds"))
    positiveBW <- file.path(testpath, "test_positive.bw")
    negativeBW <- file.path(testpath, "test_negative.bw")

    test <- truncationAssignment(test, positiveBW, negativeBW, "WT1")

    expect_is(test, "GRanges")
    expect_is(test$WT1, "numeric")
}
)
