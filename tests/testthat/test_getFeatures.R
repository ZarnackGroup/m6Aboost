context("getFeatures")
library(m6Aboost)

test_that("getFeatures works as expected",{
    testpath <- system.file("extdata", package = "m6Aboost")
    test <- readRDS(file.path(testpath, "test.rds"))
    test <- getFeatures(test, test_gff3,
        colname_reads="WTmean",
        colname_C2T="CtoTmean")

    expect_is(test, "GRanges")
}
)
