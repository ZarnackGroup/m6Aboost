context("preparingData")
library(m6Aboost)

test_that("preparingData works as expected",{
    testpath <- system.file("extdata", package = "m6Aboost")
    test <- readRDS(file.path(testpath, "test.rds"))
    test_gff3 <- file.path(testpath, "test_annotation.gff3")
    test <- preparingData(test, test_gff3,
        colname_reads="WTmean",
        colname_C2T="CtoTmean")

    expect_is(test, "GRanges")
}
)
