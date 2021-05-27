context("m6Aboost")
library(m6Aboost)

test_that("m6Aboost works as expected",{
    testpath <- system.file("extdata", package = "m6Aboost")
    test <- readRDS(file.path(testpath, "test.rds"))
    test_gff3 <- file.path(testpath, "test_annotation.gff3")
    test <- getFeatures(test, test_gff3,
                        colname_reads="WTmean",
                        colname_C2T="CtoTmean")
    test <- m6Aboost(test, "BSgenome.Mmusculus.UCSC.mm10")

    expect_is(test, "GRanges")
}
)
