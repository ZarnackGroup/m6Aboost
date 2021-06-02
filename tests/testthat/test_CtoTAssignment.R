context("CtoTAssignment")
library(m6Aboost)

if (.Platform$OS.type != "windows") {
test_that("CtoTAssignment works as expected",{
    testpath <- system.file("extdata", package = "m6Aboost")
    test <- readRDS(file.path(testpath, "test.rds"))
    ctotBw_p <- file.path(testpath, "C2T_positive.bw")
    ctotBw_n <- file.path(testpath, "C2T_negative.bw")
    test <- CtoTAssignment(test, bw_positive=ctotBw_p, bw_negative=ctotBw_n,
        sampleName = "CtoT_WT1")

    expect_is(test, "GRanges")
    expect_is(test$CtoT_WT1, "numeric")
}
)
}
