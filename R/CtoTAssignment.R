## ============================================================================
## The Methods-CtoTAssignment for the GRanges objects
## ----------------------------------------------------------------------------

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## Assign bigWig files with C-to-T transitions to GRanges
.readBW_C2T <- function(file_P, file_N){
    fp <- import.bw(con = file_P)
    fn <- import.bw(con = file_N)
    strand(fp) <- "+"
    strand(fn) <- "-"
    ## shift 1nt to the A position (easy to assign the value to peaks)
    fp <- shift(fp, shift = -1)
    fn <- shift(fn, shift = 1)
    ## Combine both strands
    bw <- c(fp, fn)
    return(bw)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "CtoTAssignment" method for GRanges objects.
##

#' @title CtoTAssignment for assigning the truncation read counts to the
#'     GRanges object with single nucleotide peaks
#'
#' @description An function for assigning the CtoT transition read counts from
#'     bigWig.
#' @author You Zhou
#'
#' @param object A GRanges object which should contains all the single
#'     nucleotide peaks of miCLIP2 experiment.
#' @param bw_positive A path to the bigWig file of C to T transition read
#'     counts at the positive strand that output from the preprocess in the
#'     m6Aboost pipeline.
#' @param bw_negative A path to the bigWig file of C to T transition read
#'     counts at the negative strand that output from the preprocess in the
#'     m6Aboost pipeline.
#' @param sampleName The column name that user would like to use for indicating
#'     the name of the sample.
#'
#' @return A GRanges object with the truncation read counts.
#'
#' @examples
#' if (.Platform$OS.type != "windows") {
#'     testpath <- system.file("extdata", package = "m6Aboost")
#'     test <- readRDS(file.path(testpath, "test.rds"))
#'     ctotBw_p <- file.path(testpath, "C2T_positive.bw")
#'     ctotBw_n <- file.path(testpath, "C2T_negative.bw")
#'     test <- CtoTAssignment(test, bw_positive=ctotBw_p, bw_negative=ctotBw_n,
#'         sampleName = "CtoT_WT1")
#' }
#'
#' @export

CtoTAssignment <-
    function(object, bw_positive, bw_negative, sampleName="")
    {
        if (nchar(sampleName) == 0)
            stop("Please add a sample name to assign the bigWig files")

        ## Import bigWig files with C-to-T transitions
        bw <- .readBW_C2T(bw_positive, bw_negative)
        ## Assign C-to-T transitions to GRanges
        object <- .assign_bw_to_grange(object, bw, name = sampleName)
        return(object)
    }
