## ============================================================================
## The Methods-truncationAssignment for the GRanges object
## ----------------------------------------------------------------------------

#' @importFrom rtracklayer import.bw
#' @importFrom S4Vectors elementMetadata

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## Import bigWig files with miCLIP2 truncation events
.readBW <- function(file_P, file_N){
    fp <- import.bw(con = file_P)
    fn <- import.bw(con = file_N)
    strand(fp) <- "+"
    strand(fn) <- "-"
    bw <- c(fp, fn)
    return(bw)
}

## Assign scores to GRanges object
.assign_bw_to_grange <- function(grange, bw, name = ""){
    o <- findOverlaps(grange,bw)
    elementMetadata(grange)[,name] <- 0
    elementMetadata(grange)[,name][queryHits(o)] <- bw$score[subjectHits(o)]
    return(grange)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "truncationAssignment" method for GRanges objects.
##

#' @title truncationAssignment for assigning the truncation read counts to the
#'     GRanges object with single nucleotide peaks
#'
#' @description An function for assigning the truncation read counts from
#'     bigWig to the GRanges peaks.
#' @author You Zhou
#'
#' @param object A GRanges object which should contains all the single
#'     nucleotide peaks of miCLIP2 experiment.
#' @param bw_positive A path to the bigWig file of truncation read counts
#'     at the positive strand that output from the preprocess in the m6Aboost
#'     pipeline.
#' @param bw_negative A path to the bigWig file of truncation read counts
#'     at the negative strand that output from the preprocess in the m6Aboost
#'     pipeline.
#' @param sampleName The column name that user would like to use for indicating
#'     the name of the sample.
#'
#' @return A GRanges object with the truncation read counts.
#'
#' @examples
#' if (.Platform$OS.type != "windows") {
#'     testpath <- system.file("extdata", package = "m6Aboost")
#'     test <- readRDS(file.path(testpath, "test.rds"))
#'     truncationBw_p <- file.path(testpath, "truncation_positive.bw")
#'     truncationBw_n <- file.path(testpath, "truncation_negative.bw")
#'     test <- truncationAssignment(test, bw_positive=truncationBw_p,
#'         bw_negative=truncationBw_n, sampleName = "WT1")
#'
#' }
#'
#' @export

truncationAssignment <-
    function(object, bw_positive, bw_negative, sampleName="")
    {
        if (nchar(sampleName) == 0)
            stop("Please add a sample name to assign the bigwig files.")
        ## Import bigWig files
        bw <- .readBW(bw_positive, bw_negative)
        ## Assign scores to GRanges
        object <- .assign_bw_to_grange(object, bw, name = sampleName)
        return(object)
    }
