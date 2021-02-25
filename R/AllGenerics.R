## All Generics

#' @title getFeatures for the miCILP2 data
#'
#' @description An function for calculating the relative signal strength and
#'     assign the features for the m6Aboost model
#'
#' @author You Zhou
#'
#' @param object A GRanges object which should contains all the single
#'     nucleotide peaks of miCLIP2 experiment.
#' @param annotation A path way to the annotation file. The format of the
#'                   annotation file should be gff3 and downloaded from
#'                   https://www.gencodegenes.org/
#' @param colname_reads The name of the meta data column which contains the
#'     mean value of the truncation reads number without C to T transition
#'     reads.
#' @param colname_C2T The name of the meta data column which contains the
#'     mean value of C to T transition read counts.
#'
#' @return A GRanges object with all the information that required by the
#'     m6Aboost model.
#'
#' @examples
#'     testpath <- system.file("extdata", package = "m6Aboost")
#'     test_gff3 <- file.path(testpath, "test_annotation.gff3")
#'     test <- readRDS(file.path(testpath, "test.rds"))
#'     test<- getFeatures(test, test_gff3, colname_reads="WTmean",
#'         colname_C2T="CtoTmean")
#'
#' @export
#'
#' @rdname getFeatures
setGeneric("getFeatures",
    def=function(object, annotation, colname_reads="", colname_C2T="") {
        standardGeneric("getFeatures")
    }
)

#' @title getFeatures for the miCILP2 data
#'
#' @description An function for calculating the relative signal strength and
#'     assign the features for the m6Aboost model
#'
#' @author You Zhou
#'
#' @param object A GRanges object which should contains all the single
#'     nucleotide peaks of miCLIP2 experiment.
#' @param genome The name of the BSgenome that you are working with. For
#'     example "BSgenome.Mmusculus.UCSC.mm10".
#' @param normalization A logical vector which indicates whether you would like
#'     normalize the RSS and C to T reads number to the mean value of our
#'     training set. This will help to reduce the false positive rate.
#'
#' @return A GRanges object with all the information that required by the
#'     m6Aboost model.
#' @examples
#'     testpath <- system.file("extdata", package = "m6Aboost")
#'     test_gff3 <- file.path(testpath, "test_annotation.gff3")
#'     test <- readRDS(file.path(testpath, "test.rds"))
#'     test<- getFeatures(test, test_gff3, colname_reads="WTmean",
#'         colname_C2T="CtoTmean")
#'
#'     ## The input of m6Aboost should be the output from getFeatures function
#'     ## Please make sure that the correct BSgenome package have installed
#'     ## before running motifProfile. For example,
#'     ## library("BSgenome.Mmusculus.UCSC.mm10")
#'
#'     test <- m6Aboost(test, "BSgenome.Mmusculus.UCSC.mm10")
#' @export
#'
#' @rdname m6Aboost
setGeneric("m6Aboost",
    def=function(object, genome="", normalization=TRUE) {
        standardGeneric("m6Aboost")
    }
)

#' @title truncationAssignment for assigning the truncation read counts
#'
#' @description An function for assigning the truncation read counts from
#'     bigwig.
#' @author You Zhou
#'
#' @param object A GRanges object which should contains all the single
#'     nucleotide peaks of miCLIP2 experiment.
#' @param bw_positive A path to the big wig file of truncation read counts
#'     at the positive strand that output from the preprocess in the m6Aboost
#'     pipeline.
#' @param bw_negative A path to the big wig file of truncation read counts
#'     at the negative strand that output from the preprocess in the m6Aboost
#'     pipeline.
#' @param sampleName The sample name of the big wig file.
#'
#' @return A GRanges object with the truncation read counts.
#'
#' @examples
#'     testpath <- system.file("extdata", package = "m6Aboost")
#'     test <- readRDS(file.path(testpath, "test.rds"))
#'     positiveBW <- file.path(testpath, "test_positive.bw")
#'     negativeBW <- file.path(testpath, "test_negative.bw")
#'
#'     test <- truncationAssignment(test, positiveBW, negativeBW, "WT1")
#'
#' @export
#'
#' @rdname truncationAssignment
setGeneric("truncationAssignment",
    def=function(object, bw_positive, bw_negative, sampleName="") {
        standardGeneric("truncationAssignment")
    }
)

#' @title CtoTAssignment for assigning the truncation read counts
#'
#' @description An function for assigning the CtoT transition read counts from
#'     bigwig.
#' @author You Zhou
#'
#' @param object A GRanges object which should contains all the single
#'     nucleotide peaks of miCLIP2 experiment.
#' @param bw_positive A path to the big wig file of C to T transition read
#'     counts at the positive strand that output from the preprocess in the
#'     m6Aboost pipeline.
#' @param bw_negative A path to the big wig file of C to T transition read
#'     counts at the negative strand that output from the preprocess in the
#'     m6Aboost pipeline.
#' @param sampleName The sample name of the big wig file.
#'
#' @return A GRanges object with the truncation read counts.
#'
#' @examples
#'     testpath <- system.file("extdata", package = "m6Aboost")
#'     test <- readRDS(file.path(testpath, "test.rds"))
#'     positiveBW <- file.path(testpath, "test_positive.bw")
#'     negativeBW <- file.path(testpath, "test_negative.bw")
#'
#'     test <- truncationAssignment(test, positiveBW, negativeBW, "CtoT1")
#' @export
#'
#' @rdname CtoTAssignment
setGeneric("CtoTAssignment",
    def=function(object, bw_positive, bw_negative, sampleName="") {
        standardGeneric("CtoTAssignment")
    }
)
