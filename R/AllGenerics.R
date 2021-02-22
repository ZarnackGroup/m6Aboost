## All Generics

#' @title getFeatures for the miCILP2 data
#'
#' @description An function for calculating the relative signal strength and
#'     assign the features for the m6Aboost model
#'
#' @author You Zhou, Kathi Zarnack
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
#' @author You Zhou, Kathi Zarnack
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
#'
#' @export
#'
#' @rdname m6Aboost
setGeneric("m6Aboost",
           def=function(object, genome="", normalization=TRUE) {
               standardGeneric("m6Aboost")
           }
)
