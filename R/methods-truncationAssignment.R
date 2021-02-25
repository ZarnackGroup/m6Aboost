## ============================================================================
## The Methods-truncationAssignment for the GRanges objects
## ----------------------------------------------------------------------------

#' @import dplyr
#' @importFrom rtracklayer import.bw
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors elementMetadata

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## Assign big wig to GRanges
.readBW <- function(file_P, file_N){
    fp <- rtracklayer::import.bw(con = file_P)
    fn <- rtracklayer::import.bw(con = file_N)
    strand(fp) <- "+"
    strand(fn) <- "-"
    bw <- c(fp, fn)
    return(bw)
}

.assign_bw_to_grange <- function(grange, bw, name = ""){
    o <- findOverlaps(grange,bw)
    n <- ncol(elementMetadata(grange))
    elementMetadata(grange)[,n+1] <- 0
    elementMetadata(grange)[,n+1][queryHits(o)] <- bw$score[subjectHits(o)]
    colnames(elementMetadata(grange))[n+1] <- name
    return(grange)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "truncationAssignment" methods for GRanges objects.
##

#' @rdname truncationAssignment
setMethod("truncationAssignment", signature(object="GRanges"),
    function(object, bw_positive, bw_negative, sampleName="")
    {
        if (nchar(sampleName) == 0)
            stop("Please add a sample name to assigned big wig file")
        ## load BW files
        bw <- .readBW(bw_positive, bw_negative)
        ## assign BW to GRanges
        object <- .assign_bw_to_grange(object, bw, name = sampleName)
        return(object)
    }
)
