## ============================================================================
## The Methods-CtoTAssignment for the GRanges objects
## ----------------------------------------------------------------------------

#' @import dplyr
#' @importFrom rtracklayer import.bw
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges findOverlaps

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## Assign C to T big wig to GRanges
.readBW_C2T <- function(file_P, file_N){
    fp <- rtracklayer::import.bw(con = file_P)
    fn <- rtracklayer::import.bw(con = file_N)
    strand(fp) <- "+"
    strand(fn) <- "-"

    start(fp) <- start(fp)-1
    end(fp) <- end(fp)-1

    start(fn) <- start(fn)+1
    end(fn) <- end(fn)+1

    bw <- c(fp, fn)
    return(bw)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "CtoTAssignment" methods for GRanges objects.
##

#' @rdname CtoTAssignment
setMethod("CtoTAssignment", signature(object="GRanges"),
    function(object, bw_positive, bw_negative, sampleName="")
    {
        if (nchar(sampleName) == 0)
            stop("Please add a sample name to assigned big wig file")

        ## load BW C to T files
        bw <- .readBW_C2T(bw_positive, bw_negative)
        ## assign BW to GRanges
        object <- .assign_bw_to_grange(object, bw, name = sampleName)
        return(object)
    }
)
