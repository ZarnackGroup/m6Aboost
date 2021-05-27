## ============================================================================
## The Methods-CtoTAssignment for the GRanges objects
## ----------------------------------------------------------------------------

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## Assign bigWig files with C-to-T transitions to GRanges
.readBW_C2T <- function(file_P, file_N){
    fp <- rtracklayer::import.bw(con = file_P)
    fn <- rtracklayer::import.bw(con = file_N)
    strand(fp) <- "+"
    strand(fn) <- "-"
    ## shift 1nt of C2T position (easy to assign the value to peaks)
    fp <- shift(fp, shift = -1)
    fn <- shift(fn, shift = 1)
    ## Combine both strands and shift to A position
    bw <- c(fp, fn)
    return(bw)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "CtoTAssignment" method for GRanges objects.
##

#' @rdname CtoTAssignment
setMethod("CtoTAssignment", signature(object="GRanges"),
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
)
