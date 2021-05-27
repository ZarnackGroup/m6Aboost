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
    fp <- rtracklayer::import.bw(con = file_P)
    fn <- rtracklayer::import.bw(con = file_N)
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

#' @rdname truncationAssignment
setMethod("truncationAssignment", signature(object="GRanges"),
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
)
