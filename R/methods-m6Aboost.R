## ============================================================================
## The Methods-m6Aboost for the GRanges objects
## ----------------------------------------------------------------------------

#' @import adabag
#' @import ExperimentHub
#' @importFrom Biostrings getSeq
#' @importFrom BSgenome getBSgenome
#' @importFrom rtracklayer import.bw

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## m6Aboost can only run agaist A sites, this step keep only the A sites
.keepA <- function(object, genome)
{
    seq <- Biostrings::getSeq(genome, object)
    object <- object[seq == "A"]
    return(object)
}

## build the data.frame for running the model
.transDF <- function(object, genome)
{
    seq <- getSeq(genome, object+10)
    df <- data.frame(seq_len(length(seq)))
    for (i in seq_len(21)) {
        k <- substr(seq,i,i)
        df[,i] <- as.factor(k)
    }
    l <- c()
    for (i in seq_len(21)) {
        l <- c(l, paste0("P",i))
    }
    colnames(df) <- l
    rownames(df) <- object$ID

    df$UTR3 <- object$UTR3
    df$UTR5 <- object$UTR5
    df$CDS <- object$CDS
    df$log2RSS <- log2(object$RSS+1)
    df$log2CtoT <- log2(object$CtoT+1)
    for (i in seq_len(24)) {
        df[,i] = as.factor(df[,i] )
    }
    return(df)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "m6Aboost" methods for GRanges objects.
##

#' @rdname m6Aboost
setMethod("m6Aboost", signature(object="GRanges"),
    function(object, genome="", normalization=TRUE)
    {
        if(!is.character(genome))
            stop("Genome should be the name of BSgenome annotation, e.g.
                BSgenome.Mmusculus.UCSC.mm10")
        if (unique(width(object)) != 1)
            stop("The input GRanges object should store the information
                of the crosslinking sites from miCLIP2 experiment. Its
                width should all equal to 1.")

        ## Keep only A sites
        genome <- BSgenome::getBSgenome(genome)
        object <- .keepA(object, genome)

        df <- .transDF(object, genome)
        ## load model
        eh <- ExperimentHub::ExperimentHub()
        model <- eh[["EH6021"]]

        df$log2SOB <- df$log2RSS
        df$log2RSS <- NULL
        if (normalization == T) {
            ## The normalization factor based on the our miLCIP2 data
            m <- mean(df$log2SOB)
            df$log2SOB <- df$log2SOB/(m/0.9122623)
            c <- mean(df$log2CtoT)
            df$log2CtoT <- df$log2CtoT/(c/0.4389567)
        }

        pred_res <- adabag::predict.boosting(model, newdata = df)
        object$class <- pred_res$class
        object$prob <- pred_res$prob[,1]
        return(object)
    }
)
