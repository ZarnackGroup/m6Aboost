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

#' @title m6Aboost for identify the m6A peaks from the miCILP2 data
#'
#' @description An function for calculating the relative signal strength and
#'     extracting all the features that required by the m6Aboost model for
#'     each peak.
#'
#' @author You Zhou
#'
#' @param object A GRanges object which should contains all the single
#'     nucleotide peaks of miCLIP2 experiment.
#' @param genome The name of the BSgenome that you are working with. For
#'     example "BSgenome.Mmusculus.UCSC.mm10".
#' @param normalization A logical vector which indicates whether you would like
#'     normalize the RSS and C to T reads number to the mean value of the
#'     training set of the model. This will help to reduce the false positive
#'     rate.
#'
#' @return A GRanges object with all the information that is required by the
#'     m6Aboost model.
#' @examples
#'     testpath <- system.file("extdata", package = "m6Aboost")
#'     test_gff3 <- file.path(testpath, "test_annotation.gff3")
#'     test <- readRDS(file.path(testpath, "test.rds"))
#'     test<- preparingData(test, test_gff3, colname_reads="WTmean",
#'         colname_C2T="CtoTmean")
#'
#'     ## The input of m6Aboost should be the output from preparingData function
#'     ## Please make sure that the correct BSgenome package have installed
#'     ## before running motifProfile. For example,
#'     ## library("BSgenome.Mmusculus.UCSC.mm10")
#'
#'     test <- m6Aboost(test, "BSgenome.Mmusculus.UCSC.mm10")
#' @export

m6Aboost <- function(object, genome="", normalization=TRUE)
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
        if (normalization == TRUE) {
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
