## ============================================================================
## The Methods-preparingData for the GRanges objects
## ----------------------------------------------------------------------------

#' @import GenomicRanges
#' @import methods
#' @rawNamespace import(IRanges, except = c(collapse, slice, desc))
#' @rawNamespace import(dplyr, except = c(union, intersect, setdiff))
#' @importFrom S4Vectors mcols
#' @importFrom utils globalVariables

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## Assign the miCLIP2 peaks to transcript regions (UTR3, UTR5, CDS)
.assignment <- function(object, UTR3, CDS, UTR5)
{
    ## Define whether the peak would overlap to at least one of the Genomic
    ## features
    object$UTR3 <- ifelse(object %over% UTR3, "YES", "NO")
    object$UTR5 <- ifelse(object %over% UTR5, "YES", "NO")
    object$CDS <- ifelse(object %over% CDS, "YES", "NO")
    return(object)
}

## To fix the global variable note
utils::globalVariables(c("elementMetadata", "queryHits" , "geneid",
    "reads", "predict.boosting", "elementMetadata<-",
    "subjectHits", "strand<-", "start<-", "end<-", "gene_id"))

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "preparingData" methods for GRanges objects.
##

#' @title preparingData for the miCILP2 data
#'
#' @description A function for calculating the relative signal strength and
#'     extract the features for running the m6Aboost.
#'
#' @author You Zhou
#'
#' @param object A GRanges object which should contain all single-
#'     nucleotide peaks from the miCLIP2 experiment.
#' @param annotation A path to the annotation file. The format of the
#'     annotation file should be a gff3 file and downloaded from
#'     https://www.gencodegenes.org/
#' @param colname_reads The name of the metadata column which contains the
#'     mean value of the truncation reads number without C to T transition
#'     reads.
#' @param colname_C2T The name of the meta data column which contains the
#'     mean value of C to T transition read counts.
#'
#' @return A GRanges object with all information that is required for running
#'     the m6Aboost model.
#'
#' @examples
#'     testpath <- system.file("extdata", package = "m6Aboost")
#'     test_gff3 <- file.path(testpath, "test_annotation.gff3")
#'     test <- readRDS(file.path(testpath, "test.rds"))
#'     test<- preparingData(test, test_gff3, colname_reads="WTmean",
#'         colname_C2T="CtoTmean")
#'
#' @export

preparingData <-
    function(object, annotation, colname_reads="", colname_C2T="")
    {
        if(!isS4(object))
            stop("The input object should be a GRanges object.")
        meta_name <- names(elementMetadata(object))
        if(!colname_reads %in% meta_name |!colname_C2T %in% meta_name)
            stop("colname_reads and colname_C2T should refer to the respective
                column names in the elementMetadata of the input object.")

        ## Load annotation file
        anno <- rtracklayer::import.gff3(con = annotation)
        anno_gene <- anno[anno$type == "gene"]
        anno_gene$level <- as.numeric(anno_gene$level)
        anno_gene <-  anno_gene[order(anno_gene$level, -width(anno_gene))]

        o <- findOverlaps(object, anno_gene, select = "first")
        object$gene_id <- anno_gene$gene_id[o]
        object$gene_id[is.na(object$gene_id)] <- "NO"

        object <- as.data.frame(object)
        ## Calculate Relative Signal Strength (RSS)
        object <- object %>% group_by(gene_id) %>%
            mutate(factor=mean(get(colname_reads)))

        object <- mutate(object, RSS=get(colname_reads)/factor)
        object$factor <- NULL
        object <- makeGRangesFromDataFrame(object, keep.extra.columns = TRUE)
        ## make the column for running the m6Aboost
        object$CtoT <- elementMetadata(object)[,colname_C2T]

        UTR5 <- anno[anno$type == "five_prime_UTR"]
        UTR3 <- anno[anno$type == "three_prime_UTR"]
        CDS  <- anno[anno$type == "CDS"]

        object <- .assignment(object, UTR3, CDS, UTR5)
        return(object)
    }

