## ============================================================================
## The Methods-getFeatures for the GRanges objects
## ----------------------------------------------------------------------------

#' @import dplyr
#' @import methods
#' @import S4Vectors
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges findOverlaps
#' @importFrom utils globalVariables

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## Assign the UTR3 CDS UTR5
.assignment <- function(object, UTR3, CDS, UTR5)
{
    object$UTR3 <- "NO"
    o <- findOverlaps(object, UTR3)
    o <- unique(queryHits(o))
    object$UTR3[o] <- "YES"
    object$UTR5 <- "NO"
    o <- findOverlaps(object, UTR5)
    o <- unique(queryHits(o))
    object$UTR5[o] <- "YES"
    object$CDS <- "NO"
    o <- findOverlaps(object, CDS)
    o <- unique(queryHits(o))
    object$CDS[o] <- "YES"
    return(object)
}

## To fix the global variable note
utils::globalVariables(c("elementMetadata", "queryHits" , "geneid",
                         "reads", "predict.boosting", "elementMetadata<-",
                         "subjectHits", "strand<-", "start<-", "end<-"))

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "getFeatures" methods for GRanges objects.
##

#' @rdname getFeatures
setMethod("getFeatures", signature(object="GRanges"),
    function(object, annotation, colname_reads="", colname_C2T="")
    {
        if(!isS4(object))
            stop("The input object should be a GRanges object.")
        meta_name <- names(elementMetadata(object))
        if(isTRUE(!colname_reads %in% meta_name) |
            isTRUE(!colname_C2T %in%meta_name))
            stop("The colname_reads and colname_C2T should be a column name
                in the elementMetadata of the input object.")

        ## Load annotation file
        anno <- rtracklayer::import.gff3(con = annotation)
        anno_gene <- anno[anno$type == "gene"]
        anno_gene$level <- as.numeric(anno_gene$level)
        anno_gene <-  anno_gene[order(anno_gene$level, -width(anno_gene))]

        o <- findOverlaps(object, anno_gene, select = "first")
        object$gene_id <- anno_gene$gene_id[o]
        object$gene_id[is.na(object$gene_id)] <- "NO"

        meta <- elementMetadata(object)
        ## Calculate Relative Signal Strength (RSS)
        df <- data.frame(geneid=object$gene_id,
            reads=meta[colnames(meta) == colname_reads])
        colnames(df)[2] <- "reads"
        df1 <- df %>% group_by(geneid) %>% summarise(reads=mean(reads))

        object$factor <- df1$reads[match(object$gene_id, df1$geneid)]
        n <- ncol(meta)
        number <- unlist(meta[colnames(meta) == colname_reads])
        object$RSS <- number/object$factor
        object$factor <- NULL

        object$CtoT <- unlist(meta[colnames(meta) == colname_C2T])

        ## Annotation assignment
        UTR5 <- anno[anno$type == "five_prime_UTR"]
        UTR3 <- anno[anno$type == "three_prime_UTR"]
        CDS  <- anno[anno$type == "CDS"]

        object <- .assignment(object, UTR3, CDS, UTR5)
        return(object)
    }
)

