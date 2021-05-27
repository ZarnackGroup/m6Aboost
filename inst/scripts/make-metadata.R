## Here makes the metadata.csv

metadata <- data.frame(Title = "The m6Aboost machine learning model",
    Description = "The machine learning model which use for identify the m6A
        signals from miCLIP2 data set",
    BiocVersion = "3.13",
    Genome = NA,
    SourceType = "RDS",
    SourceUrl = "https://github.com/Codezy99/m6Aboost/blob/main/m6Aboost.rds",
    SourceVersion = "1",
    Species = NA,
    TaxonomyId = NA,
    Coordinate_1_based = NA,
    DataProvider = "Zarnack's lab",
    Maintainer = "You Zhou <youzhoulearning@gmail.com>",
    RDataClass = "boosting",
    DispatchClass = "Rds",
    RDataPath = "m6Aboost/m6Aboost.rds",
    Tags = "")

write.csv(metadata,file = "inst/extdata/metadata.csv", row.names=FALSE)

