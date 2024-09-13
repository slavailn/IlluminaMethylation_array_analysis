library(GEOquery)
setwd(<wd>)

# Download existing RA datasets
# GEO accession: GSE121192 

# Download GSE12192
gse <- getGEO("GSE121192")
show(gse)
gse[[1]]
show(pData(phenoData(gse[[1]]))[1:4,])

# Get sample data
eset <- gse[[1]]
write.csv(pData(eset), file="GSE121192_phenodata.csv")

# Get annotation
gpl <- getGEO(eset@annotation)
class(gpl)
head(gpl@dataTable@table)
annotation <- gpl@dataTable@table
head(annotation)
# Save annotation
write.csv(annotation, file="GSE121192_annotation.csv")

# Download raw files
getGEOSuppFiles("GSE121192", makeDirectory = TRUE, baseDir = getwd(),
                fetch_files = TRUE, filter_regex = NULL)
