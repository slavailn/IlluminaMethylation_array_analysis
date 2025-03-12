library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

setwd(<working_directory>)
list.files()

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

# Location of raw files (IDAT files downloaded from GEO)
list.files("raw_files/")

# get metadata for the experiment
meta <- read.csv("docs/GSE42861_phenodata.csv", header = T, row.names = 1)
head(meta)
basename(meta$supplementary_file)
# Create a vector containing a path to IDAT files
Basename <- gsub("_Grn.idat.gz", "", basename(meta$supplementary_file))
Basename
 
# Create targets file for minfi rea.metharray.exp() function
Sample_name <- meta$geo_accession
head(meta)

# Modify variables as necessary
# Edit trait variable
Trait <- meta$disease.state.ch1
Trait[which(Trait == "rheumatoid arthritis")] <- "RA"
Trait[which(Trait == "Normal")] <- "Control"
Trait
# Gender: male and female
Gender <- meta$characteristics_ch1.4
Gender <- gsub("gender: ", "", Gender)

# Age
Age <- meta$characteristics_ch1.3
Age <- gsub("age: ", "", Age)

# Smoking status
Smoking_status <- meta$smoking.status.ch1
table(Smoking_status)

targets <- data.frame(Sample_name = Sample_name,
                      Trait = Trait,
                      Gender = Gender,
                      Age = as.numeric(Age),
                      Smoking_status = Smoking_status,
                      Basename = Basename)
head(targets)
write.csv(targets, file = "targets.csv", row.names = F)

# Read in the raw data from the IDAT files
rgSet <- read.metharray.exp(base = "raw_files", 
                            targets = targets)
rgSet
pData(rgSet)

# Give shorter column names
sampleNames(rgSet) <- targets$Sample_name
rgSet
save(rgSet, file="rgSet.RData")

######################## QUALITY CONTROL ######################################
detP <- detectionP(rgSet)
head(detP)

# examine mean detection p-values across all samples to identify any failed samples
pdf("mean_pvalues_arrays.pdf")
par(mfrow=c(1,1))
barplot(colMeans(detP), col="lightgray", las=2, 
        cex.names=0.8, ylab="Mean detection p-values",
        main = "Mean p-values across the arrays")
dev.off()

# None of samples have mean p-values over 0.05
# Print QC report
qcReport(rgSet, sampNames=targets$Sample_name, sampGroups=targets$Trait, 
         pdf="qcReport.pdf")

################# Pre-processing and normalization ##########################
# Apply preprocessFunnorm() in case of comparison groups with vastly different
# methylation profiles, for example different tissues or tumor vs normal
# Use preprocessQuantile() in the cases where conditions are not expected to be
# too different
# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet) 

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# Set color palette
pal <- brewer.pal(8,"Dark2")

# visualise what the data looks like before and after normalisation
pdf("beta_densities.pdf")
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Trait,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Trait)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Trait,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Trait)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS_PC1_PC2.pdf", height = 7, width = 10)
par(mfrow=c(1,1))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Trait)])
legend("top", legend=levels(factor(targets$Trait)), text.col=pal,
       bg="white", cex=0.7)
dev.off()

# Examine higher dimensions to look at other sources of variation
pdf("MDS_PC1_PC3.pdf", height = 7, width = 10)
par(mfrow=c(1,1))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Trait)], dim = c(1,3))
legend("top", legend=levels(factor(targets$Trait)), text.col=pal,
       bg="white", cex=0.7)
dev.off()

# The separation between RA and Controls is found at PC3 and PC4
pdf("MDS_PC3_PC4.pdf", height = 7, width = 10)
par(mfrow=c(1,1))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Trait)], dim = c(3,4))
legend("top", legend=levels(factor(targets$Trait)), text.col=pal,
       bg="white", cex=0.7)
dev.off()

################# Filtering out poor performing probes #####################
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
head(detP)

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

# exclude cross reactive probes 
# Download cross-reactive probes: 
# https://github.com/sirselim/illumina450k_filtering/blob/master/48639-non-specific-probes-Illumina450k.csv
xReactiveProbes <- read.csv("docs/48639-non-specific-probes-Illumina450k.csv", 
                            stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt

# Plot MDS using filtered and normalized data
pdf("MDS_PC1_PC2_filtered.pdf", height = 7, width = 10)
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        dim = c(1,2),
        col=pal[factor(targets$Trait)], cex=0.8)
legend("top", legend=levels(factor(targets$Trait)), text.col=pal,
       cex=0.65, bg="white")
dev.off()

# Examine higher dimensions to look at other sources of variation
pdf("MDS_PC1_PC3_filtered.pdf", height = 7, width = 10)
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Trait)], dim = c(1,3))
legend("top", legend=levels(factor(targets$Trait)), text.col=pal,
       bg="white", cex=0.7)
dev.off()

# The separation between RA and Controls is found at PC3 and PC4
pdf("MDS_PC3_PC4_filtered.pdf", height = 7, width = 10)
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Trait)], dim = c(3,4))
legend("top", legend=levels(factor(targets$Trait)), text.col=pal,
       bg="white", cex=0.7)
dev.off()

save(mSetSqFlt, file = "GSE42861_meth_norm_filt.RData")

############## Calculate M-values and beta-values ##########################
## M-values are suitable for differential methylation analysis, while beta 
## values are a direct reflection of methylation percentages and are better suited
## for the presentation of the results
# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])

## Calculate beta values
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])

# Plot M- and beta values distributions
pdf("M_Beta_densities.pdf", height = 7, width = 10)
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Trait, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Trait)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Trait, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Trait)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

save(mVals, file = "M_values.RData")
save(bVals, file = "Beta_values.RData")
write.csv(mVals, file = "M_values.csv")
write.csv(bVals, file = "Beta_values.csv")
