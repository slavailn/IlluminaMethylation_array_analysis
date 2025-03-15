library(minfi)
library(DMRcate)
library(Gviz)
library(stringr)

## Detect differentially methylated regions using DMRcate
setwd("working_directory")
list.files()

# Load M- and beta-values
load("meth_data/M_values.RData")
load("meth_data/Beta_values.RData")

# load filtered and normalized minfi GenomicRatioSet
# that was previously normalized an
load("meth_data/GSE42861_meth_norm_filt.RData")
class(mSetSqFlt)

# Create design matrix 
targets <- pData(mSetSqFlt)

load("meth_data/SVA_results_vars.RData")
surrogate <- svobj$sv
surrogate <- surrogate[,c(1,2)]
head(surrogate)
colnames(surrogate) <- c("SV1", "SV2")
targets <- cbind(targets, surrogate)
head(targets)
design <- model.matrix(~Trait + Age + Gender + Smoking_status + 
                         SV1 + SV2, 
                       data=targets)
head(design)

# Annotate CpGs and run region_level differential analysis.
# DMRcate uses limma internally

# First we create annotation object and run limma on CpG sites
Annot <- cpg.annotate(object = mVals, datatype = "array", 
                      what = "M", analysis.type = "differential", 
                      design = design, contrasts = F, 
                      coef = "TraitRA",
                      arraytype = "450K")
str(Annot)

# Delineate regions
DMRs <- dmrcate(Annot, lambda = 1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges
save(results.ranges, file = "GSE42861_GRanges.RData")
write.csv(as.data.frame(results.ranges), 
          file = "GSE42861_DMRs.csv")


# Customize visualization
# indicate which genome is being used
gen <- "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 10
# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.5*(end-start))
maxbase <- end + (0.5*(end-start))

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)
# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom, fill="darkred")

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
ann450kOrd <- ann450kSub[order(ann450kSub$chr,ann450kSub$pos),]
head(ann450kOrd)
bValsOrd <- bVals[match(ann450kOrd$Name,rownames(bVals)),]
head(bValsOrd)

# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(ann450kOrd$chr),
                   ranges=IRanges(start=ann450kOrd$pos, end=ann450kOrd$pos),
                   strand=Rle(rep("*",nrow(ann450kOrd))),
                   betas=bValsOrd)

# extract data on CpGs in DMR
cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])

# methylation data track
pal <- brewer.pal(8,"Dark2")
methTrack <- DataTrack(range=cpgData, groups=targets$Trait,
                       genome = gen,
                       chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)

# Set up the track list and indicate the relative sizes of the different tracks. 
# Finally, draw the plot using the plotTracks function.
tracks <- list(iTrack, gTrack, methTrack, dmrTrack, rTrack)
sizes <- c(2,2,5,2,3) # set up the relative sizes of the tracks
pdf("DMR10.pdf", width = 10, height = 8)
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
dev.off()
