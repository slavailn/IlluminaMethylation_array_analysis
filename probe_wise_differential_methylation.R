library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(pheatmap)

## Probe-wise methylation analysis
## We will detect differentially methylated (DM) sites
## separately within each of the cell types as well as also a
## cross all cell types.

setwd("C:/BioProjects/RA_methylation_analysis/GSE111942/")
list.files()

# load filtered and normalized minfi GenomicRatioSet
# that was previously normalized and filtered
load("meth_data/GSEGSE111942_filtered.RData")
class(mSetSqFlt)

# Create design matrix 
targets <- pData(mSetSqFlt)


# Detect DM CpGs across between RA and Control
# No adjustment variables available, adjust by predicted sex
head(targets)
design <- model.matrix(~Trait, data=targets)
head(design)

# load the matrix of M-values
load("meth_data/M_values.RData")

# fit linear model
fit <- lmFit(mVals, design)
# Fit the contrasts
fit2 <- eBayes(fit)

# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
# Get DM sites
DMPs <- topTable(fit2, num=Inf, coef="TraitRA", 
                 genelist=ann450kSub)
head(DMPs)
sum(DMPs$adj.P.Val < 0.05)

## Write the results as csv
write.table(DMPs, file="GSE111942_RA_vs_Control_DMPs.csv", sep=",", 
            row.names=FALSE)
DMPs_sig <- subset(DMPs, adj.P.Val < 0.05)
write.table(DMPs_sig, file="GSE111942_RA_vs_Control_DMPs_sig_only.csv", 
            sep=",", row.names=FALSE)

# plot the top 10 most significantly differentially methylated CpGs
# to check if the results make sense
load("meth_data/Beta_values.RData")
pdf("top10_DM_all_cell_types.pdf", width = 10, height = 8)
par(mfrow=c(2,5))
sapply(rownames(DMPs)[1:10], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Trait, 
          ylab = "Beta values")
})
dev.off()

# Create a heatmap of significantly changed CpGs based
# on M_values
# Plot a heatmap of differentially expressed probes
head(DMPs)
sig_probes <- rownames(DMPs[which(DMPs$adj.P.Val < 0.05),])
sig_cpg <- mVals[rownames(mVals) %in% sig_probes,]

# # Select metadata corresponding to the contrast of interest
meta_select <- targets 
annot <- data.frame(row.names = targets$Sample_name,
                    Trait = meta_select$Trait)
head(annot)

# Clustering and heatmap for significant probes (adj. p-values < 0.05)
pdf("GSE87095_DMP_heatmap_significant.pdf", width = 7, height = 7)
pheatmap(sig_cpg, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         scale="row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = annot, 
         show_rownames = F,
         main = "Rheumatoid Arthritis vs Controls (Females only)")
dev.off()
