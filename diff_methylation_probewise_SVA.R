library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(pheatmap)
library(sva)
library(RColorBrewer)

## Probe-wise methylation analysis
## We will detect differentially methylated (DM) sites
## separately within each of the cell types as well as also a
## cross all cell types.
setwd(<working_directory>)
list.files()

# Get annotation
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# load filtered and normalized minfi GenomicRatioSet
# that was previously normalized and filtered
load("meth_data/GSE42861_meth_norm_filt.RData")
class(mSetSqFlt)

# Create design matrix 
targets <- pData(mSetSqFlt)

# Use SVA to identify latent sources of variance
# Create null model matrix with adjustment variables only
mod0 <- model.matrix(~ Gender + Age + Smoking_status, 
                     data = targets)

# Create full model matrix that includes variable of interest and
# adjustment variables
mod <- model.matrix(~ Trait + Gender + Age + Smoking_status, 
                    data = targets)

# load the matrix of M-values
load("meth_data/M_values.RData")
svobj <- sva(mVals, mod, mod0)
save(svobj, file = "SVA_results_vars.RData")
head(svobj$sv)

# Detect DM CpGs across between RA and Control after SVA correction
surrogate <- svobj$sv[,c(1,2)]
colnames(surrogate) <- c("SV1", "SV2")
head(surrogate)
targets <- cbind(targets, surrogate) 
head(targets)

# Create model matrix with adjustment and surrogate variables
design <- model.matrix(~ Trait + Age + Gender + Smoking_status + SV1 + SV2, 
                        data = targets)
head(design)

# fit linear model
fit <- lmFit(mVals, design)

# Fit the contrasts
fit2 <- eBayes(fit)

# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[,c(1:4,12:19,24:ncol(ann450k))]
head(ann450kSub)

# Get DM sites
DMPs <- topTable(fit2, num=Inf, coef="TraitRA")
head(DMPs)
sum(DMPs$adj.P.Val < 0.05)

## Write the results as csv
write.csv(DMPs, file="RA_vs_Control_DMPs_corrected.csv")
DMPs_sig <- subset(DMPs, adj.P.Val < 0.05)
write.table(DMPs_sig, file="RA_vs_Control_DMPs_sig_only.csv", 
            sep=",", row.names=FALSE)

# Read in the results
DMPs <- read.csv("results_corrected/sites/RA_vs_Control_DMPs_corrected.csv", 
                 row.names = 1)
head(DMPs)
DMPs_annot <- merge(DMPs, ann450kSub, by = 0)
write.csv(DMPs_annot, file = "RA_vs_Control_DMPs_corrected_annot.csv")

# Select significant CpGs with annotations and write it out into a file
DMPs_annot_sig <- DMPs_annot[DMPs_annot$adj.P.Val < 0.05,]
dim(DMPs_annot_sig)
write.csv(DMPs_annot, file = "RA_vs_Control_DMPs_sig_corrected_annot.csv")

# plot the top 10 most significantly differentially methylated CpGs
# to check if the results make sense
load("meth_data/Beta_values.RData")
pdf("top10_DM_corrected.pdf", width = 10, height = 8)
par(mfrow=c(2,5))
sapply(rownames(DMPs)[1:10], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Trait, 
          ylab = "Beta values")
})
dev.off()

# Create a heatmap of top 100 DMPs by p-value
# on M_values
# Plot a heatmap of differentially expressed probes
load("meth_data/M_values.RData")
head(mVals)[,1:4]
DMPs_annot <- DMPs_annot[order(DMPs_annot$adj.P.Val),]

# Select top 100 most sigbificant probes
sig_probes <- DMPs_annot$Row.names[1:100]
sig_probes
sig_cpg <- mVals[rownames(mVals) %in% sig_probes,]
head(sig_cpg)

 # Select metadata corresponding to the contrast of interest
meta_select <- targets 
annot <- data.frame(row.names = rownames(meta_select),
                    Trait = meta_select$Trait,
                    Gender = meta_select$Gender,
                    Smoking_status = meta_select$Smoking_status
                    )
head(annot)

# Clustering and heatmap for significant probes (adj. p-values < 0.05)
pdf("GSE42861_DMP_heatmap_significant.pdf", width = 7, height = 7)
pheatmap(sig_cpg, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         scale="row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = annot, 
         show_rownames = F,
         main = "GSE42861 (PBMC); top 100 by adj. p-value")
dev.off()
