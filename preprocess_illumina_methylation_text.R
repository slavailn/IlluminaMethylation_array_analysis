library(wateRmelon)
library(reshape2)
library(lattice)
library(lumi)
library(pheatmap)
library(RColorBrewer)
library(PCAtools)

# Analyze methylation array without raw data in idat files

setwd("working_directory")
list.files()

# Read methylation data in tab delimited file
meth <- read.table("raw_files/GSE134429_signal_intensities.txt",
                   header = T, sep = "\t", quote = "\"",
                   fill = T, row.names = 1)
head(meth)[,1:4]

# Create tables of methylated and unmethylated intensities
# Methylated probe intensities
me <- meth[,grep("_Methylated_signal", names(meth))]
head(me)
dim(me)

# Unmethylated probe intensities
un <- meth[,grep("_Unmethylated_signal", names(meth))]
head(un)
dim(un)

# Create table of p-values
det_pval <- meth[,grep("_Detection_Pval", names(meth))]
head(det_pval)
dim(det_pval)
write.csv(det_pval, file = "GSE134429_det_p_values.csv")

# Filter the samples and probes based on detection p-values
pdf("mean_pvalues_arrays.pdf")
par(mfrow=c(1,1))
par(mar = c(8,8,3,3))
barplot(colMeans(det_pval), col="lightgray", las=2, 
        cex.names=0.8, ylab="",
        main = "Mean p-values across the arrays")
title(ylab="Mean detection p-values", mgp=c(5,1,0), cex.lab=1.2)
dev.off()

# Examine probe intensity distribution
me_melt <- melt(me) 
head(me_melt)
names(me_melt) <- c("Sample_name", "Intensity")

un_melt <- melt(un)
head(un_melt)
names(un_melt) <- c("Sample_name", "Intensity")

# Create density plots
pdf("Methylated_probes_raw_intensity_density.pdf")
densityplot(~Intensity, data=me_melt, 
            groups = Sample_name,
            plot.points = FALSE, ref = TRUE, 
            auto.key = F)

dev.off()

pdf("Unmethylated_probes_raw_intensity_density.pdf")
densityplot(~Intensity, data=un_melt, 
            groups = Sample_name,
            plot.points = FALSE, ref = TRUE, 
            auto.key = F)
dev.off()

# Create violin plots
pdf("Methylated_probes_raw_violin.pdf", width = 8)
bwplot(Intensity ~ Sample_name,  data = me_melt,
       panel = panel.violin, # specify panel argument to make violin plot
       xlab = "Sample_name", ylab = "Intensity",
       scales=list(x=list(rot=90)))
dev.off()

pdf("Unmethylated_probes_raw_violin.pdf", width = 8)
bwplot(Intensity ~ Sample_name,  data = un_melt,
       panel = panel.violin, # specify panel argument to make violin plot
       xlab = "Sample_name", ylab = "Intensity",
       scales=list(x=list(rot=90)))
dev.off()

# Normalize array using dasen method
# first we need to establish which probes are 'I' and which are type 'II'
manifest <- read.csv("raw_files/GPL21145_MethylationEPIC_15073387_v-1-0.csv",
                     header = T, skip = 7)
head(manifest)
manifest <- manifest[,c(1, 7)]
head(manifest)

# Merge methylation data with probe types and re-create 
# methylated and unmethylated matrices.
# We need this step to make sure that probe ids in methylated data 
# match Infinium Probe Types (I or II)
meth_merged <- merge(meth, manifest, by.x = 0, by.y = "IlmnID")
head(meth_merged)[,1:4]
sum(duplicated(meth_merged))

rownames(meth_merged) <- meth_merged$Row.names
head(meth_merged)[,1:4]

# Create vector of probe types
probe_types <- meth_merged$Infinium_Design_Type
head(probe_types)

# Re-create methylated and unmethylated matrices
meth_mat <- meth_merged[,grep("_Methylated_signal", names(meth_merged))]
head(meth_mat)[,1:4]
meth_mat <- as.matrix(meth_mat)
head(meth_mat)

unmeth_mat <- meth_merged[,grep("_Unmethylated_signal", names(meth_merged))]
head(unmeth_mat)[,1:4]
unmeth_mat <- as.matrix(unmeth_mat)
head(unmeth_mat)[,1:4]

# Normalize probe intensities
beta_norm <- nanes(meth_mat, unmeth_mat, onetwo = probe_types, fudge = 100)
head(beta_norm)[,1:10]
colnames(beta_norm) <- gsub("_Methylated_signal", "", colnames(beta_norm))
head(beta_norm)[,1:10]
# nanes function generated some 0 values, that turn into Inf
# if converted to M-values. To prevent this, let's add a small
# value to every beta.
beta_norm <- beta_norm + 0.001
write.csv(beta_norm, file = "normalized_beta_values.csv")

# Visualize beta value distributions
beta_melt <- melt(beta_norm)
names(beta_melt) <- c("Probe_id", "Sample_name", "Intensity")
head(beta_melt)

# Density plot
pdf("Norm_beta_density.pdf")
densityplot(~Intensity, data=beta_melt, 
            groups = Sample_name,
            plot.points = FALSE, ref = TRUE, 
            auto.key = list(space = "right"))
dev.off()

# Create violin plot
pdf("Norm_beta_violin.pdf", width = 8)
bwplot(Intensity ~ Sample_name,  data = beta_melt,
       panel = panel.violin, # specify panel argument to make violin plot
       xlab = "Sample_name", ylab = "Intensity",
       scales=list(x=list(rot=90)))
dev.off()

# Convert beta to M-values
m_norm <- beta2m(beta_norm)
write.csv(m_norm, file = "normalized_M_values.csv")

# Visualize beta value distributions
m_melt <- melt(m_norm)
names(m_melt) <- c("Probe_id", "Sample_name", "Intensity")
head(m_melt)

# Density plot
pdf("Norm_m_values_density.pdf")
densityplot(~Intensity, data=beta_melt, 
            groups = Sample_name,
            plot.points = FALSE, ref = TRUE, 
            auto.key = F)
dev.off()

# Create violin plot
pdf("Norm_m_violin.pdf", width = 8)
bwplot(Intensity ~ Sample_name,  data = m_melt,
       panel = panel.violin, # specify panel argument to make violin plot
       xlab = "Sample_name", ylab = "Intensity",
       scales=list(x=list(rot=90)))
dev.off()

# ------------------------------------------------------------------------ #
# Create heatmap and cluster top 500 CpG sites with the highest variability
# Get metadata
meta <- read.csv("docs/GSE134429_850k_phenodata.csv", header = T)
head(meta)
sample_name <- gsub("monocytes.", "", meta$title)
sample_name

# Match sample metadata with column names in  beta values and M-values matrixes
meta <- meta[match(colnames(m_norm), sample_name),]
sample_name <- sample_name[match(colnames(m_norm), sample_name)]
meta$sample_name <- sample_name

# Check that sample names match
data.frame(meta$sample_name, colnames(m_norm), colnames(beta_norm))

# Create metadata data frame for future use
head(meta)
Trait <- gsub("HD", "Control", meta$donor.ch1) 
targets <- data.frame(row.names = meta$sample_name,
                      Trait = Trait,
                      Age = meta$age.ch1,
                      Gender = meta$Sex.ch1)
targets
write.csv(targets, file="targets.csv", row.names = T)

# ---------------------------------------------------- #
# Cluster data before and after correction
# Select top 500 CpG sites with highest Mean Absolute Deviations (MAD)
mads_sorted <- sort(rowMads(m_norm), decreasing = T)
mads_sorted <- mads_sorted[1:500]
head(mads_sorted)

# Extract the matrix with top 500 most variable values
top500 <- m_norm[which(rownames(m_norm) %in% names(mads_sorted)),]
head(top500)[,1:4]

# Draw the heatmap based on top500 most variable CpGs
annot <- targets
# Cluster genes and samples, plot heatmap
pdf("GSE134429_heatmap_top500.pdf", width = 7, height = 10)
pheatmap(top500, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         scale="row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = annot, 
         show_rownames = F)
dev.off()

# Explore principal components
# -------------------------- PCA ---------------------------------------- #
# Plot PCA after batch correction
p <- pca(m_norm, metadata = annot, removeVar = 0.9)

# Create biplot and draw stat ellipses at 95% CI around groups
biplot(p,
       colby = 'Trait', 
       colkey = c('Control' = 'forestgreen', 'RA' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', 
       legendLabSize = 8, 
       legendIconSize = 4)
ggsave("GSE134429_PCA_biplot_after_correct.pdf", device = "pdf", units = "in", 
       width = 6, height = 6)
