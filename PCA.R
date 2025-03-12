library(ggplot2)
library(PCAtools)
library(limma)

# Conduct principal components analysis using PCAtools
setwd(<working_directory>)
list.files()

# -------------------------- PCA -------------------------------------------- #
# Run principal components analysis
# Load normalized M-values and sample meta-data
load("meth_data/GSE42861_meth_norm_filt.RData")
class(mSetSqFlt)
# Create design matrix 
targets <- pData(mSetSqFlt)
metadata <- targets 

load("meth_data/M_values.RData")

# Calculate principal components
p <- pca(mVals, metadata = metadata, removeVar = 0.9)

# Create biplot
biplot(p, showLoadings = TRUE, lab = NULL)

# Create biplot and draw stat ellipses at 95% CI around groups
biplot(p,
       colby = 'Trait', colkey = c('Control' = 'forestgreen', 'RA' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)
ggsave("GSE42861_PCA_biplot_trait.pdf", device = "pdf", units = "in", 
       width = 6, height = 6)

biplot(p,
       colby = 'Gender', colkey = c('m' = 'red', 'f' = 'blue'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)
ggsave("GSE42861_PCA_biplot_gender.pdf", device = "pdf", units = "in", 
       width = 6, height = 6)

# Create pairs plot
pairsplot(p)
ggsave("GSE42861_pairs_plot.pdf", device = "pdf", units = "in", 
       width = 6, height = 6)


# Plot loadings, fancy
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)
ggsave("GSE42861_loadings_plot_fancy.pdf", device = "pdf", units = "in", 
       width = 6, height = 7)

# Create eigencorrelation plot
pdf("GSE42861_eigencor_plot.pdf", width = 7, height = 7) 
eigencorplot(p,
             metavars = c('Trait','Gender','Age', 'Smoking_status'))
dev.off()

elbow <- findElbowPoint(p$variance)
elbow # 7

# Create screeplot
screeplot(p, components = getComponents(p, 1:10),
          vline = c(elbow)) + 
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))
ggsave("GSE42861_screeplot.pdf", device = "pdf", units = "in", 
       width = 7, height = 7)
