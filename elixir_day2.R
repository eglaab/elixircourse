#
# R code for the pathway, network and machine learning analyses in the ELIXIR course "Statistics & Machine Learning"
#
# Required input data in the working folder:
#
# - Affymetrix annotations for the HG-U133A array can be obtained via free registration on the Affymetrix website: http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu133
#
# - Pathway annotations from MSigDb can be obtained from http://software.broadinstitute.org/gsea/msigdb/ after free registration

#
# Package installations
# (Please make sure you have a recent version of R > 3.4 - otherwise install the current software as follows:
# - Install the current version of R (3.6) for your operating system from https://ftp.gwdg.de/pub/misc/cran/
# - Install the current version of R-Studio (1.2) from: https://www.rstudio.com/products/rstudio/download/
#

#
# Installation and loading of R-packages
#

## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org"
       options(repos=r)
})

#update.packages(ask = FALSE) # , dependencies = c('Suggests'))

if(!require('Biobase'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("Biobase", suppressUpdates=TRUE, ask = FALSE)
  require('Biobase')
}

# install R-package for pathway analysis
if(!require('clusterProfiler'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")

	BiocManager::install("clusterProfiler", suppressUpdates=TRUE, ask = FALSE)
	require('clusterProfiler')
}

if(!require('GSEABase'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")

	BiocManager::install("GSEABase", suppressUpdates=TRUE, ask = FALSE)
	require('GSEABase')
}

# install Limma package for statistical analyses
if(!require('limma'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")

	BiocManager::install("limma", suppressUpdates=TRUE, ask = FALSE)
	require('limma')
}

# install R-packages for clustering
if(!require('cluster'))
{
	install.packages('cluster')
	require('cluster')
}

if(!require('mclust'))
{
	install.packages('mclust')
	require('mclust')
}

# install R-packages for classification
if(!require('randomForest'))
{
	install.packages('randomForest')
	require('randomForest')
}

if(!require('e1071'))
{
	install.packages('e1071')
	require('e1071')
}

#
# End of package installations
#

# set the location of your working directory (note that there are differences between Windows & Mac concerning the use of back slash "\" vs. forward slash "/")

# format for Mac & Linux systems
# setwd('/set/your/current/directory')

# e.g.:
# setwd('/Users/enrico.glaab/Downloads')

# format for Windows
#setwd('C:/set/your/current/working/directory/here')



# load datasets from day 1
# (if you haven't saved the data from day 1, you can get the pre-processed data as Rdata-files from Moodle)
#

load(file="moran_preprocessed.Rdata") # moranvsn, moran_outcome_final
load(file="zhang_preprocessed.Rdata") # zhangvsn, zhang_outcome_final

# check identity of rownames between the two datasets
all(rownames(zhangvsn) == rownames(moranvsn))
#[1] TRUE


#
# Convert Gene IDs to Gene Symbols in order to be able to apply pathway analyses
#


#
# unzip Affymetrix annotation file (downloaded from Moodle into the working directory, see above):
#
# In Windows:
# - unzip the file "HG-U133A.na36.annot.csv.zip" manually
#
# In Mac/Linux:
# - use the following line of code:
# system('unzip HG-U133A.na36.annot.csv.zip')


# read annotations file (ignoring comments)
annot = read.csv("HG-U133A.na36.annot.csv", comment.char="#")
head(annot[,3:15])

# map probes to microarray rownames
mapids = match(rownames(zhangvsn), annot$Probe.Set.ID)

# check if all IDs were mapped successfully
any(is.na(mapids))
#[1] FALSE
# ok, no missing IDs

# extract gene symbols corresponding to microarray Probe IDs (take always the first symbol mapped)
gene_symbols = sapply( as.character(annot$Gene.Symbol[mapids]) , function(x) strsplit(x, " /// ")[[1]][1])

		
#
# Convert expression matrix with Affymetrix IDs to Gene Symbol matrix (if multiple probes match to a gene, take the max. average value probe as representative for the gene)
#

# Function to convert probe-based expression matrix to gene-based expression matrix
# Parameters:
#   matdat = input matrix with probe rownames,
#   mat_conv = vector with gene symbols corresponding to probe rownames (NA for missing conversions)
probe2genemat <- function(matdat, mat_conv)
{

	if(nrow(matdat) != length(mat_conv))
	{
	  stop("Matrix does not have the same number of rows as the gene name vector")
	}

	# take max expression vector (max = maximum of mean exp across samples), if gene occurs twice among probes
	unq <- unique(mat_conv)
	if(any(is.na(unq))){
		unq <- unq[-which(is.na(unq))]
	}
	mat <- matrix(0.0, nrow=length(unq), ncol=ncol(matdat))
	for(j in 1:nrow(mat))
	{
	  ind <- which(unq[j]==mat_conv)

	  # show conversion progress, every 1000 probes
	  if(j %% 1000 == 0){
	    print(j)
	  }

	  # 1-to-1 probe to gene symbol matching
	  if(length(ind) == 1)
	  {
	    mat[j,] = as.numeric(as.matrix(matdat[ind,]))
	  } else if(length(ind) > 1){

	    # multiple probes match to one gene symbol
	    curmat = matdat[ind,]

	    # compute average expression per row -> select row with max. avg. expression
	    avg = apply(curmat, 1, mean)
	    mat[j,] = as.numeric(as.matrix(matdat[ind[which.max(avg)],]))
	  }
	}
	rownames(mat) = unq

  return(mat)
}

# Run the conversion from probe matrix to gene matrix (Zhang data)
zhang_symb = probe2genemat(zhangvsn, gene_symbols)
# get the original column names
colnames(zhang_symb) = colnames(zhangvsn)

# show the dimensions of the new gene expression matrix
dim(zhang_symb)

# Run the conversion from probe matrix to gene matrix (Moran data)
moran_symb = probe2genemat(moranvsn, gene_symbols)
colnames(moran_symb) = colnames(moranvsn)

dim(moran_symb)



#
# Differential expression analysis at the gene level (instead of probe level)
#

# Zhang et al. DEG Analysis

zhang_label = ifelse(zhang_outcome_final == "disease state: Control","control","parkinson")
design <- model.matrix(~ -1+factor(zhang_label))
colnames(design) <- unique(zhang_label)

# compute simple linear model fit to microarray data (not robust)
fit <- lmFit(zhang_symb, design)

contrast.matrix = makeContrasts(parkinson-control, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

eb <- eBayes(fit2)

# extract the ranking table and show the top-ranked genes
ttable_zhang <- topTable(eb, n = nrow(zhang_symb))

head(ttable_zhang)

# no. of genes with adj. p-value below 0.05
length(which(ttable_zhang$adj.P.Val < 0.05))

zhang_degs = rownames(ttable_zhang)[which(ttable_zhang$adj.P.Val < 0.05)]



# Moran et al. DEG Analysis

moran_label = moran_outcome_final
design <- model.matrix(~ -1+factor(moran_label))
colnames(design) <- unique(moran_label)

# compute simple linear model fit to microarray data (not robust)
fit <- lmFit(moran_symb, design)

contrast.matrix = makeContrasts(parkinson-control, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

eb <- eBayes(fit2)

# extract the ranking table and show the top-ranked genes
ttable_moran <- topTable(eb, n = nrow(moran_symb))

head(ttable_moran)

# no. of genes with adj. p-value below 0.05
length(which(ttable_moran$adj.P.Val < 0.05))

moran_degs = rownames(ttable_moran)[which(ttable_moran$adj.P.Val < 0.05)]


#
# Pathway analyses
#


# Load pathway definitions from MSigDB:
# Decompress the file msigdb_pathway_annoations.zip obtained from Moodle (see above)
#

msigdb_go_pathways = read.gmt("c5.all.v6.2.symbols.gmt")
msigdb_kegg_pathways = read.gmt("c2.cp.kegg.v6.2.symbols.gmt")
msigdb_reactome_pathways = read.gmt("c2.cp.reactome.v6.2.symbols.gmt")
msigdb_biocarta_pathways = read.gmt("c2.cp.biocarta.v6.2.symbols.gmt")
msigdb_positional = read.gmt("c1.all.v6.2.symbols.gmt")

# Inspect top of the pathway annotation table for GO:
head(msigdb_go_pathways)


#
# Apply classical Fisher's Exact test (significance-of-overlap computation) - Zhang et al.
#

fisher_go_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_go_pathways)
head(fisher_go_zhang[,1:6])

fisher_kegg_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_kegg_pathways)
head(fisher_kegg_zhang[,1:6])

fisher_biocarta_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_biocarta_pathways)
head(fisher_biocarta_zhang[,1:6])
# no gene can be mapped

fisher_reactome_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_reactome_pathways)
head(fisher_reactome_zhang[,1:6])

fisher_positional_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_positional)
head(fisher_positional_zhang[,1:6])



#
# Apply classical Fisher's Exact test (significance-of-overlap computation) - Moran et al.
#

fisher_go_moran <- enricher(moran_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_go_pathways)
head(fisher_go_moran[,1:6])

fisher_kegg_moran <- enricher(moran_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_kegg_pathways)
head(fisher_kegg_moran[,1:6])

fisher_reactome_moran <- enricher(moran_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_reactome_pathways)
head(fisher_reactome_moran[,1:6])

fisher_positional_moran <- enricher(moran_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_positional)
head(fisher_positional_moran[,1:6])



#
# Apply GSEA - Zhang et al.
#

ranked_genelst = ttable_zhang$B
names(ranked_genelst) = rownames(ttable_zhang)

gsea_go_zhang = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
                     maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_go_pathways,
                     TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_go_zhang[,1:7])

gsea_kegg_zhang = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
                       maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_kegg_pathways,
                       TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_kegg_zhang[,1:7])

gsea_reactome_zhang = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
                           maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_reactome_pathways,
                           TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_reactome_zhang[,1:7])  

gsea_positional_zhang = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
                             maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_positional,
                             TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_positional_zhang[,1:7])


#
# Apply GSEA - Moran et al.
#

ranked_genelst = ttable_moran$B
names(ranked_genelst) = rownames(ttable_moran)

gsea_go_moran = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
                     maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_go_pathways,
                     TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_go_moran[,1:7])

gsea_kegg_moran = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
                       maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_kegg_pathways,
                       TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_kegg_moran[,1:7])

gsea_reactome_moran = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
                           maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_reactome_pathways,
                           TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_reactome_moran[,1:7])  

gsea_positional_moran = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
                             maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_positional,
                             TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_positional_moran[,1:7])


#
# Note: Since we used multiple pathway databases and multiple tests, we would normally still need to adjust the p-values across all databases and tests!
# (here omitted to save time)
#


#
# Network-based pathway analysis: EnrichNet
#

#
# - Goto www.enrichnet.org
# - Copy up top 100 genes to the EnrichNet web-interface (using the code below to move gene names to the clipboard)
# - Set the Identifier format to "HGNC symbol" and try out a few analyses
#

# Zhang et al.

# Copy to clipboard

# Mac version
clip <- pipe("pbcopy", "w")
write.table(zhang_degs, file=clip, sep = '\t', row.names = FALSE, col.names= FALSE, quote = FALSE)
close(clip)

# Windows version
write.table(zhang_degs, "clipboard", sep = '\t', row.names = FALSE, col.names= FALSE, quote = FALSE)


# Moran et al.

# Mac version
clip <- pipe("pbcopy", "w")
write.table(moran_degs[1:100], file=clip, sep = '\t', row.names = FALSE, col.names= FALSE, quote = FALSE)
close(clip)

# Windows version
write.table(moran_degs[1:100], "clipboard", sep = '\t', row.names = FALSE, col.names= FALSE, quote = FALSE)


#
# Cytoscape / Jepetto: Download Cytoscape from https://cytoscape.org/download.html and install,
# then start Cytoscape, go to "Apps" - "Apps Manager", search for "JEPETTO", install Plugin "JEPETTO"
#




# =============================================================================
# 
# Unsupervised analyses
#   Objectives:
#     - Check whether there are patterns in the microarray data.
#     - Check whether these patterns correspond to known sample annotations.
#
# =============================================================================

# Load libraries.
library("cluster")
library("mclust")

# Set a seed number for the random number generator to make the results
# reproducible.
set.seed(20221130)

# =============================================================================
#
# A. Prepare the data for the clustering analysis.
#
# =============================================================================

#' @title Perform variance filtering on a gene expression matrix.  
#' @description Function to perform variance filtering of gene expression 
#' matrix X to only retain the genes (rows) with the highest variance.
#' @param X The gene expression matrix.
#' @param filt_size The number of genes with highest variance to retain after
#' filtering. Default value to 1000.
#' @return The filtered matrix.
var_filter <- function(X, filt_size = 1000) {
  local_filt_size <- min(nrow(X), as.numeric(filt_size))
  variances       <- apply(X, 1, var)
  feat_ord        <- order(variances, decreasing = TRUE)
  X_filt          <- (X[feat_ord, ])[1:local_filt_size, ]
  return(X_filt)
}

# We check that the dimensions of the data matrices.
dim(zhangvsn)
dim(moranvsn)

# Filter the expression matrices to only retain the top 2,000 genes with the 
# highest variance.
zhang_filt <- var_filter(X = zhangvsn, filt_size = 2000)
moran_filt <- var_filter(X = moranvsn, filt_size = 2000)

# We check that the matrices have indeed been reduced to 2,000 genes.
dim(zhang_filt)
dim(moran_filt)

# Most algorithms will cluster rows (genes in our cases). We are more interested in
# clustering columns (samples) so we transpose the data.
zhang <- t(zhang_filt)
moran <- t(moran_filt)
rm(zhang_filt, moran_filt, zhangvsn, moranvsn)

# What can also be done:
#   - scaling the data (base::scale).
#   - more elegant feature selection / feature extraction (e.g., PCA).

# =============================================================================
#
# B. Perform k-means clustering.
#
# =============================================================================

# =============================================================================
#   Start with the Zhang et al. dataset.
# =============================================================================

# Compute the distance matrix, using the Euclidean distance.
zhang_distmat <- stats::dist(x = zhang, method = "euclidean")

# Two-group clustering (k = 2) with 100 random initializations to obtain 
# robust results.
zhang_kmn_k2 <- stats::kmeans(x = zhang, centers = 2, nstart = 100)

# Print clustering results (two clusters indexed 1 and 2).
zhang_kmn_k2$cluster
table(zhang_kmn_k2$cluster)

# PCA plot of the clustering results for visualization.
clusplot(zhang_distmat, zhang_kmn_k2$cluster, diss = TRUE,
         main = "Zhang et al. - PCA clustering plot")

# Now, three-group clustering (k = 3) with 100 random initializations.
zhang_kmn_k3 <- stats::kmeans(x = zhang, centers = 3, nstart = 100)

# Print clustering results (three clusters indexed 1, 2 and 3).
zhang_kmn_k3$cluster
table(zhang_kmn_k3$cluster)

# PCA plot of the clustering results for visualization.
clusplot(zhang_distmat, zhang_kmn_k3$cluster, diss = TRUE,
         main = "Zhang et al. - PCA clustering plot")

# =============================================================================
#   Continue with the Moran et al. dataset.
# =============================================================================

# Compute the distance matrix, using the Euclidean distance.
moran_distmat <- stats::dist(x = moran, method = "euclidean")

# Two-group clustering (k = 2) with 100 random initializations to obtain
# robust results.
moran_kmn_k2  <- stats::kmeans(x = moran, centers = 2, nstart = 100)

# Print clustering results (two clusters indexed 1 and 2).
moran_kmn_k2$cluster
table(moran_kmn_k2$cluster)

# PCA plot of the clustering results for visualization.
clusplot(moran_distmat, moran_kmn_k2$cluster, diss = TRUE,
         main = "Moran et al. - PCA clustering plot")

# Now, three-group clustering (k = 3) with 100 random initializations.
moran_kmn_k3 <- stats::kmeans(x = moran, centers = 3, nstart = 100)

# Show clustering results (two clusters indexed 1, 2 and 3).
moran_kmn_k3$cluster
table(moran_kmn_k3$cluster)

# PCA plot of the clustering results for visualization.
clusplot(moran_distmat, moran_kmn_k3$cluster, diss = TRUE,
         main = "Moran et al. - PCA clustering plot")

# =============================================================================
#   Questions.
# =============================================================================

# ? Does increasing the number of kmeans restarts to 1,000 change anything ?
# ? General impression from the PCA plots ?

# What can also be done:
#   - Trying more values for the parameter k.
#   - Trying variants of k-means such as k-medoids (cluster::pam or cluster::clara).

# =============================================================================
#
# C. Perform hierarchical clustering.
#
# =============================================================================

# =============================================================================
#   Start with the Zhang et al. dataset.
# =============================================================================

# Compute average linkage hierarchical clustering.
zhang_hcl <- stats::hclust(d = zhang_distmat, method = "average")

# Show cluster dendrogram.
plot(zhang_hcl)

# Turn the hierarchical clustering into a "flat" clustering by cutting
# the tree at different levels (2 and 3 again).
zhang_hcl_k2 <- stats::cutree(tree = zhang_hcl, k = 2)
zhang_hcl_k3 <- stats::cutree(tree = zhang_hcl, k = 3)

# Print clustering results.
table(zhang_hcl_k2)
table(zhang_hcl_k3)

# Plot the results on the PCA plot like for k-means.
clusplot(zhang_distmat, zhang_hcl_k2, diss = TRUE,
         main = "Zhang et al. - PCA clustering plot")
clusplot(zhang_distmat, zhang_hcl_k3, diss = TRUE,
         main = "Zhang et al. - PCA clustering plot")

# =============================================================================
#   Continue with the Moran et al. dataset.
# =============================================================================

# Compute average linkage hierarchical clustering.
moran_hcl <- stats::hclust(d = moran_distmat, method = "average")

# Show cluster dendrogram.
plot(moran_hcl)

# Turn the hierarchical clustering into a "flat" clustering by cutting
# the tree at different levels.
moran_hcl_k2 <- stats::cutree(tree = moran_hcl, k = 2)
moran_hcl_k3 <- stats::cutree(tree = moran_hcl, k = 3)

# Print clustering results.
table(moran_hcl_k2)
table(moran_hcl_k3)

# Plot the results on the PCA plot like for k-means.
clusplot(moran_distmat, moran_hcl_k2, diss = TRUE,
         main = "Moran et al. - PCA clustering plot")
clusplot(moran_distmat, moran_hcl_k3, diss = TRUE,
         main = "Moran et al. - PCA clustering plot")

# =============================================================================
#   Questions.
# =============================================================================

# ? General impression from the PCA plots wrt k-means PCA plots?
# ? What is the impact of changing the linkage method from "average" to "ward.D2" ?

# What can also be done:
#   - Trying more values for the parameter k.
#   - Trying other distances (e.g., pearson correlation based).

# =============================================================================
#
# D. Internal cluster validity assessment.
#
# =============================================================================

# Use the average silhouette width for internal cluster validity assessment.

# =============================================================================
#   Start with the Zhang et al. dataset.
# =============================================================================

# Average silhouette width (Zhang et al., k-means, k = 2).
mean((cluster::silhouette(zhang_kmn_k2$cluster, zhang_distmat))[, 3])

# Average silhouette width (Zhang et al., k-means, k = 3).
mean((cluster::silhouette(zhang_kmn_k3$cluster, zhang_distmat))[, 3])

# Average silhouette width (Zhang et al., hclust(avg), k = 2).
mean((cluster::silhouette(zhang_hcl_k2,         zhang_distmat))[, 3])

# Average silhouette width (Zhang et al., hclust(avg), k = 3).
mean((cluster::silhouette(zhang_hcl_k3,         zhang_distmat))[, 3])

# Average silhouette width (Zhang et al., hclust(wrd), k = 2).
mean((cluster::silhouette(zhang_hclw_k2,        zhang_distmat))[, 3])

# Average silhouette width (Zhang et al., hclust(wrd), k = 3).
mean((cluster::silhouette(zhang_hclw_k3,        zhang_distmat))[, 3])

# Plot of silhouette widths for the first case (k-means, k = 2).
plot(cluster::silhouette(zhang_kmn_k2$cluster, zhang_distmat),
     main = "Silhouette plot (k-means, k = 2)",
     col = c("darkred", "darkblue"))

# =============================================================================
#   Continue with the Moran et al. dataset.
# =============================================================================

# Average silhouette width (Moran et al., k-means, k = 2)
mean((cluster::silhouette(moran_kmn_k2$cluster, moran_distmat))[, 3])

# Average silhouette width (Moran et al., k-means, k = 3)
mean((cluster::silhouette(moran_kmn_k3$cluster, moran_distmat))[, 3])

# Average silhouette width (Moran et al., hclust(avg), k = 2)
mean((cluster::silhouette(moran_hcl_k2,         moran_distmat))[, 3])

# Average silhouette width (Moran et al., hclust(avg), k = 3)
mean((cluster::silhouette(moran_hcl_k3,         moran_distmat))[, 3])

# Plot of silhouette widths for the first case (k-means, k = 2).
plot(cluster::silhouette(moran_kmn_k2$cluster, moran_distmat),
     main = "Silhouette plot (k-means, k = 2)",
     col = c("darkred", "darkblue"))

# =============================================================================
#   Questions
# =============================================================================

# ? What does the silhouette width values tell us ?
# ? Which clustering for each dataset seems to be the most suitable ?
# ? Does cutting the hclust tree to get 5 clusters improve the silhouette width ?

# What can also be done:
#   - Use other metrics than silhouette (see fpc::cluster.stats)
#   - Use other clustering algorithms (e.g., gmm, DBscan).

# =============================================================================
#
# E. External cluster validity assessment.
#
# =============================================================================

# Use the adjusted rand index for external cluster validity assessment
# effectively comparing how the results agree with the disease status 
# associated with the samples.

# =============================================================================
#   Start with the Zhang et al. dataset.
# =============================================================================

# Adjusted rand index (Zhang et al., k-means, k = 2).
mclust::adjustedRandIndex(zhang_kmn_k2$cluster, zhang_outcome_final)

# Adjusted rand index (Zhang et al., k-means, k = 3).
mclust::adjustedRandIndex(zhang_kmn_k3$cluster, zhang_outcome_final)

# Adjusted rand index (Zhang et al., hclust(avg), k = 2).
mclust::adjustedRandIndex(zhang_hcl_k2,         zhang_outcome_final)

# Adjusted rand index (Zhang et al., hclust(avg), k = 3).
mclust::adjustedRandIndex(zhang_hcl_k3,         zhang_outcome_final)

# Adjusted rand index (Zhang et al., hclust(avg), k = 5).
mclust::adjustedRandIndex(zhang_hcl_k5,         zhang_outcome_final)

# Adjusted rand index (Zhang et al., hclust(wrd), k = 5).
mclust::adjustedRandIndex(zhang_hclw_k5,        zhang_outcome_final)

# =============================================================================
#   Continue with the Moran et al. dataset.
# =============================================================================

# Adjusted rand index (Moran et al., k-means, k = 2).
mclust::adjustedRandIndex(moran_kmn_k2$cluster, moran_outcome_final)

# Adjusted rand index (Moran et al., k-means, k = 3).
mclust::adjustedRandIndex(moran_kmn_k3$cluster, moran_outcome_final)

# Adjusted rand index (Moran et al., hclust(avg), k = 2).
mclust::adjustedRandIndex(moran_hcl_k2,         moran_outcome_final)

# Adjusted rand index (Moran et al., hclust(avg), k = 3).
mclust::adjustedRandIndex(moran_hcl_k3,         moran_outcome_final)

# Adjusted rand index (Moran et al., hclust(avg), k = 5).
mclust::adjustedRandIndex(moran_hcl_k5,         moran_outcome_final)

# Adjusted rand index (Moran et al., hclust(wrd), k = 2).
mclust::adjustedRandIndex(moran_hclw_k2,        moran_outcome_final)

# Adjusted rand index (Moran et al., hclust(wrd), k = 3).
mclust::adjustedRandIndex(moran_hclw_k3,        moran_outcome_final)

# Adjusted rand index (Moran et al., hclust(wrd), k = 5).
mclust::adjustedRandIndex(moran_hclw_k5,        moran_outcome_final)

# =============================================================================
#   Questions.
# =============================================================================

# ? What does the adjusted rand index values tell us ?
# ? Any agreement with the silhouette widths values ?

# What can also be done:
#   - Use other metrics than ARI (see fpc::cluster.stats)



# =============================================================================
# Sample classification analyses
# =============================================================================


# Evaluation functions

# sensitivity
sensitivity <- function(tp, fn){
  return(tp/(tp+fn))
}

# specificity
specificity <- function(tn, fp){
  return(tn/(tn+fp))
}

# Matthew's correlation coefficient (=MCC)
corcoeff <- function(tp, tn, fp, fn){
	return(  ((tp*tn)-(fp*fn))/(sqrt((tn+fn)*(tn+fp)*(tp+fn)*(tp+fp))))
}

# Huberty's proportional chance criterion - p-value calculation for classification problems
ppc <- function(tp, fp, tn, fn)
{

		total <- tp+fp+fn+tn

		c_pro <- ((tp+fn)/total)*((tp+fp)/total) + ((tn+fp)/total) *((tn+fn)/total)

		acc <- (tp+tn)/total

		cat('\nc_pro:',c_pro,'acc:',acc,'\n')


		pval <- NULL
		if(c_pro > acc)
		{
		 pval <- 1
		} else if(c_pro != 1)
		{
			z <- (acc-c_pro)/sqrt(c_pro*(1-c_pro)/total)

			pval <- pnorm(-abs(z))
			cat('\np-value: ',pval,'\n')
			cat('\np-value (rounded): ',format(pval,digits=2),'\n')

		} else {
		  pval <- 1
		}

  return (pval)
}


# make the outcome variables numeric

zhang_numout = ifelse(zhang_outcome_final=="disease state: Control",0,1)
moran_numout = ifelse(moran_outcome_final=="control",0,1)

set.seed(1234)

require('randomForest')

# Build Random Forest sample classification model for Zhang et al. data using 250 decision trees
rfmod_zhang = randomForest(t(zhangvsn), factor(zhang_numout), ntree=250, keep.forest=TRUE)

# show model evluation based on out-of-bag samples
rfmod_zhang

# compute performance statistics
sensitivity(rfmod_zhang$confusion[2,2], rfmod_zhang$confusion[2,1])
specificity(rfmod_zhang$confusion[1,1], rfmod_zhang$confusion[1,2])
corcoeff(rfmod_zhang$confusion[2,2], rfmod_zhang$confusion[1,1], rfmod_zhang$confusion[1,2], rfmod_zhang$confusion[2,1])
ppc(rfmod_zhang$confusion[2,2], rfmod_zhang$confusion[1,2], rfmod_zhang$confusion[1,1], rfmod_zhang$confusion[2,1])

# which variables were most informative for the prediction (multivariate feature selction - see day 1 lecture):
head(rfmod_zhang$importance[order(rfmod_zhang$importance, decreasing=T),])


# Random Forest model for Moran et al. data using 250 trees
rfmod_moran = randomForest(t(moranvsn), factor(moran_numout), ntree=250, keep.forest=TRUE)

# show model evluation based on out-of-bag samples
rfmod_moran

# compute performance statistics
sensitivity(rfmod_moran$confusion[2,2], rfmod_moran$confusion[2,1])
specificity(rfmod_moran$confusion[1,1], rfmod_moran$confusion[1,2])
corcoeff(rfmod_moran$confusion[2,2], rfmod_moran$confusion[1,1], rfmod_moran$confusion[1,2], rfmod_moran$confusion[2,1])
ppc(rfmod_moran$confusion[2,2], rfmod_moran$confusion[1,2], rfmod_moran$confusion[1,1], rfmod_moran$confusion[2,1])


# which variables were most informative for the prediction (multivariate feature selction - see day 1 lecture):
head(rfmod_moran$importance[order(rfmod_moran$importance, decreasing=T),])



# Support vector machine classification

# Run linear SVM - evaluate using leave-one-out cross-validation (Zhang et a.)
svmmod = svm(t(zhangvsn), factor(zhang_numout), kernel="linear", cross=ncol(zhangvsn))

# show the cross-validated accuracy (as percentage)
svmmod$tot.accuracy

# confusion matrix
conf_zhang = table(zhang_numout, ifelse(svmmod$accuracies==100,zhang_numout,1-zhang_numout))
conf_zhang

# compute performance statistics
sensitivity(conf_zhang[2,2], conf_zhang[2,1])
specificity(conf_zhang[1,1], conf_zhang[1,2])
corcoeff(conf_zhang[2,2], conf_zhang[1,1], conf_zhang[1,2], conf_zhang[2,1])
ppc(conf_zhang[2,2], conf_zhang[1,2], conf_zhang[1,1], conf_zhang[2,1])


# Run linear SVM - evaluate using leave-one-out cross-validation (Moran et a.)
svmmod = svm(t(moranvsn), factor(moran_numout), kernel="linear", cross=ncol(moranvsn))

# show the cross-validated accuracy (as percentage)
svmmod$tot.accuracy

# confusion matrix
conf_moran = table(moran_numout, ifelse(svmmod$accuracies==100,moran_numout,1-moran_numout))
conf_moran

sensitivity(conf_moran[2,2], conf_moran[2,1])
specificity(conf_moran[1,1], conf_moran[1,2])
corcoeff(conf_moran[2,2], conf_moran[1,1], conf_moran[1,2], conf_moran[2,1])
ppc(conf_moran[2,2], conf_moran[1,2], conf_moran[1,1], conf_moran[2,1])


#
# Other online resources for machine learning analysis of omics data:
# ArrayMining: www.arraymining.net
# PathVar: www.pathvar.embl.de
#

# optionally, save the session
save.image()
# reload the session data
#load(".RData")
		      
# For reproducibility: show and save information on all loaded R package versions
sessionInfo()

