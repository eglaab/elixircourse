

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
# - unzip the file "HG-U133A.na35.annot.csv.zip" manually
#
# In Mac/Linux:
# - use the following line of code:
system('unzip HG-U133A.na36.annot.csv.zip')


# read annotations file (ignoring comments)
annot = read.csv("HG-U133A.na35.annot.csv", comment.char="#")
head(annot)

# map probes to microarray rownames
mapids = match(rownames(zhangvsn), annot$Probe.Set.ID)

# check if all IDs were mapped successfully
any(is.na(mapids))
#[1] FALSE
# ok, no missing IDs

# extract gene symbols corresponding to microarray Probe IDs (take always the first symbol mapped)
mapped_symbols = sapply( as.character(annot$Gene.Symbol[mapids]) , function(x) strsplit(x, " /// ")[[1]][1])



# get annotations to convert Affymetrix gene IDs to official HGNC gene symbols (from hgu133a.db installed above)
# extract the gene symbol information
x <- hgu133aSYMBOL

# Get the probe identifiers that are mapped to a gene name
mapped_probes <- mappedkeys(x)

# Convert to a list
probes2symbol <- as.list(x[mapped_probes])

# check the top of the mapping table
head(probes2symbol)


# extract the gene description information
x <- hgu133aGENENAME

# Get the probe identifiers that are mapped to a gene name
mapped_probes <- mappedkeys(x)

# Convert to a list
probes2desc <- as.list(x[mapped_probes])

# check the top of the gene description mapping table
head(probes2desc)


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
zhang_symb = probe2genemat(zhangvsn, mapped_symbols)
# get the original column names
colnames(zhang_symb) = colnames(zhangvsn)

# show the dimensions of the new gene expression matrix
dim(zhang_symb)
# 13238    26
# the matrix has less rows than the probe expression matrix, as expected

# Run the conversion from probe matrix to gene matrix (Moran data)
moran_symb = probe2genemat(moranvsn, mapped_symbols)
colnames(moran_symb) = colnames(moranvsn)

dim(moran_symb)
# 13238    39



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

fisher_go_zhang <- enricher(zhang_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_go_pathways)
head(fisher_go_zhang)

fisher_kegg_zhang <- enricher(zhang_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_kegg_pathways)
head(fisher_kegg_zhang)

fisher_biocarta_zhang <- enricher(zhang_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_biocarta_pathways)
head(fisher_biocarta_zhang)

fisher_reactome_zhang <- enricher(zhang_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_reactome_pathways)
head(fisher_reactome_zhang)

fisher_positional_zhang <- enricher(zhang_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_positional)
head(fisher_positional_zhang)


#
# Apply classical Fisher's Exact test (significance-of-overlap computation) - Moran et al.
#

fisher_go_moran <- enricher(moran_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_go_pathways)
head(fisher_go_moran)

fisher_kegg_moran <- enricher(moran_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_kegg_pathways)
head(fisher_kegg_moran)

fisher_biocarta_moran <- enricher(moran_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_biocarta_pathways)
head(fisher_biocarta_moran)

fisher_reactome_moran <- enricher(moran_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_reactome_pathways)
head(fisher_reactome_moran)

fisher_positional_moran <- enricher(moran_degs, universe = mapped_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_positional)
head(fisher_positional_moran)




#
# Apply GSEA - Zhang et al.
#

ranked_genelst = ttable_zhang$B
names(ranked_genelst) = rownames(ttable_zhang)

gsea_go_zhang = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
  maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_go_pathways,
  TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_go_zhang)

gsea_kegg_zhang = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
  maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_kegg_pathways,
  TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_kegg_zhang)

gsea_reactome_zhang = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
  maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_reactome_pathways,
  TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_reactome_zhang)

gsea_positional_zhang = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
  maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_positional,
  TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_positional_zhang)


#
# Apply GSEA - Moran et al.
#

ranked_genelst = ttable_moran$B
names(ranked_genelst) = rownames(ttable_moran)

gsea_go_moran = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
  maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_go_pathways,
  TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_go_moran)

gsea_kegg_moran = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
  maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_kegg_pathways,
  TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_kegg_moran)

gsea_reactome_moran = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
  maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_reactome_pathways,
  TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_reactome_moran)

gsea_positional_moran = GSEA(ranked_genelst, exponent = 1, nPerm = 1000, minGSSize = 10,
  maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_positional,
  TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")
head(gsea_positional_moran)


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
# Part III : Unsupervised analyses
#   Objective: perform sample clustering on the preprocessed datasets.
#
# =============================================================================

# =============================================================================
# A. Prepare the data for the clustering analysis.
# =============================================================================

#' @title Perform variance filtering on a gene expression matrix.
#' @description Function to perform variance filtering of gene expression
#' matrix X to only retain the genes (rows) with the highest variance.
#' @param X The gene expression matrix.
#' @param filt_size The number of genes with highest variance to retain after
#' filtering. Default value to 1000.
#' @return The filtered matrix.
var_filter = function(X, filt_size = 1000) {
	local_filt_size <- min(nrow(X), as.numeric(filt_size))
	variances <- apply(X, 1, var)
	feat_ord  <- order(variances, decreasing = TRUE)
	X_filt    <- (X[feat_ord, ])[1:local_filt_size, ]
	return(X_filt)
}

# Filter the expression matrices to only retain the top 2,000 genes with the
# highest variance.
zhang_filt <- var_filter(zhangvsn, filt_size = 2000)
moran_filt <- var_filter(moranvsn, filt_size = 2000)

# We carrefully check that the dimensions are reduced to 2,000 genes.
dim(zhang_filt)
#[1] 2000   26
dim(moran_filt)
#[1] 2000   39

# =============================================================================
# B. Perform k-means clustering.
# =============================================================================

# Set a seed number for the random number generator to make the results
# reproducible.
set.seed(1234)

# Start with the Zhang et al. dataset.

# Compute Euclidean distance matrix.
distmat_zhang <- dist(t(zhang_filt), method = "euclidean")

# Two-group clustering (k = 2) with 100 random initializations to obtain
# robust results.
kclust2_zhang <- kmeans(t(zhang_filt), centers = 2, nstart = 100)

# Show clustering results (two clusters indexed 1 and 2).
kclust2_zhang$cluster

# PCA plot of the clustering results for visualization.
clusplot(distmat_zhang, kclust2_zhang$cluster, diss = TRUE,
         main = "Zhang et al. - PCA clustering plot")

# Now, three-group clustering (k = 3) with 100 random initializations to
# obtain robust results.
kclust3_zhang <- kmeans(t(zhang_filt), centers = 3, nstart = 100)

# Show clustering results (three clusters indexed 1, 2 and 3).
kclust3_zhang$cluster

# PCA plot of the clustering results for visualization.
clusplot(distmat_zhang, kclust3_zhang$cluster, diss = TRUE,
         main = "Zhang et al. - PCA clustering plot")

# Continue with the Moran et al. dataset.

# Compute Euclidean distance matrix
distmat_moran <- dist(t(moran_filt), method = "euclidean")

# Two-group clustering (k = 2) with 100 random initializations to obtain
# robust results.
kclust2_moran <- kmeans(t(moran_filt), centers = 2, nstart = 100)

# Show clustering results (two clusters indexed 1 and 2):
kclust2_moran$cluster

# PCA plot of the clustering results for visualization.
clusplot(distmat_moran, kclust2_moran$cluster, diss = TRUE,
         main = "Moran et al. - PCA clustering plot")

# Now, three-group clustering (k = 3) with 100 random initializations
# to obtain robust results.
kclust3_moran <- kmeans(t(moran_filt), centers = 3, nstart = 100)

# Show clustering results (two clusters indexed 1, 2 and 3).
kclust3_moran$cluster

# PCA plot of the clustering results for visualization.
clusplot(distmat_moran, kclust3_moran$cluster, diss = TRUE,
         main = "Moran et al. - PCA clustering plot")

# =============================================================================
# C. Perform hierarchical clustering.
# =============================================================================

# Start with the Zhang et al. dataset.

# Compute average linkage hierarchical clustering.
hcldat_zhang <- hclust(distmat_zhang, method = "average")

# Show cluster dendrogram.
plot(hcldat_zhang)

# Turn the hierarchical clustering into a "flat" clustering by cutting
# the tree at different levels.
hcl2_zhang <- cutree(hcldat_zhang, k = 2)
hcl3_zhang <- cutree(hcldat_zhang, k = 3)

# Continue with the Moran et al. dataset.

# Compute average linkage hierarchical clustering.
hcldat_moran <- hclust(distmat_moran, method = "average")

# Show cluster dendrogram
plot(hcldat_moran)

# Turn the hierarchical clustering into a "flat" clustering by cutting
# the tree at different levels.
hcl2_moran <- cutree(hcldat_moran, k = 2)
hcl3_moran <- cutree(hcldat_moran, k = 3)

# =============================================================================
# D. Internal cluster validity assessment.
# =============================================================================

# Use the average silhouette width for internal cluster validity assessment.

# Start with the Zhang et al. dataset.

# Average silhouette width (Zhang et al., k-means, k = 2).
kclust2_zhang_score <- mean((silhouette(kclust2_zhang$cluster,
                                        distmat_zhang))[, 3])
kclust2_zhang_score

# Average silhouette width (Zhang et al., k-means, k = 3).
kclust3_zhang_score <- mean((silhouette(kclust3_zhang$cluster,
                                        distmat_zhang))[, 3])
kclust3_zhang_score

# Average silhouette width (Zhang et al., hclust, k = 2).
hcl2_zhang_score    <- mean((silhouette(hcl2_zhang,
                                        distmat_zhang))[, 3])
hcl2_zhang_score

# Average silhouette width (Zhang et al., hclust, k = 3).
hcl3_zhang_score    <- mean((silhouette(hcl3_zhang,
                                        distmat_zhang))[, 3])
hcl3_zhang_score

# Plot of silhouette widths.
plot(silhouette(kclust3_zhang$cluster, distmat_zhang),
     main = "Silhouette plot of best clustering result",
     col = c("darkred", "darkgreen", "darkblue"))

# Continue with the Moran et al. dataset.

# Average silhouette width (Moran et al., k-means, k = 2)
kclust2_moran_score <- mean((silhouette(kclust2_moran$cluster,
                                        distmat_moran))[, 3])
kclust2_moran_score

# Average silhouette width (Moran et al., k-means, k = 3)
kclust3_moran_score <- mean((silhouette(kclust3_moran$cluster,
                                        distmat_moran))[, 3])
kclust3_moran_score

# Average silhouette width (Moran et al., hclust, k = 2)
hcl2_moran_score    <- mean((silhouette(hcl2_moran,
                                        distmat_moran))[, 3])
hcl2_moran_score

# Average silhouette width (Moran et al., hclust, k = 3)
hcl3_moran_score    <- mean((silhouette(hcl3_moran,
                                        distmat_moran))[, 3])
hcl3_moran_score

# Plot of silhouette widths.
plot(silhouette(kclust2_moran$cluster, distmat_moran),
     main = "Silhouette plot of best clustering result",
     col = c("darkorange", "darkgrey"))

# =============================================================================
# D. External cluster validity assessment.
# =============================================================================

# Use the adjusted rand index for external cluster validity assessment.

# Start with the Zhang et al. dataset.

# Adjusted rand index (Zhang et al., k-means, k = 2).
adjustedRandIndex(kclust2_zhang$cluster, zhang_outcome_final)

# Adjusted rand index (Zhang et al., k-means, k = 3).
adjustedRandIndex(kclust3_zhang$cluster, zhang_outcome_final)

# Adjusted rand index (Zhang et al., hclust, k = 2).
adjustedRandIndex(hcl2_zhang, zhang_outcome_final)

# Adjusted rand index (Zhang et al., hclust, k = 3).
adjustedRandIndex(hcl3_zhang, zhang_outcome_final)

# Adjusted rand index (Moran et al., k-means, k = 2).
adjustedRandIndex(kclust2_moran$cluster, moran_outcome_final)

# Adjusted rand index (Moran et al., k-means, k = 3).
adjustedRandIndex(kclust3_moran$cluster, moran_outcome_final)

# Adjusted rand index (Moran et al., hclust, k = 2).
adjustedRandIndex(hcl2_moran, moran_outcome_final)

# Adjusted rand index (Moran et al., hclust, k = 3).
adjustedRandIndex(hcl3_moran, moran_outcome_final)


#
# Other online resources for machine learning analysis of omics data:
# ArrayMining: www.arraymining.net
# PathVar: www.pathvar.embl.de
#


# For reproducibility: show and save information on all loaded R package versions
sessionInfo()

