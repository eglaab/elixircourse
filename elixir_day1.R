#
# R code for the data quality control and differential expression meta-analysis in the ELIXIR course "Statistics & Machine Learning"
#

#
# Installation instructions:
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

update.packages(ask = FALSE, dependencies = c('Suggests'))

# load annotation package for gene ID conversion

# for old R version:
# source("http://bioconductor.org/biocLite.R")
# biocLite("hgu133a.db")

if(!require('hgu133a.db'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")
	BiocManager::install("hgu133a.db", suppressUpdates=TRUE, ask = FALSE)
	require('hgu133a.db')
}

# load R-packages for quality control

if(!require('arrayQualityMetrics'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("arrayQualityMetrics", suppressUpdates=TRUE, ask = FALSE)
  install.packages("gridSVG")
  # install.packages("https://cran.r-project.org/src/contrib/Archive/gridSVG/gridSVG_1.4-3.tar.gz", repos=NULL)
  require('arrayQualityMetrics')
}
if(!require('Biobase'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("Biobase", suppressUpdates=TRUE, ask = FALSE)
  require('Biobase')
}

# load R-packages for power calculation
if(!require('impute'))
{
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 BiocManager::install("impute", suppressUpdates=TRUE, ask = FALSE)
 require('impute')
}

if(!require('samr'))
{
  install.packages("samr")
  require('samr')
}
options(error=NULL)

# load R-package for Variance stabilizing normalization
if(!require('vsn'))
{
	if (!requireNamespace("vsn", quietly = TRUE))
	    install.packages("BiocManager")
	BiocManager::install("vsn", suppressUpdates=TRUE, ask = FALSE)
	require('vsn')
}

# load R-package for meta-analysis
if(!require('metaMA'))
{
	install.packages('metaMA')
	require('metaMA')
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
setwd('/Users/set/your/current/working/directory/here')

# format for Windows
setwd('C:/set/your/current/working/directory/here')


#
# 1) Dataset: GEO-ID: GSE20295, Y. Zhang et al., Am J Med Genet B Neuropsychiatr Genet, 2005 multiple brain regions, post mortem, PD (40), healthy (53)
#    Array platform: Affymetrix HG-U133A
#

# Download the data into the current working directory

# for Windows - manually via the web-browser using this url:
#
#ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20295/matrix/GSE20295_series_matrix.txt.gz
#

# for Mac/Linux - automatically via R command line
system('wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20295/matrix/GSE20295_series_matrix.txt.gz')

# Read the data into R
zhangdatgeo = read.table(gzfile("GSE20295_series_matrix.txt.gz"), header=T, comment.char="!", sep="\t")

# Use the labels in the first column as row names
zhangdat = zhangdatgeo[,2:ncol(zhangdatgeo)]
rownames(zhangdat) = zhangdatgeo[,1]

# Filter out tissue samples which are not from the midbrain / substantia nigra region
zhang_tissues = as.matrix(read.table(gzfile("GSE20295_series_matrix.txt.gz"), header=F, nrows=1, skip=39, sep="\t"))
zhang_tissues = zhang_tissues[2:length(zhang_tissues)]

table(zhang_tissues)
# zhang_tissues
#     Postmortem brain prefrontal cortex                Postmortem brain putamen Postmortem brain whole substantia nigra 
#                                     29                                      35                                      29 

# select only substantia nigra samples
nigra_ind = which(zhang_tissues == "Postmortem brain whole substantia nigra")

zhang_outcome = as.matrix(read.table(gzfile("GSE20295_series_matrix.txt.gz"), header=F, nrows=1, skip=41, sep="\t"))
zhang_outcome = zhang_outcome[2:length(zhang_outcome)]
table(zhang_outcome)
#zhang_outcome
#            disease state: control             disease state: Control disease state: Parkinson's disease 
#                                15                                 38                                 14 
# disease state: Parkinsons disease 
#                                26 

zhangfilt = zhangdat[,nigra_ind]
dim(zhangfilt)
#[1] 22283    29

zhang_outcomefilt = zhang_outcome[nigra_ind]
table(zhang_outcomefilt)
#zhang_outcomefilt
#           disease state: Control disease state: Parkinsons disease 
#                               18                                11

# convert Affymetrix probe set IDs to gene symbols
conv_ids <- mapIds(hgu133a.db, keys=as.character(rownames(zhangfilt)), c("SYMBOL"), keytype="PROBEID")
head(conv_ids)


#
# 2) Dataset GSE8397: L. B. Moran et al., Neurogenetics, 2006, SN + frontal gyrus, post mortem,	PD (29), healthy (18)
#    Array platform: Affymetrix HG-U133A
#

# Download the data into the current working directory

# for Windows - manually via the web-browser using this url:
#
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8397/matrix/GSE8397-GPL96_series_matrix.txt.gz
#

# for Mac/Linux - automatically via R command line
system('wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8397/matrix/GSE8397-GPL96_series_matrix.txt.gz')

# Read the data into R
morandatgeo = read.table(gzfile("GSE8397-GPL96_series_matrix.txt.gz"), header=T, comment.char="!", sep="\t")

# Use the labels in the first column as row names
morandat = morandatgeo[,2:ncol(morandatgeo)]
rownames(morandat) = morandatgeo[,1]

# Filter out tissue samples which are not from the midbrain / substantia nigra region
moran_tissues = as.matrix(read.table(gzfile("GSE8397-GPL96_series_matrix.txt.gz"), header=F, nrows=1, skip=36, sep="\t"))
moran_tissues = moran_tissues[2:length(moran_tissues)]

nigra_ind = grep("substantia nigra",moran_tissues)

moran_outcome = as.matrix(read.table(gzfile("GSE8397-GPL96_series_matrix.txt.gz"), header=F, nrows=1, skip=28, sep="\t"))
moran_outcome = moran_outcome[2:length(moran_outcome)]
moran_outcome[grep("control",moran_outcome)] = rep("control",length(grep("control",moran_outcome)))
moran_outcome[grep("Parkinson",moran_outcome)] = rep("parkinson",length(grep("Parkinson",moran_outcome)))

moranfilt = morandat[,nigra_ind]
dim(moranfilt)
#[1] 22283    39

moran_outcomefilt = moran_outcome[nigra_ind]
table(moran_outcomefilt)
#moran_outcomefilt
#  control Parkinson 
#       15        24

all(rownames(zhangfilt) == rownames(moranfilt))
# TRUE


# Quality control

# Create quality report for Zhang et al. dataset
minimalSet <- ExpressionSet(assayData=as.matrix(zhangfilt))
arrayQualityMetrics(expressionset = minimalSet, outdir = "Quality_Report_Zhang", force = TRUE, do.logtransform = FALSE)


# Create quality report for Moran et al. dataset
minimalSet <- ExpressionSet(assayData=as.matrix(moranfilt))
arrayQualityMetrics(expressionset = minimalSet, outdir = "Quality_Report_Moran", force = TRUE, do.logtransform = FALSE)


# remove all samples failing at least two quality tests

remove_samples_zhang = match(c("GSM606624","GSM606625","GSM606626"),colnames(zhangfilt))

zhangfilt2 = zhangfilt[,-remove_samples_zhang]
zhang_outcome_final = zhang_outcomefilt[-remove_samples_zhang]

# Moran dataset: no failing samples to remove
moran_outcome_final = moran_outcomefilt


# Manual outlier check

# detect outliers via hierarchical clustering

distmat = dist(t(zhangfilt2))

hcldat = hclust(distmat, method="average")

plot(hcldat)

# GSM508732 potential outlier


# detect outliers via PCoA

medianscale <- cmdscale(dist(t(zhangfilt2)), k = 2)

plot(medianscale[,1], medianscale[,2], col=rainbow(2)[match(zhang_outcome_final, unique(zhang_outcome_final))], pch=20, main="PCoA plot", labels=NULL, cex=2, cex.axis=0.1, tck=0, xlab="Dimension 1", ylab="Dimension 2")

# no apparent outliers



distmat = dist(t(moranfilt))

hcldat = hclust(distmat, method="average")

plot(hcldat)

# no apparent outliers


# detect outliers via PCoA

medianscale <- cmdscale(dist(t(moranfilt)), k = 2)

plot(medianscale[,1], medianscale[,2], col=rainbow(2)[match(moran_outcome_final, unique(moran_outcome_final))], pch=20, main="PCoA plot", labels=NULL, cex=2, cex.axis=0.1, tck=0, xlab="Dimension 1", ylab="Dimension 2")


#
# Data transformation using Variance stabilising normalization (VSN)
# 

# check for intensity-dependent variance
meanSdPlot(as.matrix(zhangfilt2))
# yes, variance dependence on average intensity -- apply VSN transformation

zhangvsn = exprs(vsn2(as.matrix(zhangfilt2)))

# verify the fit
meanSdPlot(zhangvsn)


# check for intensity-dependent variance: Moran dataset
meanSdPlot(as.matrix(moranfilt))
# yes, variance dependence on average intensity -- apply VSN transformation

moranvsn = exprs(vsn2(as.matrix(moranfilt)))

# verify the fit
meanSdPlot(moranvsn)


#
# Power calculation
#

#
# Zhang et al. data - power calculation
#

# outcome "y" for two unpaired classes must be numeric labels 1, 2
data = list(x=zhangvsn, y=ifelse(zhang_outcome_final=="disease state: Control",1,2), geneid=as.character(1:nrow(zhangvsn)),genenames=paste("g",as.character(1:nrow(zhangvsn)),sep=""), logged2=TRUE)


# run the simulation with 1000 permutations
samr.obj <- samr(data,  resp.type="Two class unpaired", nperms=1000, random.seed=1234)


# investigate the following sample sizes of interest: 10, 20, 30, 50
colnum = ncol(data$x)
sfactors = c(10/colnum, 20/colnum, 30/colnum, 50/colnum)

# set seed value for the random number generator
set.seed(1234)


# determine power to detect 1.5-fold changes
samr.assess15 <- samr.assess.samplesize(samr.obj, data, log2(1.5), samplesize.factors=sfactors)
samr.assess.samplesize.plot(samr.assess15)

# determine power to detect 1.1-fold changes
samr.assess11 <- samr.assess.samplesize(samr.obj, data, log2(1.1), samplesize.factors=sfactors)
samr.assess.samplesize.plot(samr.assess11)


#
# Moran et al. data - power calculation
#

# outcome "y" for two unpaired classes must be numeric labels 1, 2
data = list(x=moranvsn, y=ifelse(moran_outcome_final=="control",1,2), geneid=as.character(1:nrow(zhangvsn)),genenames=paste("g",as.character(1:nrow(zhangvsn)),sep=""), logged2=TRUE)


# run the simulation with 1000 permutations
samr.obj <- samr(data,  resp.type="Two class unpaired", nperms=1000, random.seed=1234)


# investigate the following sample sizes of interest: 10, 20, 30, 50
colnum = ncol(data$x)
sfactors = c(10/colnum, 20/colnum, 30/colnum, 50/colnum)

# set seed value for the random number generator
set.seed(1234)


# determine power to detect 1.5-fold changes
samr.assess15 <- samr.assess.samplesize(samr.obj, data, log2(1.5), samplesize.factors=sfactors)
samr.assess.samplesize.plot(samr.assess15)

# determine power to detect 1.1-fold changes
samr.assess11 <- samr.assess.samplesize(samr.obj, data, log2(1.1), samplesize.factors=sfactors)
samr.assess.samplesize.plot(samr.assess11)


#
# DEG Analysis of individual datasets
#

# Limma analysis of Zhang dataset
zhang_label = ifelse(zhang_outcome_final == "disease state: Control","control","parkinson")
design <- model.matrix(~ -1+factor(zhang_label))
colnames(design) <- unique(zhang_label)

# compute simple linear model fit to microarray data (not robust)
fit <- lmFit(zhangvsn, design)

contrast.matrix = makeContrasts(parkinson-control, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

eb <- eBayes(fit2)

# extract the ranking table and show the top-ranked genes
ttable_zhang <- topTable(eb, n = nrow(zhangfilt2)) 

head(ttable_zhang)


# Limma analysis of Moran dataset
design <- model.matrix(~ -1+factor(moran_outcome_final))
colnames(design) <- unique(moran_outcome_final)

# compute simple linear model fit to microarray data (not robust)
fit <- lmFit(moranvsn, design)

contrast.matrix = makeContrasts(parkinson-control, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

eb <- eBayes(fit2)

# extract the ranking table and show the top-ranked genes
ttable_moran <- topTable(eb, n = nrow(moranfilt)) 

head(ttable_moran)


#
# Meta-analysis
#

metarank = pvalcombination(esets=list(zhangvsn, moranvsn), classes=list(zhang_label, moran_outcome_final), moderated = "limma", BHth = 0.05)


metacomb <- cbind(rownames(zhangvsn)[metarank$Meta], metarank$TestStatistic[metarank$Meta])
metaord <- order(abs(as.numeric(metacomb[,2])), decreasing=TRUE)

logfcmat = cbind(ttable_zhang[match(rownames(zhangvsn), rownames(ttable_zhang)),]$logFC, ttable_moran[match(rownames(zhangvsn), rownames(ttable_moran)),]$logFC)
rownames(logfcmat)= rownames(zhangfilt2)

pmat = cbind(ttable_zhang[match(rownames(zhangvsn), rownames(ttable_zhang)),]$P, ttable_moran[match(rownames(zhangvsn), rownames(ttable_moran)),]$P)
rownames(pmat)= rownames(zhangvsn)


compresmeta = cbind(metacomb[metaord,1], conv_ids[match(metacomb[metaord,1],rownames(logfcmat))], logfcmat[match(metacomb[metaord,1],rownames(logfcmat)),], pmat[match(metacomb[metaord,1],rownames(pmat)),], metacomb[metaord,2])
colnames(compresmeta) = c("ID", "Gene Symbol", paste(c("Zhang","Moran"),"logFC"), paste(c("Zhang","Moran"),"P"), "Comb. Z")

head(compresmeta)


dim(compresmeta)


# save datasets to continue with pathway analysis tomorrow
save(moranvsn, moran_outcome_final, file="moran_preprocessed.Rdata")
save(zhangvsn, zhang_outcome_final, file="zhang_preprocessed.Rdata")

# For reproducibility: show and save information on all loaded R package versions
sessionInfo()
