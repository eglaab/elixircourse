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

update.packages(ask = FALSE)

# load annotation package for gene ID conversion

# for old R version:
# source("http://bioconductor.org/biocLite.R")
# biocLite("hgu133a.db")

if(!require('hgu133a.db'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")
	BiocManager::install("hgu133a.db", suppressUpdates=TRUE)
	require('hgu133a.db')
}

# load R-packages for quality control

if(!require('arrayQualityMetrics'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("arrayQualityMetrics", suppressUpdates=TRUE)
  install.packages("gridSVG")
  # install.packages("https://cran.r-project.org/src/contrib/Archive/gridSVG/gridSVG_1.4-3.tar.gz", repos=NULL)
  require('arrayQualityMetrics')
}
if(!require('Biobase'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("Biobase", suppressUpdates=TRUE)
  require('Biobase')
}

# load R-packages for power calculation
if(!require('impute'))
{
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 BiocManager::install("impute", suppressUpdates=TRUE)
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
	BiocManager::install("vsn", suppressUpdates=TRUE)
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
	    install.packages("BiocManager", suppressUpdates=TRUE)

	BiocManager::install("clusterProfiler")
	require('clusterProfiler')
}

if(!require('GSEABase'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")

	BiocManager::install("GSEABase", suppressUpdates=TRUE)
	require('GSEABase')
}

# install Limma package for statistical analyses
if(!require('limma'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")

	BiocManager::install("limma", suppressUpdates=TRUE)
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
fit <- lmFit(zhangfilt2, design)

contrast.matrix = makeContrasts(parkinson-control, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

eb <- eBayes(fit2)

# extract the ranking table and show the top-ranked genes
ttable_zhang <- topTable(eb, n = nrow(zhangfilt2)) 

head(ttable_zhang)
#                 logFC  AveExpr         t      P.Value  adj.P.Val        B
#215812_s_at  0.4977913 7.653852  6.433649 8.067726e-07 0.01513624 5.470990
#210854_x_at  0.5029501 8.069007  6.229166 1.358546e-06 0.01513624 5.030257
#201658_at   -0.9458851 9.560856 -6.044963 2.180163e-06 0.01619352 4.628313
#219718_at   -0.4486547 7.372100 -5.670098 5.761279e-06 0.01907364 3.797219
#213843_x_at  0.4714778 8.020686  5.643794 6.170385e-06 0.01907364 3.738298
#221806_s_at  0.9908612 8.231832  5.612125 6.702084e-06 0.01907364 3.667263


# Limma analysis of Moran dataset
design <- model.matrix(~ -1+factor(moran_outcome_final))
colnames(design) <- unique(moran_outcome_final)

# compute simple linear model fit to microarray data (not robust)
fit <- lmFit(moranfilt, design)

contrast.matrix = makeContrasts(parkinson-control, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

eb <- eBayes(fit2)

# extract the ranking table and show the top-ranked genes
ttable_moran <- topTable(eb, n = nrow(moranfilt)) 

head(ttable_moran)
#                logFC  AveExpr          t      P.Value    adj.P.Val        B
#213920_at   -1.568310 4.203679 -11.013076 1.799182e-13 2.437296e-09 20.29504
#207087_x_at -1.702229 4.980834 -10.940202 2.187583e-13 2.437296e-09 20.11309
#209797_at   -1.023059 7.413173 -10.457428 8.126492e-13 6.036087e-09 18.88803
#209560_s_at -3.335964 4.377195  -9.804797 5.026856e-12 2.800336e-08 17.17761
#205391_x_at -2.052499 5.183501  -9.706350 6.648986e-12 2.963187e-08 16.91421
#205390_s_at -1.427009 7.179817  -9.612307 8.695324e-12 3.229298e-08 16.66129

#
# Meta-analysis
#

metarank = pvalcombination(esets=list(zhangfilt2, moranfilt), classes=list(zhang_label, moran_outcome_final), moderated = "limma", BHth = 0.05)
#     DE     IDD    Loss     IDR     IRR 
#4486.00 1118.00 2254.00   24.92   40.09

metacomb <- cbind(rownames(zhangfilt2)[metarank$Meta], metarank$TestStatistic[metarank$Meta])
metaord <- order(abs(as.numeric(metacomb[,2])), decreasing=TRUE)

logfcmat = cbind(ttable_zhang[match(rownames(zhangfilt2), rownames(ttable_zhang)),]$logFC, ttable_moran[match(rownames(zhangfilt2), rownames(ttable_moran)),]$logFC)
rownames(logfcmat)= rownames(zhangfilt2)

pmat = cbind(ttable_zhang[match(rownames(zhangfilt2), rownames(ttable_zhang)),]$P, ttable_moran[match(rownames(zhangfilt2), rownames(ttable_moran)),]$P)
rownames(pmat)= rownames(zhangfilt2)


compresmeta = cbind(metacomb[metaord,1], conv_ids[match(metacomb[metaord,1],rownames(logfcmat))], logfcmat[match(metacomb[metaord,1],rownames(logfcmat)),], pmat[match(metacomb[metaord,1],rownames(pmat)),], metacomb[metaord,2])
colnames(compresmeta) = c("ID", "Gene Symbol", paste(c("Zhang","Moran"),"logFC"), paste(c("Zhang","Moran"),"P"), "Comb. Z")

head(compresmeta)
#            ID            Gene Symbol Zhang logFC         Moran logFC         Zhang P               
#213920_at   "213920_at"   "CUX2"      "-1.5683104206"     "-1.5683104206"     "1.7991820716903e-13" 
#207087_x_at "207087_x_at" "ANK1"      "-1.70222868664166" "-1.70222868664166" "2.18758321067476e-13"
#203282_at   "203282_at"   "GBE1"      "-2.13122487826666" "-2.13122487826666" "1.78042429089719e-10"
#209797_at   "209797_at"   "CNPY2"     "-1.02305868425"    "-1.02305868425"    "8.12649167128433e-13"
#205391_x_at "205391_x_at" "ANK1"      "-2.05249918844167" "-2.05249918844167" "6.64898637205384e-12"
#221509_at   "221509_at"   "DENR"      "-1.16352955708333" "-1.16352955708333" "3.05431887233395e-11"
#            Moran P                Comb. Z            
#213920_at   "1.7991820716903e-13"  "-7.29642343429427"
#207087_x_at "2.18758321067476e-13" "-7.2172388962685" 
#203282_at   "1.78042429089719e-10" "-6.9554449934437" 
#209797_at   "8.12649167128433e-13" "-6.95527826366323"
#205391_x_at "6.64898637205384e-12" "-6.93116982534301"
#221509_at   "3.05431887233395e-11" "-6.91812814470966"

dim(compresmeta)
#[1] 4486    6


# save datasets to continue with pathway analysis tomorrow
save(moranvsn, moran_outcome_final, file="moran_preprocessed.Rdata")
save(zhangvsn, zhang_outcome_final, file="zhang_preprocessed.Rdata")

# For reproducibility: show and save information on all loaded R package versions
sessionInfo()
