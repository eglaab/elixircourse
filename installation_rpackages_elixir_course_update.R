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
