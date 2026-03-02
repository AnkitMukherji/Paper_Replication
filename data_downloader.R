# Data Downloader for Baron et al. (2016) Replication

# This script downloads the necessary GEO datasets for the replication.
# It requires the 'GEOquery' Bioconductor package.

library(GEOquery)

# Set output directory
output_path <- "./data"
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}
setwd(output_path)

# 1. Download Single-Cell RNA-seq data from Baron et al. (GSE84133)
print("Downloading GSE84133 (Single-Cell RNA-seq)...")
gse84133 <- getGEO("GSE84133", destdir = ".")
# The supplemental files contain the count matrices
getGEOSuppFiles("GSE84133")

# 2. Download Bulk RNA-seq data from Fadista et al. (GSE50244)
# This is used for the deconvolution part of the paper
print("Downloading GSE50244 (Bulk RNA-seq)...")
gse50244 <- getGEO("GSE50244", destdir = ".")
getGEOSuppFiles("GSE50244")

print("Downloads complete. Please check the 'GSE84133' and 'GSE50244' folders.")
