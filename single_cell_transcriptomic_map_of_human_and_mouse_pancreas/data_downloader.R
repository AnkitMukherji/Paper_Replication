library(GEOquery)

# Set output directory
output_path <- "C:/Users/hp/OneDrive - Savitribai Phule Pune University/Paper_Replication/single_cell_transcriptomic_map_of_human_and_mouse_pancreas"
setwd(output_path)

# 1. Downloading Single-Cell RNA-seq data from Baron et al. (GSE84133)
gse84133 <- getGEO("GSE84133", destdir = ".")
# The supplemental files contain the count matrices
getGEOSuppFiles("GSE84133")

# 2. Downloading Bulk RNA-seq data from Fadista et al. (GSE50244)
# This is used for the deconvolution part of the paper
gse50244 <- getGEO("GSE50244", destdir = ".")
getGEOSuppFiles("GSE50244")
