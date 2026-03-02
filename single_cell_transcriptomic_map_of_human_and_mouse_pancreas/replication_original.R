# Replication of Baron et al. (2016) - Original Methodology
# Methodology focused on inDrop pipeline and recursive hierarchical clustering.

library(bseqsc)
library(edgeR)
library(cluster)

# Set working directory (same as downloader)
output_path <- "./data"
if (dir.exists(output_path)) {
    setwd(output_path)
}

# --- 1. Data Loading ---
# Note: The user will need to unzip the files downloaded by data_downloader.R
# Usually: GSE84133_human_islets_all_cells.txt.gz
human_data_file <- file.path("GSE84133", "GSE84133_human_islets_all_cells.csv.gz")

if (file.exists(human_data_file)) {
    # Use fread for speed and reliability, then convert to matrix
    library(data.table)
    counts_dt <- fread(human_data_file)
    genes <- counts_dt$gene
    counts_matrix <- as.matrix(counts_dt[, -1])
    rownames(counts_matrix) <- genes
    counts <- as.data.frame(counts_matrix)
} else {
    stop("Data file not found. Please run merge_donors.R and ensure files are in the expected path (./data/GSE84133).")
}

# --- 2. Preprocessing (Original Style) ---
# Conversion to TPM (Transcripts Per Million)
# Baron et al. quantified expression as TPM
counts_matrix <- as.matrix(counts)
tpm <- apply(counts_matrix, 2, function(x) (x / sum(x)) * 1e6)

# --- 3. Recursive Hierarchical Clustering ---
# The paper used a custom recursive approach. Here we implement the logic:
# 1. Filter genes: High expression (fraction of total TPM > 0.5) and high variance (Fano factor).
# 2. Hierarchical clustering with Ward's criterion.

calculate_fano <- function(x) {
    var(x) / mean(x)
}

# Simple implementation of the first level of clustering
gene_means <- rowMeans(tpm)
gene_fano <- apply(tpm, 1, calculate_fano)

# Thresholds as per paper: Fano factor above mean-dependent threshold
# Here we take top variable genes for demonstration
variable_genes <- names(sort(gene_fano, decreasing = TRUE)[1:2000])

# log TPM for clustering
log_tpm_subset <- log2(tpm[variable_genes, ] + 1)

# Hierarchical clustering (Ward's D2)
d <- dist(t(log_tpm_subset))
hc <- hclust(d, method = "ward.D2")
clusters <- cutree(hc, k = 14) # Paper identified ~14-15 types

# --- 4. Marker Identification ---
# Find genes exclusively expressed in clusters
# Using KS test as mentioned in the paper
markers <- list()
for (i in unique(clusters)) {
    # Logic: KS test comparing cluster i vs others
    # Simplified for the script
    cluster_cells <- which(clusters == i)
    other_cells <- which(clusters != i)
    # (Placeholder for full KS test implementation)
}

# --- 5. Deconvolution with bseqsc ---
# Using the bseqsc package developed by the authors
# Requires a basis matrix of markers

# Example of using bseqsc to deconvolve bulk data (GSE50244)
# (Assuming bulk_tpm is loaded from GSE50244)
# bseq_sc_results <- bseqsc_deconvolve(bulk_tpm, basis_matrix)

print("Original replication script skeleton complete.")
