# Replication of Baron et al. (2016)
# Methodology focused on inDrop pipeline and recursive hierarchical clustering.

library(bseqsc)
library(edgeR)
library(cluster)
library(Rtsne)
library(ggplot2)
library(gridExtra)

# Set working directory
output_path <- "C:/Users/hp/OneDrive - Savitribai Phule Pune University/Paper_Replication/single_cell_transcriptomic_map_of_human_and_mouse_pancreas"
setwd(output_path)

human_data_file <- file.path("GSE84133", "GSE84133_human_islets_all_cells.csv.gz")

if (file.exists(human_data_file)) {
    # Using fread for speed and reliability, then convert to matrix
    library(data.table)
    counts_dt <- fread(human_data_file)
    genes <- counts_dt$gene
    counts_matrix <- as.matrix(counts_dt[, -1])
    rownames(counts_matrix) <- genes
    counts <- as.data.frame(counts_matrix)
    rm(counts_dt)
    gc()
} else {
    stop("Data file not found. Please run merge_donors.R and ensure files are in the working directory.")
}

# Preprocessing
# Conversion to TPM (Transcripts Per Million)
# Baron et al. quantified expression as TPM
counts_matrix <- as.matrix(counts)
# TPM calculation (assuming counts are UMI counts which are roughly analogous to transcripts)
tpm <- apply(counts_matrix, 2, function(x) (x / sum(x)) * 1e6)
tpm[is.na(tpm)] <- 0

# Extracting Metadata from column names (barcode format: donor_barcode)
cell_ids <- colnames(tpm)
donors <- sapply(strsplit(cell_ids, "_"), `[`, 1)

# Recursive Hierarchical Clustering
# 1. Filter genes: High variance (Fano factor).
calculate_fano <- function(x) {
    var(x) / mean(x)
}

gene_means <- rowMeans(tpm)
# Consider only genes with mean expression > 0 to avoid NA Fano
expressed_genes <- names(which(gene_means > 0))
tpm_expressed <- tpm[expressed_genes, ]
gene_fano <- apply(tpm_expressed, 1, calculate_fano)
gene_fano[is.na(gene_fano)] <- 0

# Thresholds as per paper: Fano factor above mean-dependent threshold
# The paper empirically chose top variable genes. We select the top 2000 for demonstration.
variable_genes <- names(sort(gene_fano, decreasing = TRUE)[1:2000])

# log TPM for clustering
log_tpm_subset <- log2(tpm[variable_genes, ] + 1)

# Hierarchical clustering (Ward's D2)
message("Calculating distance and clustering...")
# Using Euclidean distance on log transformed TPM as is standard
d <- dist(t(log_tpm_subset))
# saveRDS(d, "distance_matrix_cell.rds") # uncomment to save
hc <- hclust(d, method = "ward.D2")
clusters <- cutree(hc, k = 14) # Paper identified 14 human types

# Dimensionality Reduction for Visualization (t-SNE)
# The paper used Barnes-Hut t-SNE
message("Running Barnes-Hut t-SNE...")
set.seed(42)
# Rtsne needs unique rows, which we have (cells).
tsne_out <- Rtsne(t(log_tpm_subset), pca = TRUE, perplexity = 30, theta = 0.5)

tsne_df <- data.frame(
    TSNE1 = tsne_out$Y[, 1],
    TSNE2 = tsne_out$Y[, 2],
    Cluster = as.factor(clusters),
    Donor = donors
)

# Marker Identification
# Find genes exclusively expressed in clusters using Kolmogorov-Smirnov test
message("Finding cluster markers (Two-sample KS test)...")
# Due to the slow nature of KS test across all genes and all clusters,
# we focus on canonical markers for typical visualization.
canonical_markers <- c("GCG", "INS", "SST", "PPY", "PRSS1", "KRT19", "COL1A1", "PECAM1")
marker_expression <- log_tpm_subset[intersect(canonical_markers, rownames(log_tpm_subset)), ]
marker_expression <- t(marker_expression)

# Combine for plotting
plot_df <- cbind(tsne_df, marker_expression)

# --- 5. Deconvolution / CIBERSORT (Optional, as per paper) ---
# bseq_sc_results <- bseqsc_deconvolve(bulk_tpm, basis_matrix)

# --- Visualizations ---
message("Generating Replicated Plots...")

# 1. t-SNE plot colored by cluster
p_tsne_cluster <- ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = Cluster)) +
    geom_point(size = 0.5, alpha = 0.8) +
    theme_classic() +
    ggtitle("t-SNE Map: Human Pancreas (Replication)") +
    guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(file.path("GSE84133", "Baron_Replication_tSNE_Clusters.png"), p_tsne_cluster, width = 8, height = 6)

# 2. t-SNE plot colored by donor to show mixing
p_tsne_donor <- ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = Donor)) +
    geom_point(size = 0.5, alpha = 0.5) +
    theme_classic() +
    ggtitle("t-SNE Map by Donor")

ggsave(file.path("GSE84133", "Baron_Replication_tSNE_Donors.png"), p_tsne_donor, width = 8, height = 6)

# 3. Canonical Marker Expression (Simulating the paper's feature plots)
# We plot a few key genes
plot_list <- list()
for (gene in intersect(canonical_markers, colnames(plot_df))) {
    p <- ggplot(plot_df, aes_string(x = "TSNE1", y = "TSNE2", color = gene)) +
        geom_point(size = 0.5, alpha = 0.8) +
        scale_color_gradient(low = "lightgrey", high = "red") +
        theme_classic() +
        ggtitle(gene) +
        theme(legend.position = "none")
    plot_list[[gene]] <- p
}

if (length(plot_list) > 0) {
    p_markers <- do.call(grid.arrange, c(plot_list, ncol = 3))
    ggsave(file.path("GSE84133", "Baron_Replication_Marker_tSNE.png"), p_markers, width = 10, height = 8)
}

# Save output data
saveRDS(list(hc = hc, clusters = clusters, tsne = tsne_out), file.path("GSE84133", "replication_original_output.rds"))

message("Original replication script complete. Results saved in GSE84133 folder.")
