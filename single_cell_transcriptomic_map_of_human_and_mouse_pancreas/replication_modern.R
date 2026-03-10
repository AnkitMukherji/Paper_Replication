# Replication of Baron et al. (2016) - Modern Methodology
# Using Seurat v5 for state-of-the-art analysis.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set working directory
output_path <- "C:/Users/hp/OneDrive - Savitribai Phule Pune University/Paper_Replication/single_cell_transcriptomic_map_of_human_and_mouse_pancreas"
data_dir <- file.path(output_path, "GSE84133")

# --- 1. Data Loading ---
human_data_file <- file.path(data_dir, "GSE84133_human_islets_all_cells.csv.gz")

if (file.exists(human_data_file)) {
    message("Loading merged human data from ", human_data_file)
    # Using data.table::fread for significantly faster loading of large gzipped CSVs
    counts_dt <- data.table::fread(human_data_file)
    # Convert to matrix with gene names as rownames
    counts <- as.matrix(counts_dt[, -1])
    rownames(counts) <- counts_dt[[1]]
    rm(counts_dt)
} else {
    stop("Merged data file not found. Please ensure merge_donors.R has been run successfully.")
}

# Create Seurat Object
pancreas <- CreateSeuratObject(counts = counts, project = "PancreasReplication", min.cells = 3, min.features = 200)
rm(counts) # Free memory
gc()

# --- 2. Quality Control (Modern) ---
# Calculate mitochondrial percentage (human genes usually start with MT-)
pancreas[["percent.mt"]] <- PercentageFeatureSet(pancreas, pattern = "^MT-")

# Extract donor info
pancreas$donor <- sapply(strsplit(colnames(pancreas), "_"), `[`, 1)

# Visualize QC metrics
# VlnPlot(pancreas, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells (Adjust thresholds based on VlnPlot results)
# Standard thresholds: nFeature [200, 2500], percent.mt < 5%
pancreas <- subset(pancreas, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# --- 3. Normalization and Scaling (SCTransform) ---
# SCTransform is the modern replacement for NormalizeData, ScaleData, and FindVariableFeatures
pancreas <- SCTransform(pancreas, vars.to.regress = "percent.mt", verbose = FALSE)

# --- 4. Doublet Removal (DoubletFinder) ---
message("Running DoubletFinder...")
if (requireNamespace("DoubletFinder", quietly = TRUE)) {
    library(DoubletFinder)
    pancreas <- RunPCA(pancreas, verbose = FALSE)
    pancreas <- RunUMAP(pancreas, dims = 1:20, verbose = FALSE)

    # Assuming ~7.5% doublet formation rate for 10k cells
    nExp_poi <- round(0.075 * ncol(pancreas))
    pancreas <- doubletFinder_v3(pancreas, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

    # Remove doublets
    DF.name <- grep("DF.classification", colnames(pancreas@meta.data), value = TRUE)
    pancreas <- subset(pancreas, subset = get(DF.name) == "Singlet")
    message(paste("Removed doublets. Remaining cells:", ncol(pancreas)))
} else {
    message("DoubletFinder package not found. Skipping doublet removal.")
}

# --- 5. Dimensionality Reduction & Integration (Harmony) ---
message("Running Dimensionality Reduction and Harmony Integration...")
pancreas <- RunPCA(pancreas, verbose = FALSE)

if (requireNamespace("harmony", quietly = TRUE)) {
    library(harmony)
    pancreas <- RunHarmony(pancreas, group.by.vars = "donor", assay.use = "SCT")
    reduction_use <- "harmony"
} else {
    message("harmony package not found. Proceeding without integration.")
    reduction_use <- "pca"
}

pancreas <- RunUMAP(pancreas, reduction = reduction_use, dims = 1:20, verbose = FALSE)
pancreas <- FindNeighbors(pancreas, reduction = reduction_use, dims = 1:20, verbose = FALSE)
pancreas <- FindClusters(pancreas, resolution = 0.5, verbose = FALSE)

# --- 6. Visualization ---
p_umap_cluster <- DimPlot(pancreas, reduction = "umap", label = TRUE) + NoLegend() + ggtitle("Seurat UMAP (Modern)")
p_umap_donor <- DimPlot(pancreas, reduction = "umap", group.by = "donor") + ggtitle("Seurat UMAP by Donor")

ggsave(file.path(data_dir, "Baron_Modern_UMAP_Clusters.png"), p_umap_cluster, width = 8, height = 6)
ggsave(file.path(data_dir, "Baron_Modern_UMAP_Donors.png"), p_umap_donor, width = 8, height = 6)

# --- 7. Marker Identification and Cluster Annotation ---
# Identify markers for every cluster compared to all remaining cells
all_markers <- FindAllMarkers(pancreas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all_markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)

# Canonical Markers for Pancreas:
# Alpha: GCG
# Beta: INS
# Delta: SST
# Gamma (PP): PPY
# Acinar: PRSS1
# Ductal: KRT19
# Stellate: COL1A1
# Endothelial: PECAM1
# Immune: PTPRC (CD45), CD68

p_canonical <- FeaturePlot(pancreas, features = c("GCG", "INS", "SST", "PPY", "PRSS1", "KRT19"))
ggsave(file.path(data_dir, "Baron_Modern_Canonical_Markers.png"), p_canonical, width = 10, height = 8)

# --- 8. Cell Type Assignment ---
# (Manual assignment based on markers)
new_cluster_ids <- c("Alpha", "Beta", "Ductal", "Acinar", "Delta", "Gamma", "Stellate", "Endothelial", "Macrophage", "Mast", "Schwann", "Epsilon")
# Note: Ensure the vector matches the number of clusters found
# names(new_cluster_ids) <- levels(pancreas)
# pancreas <- RenameIdents(pancreas, new_cluster_ids)

# --- 8. Comparison with Original Analysis ---
# Export metadata for comparison with Baron et al. assignments if available
# write.csv(pancreas@meta.data, "modern_clustering_metadata.csv")

# --- 9. Save Object for Extension Analysis ---
saveRDS(pancreas, file.path(data_dir, "pancreas_processed.rds"))
message("Processed Seurat object saved to pancreas_processed.rds")

print("Modern replication script (Seurat v5) complete.")
