# Replication of Baron et al. (2016) - Modern Methodology
# Using Seurat v5 for analysis.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set working directory
output_path <- "C:/Users/hp/OneDrive/Paper_Replication/single_cell_transcriptomic_map_of_human_and_mouse_pancreas"
data_dir <- file.path(output_path, "GSE84133")

# Data Loading
human_data_file <- file.path(data_dir, "GSE84133_human_islets_all_cells.csv.gz")

if (file.exists(human_data_file)) {
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
pancreas <- CreateSeuratObject(counts = counts, project = "Pancreas_sc", min.cells = 3, min.features = 200)
rm(counts) # Free memory
gc()

# Quality Control
# Calculate mitochondrial percentage (human genes usually start with MT-)
# no information about MT present in the dataset
# pancreas[["percent.mt"]] <- PercentageFeatureSet(pancreas, pattern = "^MT-")

# Extract donor info
pancreas$donor <- sapply(strsplit(colnames(pancreas), "_"), `[`, 1)

# Visualize QC metrics
VlnPlot(pancreas, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Filter cells (Adjust thresholds based on VlnPlot results)
# Standard thresholds: nFeature [200, 2500], percent.mt < 5%
pancreas <- subset(pancreas, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# Normalization and Scaling (SCTransform)
# SCTransform is the modern replacement for NormalizeData, ScaleData, and FindVariableFeatures
pancreas <- SCTransform(pancreas)

# Doublet Removal (DoubletFinder)
library(DoubletFinder)
pancreas <- RunPCA(pancreas, verbose = FALSE)
pancreas <- RunUMAP(pancreas, dims = 1:20, verbose = FALSE)

# Assuming ~7.5% doublet formation rate for 10k cells
nExp_poi <- round(0.075 * ncol(pancreas))
pancreas <- doubletFinder(pancreas, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)

# Remove doublets
DF.name <- grep("DF.classification", colnames(pancreas@meta.data), value = TRUE)
pancreas <- subset(pancreas, subset = get(DF.name) == "Singlet")
message(paste("Removed doublets. Remaining cells:", ncol(pancreas)))

# Dimensionality Reduction & Integration (Harmony)
pancreas <- RunPCA(pancreas, verbose = FALSE)
library(harmony)
pancreas <- RunHarmony(pancreas, group.by.vars = "donor", assay.use = "SCT")
reduction_use <- "harmony"

pancreas <- RunUMAP(pancreas, reduction = reduction_use, dims = 1:20, verbose = FALSE)
pancreas <- FindNeighbors(pancreas, reduction = reduction_use, dims = 1:20, verbose = FALSE)
pancreas <- FindClusters(pancreas, resolution = 0.5, verbose = FALSE)

# Visualization
p_umap_cluster <- DimPlot(pancreas, reduction = "umap", label = TRUE) + NoLegend() + ggtitle("Seurat UMAP")
p_umap_donor <- DimPlot(pancreas, reduction = "umap", group.by = "donor") + ggtitle("Seurat UMAP by Donor")

# Marker Identification and Cluster Annotation
# Identify markers for every cluster compared to all remaining cells
all_markers <- FindAllMarkers(pancreas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all_markers |> 
    group_by(cluster) |> 
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

# Cell Type Assignment
# (Manual assignment based on markers)
# new_cluster_ids <- c("Alpha", "Beta", "Ductal", "Acinar", "Delta", "Gamma", "Stellate", "Endothelial", "Macrophage", "Mast", "Schwann", "Epsilon")
# Note: Ensure the vector matches the number of clusters found
# names(new_cluster_ids) <- levels(pancreas)
# pancreas <- RenameIdents(pancreas, new_cluster_ids)

# Save Object for Extension Analysis
saveRDS(pancreas, file.path(data_dir, "pancreas_processed.rds"))
