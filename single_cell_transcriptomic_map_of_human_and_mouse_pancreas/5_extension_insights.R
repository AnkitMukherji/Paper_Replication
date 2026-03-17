# Pancreatic Single-Cell Transcriptomics: Extension Analysis
# Aim: Gaining new biological insights beyond cell-type mapping.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set paths
output_path <- "C:/Users/hp/OneDrive/Paper_Replication/single_cell_transcriptomic_map_of_human_and_mouse_pancreas"
data_dir <- file.path(output_path, "GSE84133")
rds_file <- file.path(data_dir, "pancreas_processed.rds")

pancreas <- readRDS(rds_file)

# 1: Functional Regulatory Networks (DoRothEA)
# Infer Transcription Factor activity.
library(dorothea)
library(viper)

# Get regulons
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs |> filter(confidence %in% c("A", "B", "C"))
# Convert to VIPER regulon object
regulon_list <- df2regulon(regulons)

# Run Viper to infer TF activity
expr <- GetAssayData(pancreas, assay = "SCT", layer = "data")
expr <- as.matrix(expr)
viper_res <- viper(expr, regulon_list, minsize = 4, verbose = FALSE)

# Add TF activity to Seurat
pancreas[["dorothea"]] <- CreateAssayObject(data = viper_res)

# Identify top TFs per cell type
DefaultAssay(pancreas) <- "dorothea"
# Scale TF activity for visualization
pancreas <- ScaleData(pancreas, assay = "dorothea")

tf_markers <- FindAllMarkers(pancreas, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

top5_tfs <- tf_markers |> 
    group_by(cluster) |> 
    top_n(n = 5, wt = avg_log2FC)

p1 <- DoHeatmap(pancreas, features = top5_tfs$gene, slot = "scale.data") +
    ggtitle("Top Inferred TF Activities by Cell Type")

# 2: Cell-Cell Communication (CellChat)
library(CellChat)

# Extract normalized data
data.input <- GetAssayData(pancreas, assay = "SCT", layer = "data")

meta <- pancreas@meta.data
meta$ident <- paste0("Cluster_", Idents(pancreas))
# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "ident")

# Set database
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# Processing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability
cellchat <- computeCommunProb(cellchat)
cellchat <- aggregateNet(cellchat)

# Visualization
netVisual_circle(cellchat@net$count,
    vertex.weight = as.numeric(table(pancreas@active.ident)),
    weight.scale = T, label.edge = F, title.name = "Number of interactions"
)

# 3: Metabolic Profiling (scMetabolism)
library(scMetabolism)

counts_mat <- GetAssayData(
  pancreas,
  assay = "RNA",
  layer = "counts"
)
rna_v4 <- CreateAssayObject(counts = counts_mat)
Key(rna_v4) <- "RNA_"
pancreas_clean <- new(
  Class = "Seurat",
  assays = list(RNA = rna_v4),
  meta.data = pancreas@meta.data
)
pancreas_clean <- new(
  Class = "Seurat",
  assays = list(RNA = rna_v4),
  meta.data = pancreas@meta.data
)
DefaultAssay(pancreas_clean) <- "RNA"
pancreas_clean <- NormalizeData(pancreas_clean)

pancreas_clean <- sc.metabolism.Seurat(
  obj = pancreas_clean,
  method = "AUCell",
  imputation = FALSE,
  ncores = 1,
  metabolism.type = "KEGG"
)

# 4: Trajectory Inference (Monocle 3)
# Install grr dependency
# install.packages(
#   "https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz",
#   repos = NULL,
#   type = "source"
# )
library(monocle3)
library(SeuratWrappers)

# Convert to Monocle3
cds <- as.cell_data_set(pancreas)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

p3 <- plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = FALSE) +
    ggtitle("Monocle 3 Trajectory Graph")
