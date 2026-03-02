# Pancreatic Single-Cell Transcriptomics: Exhaustive Extension Analysis
# Aim: Gain new biological insights beyond cell-type mapping.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set paths
data_dir <- "e:/Bioinformatics/Paper_Replication/data/GSE84133"
rds_file <- file.path(data_dir, "pancreas_processed.rds")

if (!file.exists(rds_file)) {
    stop("Processed Seurat object not found. Please run replication_modern.R first.")
}

message("Loading processed Seurat object...")
pancreas <- readRDS(rds_file)

# --- 1. Module 1: Functional Regulatory Networks (DoRothEA) ---
# Infer Transcription Factor activity.
message("Running DoRothEA insight module...")

# Load DoRothEA regulons (Human, Confidence A, B, C)
if (requireNamespace("dorothea", quietly = TRUE) && requireNamespace("viper", quietly = TRUE)) {
    library(dorothea)
    library(viper)

    # Get regulons
    data(dorothea_hs, package = "dorothea")
    regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))

    # Run Viper to infer TF activity
    # We use the SCT assay for counts
    # Standard dorothea::run_viper signature for Seurat:
    pancreas <- run_viper(pancreas, regulons,
        assay = "SCT", slot = "data",
        minsize = 4,
        verbose = FALSE
    )

    # Identify top TFs per cell type
    DefaultAssay(pancreas) <- "dorothea"
    # Ensure scale.data is present for heatmap
    pancreas <- ScaleData(pancreas, assay = "dorothea")

    tf_markers <- FindAllMarkers(pancreas, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

    top5_tfs <- tf_markers %>%
        group_by(cluster) %>%
        top_n(n = 5, wt = avg_log2FC)

    p1 <- DoHeatmap(pancreas, features = top5_tfs$gene, slot = "scale.data") +
        ggtitle("Top Inferred TF Activities by Cell Type")
    ggsave(file.path(data_dir, "TF_activity_heatmap.png"), p1, width = 12, height = 10)
} else {
    message("dorothea or viper not installed. Skipping TF module.")
}

# --- 2. Module 2: Cell-Cell Communication (CellChat) ---
message("Running CellChat insight module...")

if (requireNamespace("CellChat", quietly = TRUE)) {
    library(CellChat)

    # Create CellChat object
    cellchat <- createCellChat(object = pancreas, group.by = "ident", assay = "SCT")

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
    png(file.path(data_dir, "cell_communication_network.png"), width = 800, height = 800)
    netVisual_circle(cellchat@net$count,
        vertex.weight = as.numeric(table(pancreas@active.ident)),
        weight.scale = T, label.edge = F, title.name = "Number of interactions"
    )
    dev.off()
} else {
    message("CellChat not installed. Skipping communication module.")
}

# --- 3. Module 3: Metabolic Profiling (scMetabolism) ---
message("Running scMetabolism insight module...")

if (requireNamespace("scMetabolism", quietly = TRUE)) {
    library(scMetabolism)

    # Run metabolic scoring
    pancreas <- sc.metabolism.Seurat(
        obj = pancreas, method = "AUCell",
        imputation = FALSE, ncores = 1,
        metabolism.type = "KEGG"
    )

    # Plot top metabolic pathways
    # Seurat automatically adds the scores to the metadata
    metabolic_cols <- grep("KEGG.", colnames(pancreas@meta.data), value = TRUE)
    if (length(metabolic_cols) > 0) {
        p2 <- DotPlot(pancreas, features = metabolic_cols[1:10]) +
            RotatedAxis() + ggtitle("Top Metabolic Pathway Scores")
        ggsave(file.path(data_dir, "metabolic_pathways_dotplot.png"), p2, width = 12, height = 6)
    }
} else {
    message("scMetabolism not installed. Skipping metabolism module.")
}

# --- 4. Module 4: Trajectory Inference (Monocle 3) ---
message("Running Monocle 3 insight module...")

if (requireNamespace("monocle3", quietly = TRUE) && requireNamespace("SeuratWrappers", quietly = TRUE)) {
    library(monocle3)
    library(SeuratWrappers)

    # Convert to Monocle3
    cds <- as.cell_data_set(pancreas)
    cds <- cluster_cells(cds)
    cds <- learn_graph(cds)

    p3 <- plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = FALSE) +
        ggtitle("Monocle 3 Trajectory Graph")
    ggsave(file.path(data_dir, "monocle3_trajectory.png"), p3, width = 8, height = 6)
} else {
    message("monocle3 or SeuratWrappers not installed. Skipping trajectory module.")
}

message("Exhaustive Extension Analysis Modules Complete.")
