# Merge Baron et al. (2016) Human Donor Data

library(data.table)

output_path <- "C:/Users/hp/OneDrive - Savitribai Phule Pune University/Paper_Replication/single_cell_transcriptomic_map_of_human_and_mouse_pancreas"
data_dir <- paste0(output_path, "/GSE84133")
setwd(data_dir)
human_files <- list.files(data_dir, pattern = "human.*_counts.csv.gz", full.names = TRUE)

all_data <- list()

for (f in human_files) {
    donor_id <- sub(".*(human[1-4]).*", "\\1", basename(f))
    message("Processing ", donor_id, "...")

    # Read CSV (rows = cells, cols = metadata + genes)
    dt <- fread(f)

    # Finding barcode column
    barcode_col <- intersect(names(dt), "barcode")
    if (is.na(barcode_col)) barcode_col <- names(dt)[1] # Fallback to first column

    # Ensuring cell names are unique across donors
    dt[[barcode_col]] <- paste0(donor_id, "_", dt[[barcode_col]])

    all_data[[donor_id]] <- dt
}

# 1. Combining all donors
merged_dt <- rbindlist(all_data, fill = TRUE)

# 2. Extracting metadata
metadata <- merged_dt[, .(barcode, assigned_cluster)]
metadata[, donor := sub("_.*", "", barcode)]

# 3. Processing count matrix
gene_cols <- setdiff(names(merged_dt), c("V1", "barcode", "assigned_cluster"))
counts_dt <- merged_dt[, ..gene_cols]

# Converting to matrix
counts_matrix <- as.matrix(counts_dt)
rownames(counts_matrix) <- merged_dt$barcode
counts_matrix <- t(counts_matrix)

# Replacing NAs with 0
counts_matrix[is.na(counts_matrix)] <- 0

# 4. Saving merged counts
output_file <- file.path(data_dir, "GSE84133_human_islets_all_cells.csv.gz")

# Convert back to data.table for fwrite
final_dt <- as.data.table(counts_matrix, keep.rownames = "gene")
fwrite(final_dt, output_file)
