# Merge Baron et al. (2016) Human Donor Data
# Properly aggregates genes with synonymous names (e.g., 7-Sep vs 7_Sep).

library(data.table)

data_dir <- "e:/Bioinformatics/Paper_Replication/data/GSE84133"
human_files <- list.files(data_dir, pattern = "human.*_counts.csv.gz", full.names = TRUE)

message("Merging files:")
print(human_files)

all_data <- list()

for (f in human_files) {
    donor_id <- sub(".*(human[1-4]).*", "\\1", basename(f))
    message("Processing ", donor_id, "...")

    # Read CSV (rows = cells, cols = metadata + genes)
    dt <- fread(f)

    # Standardize index/barcode columns
    # Find barcode column (sometimes called 'barcode', sometimes 'V1' or 'X' in headers)
    barcode_col <- intersect(names(dt), c("barcode", "V1", "X"))[1]
    if (is.na(barcode_col)) barcode_col <- names(dt)[1] # Fallback to first column

    # Ensure cell names are unique across donors
    dt[[barcode_col]] <- paste0(donor_id, "_", dt[[barcode_col]])

    all_data[[donor_id]] <- dt
}

# 1. Combine all donors
message("Combining all donor data...")
merged_dt <- rbindlist(all_data, fill = TRUE)

# 2. Extract metadata
message("Extracting metadata...")
barcode_col <- intersect(names(merged_dt), c("barcode", "V1", "X"))[1]
metadata_cols <- intersect(names(merged_dt), c(barcode_col, "assigned_cluster"))
metadata <- merged_dt[, ..metadata_cols]
setnames(metadata, barcode_col, "barcode")
metadata[, donor := sub("_.*", "", barcode)]

# 3. Process count matrix
message("Processing count matrix and aggregating duplicate genes...")
# Exclude metadata for count part
gene_cols <- setdiff(names(merged_dt), c("V1", "barcode", "X", "assigned_cluster"))
counts_dt <- merged_dt[, ..gene_cols]


# Convert to matrix
counts_matrix <- as.matrix(counts_dt)
rownames(counts_matrix) <- merged_dt$cell_id
counts_matrix <- t(counts_matrix)

# Replace NAs with 0
counts_matrix[is.na(counts_matrix)] <- 0

# 4. Save merged counts
output_file <- file.path(data_dir, "GSE84133_human_islets_all_cells.csv.gz")
message("Saving merged matrix to ", output_file, "...")

# Convert back to data.table for fwrite
final_dt <- as.data.table(counts_matrix, keep.rownames = "gene")
fwrite(final_dt, output_file)

message("Merge complete.")
