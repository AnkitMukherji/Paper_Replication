# Pancreas single-cell paper replication

Use the root renv record, review study-specific paths, and read the original study alongside the code. Run acquisition and donor merging before the replication paths. Check the intermediate object named in the extension script before running it.

## Script and notebook guide

| File | Purpose |
| --- | --- |
| [`1_data_downloader.R`](1_data_downloader.R) | Uses GEOquery to obtain GSE84133 single-cell data and GSE50244 bulk RNA-seq data for the study context. |
| [`2_merge_donors.R`](2_merge_donors.R) | Loads donor count tables and prepares merged donor inputs for downstream analyses. |
| [`3_replication_original.R`](3_replication_original.R) | Implements the project’s reconstruction of the original analysis route using the merged data. |
| [`4_replication_modern.R`](4_replication_modern.R) | Applies a modern single-cell analysis route for comparison and saves its outputs for downstream inspection. |
| [`5_extension_insights.R`](5_extension_insights.R) | Explores additional biology with cell communication, metabolism and trajectory methods. |

[Repository overview](../README.md)
