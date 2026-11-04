# Single-cell Paper Replication and Methodological Extensions

R scripts for revisiting a published single-cell pancreas study, comparing an original-method replication with a modern workflow, and exploring additional biological analyses.

## Study and data

The project directory [single_cell_transcriptomic_map_of_human_and_mouse_pancreas/](single_cell_transcriptomic_map_of_human_and_mouse_pancreas/) contains a Baron et al. pancreas analysis. The downloader requests **GSE84133** for single-cell data and **GSE50244** for bulk RNA-seq used in the deconvolution context.

## Reading order

1. [1_data_downloader.R](single_cell_transcriptomic_map_of_human_and_mouse_pancreas/1_data_downloader.R) — retrieve GEO data.
2. [2_merge_donors.R](single_cell_transcriptomic_map_of_human_and_mouse_pancreas/2_merge_donors.R) — assemble donor inputs.
3. [3_replication_original.R](single_cell_transcriptomic_map_of_human_and_mouse_pancreas/3_replication_original.R) — revisit the original analysis route.
4. [4_replication_modern.R](single_cell_transcriptomic_map_of_human_and_mouse_pancreas/4_replication_modern.R) — compare a modern analysis route.
5. [5_extension_insights.R](single_cell_transcriptomic_map_of_human_and_mouse_pancreas/5_extension_insights.R) — explore extensions including cell communication, metabolism, and trajectory analysis.

## Setup and execution

The repository records an R environment in [renv.lock](renv.lock). Start at the repository root:

```r
renv::restore()
```

Read scripts interactively before executing them. The downloader contains a machine-specific Windows working directory that must be changed. Review paths, package availability, downloaded archive layout, and object dependencies throughout the sequence. Restoring packages alone does not make every analysis portable.

## What a replication should establish

Compare cell filtering, normalization, feature selection, clustering, marker expression, donor composition, and annotation against the source study. Record where package changes or deliberate methodological choices affect agreement.

The presence of replication scripts is not evidence that every published result has been reproduced. Modern-workflow differences should be distinguished from biological discoveries. Communication, metabolic, and trajectory extensions are exploratory and require checks for donor effects, cell composition, and model assumptions.

Related: [single-cell workflows in Python and R](https://github.com/AnkitMukherji/single_cell_analysis).
