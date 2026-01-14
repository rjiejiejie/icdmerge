# icdmerge

`icdmerge` is an R package for merging fine-grained ICD diagnosis nodes into clinically interpretable merged nodes using **block-wise consensus k-means** with **PAC-based K selection**.

## Overview

Key design choices:

- **Global structural profiles**: each disease is represented by a co-occurrence–based structural profile defined against *all diseases in the network* (global roles), not only within its own block.
- **Clinical interpretability constraint**: merging is only allowed within the same clinical `block` (e.g., ICD-10 subchapter / sub-block), preventing cross-block merges.
- **Consensus clustering + PAC**: for each block, we resample patients, run k-means, build a consensus co-clustering matrix, and select the block-specific `K*` by minimizing PAC.
- **Optional homogeneity gate (K=1 allowed)**: if the block shows little separable structure (e.g., the `K=2` consensus matrix is nearly uniform), we treat it as a homogeneous block and allow a full-block merge under `K=1`.

## Installation

Install from local source during development:

```r
devtools::install()
```

## Input data format
 - raw_data
 A table with at least these columns:
  - eid: patient identifier
  - Category_No: disease identifier (mapped ICD node ID)
 Each row means patient eid has diagnosis Category_No.
 
 - codebook
 A table describing each disease node, with at least:
  - No: disease identifier (must match Category_No)
  - Frequency: number of patients with the disease
  - Codes: original ICD codes (or a string)
  - Names: disease names (or a string)
  - block: clinical block label used to restrict merging

## Quick start
```r
library(data.table)
library(icdmerge)

set.seed(1)

# ---- toy raw_data: patient-disease table ----
n_patients <- 200
diseases <- sprintf("D%02d", 1:12)

raw_data <- data.table(
  eid = sample(sprintf("P%03d", 1:n_patients), size = 900, replace = TRUE),
  Category_No = sample(diseases, size = 900, replace = TRUE)
)
raw_data <- unique(raw_data)

# ---- toy codebook ----
codebook <- data.table(
  No = diseases,
  Frequency = as.integer(tabulate(match(raw_data$Category_No, diseases), nbins = length(diseases))),
  Codes = diseases,
  Names = paste("Disease", diseases),
  block = rep(c("Block_A", "Block_B", "Block_C"), length.out = length(diseases))
)

# ---- run ----
res <- run_block_merge(
  raw_data = raw_data,
  codebook = codebook,
  prevalence_cutoff = 0.02,  # toy data needs a higher cutoff
  n_resamples = 10,          # keep it fast for README
  k_max_per_block = 3,
  write_files = FALSE,
  verbose = TRUE
)

print(res)

```

## Main function

### `run_block_merge()`

This is the main entry point. It performs:

1. Prevalence filtering (remove diseases with frequency below `prevalence_cutoff`)
2. Build a patient × disease sparse matrix
3. For each block:
   - resample patients `n_resamples` times (fraction `sample_frac`)
   - compute global structural profiles (block rows × all-diseases columns)
   - run k-means for candidate K values
   - build consensus co-clustering matrices
   - compute PAC(K) and choose `K*`
   - optionally trigger the homogeneity gate to allow K=1 full-block merge
4. Apply merges and write outputs (optional)


## Outputs

When write_files = TRUE, the following files are written to output_dir:
- PAC_Curves_ByBlock.csv: PAC(K) per block
- PAC_Summary_Kstar.csv: selected K* per block (plus homogeneity diagnostics if enabled)
- Final_Block_Memberships.csv: final cluster membership within each block
- Merge_Update_Map.csv: mapping from old disease IDs to merged disease IDs
- Final_Merged_Data.csv: updated patient–disease table after merging
- Final_Merged_Dict.csv: updated lineage dictionary for merged nodes
The returned object (class icdmerge_result) contains the same outputs in memory:
pac_curves, kstar, membership, update_map, merged_data, lineage, runtime_sec

## Notes

- K_star = 1 does not necessarily mean the block is “homogeneous”; it means the algorithm did not find stable, low-PAC splits that justify within-block merging under the current settings.
- If many blocks yield K_star = 2, consider revisiting the PAC interval (pac_x1, pac_x2), pac_no_merge, or the structural-profile design.


## License

MIT.

