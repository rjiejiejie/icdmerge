# icdmerge 0.1.0

## Initial release

- Implemented block-wise ICD node merging using:
  - global structural profiles (block diseases × all diseases),
  - patient resampling and consensus co-clustering matrices,
  - PAC-based selection of block-specific `K*`,
  - optional homogeneity gate allowing `K=1` full-block merge when `K=2` consensus is nearly uniform.
- Added the main user-facing function:
  - `run_block_merge()`
- Added outputs:
  - PAC curves and K* summaries by block,
  - final memberships,
  - merge mapping table,
  - merged patient–disease data and updated lineage dictionary.

