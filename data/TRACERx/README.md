# TRACERx cohort

`graphs_tracerx.txt` and `owner_tracerx.txt` are the TRACERx non-small-cell lung
cancer cohort used in the paper: 394 patients and 469 candidate trees (at most 24
per patient), one tree per line of `graphs_tracerx.txt` in the input format of
the main README, and the patient of each line in `owner_tracerx.txt`.
They are ready to use:

```
python3 run_pipeline.py --graphs ../data/TRACERx/graphs_tracerx.txt \
    --owner ../data/TRACERx/owner_tracerx.txt --sigma 5 --theta_list 0.5,1.0 \
    --significance --outdir results/demo_tracerx
```

## Source

The trees come from the TRACERx 421 study and were inferred with CONIPHER:

> A. M. Frankell et al. The evolution of lung cancer and impact of subclonal
> selection in TRACERx. *Nature* 616, 525-533 (2023).
> https://doi.org/10.1038/s41586-023-05783-5

The release files are distributed on Zenodo (https://doi.org/10.5281/zenodo.7822002,
`figurecode.zip`) under an academic-use-only data licence, so they are not
included in this repository.

## Rebuilding the two files

The preprocessing follows POTTR (https://github.com/AlBi-HHU/POTTR), with the
same filters: driver mutations only, labels `<gene>_<type>`, and identical trees
of a patient kept once. `create_tree_files.r` and `create_tracerx_trees.py` are
adapted from POTTR (MIT License, (c) 2025 Sara Schulte).

1. Download and unpack `figurecode.zip` from Zenodo, and copy
   `20221109_TRACERx421_phylogenetic_trees.rds` and
   `20221109_TRACERx421_mutation_table.rds` into `original_data/`.
2. Export the candidate trees of every patient (R with the `openxlsx` package):
   ```
   Rscript create_tree_files.r original_data processed_data
   ```
3. Convert them to the MASTRO format (needs `pyreadr` and `openpyxl`):
   ```
   python3 create_tracerx_trees.py -i processed_data -o data_tracerx \
       -mut original_data/20221109_TRACERx421_mutation_table.rds -s Hugo+Type -drivers
   ```
4. Merge the per-tree files into the cohort files:
   ```
   python3 merge_tracerx.py --indir data_tracerx --outdir .
   ```
