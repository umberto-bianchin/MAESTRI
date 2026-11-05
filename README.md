# MAESTRI: Multi-phylogeny-Aware Extraction of Significant TRajectorIes

This repository contains MAESTRI, our algorithm to discover statistically significant conserved evolutionary trajectories of alterations when each patient is described by a weighted *set* of candidate phylogenetic trees rather than by a single tree.
Every candidate tree of a patient contributes with its weight, so the evidence of a patient for a trajectory is a quantity in [0,1], and no patient is counted more than once.
Significance is assessed under a null model that relabels the alterations of a patient coherently across all of its trees, and multiple testing is corrected with resampling-based procedures (Westfall-Young FWER and an empirical FDR).

MAESTRI extends MASTRO (L. Pellegrina and F. Vandin, Bioinformatics 2022, https://doi.org/10.1093/bioinformatics/btac467), which solves the single-tree version of the problem.
The single-tree baseline of the experiments is MAESTRI run on one sampled tree per patient (`run_pipeline.py --single_tree_only`); the original MASTRO implementation is at https://github.com/VandinLab/MASTRO.

Please address questions and bug reports to Umberto Bianchin (umberto.bianchin@phd.unipd.it).
For the original MASTRO, to Leonardo Pellegrina (leonardo.pellegrina@unipd.it) and Fabio Vandin (fabio.vandin@unipd.it).

## Contents

- `MAESTRI/`: the implementation of MAESTRI and the scripts that reproduce the experiments of the paper.
- `MAESTRI/scripts/`: one script per experiment of the paper, each documented in its header (`scripts/dei/` holds the SLURM wrappers used to run them on a cluster).
- `MAESTRI/lcm53/`: LCM version 5.3 by Takeaki Uno (http://research.nii.ac.jp/~uno/codes.htm), used to mine frequent itemsets.
- `data/breastCancer.npy`: the breast cancer cohort of the paper, 1315 patients and 37,809 candidate trees reconstructed with SPRUCE and analyzed in CloMu (S. Ivanovic and M. El-Kebir, Genome Research 2023).
- `data/TRACERx/`: the TRACERx non-small-cell lung cancer cohort of the paper, 394 patients and 469 candidate trees, and the scripts that rebuild it from the TRACERx 421 release (see `data/TRACERx/README.md`).

## Requirements

Python 3.10 or later (tested with 3.10 and 3.12) and the packages in `requirements.txt`:

```
pip install -r requirements.txt
```

plus a C compiler, `make` and `perl`. LCM is compiled from `MAESTRI/lcm53/` automatically the first time it is needed (or by hand with `make -C MAESTRI/lcm53`), and its `transnum.pl` maps edge labels to integer item ids.

## Running MAESTRI

All commands are run from the `MAESTRI/` directory. The entry point is `run_pipeline.py`; `-h` lists its parameters:

```
usage: run_pipeline.py [-h] (--npy NPY | --graphs GRAPHS) [--owner OWNER]
                       [--weights WEIGHTS] [--sigma SIGMA] [--seed SEED]
                       [--theta_list THETA_LIST]
                       [--min_mine_sigma MIN_MINE_SIGMA] [--lcmdir LCMDIR]
                       [--outdir OUTDIR] [--keep_gl] [--significance]
                       [--single_tree_only] [--sig_null {indep,perm}]
                       [--sig_mc_cutoff SIG_MC_CUTOFF]
                       [--sig_mc_samples SIG_MC_SAMPLES]
                       [--sig_exact_when_cheaper] [--sig_n_jobs SIG_N_JOBS]
```

A typical run on the breast cancer cohort:

```
python3 run_pipeline.py --npy ../data/breastCancer.npy --sigma 5 --theta_list 0.5,1.0 \
    --min_mine_sigma 2 --significance --sig_exact_when_cheaper --sig_mc_samples 20000 \
    --outdir results/demo
```

`--npy` gives a dataset as an object array indexed `data[patient][tree][edge] = [parent, child]`.
Alternatively, `--graphs` takes a file of transactions, one tree per line, in the input format below, and `--owner` a file with the patient of each line (any label, one per line); without `--owner` every line is a separate patient, which is the single-tree setting.
By default the trees of a patient with `M_i` trees have weight `1/M_i`; `--weights` gives other per-tree weights, which must sum to 1 within each patient.

`--sigma` is the minimum support of the reported trajectories.
Support is measured in patients, and with weighted trees it is a sum of per-patient quantities in [0,1].

`--theta_list` sets the consensus thresholds: a trajectory is theta-frequent for a patient when the trees containing it carry at least a fraction `theta` of that patient's weight.

One run mines four families of trajectories, whose output files carry a short prefix:

1. `alg0_`, the single-tree baseline: one tree per patient is sampled at random (controlled by `--seed`) and mined without weights.
2. `alg1_`, the expected-support family: all candidate trees are pooled, each with its weight, and mined with weighted frequent-itemset mining.
3. `alg2_`, the theta-frequent family: trajectories that are theta-frequent in at least `sigma` patients.
4. `alg3_`, the theta-maximal family: the maximal elements of the theta-frequent family.

## Significance

Two notions of support are tested, each on its own mined family, so that observed and null distributions refer to the same set of hypotheses.
The *expected support* `s_exp(P)` is the expected number of patients carrying `P` if one tree were drawn per patient according to the weights.
The *theta-consensus support* `s_theta(P)` is the number of patients for which `P` is theta-frequent.
The expected-support test is applied to the `alg0_` and `alg1_` families, the theta-consensus test to the `alg2_` and `alg3_` families.

`--sig_null` selects the null model.
`perm` (the default, used in the paper) draws one random permutation of the alteration set of each patient and applies it to all of that patient's candidate trees, so the trees of a patient stay coupled under the null.
`indep` places each alteration of the trajectory independently and uniformly at random, as the null of the original MASTRO.

The per-patient null distributions are computed exactly when a patient has at most `--sig_mc_cutoff` candidate trees.
Above that they are estimated by Monte Carlo with `--sig_mc_samples` draws, unless `--sig_exact_when_cheaper` is set and the number of distinct placements is at most `--sig_mc_samples`, in which case they are enumerated exactly.
Monte-Carlo p-values use the add-one estimator `(B*p + 1)/(B + 1)`, so none is exactly zero, and none can fall below `1/(B+1)`: when a corrected threshold lands near that floor, raise `--sig_mc_samples`.

Multiple testing is corrected by `run_wy_correction_multitree.py`, which generates resampled datasets under the null, re-mines each test's family on every resample, and writes the Westfall-Young thresholds and the empirical FDR curve.
`plot_fdr_v2.py` draws the FDR curves; the reported max-k is the largest k before the curve first crosses the level, since the estimated FDR is not monotone in k.
`empirical_fwer_multitree.py` measures the empirical FWER on null datasets derived from a cohort, to verify that the correction achieves its nominal level.

## Input format

MAESTRI takes as input a file containing, in each line, the edges of the complete expanded tumor graph of one tumor phylogenetic tree, separated by a space.
Alterations can be of any type (e.g., SNV or CNA) but must have a distinct id in each tree (e.g., `TP53_SNV` and `TP53_DEL` are allowed in the same tree).
The complete expanded tumor graph has three types of edges:

1. a directed edge `A->-B` denotes that alteration `A` is an ancestor of alteration `B`;
2. an undirected edge `A-/-B` denotes that `A` and `B` belong to different branches of the tree;
3. an undirected edge `A-?-B` denotes that `A` and `B` belong to the same node, so their order is not resolved.

The last two edges are undirected, so their alterations are sorted alphabetically (`A-?-B` is correct, `B-?-A` is not).
Edges incident to the germline node are not needed to find trajectories of at least two alterations; otherwise, use `g` to denote the germline node.

For example, the following file describes three trees with edges `{g->-A, A->-B, B->-C}`, `{g->-A, A->-B, A->-C}` and `{g->-A, g->-B, A->-C}`:

```
A->-B A->-C B->-C
A->-B A->-C B-/-C
A-/-B A->-C B-/-C
```

When a patient has several candidate trees, its trees occupy several lines, and the owner file (`--owner`) gives the patient of each line.

## Output format

`run_pipeline.py` writes under `--outdir`:

```
inputs/              graphs_all.txt, owner.txt (0-based patient index of each line),
                     weights_uniform.txt, weights.txt (the weights used),
                     graphs_sampled.txt, owner_sampled.txt, weights_sampled.txt (baseline)
runs/                per-family mining output, including the filtered families
significance/        per-trajectory p-values, one CSV per family
analysis/            pairwise Jaccard similarity and per-family unique patterns
dataset_stats.json   size of the cohort
```

The mined families list frequent trajectories in the same format as the input trees; each trajectory line is followed by a line listing the transactions in which it occurs.
The significance CSVs give, for each trajectory, its number of nodes and edges, its expected support `s_exp`, its theta-consensus support `s_theta` and the p-value of the test that applies to the family.
`run_wy_correction_multitree.py` writes `wy_thresholds.txt` (corrected thresholds at alpha = 0.01, 0.05 and 0.10), `min_pvalues.csv` (the minimum p-value of each resample) and `fdr_exp.csv` / `fdr_theta.csv` (the observed p-values in rank order with the estimated FDR).

## Reproducing the experiments

`MAESTRI/scripts/` holds one script per experiment of the paper: calibration of the error control (`run_exp_calibration_*.sh`), power on implanted trajectories (`run_exp_power.sh`), discovery on the two cohorts (`run_exp_discovery_*.sh`), the single-tree seed sweep (`run_exp_singletree_seeds.sh`), the ablation of the design choices (`run_ablation.sh`), and the comparison with POTTR (https://github.com/AlBi-HHU/POTTR), including running times (`run_pottr_*.sh`, `run_timing_pottr_vs_mastro.sh`).
Run them from `MAESTRI/` as `bash scripts/<name>.sh`; their defaults are the configuration used in the paper, their settings can be changed through the environment variables defined at the top of each script.
As a check that the toolchain works, on `data/breastCancer.npy` the expected-support family contains 2354 trajectories at `sigma = 2` and 372 at `sigma = 5`.

## License

MIT, see `LICENSE`.
LCM 5.3 (`MAESTRI/lcm53/`) is redistributed under its own terms, stated in `MAESTRI/lcm53/readme.txt`: free for academic use, commercial use requires contacting its author.
`data/TRACERx/create_tree_files.r` and `data/TRACERx/create_tracerx_trees.py` are adapted from POTTR (MIT License, (c) 2025 Sara Schulte).
