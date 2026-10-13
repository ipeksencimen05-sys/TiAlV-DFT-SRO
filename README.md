# Ti-Al-V DFT Energy and Chemical Short-Range Order

## Research Objective

This project investigates the relationship between chemical
short-range order (SRO), local atomic geometry, and DFT
energy in Ti-Al-V alloy configurations.

## Data Source

Public Ti-Al-V DFT dataset:
https://zenodo.org/records/15114121

Original file: ti64_combined_data_Apr08.extxyz

## Dataset Preparation

- Original configurations: 8,507
- Selected configurations: 50
- Composition: Ti13Al2V
- Structure family: HCP-labeled configurations
- Selection: Nearly identical cell volumes
- Sampling: Random state 42

## Descriptors

Chemical descriptors:
- SRO_Ti_Al
- SRO_Ti_V
- SRO_Al_V

Geometric descriptors:
- mean_neighbor_distance
- std_neighbor_distance
- min_shell_gap

Target:
- DFT_energy_per_atom (eV/atom)

## Preliminary Findings

The standard deviation of nearest-neighbor distances
showed a Pearson correlation of approximately 0.690
with DFT energy.

The individual SRO descriptors showed weaker linear
correlations in this initial sample.

These correlations do not establish causation.

## Reproducibility

Required Python packages:
- ase
- numpy
- pandas
- matplotlib

Place the original EXTXYZ file inside the raw/ folder.

Run:

    python scripts/prepare_tialv.py

The processed CSV file will be saved inside processed/.

## Limitations

The dataset contains perturbed atomic configurations.
Chemical-ordering effects and geometric distortions
have not yet been isolated.

SRO descriptors are calculated using the 12 nearest
periodic neighbors in small simulation cells.
They should not be interpreted as converged bulk SRO.
