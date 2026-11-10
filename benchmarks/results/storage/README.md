
# Storage benchmark

Comparison of tabulated DFTB reference data and four
Fourier–Morse potential representations.

## Reproduction

Run:

python benchmark_storage.py

## Reference data

- Interaction classes: EP, EA, OP, OA
- Selected axial displacement: zeta = 0 Å
- Selection tolerance: 1e-8 Å
- Reference dataset: DOI 10.5281/zenodo.21904447

## Analytical models

- Compact: constant alpha
- Full: Fourier-dependent alpha
- Pruned: selected coefficients retained, with zero padding
- Unpruned: all symmetry-allowed coefficients retained

## Results

- tabulated_data.csv: original files and selected zeta slice
- chimorse_models.csv: model sizes and coefficient counts
- chimorse_by_class.csv: per-interaction coefficient counts

## Interpretation

The results describe serialized storage and approximate
model-data memory requirements.

The tabulated zeta=0 text size is calculated from the
selected rows, not from a separately generated file.

Estimated C++ memory values are not measurements of
actual allocated memory.
