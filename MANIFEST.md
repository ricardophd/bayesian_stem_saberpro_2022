# Manifest

## Documentation

- `README.md`: project overview, reproduction commands, analysis summary.
- `CODEBOOK.md`: variable definitions and STEM/quintile construction notes.
- `data/raw/README.md`: optional raw-data placement instructions.

## Scripts

- `scripts/00_install_packages.R`: installs required R packages into `r_libs/`.
- `scripts/01_build_from_raw_and_analyze.R`: rebuilds the linked dataset and
  full analysis from authorized raw ICFES/SNIES files.
- `scripts/02_reproduce_from_processed.R`: reproduces the main Bayesian
  analysis from the included processed data.

## Data

- `data/processed/analytic_model_data.csv`: de-identified analytic data used by
  the processed-data reproduction script.
- `data/processed/analytic_model_cell_counts.csv`: cell counts by gender and
  math/science/language quintile.

## Outputs

- `outputs/tables/`: CSV result tables.
- `outputs/figures/`: publication-ready PNG and PDF figures.
- `outputs/models/bayesian_logistic_laplace_fit.rds`: saved fitted model object.
