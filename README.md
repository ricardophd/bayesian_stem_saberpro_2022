# Gender, Saber 11 Performance, and STEM Selection in Saber Pro 2022

This repository contains the reproducible analysis bundle for estimating whether
Saber 11 performance quintiles are associated with choosing a STEM field in
higher education, with special attention to whether women and men with the same
high-school performance quintile have different probabilities of selecting STEM.

The main reproducible dataset included here is a de-identified, model-ready CSV
derived from the linked Saber 11 -> Saber Pro 2022 data. It contains only the
variables needed for the Bayesian model: gender, socioeconomic level, Saber 11
cohort, math/science/language quintiles, and the binary STEM outcome.

## Repository Structure

```text
.
|-- README.md
|-- CODEBOOK.md
|-- data
|   |-- processed
|   |   |-- analytic_model_data.csv
|   |   `-- analytic_model_cell_counts.csv
|   `-- raw
|       `-- README.md
|-- outputs
|   |-- figures
|   |-- models
|   `-- tables
`-- scripts
    |-- 00_install_packages.R
    |-- 01_build_from_raw_and_analyze.R
    `-- 02_reproduce_from_processed.R
```

## Quick Reproduction

From the repository root:

```r
source("scripts/00_install_packages.R")
source("scripts/02_reproduce_from_processed.R")
```

Or from a terminal:

```bash
Rscript scripts/00_install_packages.R
Rscript scripts/02_reproduce_from_processed.R
```

The reproduction script reads:

```text
data/processed/analytic_model_data.csv
```

and regenerates:

```text
outputs/tables/04_observed_beta_binomial_stem_rates.csv
outputs/tables/05_observed_beta_binomial_gender_contrasts.csv
outputs/tables/06_bayesian_logistic_coefficients.csv
outputs/tables/07_bayesian_predicted_probabilities.csv
outputs/tables/08_bayesian_gender_contrasts_same_quintile.csv
outputs/tables/09_bayesian_one_quintile_offset_contrasts.csv
outputs/figures/*.png
outputs/figures/*.pdf
outputs/models/bayesian_logistic_laplace_fit.rds
```

## Analysis Summary

Outcome:

- `STEM_AREA = 1` if the Saber Pro program is classified as STEM.
- STEM fields are engineering/architecture/urbanism, mathematics and natural
  sciences, agronomy/veterinary, and health sciences.

Predictors:

- Gender: men and women.
- Saber 11 math quintile.
- Saber 11 natural science quintile.
- Saber 11 language/critical reading quintile.
- Socioeconomic level.
- Saber 11 cohort/period.

Model:

- Bayesian logistic regression estimated with a Laplace approximation.
- Weakly informative normal priors:
  - Intercept: Normal(0, 2.5)
  - Main effects: Normal(0, 1.5)
  - Gender-by-quintile interactions: Normal(0, 1.0)
- Posterior draws default to 8,000.

The primary inferential estimand is the posterior probability difference:

```text
Pr(STEM | Women, quintile q) - Pr(STEM | Men, quintile q)
```

for math, science, and language quintiles.

## Data Notes

The raw ICFES files are not included in this GitHub bundle because they are very
large and include student-level link keys. The included processed data is enough
to reproduce the published descriptive and Bayesian inferential results without
exposing individual identifiers.

To rebuild the processed data from raw files, see:

```text
scripts/01_build_from_raw_and_analyze.R
data/raw/README.md
```

That raw rebuild script expects the same folder layout used in the original
project workspace.

## Main Outputs

The most important table is:

```text
outputs/tables/08_bayesian_gender_contrasts_same_quintile.csv
```

The two main publication figures are:

```text
outputs/figures/02_bayesian_predicted_probability_by_quintile_gender.png
outputs/figures/03_bayesian_gender_difference_same_quintile.png
```

Figures are exported without embedded titles; titles and captions are added in
the manuscript. Facet labels are in English.

## Reproducibility Controls

Optional environment variables:

- `BAYES_STEM_DRAWS`: number of posterior draws. Default: `8000`.
- `BAYES_STEM_SEED`: random seed. Default: `20260616`.
- `BAYES_STEM_ROPE`: equivalence threshold in probability points. Default: `0.01`.

Example:

```bash
BAYES_STEM_DRAWS=2000 Rscript scripts/02_reproduce_from_processed.R
```

## Important Scope Limitation

The analytic population is students linked to Saber Pro 2022. The analysis
therefore estimates STEM selection among students observed in higher education
and taking Saber Pro, not among all Saber 11 test takers.
