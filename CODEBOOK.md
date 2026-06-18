# Codebook

## `data/processed/analytic_model_data.csv`

This is the de-identified analytic dataset used by
`scripts/02_reproduce_from_processed.R`.

Each row represents one linked Saber 11 -> Saber Pro 2022 observation after
filtering to valid gender, valid STEM classification, and non-missing math,
science, and language quintiles.

| Variable | Type | Description |
|---|---:|---|
| `gender` | factor/string | Student gender, harmonized to `Men` and `Women`. |
| `nse_model` | factor/string | Saber 11 individual socioeconomic level used as a model covariate. Missing values are coded as `NO INFO` in the raw-build script. |
| `sb11_period` | factor/string | Saber 11 period/cohort used as a model covariate. |
| `math_q` | ordered category | Saber 11 mathematics performance quintile, `Q1` lowest to `Q5` highest. |
| `science_q` | ordered category | Saber 11 natural sciences performance quintile, `Q1` lowest to `Q5` highest. |
| `language_q` | ordered category | Saber 11 language/critical reading performance quintile, `Q1` lowest to `Q5` highest. |
| `STEM_AREA` | integer | Binary outcome. `1` if the Saber Pro higher education program is classified as STEM; `0` otherwise. |

## Quintile Construction

Quintiles are constructed within Saber 11 cohort/period in the raw-data build.
This keeps performance ranks comparable inside each Saber 11 test administration.

Included quintiles:

- `math_q`: based on `PUNT_MATEMATICAS`.
- `science_q`: based on `PUNT_C_NATURALES`.
- `language_q`: based on `PUNT_LECTURA_CRITICA`, renamed to language for
  continuity with older scripts.

## STEM Classification

Programs are classified as STEM when the SNIES area of knowledge is one of:

- `Ingenieria, arquitectura, urbanismo y afines`
- `Matematicas y ciencias naturales`
- `Agronomia, veterinaria y afines`
- `Ciencias de la salud`

All other classified areas are coded `0`.

## `data/processed/analytic_model_cell_counts.csv`

Aggregated cell counts by gender and the three performance quintiles. This file
is not required by the reproduction script, but it is included for transparent
checking of how the model-ready data are distributed across performance cells.
