# Raw Data Placement

Raw ICFES files are not included in this GitHub bundle because they are large
and contain student-level link keys. The processed, de-identified dataset in
`data/processed/analytic_model_data.csv` is sufficient to reproduce the Bayesian
analysis.

If you have authorized access to the raw files and want to rebuild the processed
dataset, place the raw folders in the same layout expected by
`scripts/01_build_from_raw_and_analyze.R`:

```text
cruces-SB11-SBPRO/
|-- llaveSB112006_SBPro2019.txt
`-- Llave_Saber11_SaberPro.txt

dataSaberPro/
|-- SaberPro_Genericas_20221.TXT
|-- SaberPro_Genericas_20222.TXT
`-- programas_SNIES.xlsx

dataSaber11/
|-- SB11_20142.txt
|-- SB11_20151.txt
|-- SB11_20152.txt
|-- SB11_20161.txt
|-- SB11_20162.txt
|-- SB11_20171.TXT
|-- SB11_20172.TXT
|-- SB11_20181.TXT
`-- SB11_20182.TXT
```

The raw-build script should be run from a project root containing those folders.
It creates linked analytic tables and the same final figures/tables as the
processed-data reproduction script.
