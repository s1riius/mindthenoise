# Mind the Noise — data and analysis package

Data and code accompanying *Mind the Noise: Thermodynamic Non-Ideality as a Possible
Contributor to the Excessive Noise in Binding SAR Data*.

Tag `v1.0` is the version corresponding to the submitted manuscript.

## What is here

```
data/raw/          sources as extracted or as published, not modified by any script
data/processed/    the analysis inputs
data/workbooks/    activity-coefficient spreadsheets (source data, formulas, fits, corrections log)
scripts/           extraction and analysis code
results/           generated tables and figures (regenerable; not edited by hand)
```

Nothing in `data/raw/` is written to by any script, so the raw-to-processed step is
itself inspectable.

## Map from manuscript to code

| Manuscript item | Produced by | Inputs |
|---|---|---|
| Tables S8–S10 | `scripts/affinity_agreement.py vegfr` | `data/processed/vegfr.csv` |
| Tables S11–S12 | `scripts/affinity_agreement.py lindberg` | `data/processed/lindberg.csv` |
| Tables S5, S7 | `scripts/affinity_agreement.py ratios` | both CSVs |
| Figures S14, S17, S18 | `scripts/figures.py` | both CSVs |
| Figures S3D, S9 | `scripts/bsa_gamma.py` | `data/workbooks/BSA_NaCl_corrected.xlsx` |
| Figures S1, S16 | `scripts/fig_s1.py`, `scripts/fig_s16.py` | none (computed) |
| Figures S19–S22 | `scripts/fig_s9supp.py` | `data/processed/vegfr1.csv`, `data/processed/vegfr2.csv` |
| Tables S18 and S20, Supplement S9 statistics | `scripts/affinity_agreement.py` (`A.vegfr(path)`, see `run_all.sh`) | `data/processed/vegfr1.csv`, `vegfr2.csv`, `vegfr1_2023_uncensored.csv` |
| `data/raw/vegfr1_raw.csv`, `vegfr2_raw.csv` and the processed files from them | `scripts/chembl_extract.py` (retrieval), `scripts/chembl_clean.py` (selection rules) | ChEMBL web services; run once, not part of `run_all.sh` |
| Supplement S9.6 statistics | `scripts/landrum_vegfr2.py` | `data/raw/landrum_vegfr2_ic50.csv.gz` |
| all of the above | `scripts/run_all.sh` | everything in `data/` |
| Figures S2–S9 | `data/workbooks/*.xlsx` | literature data cited in each workbook |
| Table S13 | `data/processed/vegfr.csv`, `data/raw/vegfr_compound_ids.csv`, `data/raw/vegfr_records.csv` | `data/raw/table_s11_source/` |
| Table S14 | `data/processed/lindberg.csv` | Lindberg et al., Supporting Information |

## Reproducing

```
conda env create -f environment.yml && conda activate mindthenoise
cd scripts && ./run_all.sh
```

which writes every table and figure into `results/`. The individual steps are:

```
python affinity_agreement.py vegfr | lindberg | ratios | extras
python figures.py          # Figures S14, S17, S18
python bsa_gamma.py [fig]  # table; Figures S3D and S9
python fig_s1.py           # Figure S1
python fig_s16.py          # Figure S16
python fig_s9supp.py vegfr1.csv vegfr2.csv   # Figures S19-S22 (~5 min)
```

Random seeds are fixed (2026 for VEGFR, 2027 for Lindberg, 11 for the figures), so the
reported numbers are reproduced exactly. Bootstrap and null-model ranges may differ in the
second decimal place under other seeds.

## Data provenance

- **VEGFR2 (CHEMBL279).** Taken from the published ChEMBL32 maximal-curation extraction of
  Landrum and Riniker (https://github.com/rinikerlab/overlapping_assays). The curation was
  performed by those authors; the SQL queries and the maximal-curation rules (unqualified
  values in nM with a non-null pChEMBL value and no data validity comment; assay confidence
  score 9; no mutant or variant assays; one assay per document and target) are theirs and are
  documented in their repository. From their released data we selected the records for
  CHEMBL279 with standard types Ki and IC50, and kept the compounds having at least one of
  each. Our selection and processing steps are described in Supplementary S7.1.
  `data/raw/vegfr_compound_ids.csv` maps each compound to its ChEMBL identifier.
- **Lindberg et al.** Values as published in the paper and its Supporting Information, with
  censoring flags preserved.
- **Activity coefficients.** Literature osmometry and light-scattering data; the source DOI
  for each workbook is given in `data/raw/sources.txt` and inside each workbook.
- **Values read from published figures** are listed in Supplementary S8.1, with the reading
  precision stated in the caption of the corresponding table.

## Licence

Data: CC-BY-4.0. Code: MIT.
