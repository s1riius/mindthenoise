# Data

## raw/ — as extracted or as published; no script writes here

| File | Contents | Source |
|---|---|---|
| `vegfr_compound_ids.csv` | the 24 VEGFR1 compounds of Table S13 with their ChEMBL identifiers | ChEMBL |
| `vegfr_records.csv` | the same 24 compounds with the Ki value, the number of IC50 reports and every individual IC50 value (Table S13) | ChEMBL |
| `vegfr1_raw.csv`, `vegfr2_raw.csv` | every retrieved record of the current extractions (October 2026), with compound, document and assay identifier (Supplement S9) | ChEMBL web services |
| `table_s11_source/vegfr1_chembl1868_records_2024.txt` | the record table for CHEMBL1868 from which Table S13 was compiled (4027 records: compound, activity type, units, value), March 2024 | see Supplementary S7.1 |
| `table_s11_source/vegfr1_table_s11_workbook.xlsx` | the workbook that reduces those records to the 24 compounds of Table S13 (selection, unit conversion, IC50 averaging) | this work |
| `landrum_vegfr2_ic50.csv.gz` | the maximally curated VEGFR2 (CHEMBL279) IC50 dataset of Landrum and Riniker, as distributed: assay and compound identifier, SMILES and pchembl value (Supplement S9.6). Their release provides no K_I data for this target | Landrum and Riniker, ChEMBL32 |
| `sources.txt` | DOI of the source publication for each activity-coefficient workbook | — |
| `literature/tableS0_fig_s2a_t4l_alkylbenzenes.csv` | alkyl benzenes binding T4L L99A: polarisability, K_D, dG (Table S2, Figure S2A) | refs 28–30; values read from the published figures |
| `literature/tableS1_host_guest_salts.csv` | adamantane guest in a macromolecular host, association constants per salt (Table S3) | literature |
| `literature/tableS13_thermolysin_itc_krimmer_klebe.csv` | ITC of ligand binding to thermolysin at 50–300 uM protein (Table S15, Figure S10) | Krimmer & Klebe |
| `literature/tableS14_calbindin_d9k_n56a_linse.csv` | Ca2+ binding to site II of calbindin D9k N56A/P43M (Table S16, Figure S11) | Linse et al., Table 1, reproduced in full |
| `literature/tableS15_benzene_activity_coefficients_mcdevit_long.csv` | benzene activity coefficients in NaCl, Na2SO4 and N(Me)4Br at 25 C (Table S17, Figure S2C) | McDevit & Long; values read from the published figure |
| `structures/*.cif` | crystal structures used for Figure 2 and the ion-site analysis | RCSB PDB |
| `structures/data_summary.xlsx` | structure-by-structure summary of the modelled ions and occupancies | this work |

Values marked as read from published figures are listed in Supplementary S8.1 with the
reading precision; they are used only in the illustrative figures of Supplement S2 and in no
statistical analysis.

## processed/ — the analysis inputs

| File | Contents |
|---|---|
| `vegfr.csv` | the dataset of Table S13 (VEGFR1, ChEMBL32 extraction of Landrum and Riniker), 24 compounds: Ki (nM), mean IC50 (nM), SD of the IC50 reports (0 = single report) |
| `lindberg.csv` | 37 compounds: Kd for DYRK1A, CLK1, GSK-3b (nM); radiometric and ADP-Glo IC50 (uM); empty = not available |
| `vegfr1.csv`, `vegfr2.csv` | the current VEGFR1 (22 compounds) and VEGFR2 (128 compounds) extractions, same columns as `vegfr.csv` (Supplement S9) |
| `vegfr1_2023_uncensored.csv` | `vegfr.csv` without the four compounds carrying values that ChEMBL reports as bounds (20 compounds; Supplement S9.2.2) |
| `vegfr1_rule_applied.csv` | `vegfr.csv` with the two corrections of Supplement S9.2.2: CHEMBL262433 restored and the merged CHEMBL424672/CHEMBL424872 row removed (24 compounds) |

Censoring and range rules are applied inside the scripts, as described in Supplementary S7.1.

## workbooks/ — activity-coefficient calculations

Nine workbooks (alpha-chymotrypsin, HEWL, BSA). Each contains the source data with its
reference, all formulas, the fits, and a "Corrections log" sheet listing the changes made in
revision 1.
