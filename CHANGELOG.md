# Changelog

## v1.0 — version accompanying the submitted manuscript (revision 1)

Corrections made during revision 1. Each activity-coefficient workbook also carries a
"Corrections log" sheet listing its own changes.

- BSA activity coefficients recomputed; an error in the earlier virial calculation corrected
  (Figures S3D and S9).
- Figure S1 redrawn: the reference line is the complete-cancellation case, K_D^gamma = 1.
- Figure S16 axis corrected to micromolar (Figure S17 before the renumbering of this revision).
- Figure S2A legend corrected (reference number).
- Tables 1 and S4 corrected against the original publications and the Supporting Information
  of Lindberg et al.: ADP-Glo kinase amounts, DMSO concentration, and the radiometric rows.
- Protein identity in Figure S11 corrected to calbindin D9k N56A.
- Statistics redone with errors-in-variables regression, matched null models and a direct
  estimate of rank inversions.
- New Supplement S9: the analysis repeated on a current ChEMBL extraction of VEGFR1 and
  extended to VEGFR2, with the same code, selection rules, null models and seeds
  (`scripts/chembl_extract.py`, `scripts/chembl_clean.py`, `scripts/fig_s9supp.py`;
  Figures S19-S22, Tables S17 and S18).
- Figure files renamed to match the figure numbers of the revised supplement:
  FigS1.png, FigS14.png, FigS16.png, FigS17.png, FigS18.png and FigS19-S22 for the
  Supplement S9 figures.
- New Supplement S9.6: assay-to-assay reproducibility of IC50 for VEGFR2 in the maximally
  curated extraction of Landrum and Riniker (`scripts/landrum_vegfr2.py`).
