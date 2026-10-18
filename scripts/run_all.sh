#!/bin/sh
# Regenerates every table and figure in the paper that is produced by code.
# Run from this folder; writes figures into ../results.
set -e
cp ../data/processed/vegfr.csv ../data/processed/lindberg.csv .
python affinity_agreement.py vegfr    > ../results/tables_S8_S10.txt
python affinity_agreement.py lindberg > ../results/tables_S11_S12.txt
python affinity_agreement.py ratios   > ../results/tables_S5_S7.txt
python affinity_agreement.py extras   > ../results/extras.txt
python figures.py                                  # Figures S14, S17, S18
python bsa_gamma.py                   > ../results/bsa_gamma_table.txt
python bsa_gamma.py fig                            # Figures S3D, S9
python fig_s1.py                                   # Figure S1
python fig_s16.py                                  # Figure S16
cp ../data/processed/vegfr1.csv ../data/processed/vegfr2.csv ../data/processed/vegfr1_2023_uncensored.csv .
# Supplement S9: current extractions, and the 2023 dataset without the censored records
python -c "import affinity_agreement as A; A.vegfr('vegfr1.csv')"                 > ../results/vegfr1_agreement.txt
python -c "import affinity_agreement as A; A.vegfr('vegfr2.csv')"                 > ../results/vegfr2_agreement.txt
python -c "import affinity_agreement as A; A.vegfr('vegfr1_2023_uncensored.csv')" > ../results/vegfr1_2023_uncensored_agreement.txt
python fig_s9supp.py vegfr1.csv vegfr2.csv         # Figures S19-S22 (~5 min)
python landrum_vegfr2.py ../data/raw/landrum_vegfr2_ic50.csv.gz > ../results/landrum_vegfr2.txt   # Supplement S9.6
mv -f *.png ../results/ 2>/dev/null || true
rm -f vegfr.csv lindberg.csv vegfr1.csv vegfr2.csv vegfr1_2023_uncensored.csv
echo "done; see ../results"
