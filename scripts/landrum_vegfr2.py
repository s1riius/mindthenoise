"""landrum_vegfr2.py - reproducibility of a single metric across assays, from the
maximally curated VEGFR2 (CHEMBL279) IC50 data of Landrum and Riniker (Supplementary S9.6).
Their release provides IC50 but no K_I data for this target, so the cross-metric
comparison of Supplement S9.3 cannot be repeated on it; what their data measure is the
assay-to-assay reproducibility of one metric, which is the quantity the null models assume.
Their curation keeps one assay per document, so distinct assays are distinct publications.

Usage:  python landrum_vegfr2.py ../data/raw/landrum_vegfr2_ic50.csv.gz
"""
import sys, gzip, csv, math, itertools, collections
import numpy as np
import affinity_agreement as A                      # S_REP, S_IC: the null-model noise levels

PATH = sys.argv[1] if len(sys.argv) > 1 else 'landrum_vegfr2_ic50.csv.gz'
NBOOT, SEED = 2000, 2026
L = math.log(10)                                    # pchembl values are base-10

rows = list(csv.DictReader(gzip.open(PATH, 'rt')))
rec = collections.defaultdict(lambda: collections.defaultdict(list))
for r in rows:
    rec[r['compound_chembl_id']][r['assay_chembl_id']].append(float(r['pchembl_value']))
multi = {c: {a: sum(v) / len(v) for a, v in m.items()} for c, m in rec.items() if len(m) > 1}


def deltas(compounds):
    """Signed differences in ln units between assay means of the same compound."""
    out = []
    for c in compounds:
        m = multi[c]
        out += [(m[a] - m[b]) * L for a, b in itertools.combinations(sorted(m), 2)]
    return np.array(out)


D = deltas(multi); AD = np.abs(D)
sd = math.sqrt((D ** 2).mean() / 2)                 # two measurements, each with SD sd
rng = np.random.default_rng(SEED); keys = list(multi)
boot = [math.sqrt((deltas([keys[i] for i in rng.integers(0, len(keys), len(keys))]) ** 2).mean() / 2)
        for _ in range(NBOOT)]
lo, hi = np.percentile(boot, [2.5, 97.5])

print(f'records {len(rows)}; compounds {len(rec)}; assays {len({r["assay_chembl_id"] for r in rows})}')
print(f'compounds measured in >1 assay: {len(multi)} across '
      f'{len({a for c in multi for a in multi[c]})} assays; {len(D)} pairwise comparisons')
print(f'SD of a single measurement: {sd:.2f} ln units [{lo:.2f}, {hi:.2f}] '
      f'(bootstrap over compounds, {NBOOT} resamples, seed {SEED})')
print(f'median |difference| {np.median(AD):.2f} ln ({math.exp(np.median(AD)):.1f}-fold); '
      f'90th percentile {np.percentile(AD, 90):.2f} ln ({math.exp(np.percentile(AD, 90)):.0f}-fold); '
      f'maximum {AD.max():.2f} ln ({math.exp(AD.max()):.0f}-fold)')
for t in (2, 5, 10):
    print(f'pairs differing by more than {t:2d}-fold: {(AD > math.log(t)).mean():.2f}')
print(f'observed RMS difference: {math.sqrt((D ** 2).mean()):.2f} ln units')
for lab, s in (('within-assay repeatability', A.S_REP), ('inter-laboratory reproducibility', A.S_IC)):
    e = s * math.sqrt(2)
    frac = math.erfc(math.log(10) / (e * math.sqrt(2)))
    print(f'null, {lab} ({s:.2f} ln per metric): expected RMS difference {e:.2f} ln, '
          f'expected fraction beyond 10-fold {frac:.2f}')
