"""chembl_clean.py - build the analysis inputs for the current ChEMBL extractions (S9.1).

Selection, identical for both targets: unqualified records (standard relation '='),
standard units nM, positive value; a compound is kept only if it has at least one Ki
and at least one IC50 record after filtering.

Re-tabulated K_I records are removed by SOURCE PROVENANCE, not by value equality:
RETABULATED lists the (document, assay) sources whose K_I records reproduce, to the
digit, IC50 values already recorded for the same compounds in earlier documents. Two
alternative rules are provided for the sensitivity analysis of S9.1:

    python3 chembl_clean.py              # provenance rule (primary)
    python3 chembl_clean.py equality     # drop any Ki numerically equal to an IC50
    python3 chembl_clean.py none         # no exclusion

Outputs <tag>_clean.csv, <tag>_equality.csv or <tag>_none.csv respectively.
"""
import pandas as pd, statistics as st, sys
from collections import defaultdict

# K_I sources identified as re-tabulations of earlier IC50 records (S9.1).
RETABULATED = {'CHEMBL5726637'}


def build(tag, rule='provenance'):
    d = pd.read_csv(f'{tag}_raw.csv')
    k = d[(d.standard_relation == '=') & (d.standard_units == 'nM') &
          d.standard_value.notna() & (d.standard_value > 0)]
    name = {r.molecule_chembl_id: r.molecule_pref_name for r in d.itertuples()
            if isinstance(r.molecule_pref_name, str)}
    rec = defaultdict(lambda: {'Ki': [], 'IC50': []})
    for r in k.itertuples():
        if r.standard_type == 'Ki' and rule == 'provenance' and r.document_chembl_id in RETABULATED:
            continue
        rec[r.molecule_chembl_id][r.standard_type].append(float(r.standard_value))
    excluded = (len(k[(k.standard_type == 'Ki') & k.document_chembl_id.isin(RETABULATED)])
                if rule == 'provenance' else 0)
    rows, dup = [], 0
    for cid, v in rec.items():
        ic = v['IC50']
        ki = [x for x in v['Ki'] if x not in ic] if rule == 'equality' else v['Ki']
        dup += len(v['Ki']) - len(ki)
        if not (ki and ic):
            continue
        rows.append([name.get(cid, cid), round(sum(ki)/len(ki), 3),
                     round(sum(ic)/len(ic), 3),
                     round(st.stdev(ic), 3) if len(ic) > 1 else 0.0])
    rows.sort()
    suffix = {'provenance': '_clean', 'equality': '_equality', 'none': '_none'}[rule]
    pd.DataFrame(rows, columns=['compound', 'Ki', 'IC50', 'stdev']).to_csv(f'{tag}{suffix}.csv', index=False)
    print(f'{tag} [{rule}]: {len(k)} unqualified nM records; '
          f'{excluded} K_I records excluded by source, {dup} by value equality; '
          f'{len(rows)} compounds with both metrics')


if __name__ == '__main__':
    rule = sys.argv[1] if len(sys.argv) > 1 else 'provenance'
    for t in ('vegfr1', 'vegfr2'):
        build(t, rule)
