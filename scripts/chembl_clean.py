import pandas as pd, statistics as st, sys
from collections import defaultdict
def build(tag):
    d = pd.read_csv(f'{tag}_raw.csv')
    k = d[(d.standard_relation == '=') & (d.standard_units == 'nM') &
          d.standard_value.notna() & (d.standard_value > 0)]
    name = {r.molecule_chembl_id: r.molecule_pref_name for r in d.itertuples()
            if isinstance(r.molecule_pref_name, str)}
    rec = defaultdict(lambda: {'Ki': [], 'IC50': []})
    for r in k.itertuples():
        rec[r.molecule_chembl_id][r.standard_type].append(float(r.standard_value))
    rows, dup = [], 0
    for cid, v in rec.items():
        ic = v['IC50']
        ki = [x for x in v['Ki'] if x not in ic]      # drop re-tabulated duplicates
        dup += len(v['Ki']) - len(ki)
        if not (ki and ic): continue
        rows.append([name.get(cid, cid), round(sum(ki)/len(ki), 3),
                     round(sum(ic)/len(ic), 3),
                     round(st.stdev(ic), 3) if len(ic) > 1 else 0.0])
    rows.sort()
    pd.DataFrame(rows, columns=['compound','Ki','IC50','stdev']).to_csv(f'{tag}_clean.csv', index=False)
    print(f'{tag}: {len(k)} unqualified nM records, {dup} Ki records dropped as '
          f're-tabulated duplicates, {len(rows)} compounds with both metrics')
for t in ('vegfr1','vegfr2'): build(t)
