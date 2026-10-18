#!/usr/bin/env python3
"""chembl_extract.py - extract and clean ChEMBL affinity data for one target.

Usage:
    python3 chembl_extract.py CHEMBL279  vegfr2     # VEGFR2 (KDR)
    python3 chembl_extract.py CHEMBL1868 vegfr1     # VEGFR1 (FLT1)

Writes, in the current folder:
    <name>_raw.csv   every record retrieved, unmodified
    <name>.csv       analysis input: compound, Ki, IC50, stdev   (nM)

Needs only the Python standard library (urllib). No pip install required.
Prints progress as it goes; the whole run takes a few minutes per target.

Selection rule, identical for every target:
  - activity types Ki, pKi, IC50, pIC50 for the given target;
  - pX converted as value_nM = 10**(9 - pX); uM x1000; M x1e9;
    records in other units are dropped and counted;
  - a compound is kept only if it has at least one Ki-type AND one IC50-type record;
  - Ki and IC50 are arithmetic means in nM; stdev is the sample SD of the IC50 records.
"""
import sys, csv, json, time, statistics as st
from collections import defaultdict
from urllib.request import urlopen
from urllib.error import URLError, HTTPError

BASE = 'https://www.ebi.ac.uk/chembl/api/data/activity.json'
FIELDS = ['molecule_chembl_id', 'molecule_pref_name', 'target_chembl_id',
          'document_chembl_id', 'assay_chembl_id', 'standard_type',
          'standard_units', 'standard_value', 'standard_relation']

def fetch(target):
    rows = []
    for typ in ('Ki', 'pKi', 'IC50', 'pIC50'):
        offset, n = 0, 0
        while True:
            url = (f'{BASE}?target_chembl_id={target}&standard_type={typ}'
                   f'&limit=1000&offset={offset}')
            for attempt in range(5):
                try:
                    with urlopen(url, timeout=120) as r:
                        data = json.loads(r.read().decode())
                    break
                except (URLError, HTTPError) as e:
                    if attempt == 4:
                        print(f'  ! giving up on {typ} offset {offset}: {e}', flush=True)
                        data = {'activities': [], 'page_meta': {}}
                        break
                    print(f'  . retrying {typ} offset {offset} ({e})', flush=True)
                    time.sleep(5 * (attempt + 1))
            acts = data.get('activities', [])
            rows += [{k: a.get(k) for k in FIELDS} for a in acts]
            n += len(acts)
            print(f'  {typ}: {n} records', flush=True)
            if not data.get('page_meta', {}).get('next'):
                break
            offset += 1000
    return rows

def to_nM(typ, units, value):
    if value in (None, ''):
        return None
    try:
        v = float(value)
    except ValueError:
        return None
    if typ in ('pKi', 'pIC50'):
        return 10 ** (9 - v)
    return {'nM': v, 'uM': v * 1000.0, 'M': v * 1e9}.get(units)

def main():
    if len(sys.argv) < 3:
        sys.exit('usage: python3 chembl_extract.py <TARGET_CHEMBL_ID> <name>')
    target, name = sys.argv[1], sys.argv[2]
    print(f'fetching {target} ...', flush=True)
    rows = fetch(target)
    if not rows:
        sys.exit('no records retrieved - check the target id and the network')

    with open(f'{name}_raw.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    print(f'wrote {name}_raw.csv ({len(rows)} records)', flush=True)

    d = defaultdict(lambda: {'K': [], 'I': [], 'name': ''})
    dropped = 0
    for r in rows:
        nm = to_nM(r['standard_type'], r['standard_units'], r['standard_value'])
        if nm is None:
            dropped += 1; continue
        cid = r['molecule_chembl_id']
        if r['molecule_pref_name']:
            d[cid]['name'] = r['molecule_pref_name']
        d[cid]['K' if r['standard_type'] in ('Ki', 'pKi') else 'I'].append(nm)

    kept = {c: v for c, v in d.items() if v['K'] and v['I']}
    with open(f'{name}.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['compound', 'Ki', 'IC50', 'stdev'])
        for cid, v in sorted(kept.items()):
            ic = v['I']
            w.writerow([v['name'] or cid,
                        round(sum(v['K']) / len(v['K']), 3),
                        round(sum(ic) / len(ic), 3),
                        round(st.stdev(ic), 3) if len(ic) > 1 else 0.0])
    print(f'{len(d)} compounds, {dropped} records dropped for unusable units, '
          f'{len(kept)} compounds with both metrics -> {name}.csv', flush=True)

if __name__ == '__main__':
    main()
