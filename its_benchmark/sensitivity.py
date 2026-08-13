import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def main():
    ap = argparse.ArgumentParser(description='Summarize assay rankings under altered strict-call thresholds.')
    ap.add_argument('--calls', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    calls = pd.read_csv(args.calls)
    scenarios = [('default_mm2_no3p', 2, False), ('mm3_no3p', 3, False), ('mm2_allow3p', 2, True)]
    rows = []
    for name, mm_thr, allow3p in scenarios:
        df = calls.copy()
        eligible = df['eligible'].fillna(0).astype(int) == 1
        total_mm = pd.to_numeric(df['total_mm'], errors='coerce')
        any3p = df['any_3p'].fillna(0).astype(int) == 1
        df['strict_scenario'] = (eligible & (total_mm <= mm_thr) & (allow3p | ~any3p)).astype(int)
        clean = df.groupby(['species','pair_name','accession'], as_index=False).agg(N_raw_rows=('accession','size'), eligible=('eligible','max'), strict_hit=('strict_scenario','max'))
        summ = clean.groupby(['species','pair_name'], as_index=False).agg(N=('accession','nunique'), eligible=('eligible','sum'), covered=('strict_hit','sum'))
        summ['coverage_pct_eligible'] = np.where(summ['eligible']>0, 100*summ['covered']/summ['eligible'], np.nan)
        summ['scenario'] = name
        rows.append(summ)
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(out, index=False)
    print(f'Wrote {out}')

if __name__ == '__main__':
    main()
