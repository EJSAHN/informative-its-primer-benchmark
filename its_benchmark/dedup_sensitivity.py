import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def summarize_unit(calls: pd.DataFrame, unit: str, key_col: str) -> pd.DataFrame:
    clean = (calls.groupby(['species','pair_name',key_col], as_index=False)
             .agg(informative=('informative_pair','max'), eligible=('eligible','max'),
                  strict_hit=('strict_hit','max')))
    out = (clean.groupby(['species','pair_name'], as_index=False)
           .agg(N=(key_col,'nunique'), informative=('informative','sum'),
                eligible=('eligible','sum'), covered=('strict_hit','sum')))
    out['coverage_pct_eligible'] = np.where(out.eligible>0, 100*out.covered/out.eligible, np.nan)
    out['deduplication_unit'] = unit
    return out


def main():
    ap = argparse.ArgumentParser(description='Compare accession-level and exact-sequence-level deduplication.')
    ap.add_argument('--calls', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    calls = pd.read_csv(args.calls)
    if 'sequence_hash' not in calls.columns:
        raise SystemExit('calls file lacks sequence_hash; rerun scan_pairs with the current code.')
    a = summarize_unit(calls, 'accession', 'accession')
    s = summarize_unit(calls, 'exact_full_sequence', 'sequence_hash')
    out = pd.concat([a,s], ignore_index=True)
    path = Path(args.out); path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    print(f'Wrote {path}')

if __name__ == '__main__':
    main()
