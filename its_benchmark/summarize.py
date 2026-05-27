import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def summarize(calls_path: Path, outdir: Path):
    calls = pd.read_csv(calls_path)
    outdir.mkdir(parents=True, exist_ok=True)
    clean = (calls.groupby(['species','pair_name','accession'], as_index=False)
             .agg(N_raw_rows=('accession','size'), informative_pair=('informative_pair','max'),
                  eligible=('eligible','max'), strict_hit=('strict_hit','max'),
                  fail_MM=('fail_MM','max'), fail_3p=('fail_3p','max'),
                  product_len_bp=('product_len_bp','median'),
                  elig_reason=('elig_reason', lambda x: ';'.join(sorted(set(str(v) for v in x if str(v) and str(v)!='nan')))),
                  pair_status=('pair_status', lambda x: ';'.join(sorted(set(str(v) for v in x if str(v) and str(v)!='nan'))))))
    clean.to_csv(outdir / 'calls_deduplicated.csv', index=False)
    summary = (clean.groupby(['species','pair_name'], as_index=False)
               .agg(N=('accession','nunique'), informative=('informative_pair','sum'),
                    eligible=('eligible','sum'), covered=('strict_hit','sum'),
                    fail_MM=('fail_MM','sum'), fail_3p=('fail_3p','sum')))
    summary['eligible_not_strict'] = summary['eligible'] - summary['covered']
    summary['noninformative_or_ineligible'] = summary['N'] - summary['eligible']
    summary['coverage_pct_eligible'] = np.where(summary['eligible']>0, 100*summary['covered']/summary['eligible'], np.nan)
    summary['coverage_pct_all'] = 100*summary['covered']/summary['N']
    summary['delta_coverage_pp'] = summary['coverage_pct_eligible'] - summary['coverage_pct_all']
    summary.to_csv(outdir / 'Table_all_pairs_by_species_informative.csv', index=False)
    bad = clean[clean['eligible'] == 0].copy()
    if not bad.empty:
        reason = bad.groupby(['species','pair_name','elig_reason'], as_index=False).size().rename(columns={'size':'count'})
    else:
        reason = pd.DataFrame(columns=['species','pair_name','elig_reason','count'])
    reason.to_csv(outdir / 'Ineligibility_breakdown_informative.csv', index=False)
    dup = calls.groupby(['species','pair_name','accession'], as_index=False).size().rename(columns={'size':'rows_per_accession'})
    dup_summary = (dup.groupby(['species','pair_name'], as_index=False)
                   .agg(n_accessions=('accession','nunique'), total_rows=('rows_per_accession','sum'), max_rows_per_accession=('rows_per_accession','max')))
    dup_summary['has_duplicates'] = dup_summary['total_rows'] > dup_summary['n_accessions']
    dup_summary['duplication_ratio'] = dup_summary['total_rows'] / dup_summary['n_accessions']
    dup_summary.to_csv(outdir / 'duplicate_report_by_pair.csv', index=False)
    return summary

def main():
    ap = argparse.ArgumentParser(description='Summarize informative-site primer-pair calls.')
    ap.add_argument('--calls', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    summarize(Path(args.calls), Path(args.outdir))
    print('Wrote summary tables to', args.outdir)

if __name__ == '__main__':
    main()
