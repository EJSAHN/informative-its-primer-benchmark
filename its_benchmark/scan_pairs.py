import argparse
from pathlib import Path
import pandas as pd
from .common import iter_species_fastas, scan_pair

def main():
    ap = argparse.ArgumentParser(description='Score primer pairs using informative-site flags and strict mismatch rules.')
    ap.add_argument('--fasta-root', required=True)
    ap.add_argument('--pair-catalog', required=True)
    ap.add_argument('--site-flags', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--max-mm-per-primer', type=int, default=4)
    ap.add_argument('--strict-total-mm', type=int, default=2)
    ap.add_argument('--allow-3p-mismatch', action='store_true')
    args = ap.parse_args()

    pairs = pd.read_csv(args.pair_catalog)
    flags = pd.read_csv(args.site_flags)
    flag_map = {(r.species, r.accession, r.pair_name): r for r in flags.itertuples(index=False)}
    rows = []
    for species, acc, fp, header, seq in iter_species_fastas(Path(args.fasta_root)):
        for _, p in pairs.iterrows():
            fl = flag_map.get((species, acc, p['pair_name']))
            f, r = scan_pair(seq, p['fwd'], p['rev'], max_mm=args.max_mm_per_primer)
            base = {
                'species': species, 'accession': acc, 'pair_name': p['pair_name'],
                'tag': p.get('tag', ''), 'fwd_name': p['fwd_name'], 'rev_name': p['rev_name'],
                'tmin': int(p['tmin']), 'tmax': int(p['tmax']),
                'informative_pair': False if fl is None else bool(fl.informative_pair),
                'pair_status': 'missing_flag' if fl is None else fl.pair_status,
                'fwd_site_status': '' if fl is None else fl.fwd_site_status,
                'rev_site_status': '' if fl is None else fl.rev_site_status,
            }
            out = dict(base)
            out.update({'eligible':0, 'strict_hit':0, 'product_len_bp':None, 'total_mm':None, 'any_3p':0,
                        'fail_MM':0, 'fail_3p':0, 'fail_size':0, 'elig_reason':''})
            if fl is None:
                out['elig_reason'] = 'missing_flag'; rows.append(out); continue
            if not bool(fl.informative_pair):
                out['elig_reason'] = fl.pair_status; rows.append(out); continue
            if not f and not r:
                out['elig_reason'] = 'no_fwd_no_rev'; rows.append(out); continue
            if not f:
                out['elig_reason'] = 'no_fwd'; rows.append(out); continue
            if not r:
                out['elig_reason'] = 'no_rev'; rows.append(out); continue
            product_len = r['end_plus'] - f['start_plus']
            if product_len <= 0:
                out['elig_reason'] = 'wrong_orientation'; rows.append(out); continue
            out['product_len_bp'] = product_len
            if not (int(p['tmin']) <= product_len <= int(p['tmax'])):
                out['elig_reason'] = 'outside_window'; out['fail_size'] = 1; rows.append(out); continue
            out['eligible'] = 1
            total_mm = f['mm'] + r['mm']
            any_3p = bool(f['any_3p'] or r['any_3p'])
            out['total_mm'] = total_mm
            out['any_3p'] = int(any_3p)
            strict = (total_mm <= args.strict_total_mm) and (args.allow_3p_mismatch or not any_3p)
            out['strict_hit'] = int(strict)
            if not strict:
                out['fail_MM'] = int(total_mm > args.strict_total_mm)
                out['fail_3p'] = int(any_3p and not args.allow_3p_mismatch)
            rows.append(out)
    out_path = Path(args.out); out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f'Wrote {len(rows)} rows to {out_path}')

if __name__ == '__main__':
    main()
