import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from .common import iter_species_fastas, scan_pair, is_terminal_hit, primer_appears_in_metadata

def main():
    ap = argparse.ArgumentParser(description='Flag terminal, primer-derived, and non-informative primer sites.')
    ap.add_argument('--fasta-root', required=True)
    ap.add_argument('--pair-catalog', required=True)
    ap.add_argument('--genbank-dir', default='', help='Optional folder containing GenBank flatfiles for metadata flagging.')
    ap.add_argument('--out', required=True)
    ap.add_argument('--max-mm', type=int, default=4)
    ap.add_argument('--terminal-mode', choices=['fixed', 'primer_length'], default='primer_length')
    ap.add_argument('--terminal-bp', type=int, default=30, help='Used only with --terminal-mode fixed.')
    ap.add_argument('--terminal-extra-bp', type=int, default=0)
    args = ap.parse_args()

    pairs = pd.read_csv(args.pair_catalog)
    required = {'pair_name','fwd_name','rev_name','fwd','rev','tmin','tmax'}
    missing = required - set(pairs.columns)
    if missing:
        raise SystemExit(f'Pair catalog missing columns: {missing}')

    def terminal_window(primer_seq: str) -> int:
        if args.terminal_mode == 'fixed':
            return int(args.terminal_bp)
        return int(len(str(primer_seq).replace(' ', '')) + args.terminal_extra_bp)

    metadata = {}
    gb_root = Path(args.genbank_dir) if args.genbank_dir else None
    if gb_root and gb_root.exists():
        for gb in gb_root.rglob('*.gb'):
            try:
                rec = next(SeqIO.parse(str(gb), 'genbank'))
                text = ' '.join([rec.description or '', str(rec.annotations), str(rec.features)])
                for key in {rec.id, rec.id.split('.')[0], getattr(rec, 'name', rec.id), getattr(rec, 'name', rec.id).split('.')[0]}:
                    metadata[key] = text
            except Exception as e:
                print(f'WARN: could not parse {gb}: {e}')

    rows = []
    for species, acc, fp, header, seq in iter_species_fastas(Path(args.fasta_root)):
        meta_text = metadata.get(acc, metadata.get(acc.split('.')[0], ''))
        seq_len = len(seq)
        for _, p in pairs.iterrows():
            f, r = scan_pair(seq, p['fwd'], p['rev'], max_mm=args.max_mm)
            f_found = f is not None
            r_found = r is not None
            f_win = terminal_window(p['fwd'])
            r_win = terminal_window(p['rev'])
            f_terminal = bool(f_found and is_terminal_hit(f['start_plus'], f['end_plus'], seq_len, f_win))
            r_terminal = bool(r_found and is_terminal_hit(r['start_plus'], r['end_plus'], seq_len, r_win))
            f_meta = primer_appears_in_metadata(meta_text, p['fwd_name'], p['fwd'])
            r_meta = primer_appears_in_metadata(meta_text, p['rev_name'], p['rev'])
            f_info = bool(f_found and not f_terminal and not f_meta)
            r_info = bool(r_found and not r_terminal and not r_meta)

            def status(found, terminal, meta, informative):
                if not found:
                    return 'no_site'
                if meta:
                    return 'metadata_same_primer'
                if terminal:
                    return 'terminal_primer_derived_or_uninformative'
                if informative:
                    return 'informative_internal'
                return 'noninformative_unknown'

            if f_info and r_info:
                pair_status = 'informative_pair'
            elif not f_found and not r_found:
                pair_status = 'no_fwd_no_rev'
            elif not f_found:
                pair_status = 'no_fwd'
            elif not r_found:
                pair_status = 'no_rev'
            elif not f_info and not r_info:
                pair_status = 'noninformative_fwd_rev'
            elif not f_info:
                pair_status = 'noninformative_fwd'
            else:
                pair_status = 'noninformative_rev'

            rows.append({
                'species': species, 'accession': acc, 'pair_name': p['pair_name'],
                'fwd_name': p['fwd_name'], 'rev_name': p['rev_name'], 'seq_len': seq_len,
                'terminal_mode': args.terminal_mode, 'terminal_extra_bp': args.terminal_extra_bp,
                'fwd_terminal_bp_used': f_win, 'rev_terminal_bp_used': r_win,
                'fwd_found': f_found, 'rev_found': r_found,
                'fwd_start': f['start_plus'] if f_found else None, 'fwd_end': f['end_plus'] if f_found else None,
                'rev_start': r['start_plus'] if r_found else None, 'rev_end': r['end_plus'] if r_found else None,
                'fwd_mm': f['mm'] if f_found else None, 'rev_mm': r['mm'] if r_found else None,
                'fwd_terminal': f_terminal, 'rev_terminal': r_terminal,
                'fwd_metadata_same_primer': f_meta, 'rev_metadata_same_primer': r_meta,
                'fwd_site_status': status(f_found, f_terminal, f_meta, f_info),
                'rev_site_status': status(r_found, r_terminal, r_meta, r_info),
                'informative_fwd': f_info, 'informative_rev': r_info,
                'informative_pair': bool(f_info and r_info), 'pair_status': pair_status,
            })
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f'Wrote {len(rows)} rows to {out}')

if __name__ == '__main__':
    main()
