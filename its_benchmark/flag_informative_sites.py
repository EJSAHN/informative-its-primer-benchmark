import argparse
from pathlib import Path
import pandas as pd
from .common import iter_species_fastas, scan_pair, is_terminal_hit, sequence_hash
from .metadata import load_genbank_metadata, match_primer_metadata


def main():
    ap = argparse.ArgumentParser(description='Flag terminal, metadata-derived, ambiguous, and non-informative primer sites.')
    ap.add_argument('--fasta-root', required=True)
    ap.add_argument('--pair-catalog', required=True)
    ap.add_argument('--genbank-dir', default='', help='Optional folder containing GenBank flatfiles for metadata flagging.')
    ap.add_argument('--out', required=True)
    ap.add_argument('--max-mm', type=int, default=4)
    ap.add_argument('--terminal-mode', choices=['fixed', 'primer_length'], default='primer_length')
    ap.add_argument('--terminal-bp', type=int, default=30, help='Used only with --terminal-mode fixed.')
    ap.add_argument('--terminal-extra-bp', type=int, default=0)
    ap.add_argument('--ambiguity-followup-threshold', type=float, default=0.05,
                    help='Flag, but do not automatically exclude, binding windows above this ambiguous-base fraction.')
    args = ap.parse_args()

    pairs = pd.read_csv(args.pair_catalog)
    required = {'pair_name', 'fwd_name', 'rev_name', 'fwd', 'rev', 'tmin', 'tmax'}
    missing = required - set(pairs.columns)
    if missing:
        raise SystemExit(f'Pair catalog missing columns: {missing}')

    def terminal_window(primer_seq: str) -> int:
        if args.terminal_mode == 'fixed':
            return int(args.terminal_bp)
        return int(len(str(primer_seq).replace(' ', '')) + args.terminal_extra_bp)

    metadata = load_genbank_metadata(args.genbank_dir) if args.genbank_dir else {}

    rows = []
    for species, acc, fp, header, seq in iter_species_fastas(Path(args.fasta_root)):
        meta_fields = metadata.get(acc, metadata.get(acc.split('.')[0], {}))
        seq_len = len(seq)
        seq_hash = sequence_hash(seq)
        for _, p in pairs.iterrows():
            f, r = scan_pair(seq, p['fwd'], p['rev'], max_mm=args.max_mm)
            f_found = f is not None
            r_found = r is not None
            f_win = terminal_window(p['fwd'])
            r_win = terminal_window(p['rev'])
            f_terminal = bool(f_found and is_terminal_hit(f['start_plus'], f['end_plus'], seq_len, f_win))
            r_terminal = bool(r_found and is_terminal_hit(r['start_plus'], r['end_plus'], seq_len, r_win))
            f_meta_rec = match_primer_metadata(meta_fields, p['fwd_name'], p['fwd'])
            r_meta_rec = match_primer_metadata(meta_fields, p['rev_name'], p['rev'])
            f_meta = bool(f_meta_rec['matched'])
            r_meta = bool(r_meta_rec['matched'])
            f_info = bool(f_found and not f_terminal and not f_meta)
            r_info = bool(r_found and not r_terminal and not r_meta)
            f_ambig = float(f['ambiguous_fraction']) if f_found else None
            r_ambig = float(r['ambiguous_fraction']) if r_found else None
            f_follow = bool(f_found and f_ambig > args.ambiguity_followup_threshold)
            r_follow = bool(r_found and r_ambig > args.ambiguity_followup_threshold)

            def status(found, terminal, meta, informative):
                if not found:
                    return 'no_site'
                if meta and terminal:
                    return 'metadata_and_terminal'
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
                'species': species,
                'accession': acc,
                'sequence_hash': seq_hash,
                'source_fasta': str(Path(fp).relative_to(Path(args.fasta_root))),
                'pair_name': p['pair_name'],
                'fwd_name': p['fwd_name'],
                'rev_name': p['rev_name'],
                'seq_len': seq_len,
                'terminal_mode': args.terminal_mode,
                'terminal_extra_bp': args.terminal_extra_bp,
                'metadata_match_rule': 'trusted_primer_fields',
                'ambiguity_followup_threshold': args.ambiguity_followup_threshold,
                'fwd_terminal_bp_used': f_win,
                'rev_terminal_bp_used': r_win,
                'fwd_found': f_found,
                'rev_found': r_found,
                'fwd_start': f['start_plus'] if f_found else None,
                'fwd_end': f['end_plus'] if f_found else None,
                'rev_start': r['start_plus'] if r_found else None,
                'rev_end': r['end_plus'] if r_found else None,
                'fwd_window': f['window'] if f_found else '',
                'rev_window': r['window'] if r_found else '',
                'fwd_mm': f['mm'] if f_found else None,
                'rev_mm': r['mm'] if r_found else None,
                'fwd_ambiguous_count': f['ambiguous_count'] if f_found else None,
                'rev_ambiguous_count': r['ambiguous_count'] if r_found else None,
                'fwd_ambiguous_fraction': f_ambig,
                'rev_ambiguous_fraction': r_ambig,
                'fwd_ambiguity_followup': f_follow,
                'rev_ambiguity_followup': r_follow,
                'fwd_terminal': f_terminal,
                'rev_terminal': r_terminal,
                'fwd_metadata_available': f_meta_rec['metadata_available'],
                'rev_metadata_available': r_meta_rec['metadata_available'],
                'fwd_metadata_same_primer': f_meta,
                'rev_metadata_same_primer': r_meta,
                'fwd_metadata_name_match': f_meta_rec['name_match'],
                'rev_metadata_name_match': r_meta_rec['name_match'],
                'fwd_metadata_sequence_match': f_meta_rec['sequence_match'],
                'rev_metadata_sequence_match': r_meta_rec['sequence_match'],
                'fwd_metadata_match_fields': f_meta_rec['match_fields'],
                'rev_metadata_match_fields': r_meta_rec['match_fields'],
                'fwd_metadata_match_terms': f_meta_rec['match_terms'],
                'rev_metadata_match_terms': r_meta_rec['match_terms'],
                'fwd_metadata_context_name_match': f_meta_rec['context_name_match'],
                'rev_metadata_context_name_match': r_meta_rec['context_name_match'],
                'fwd_metadata_context_sequence_match': f_meta_rec['context_sequence_match'],
                'rev_metadata_context_sequence_match': r_meta_rec['context_sequence_match'],
                'fwd_metadata_context_match_fields': f_meta_rec['context_match_fields'],
                'rev_metadata_context_match_fields': r_meta_rec['context_match_fields'],
                'fwd_metadata_context_match_terms': f_meta_rec['context_match_terms'],
                'rev_metadata_context_match_terms': r_meta_rec['context_match_terms'],
                'fwd_site_status': status(f_found, f_terminal, f_meta, f_info),
                'rev_site_status': status(r_found, r_terminal, r_meta, r_info),
                'informative_fwd': f_info,
                'informative_rev': r_info,
                'informative_pair': bool(f_info and r_info),
                'pair_status': pair_status,
            })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f'Wrote {len(rows)} rows to {out}')


if __name__ == '__main__':
    main()
