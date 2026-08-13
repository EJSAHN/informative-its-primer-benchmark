import argparse, random
from pathlib import Path
import pandas as pd
import numpy as np
from .common import iter_species_fastas, scan_best, rc, is_terminal_hit

def main():
    ap = argparse.ArgumentParser(description='Extract informative internal binding-site windows and rarefy observed variants.')
    ap.add_argument('--fasta-root', required=True)
    ap.add_argument('--pair-catalog', required=True)
    ap.add_argument('--species', required=True, help='Example: "P. palmivora"')
    ap.add_argument('--primer-name', required=True, help='Example: ITS6')
    ap.add_argument('--direction', choices=['forward','reverse'], required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--site-flags', default='', help='Optional site_flags CSV; recommended for final analyses.')
    ap.add_argument('--max-mm', type=int, default=4)
    ap.add_argument('--terminal-mode', choices=['fixed', 'primer_length'], default='primer_length')
    ap.add_argument('--terminal-bp', type=int, default=30)
    ap.add_argument('--terminal-extra-bp', type=int, default=0)
    ap.add_argument('--reps', type=int, default=50)
    ap.add_argument('--seed', type=int, default=17)
    args = ap.parse_args()

    pairs = pd.read_csv(args.pair_catalog)
    seqs = []
    for _, r in pairs.iterrows():
        if r['fwd_name'] == args.primer_name:
            seqs.append(r['fwd'])
        if r['rev_name'] == args.primer_name:
            seqs.append(r['rev'])
    seqs = list(dict.fromkeys(seqs))
    if not seqs:
        raise SystemExit(f'Primer {args.primer_name} not found in pair catalog')
    primer_seq = seqs[0]

    def terminal_window() -> int:
        if args.terminal_mode == 'fixed':
            return int(args.terminal_bp)
        return int(len(primer_seq.replace(' ', '')) + args.terminal_extra_bp)

    seq_map = {}
    for species, acc, fp, header, seq in iter_species_fastas(Path(args.fasta_root)):
        if species == args.species:
            seq_map[acc] = seq
            seq_map[acc.split('.')[0]] = seq

    rows = []
    windows = []
    if args.site_flags:
        flags = pd.read_csv(args.site_flags)
        flags = flags[flags['species'].eq(args.species)].copy()
        if args.direction == 'forward':
            flags = flags[(flags['fwd_name'].eq(args.primer_name)) & (flags['fwd_site_status'].eq('informative_internal'))]
            coord_cols = ('fwd_start', 'fwd_end', 'fwd_mm')
        else:
            flags = flags[(flags['rev_name'].eq(args.primer_name)) & (flags['rev_site_status'].eq('informative_internal'))]
            coord_cols = ('rev_start', 'rev_end', 'rev_mm')
        flags = flags.drop_duplicates(subset=['species', 'accession', coord_cols[0], coord_cols[1]])
        for _, row in flags.iterrows():
            acc = str(row['accession'])
            seq = seq_map.get(acc, seq_map.get(acc.split('.')[0]))
            if seq is None:
                rows.append({'species': args.species, 'accession': acc, 'primer_name': args.primer_name, 'status': 'sequence_not_found'})
                continue
            start = int(row[coord_cols[0]])
            end = int(row[coord_cols[1]])
            mm = int(row[coord_cols[2]]) if pd.notna(row[coord_cols[2]]) else None
            window = seq[start:end] if args.direction == 'forward' else rc(seq[start:end])
            rows.append({'species': args.species, 'accession': acc, 'primer_name': args.primer_name,
                         'status': 'informative_internal', 'start': start, 'end': end, 'mm': mm, 'window': window})
            windows.append(window)
    else:
        tbp = terminal_window()
        for species, acc, fp, header, seq in iter_species_fastas(Path(args.fasta_root)):
            if species != args.species:
                continue
            hit = None
            if args.direction == 'forward':
                hit = scan_best(seq, primer_seq, max_mm=args.max_mm)
                if hit:
                    start, end = hit['pos'], hit['pos'] + len(primer_seq)
                    window = hit['window']
            else:
                hit = scan_best(rc(seq), primer_seq, max_mm=args.max_mm)
                if hit:
                    start = len(seq) - (hit['pos'] + len(primer_seq))
                    end = len(seq) - hit['pos']
                    window = hit['window']
            if not hit:
                rows.append({'species': species, 'accession': acc, 'primer_name': args.primer_name, 'status':'no_site'})
                continue
            terminal = is_terminal_hit(start, end, len(seq), tbp)
            status = 'terminal_noninformative' if terminal else 'informative_internal'
            rows.append({'species': species, 'accession': acc, 'primer_name': args.primer_name,
                         'status': status, 'start': start, 'end': end, 'mm': hit['mm'], 'window': window})
            if status == 'informative_internal':
                windows.append(window)

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    tag = args.species.replace(' ','').replace('.','')
    pd.DataFrame(rows).to_csv(outdir / f'{tag}_{args.primer_name}_site_windows.csv', index=False)

    if windows:
        random.seed(args.seed)
        n = len(windows)
        arr = np.zeros((args.reps, n))
        for r in range(args.reps):
            order = windows[:]
            random.shuffle(order)
            seen = set()
            for i, w in enumerate(order):
                seen.add(w)
                arr[r, i] = len(seen)
        mean = arr.mean(axis=0)
        total = len(set(windows))
        n95 = next((i+1 for i, v in enumerate(mean) if v >= 0.95 * total), n)
        pd.DataFrame({'n_sequences': range(1, n+1), 'mean_unique_observed_variants': mean}).to_csv(outdir / f'{tag}_{args.primer_name}_rarefaction.csv', index=False)
        pd.DataFrame([{'target': f'{tag}_{args.primer_name}', 'informative_windows': n, 'unique_observed_variants': total, 'n95': n95}]).to_csv(outdir / f'{tag}_{args.primer_name}_rarefaction_summary.csv', index=False)
        print(f'Informative windows={n}; unique observed variants={total}; n95={n95}')
    else:
        print('No informative internal windows found')

if __name__ == '__main__':
    main()
