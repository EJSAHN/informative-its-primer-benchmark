import argparse
from pathlib import Path
import pandas as pd
from .common import iter_species_fastas, rc


def context(seq, start, end, flank=18):
    if pd.isna(start) or pd.isna(end):
        return ''
    start, end = int(start), int(end)
    return seq[max(0,start-flank):min(len(seq),end+flank)]


def main():
    ap = argparse.ArgumentParser(description='Select deterministic representative records for filtering audits.')
    ap.add_argument('--fasta-root', required=True)
    ap.add_argument('--site-flags', required=True)
    ap.add_argument('--pair-catalog', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    flags = pd.read_csv(args.site_flags).sort_values(['species','accession','pair_name'])
    pairs = pd.read_csv(args.pair_catalog).set_index('pair_name')
    seqs = {}
    for species, acc, fp, header, seq in iter_species_fastas(Path(args.fasta_root)):
        seqs[(species,acc)] = seq
    flags['any_terminal'] = flags[['fwd_terminal','rev_terminal']].fillna(False).astype(bool).any(axis=1)
    flags['fwd_metadata_contributed'] = flags['fwd_found'].fillna(False).astype(bool) & flags['fwd_metadata_same_primer'].fillna(False).astype(bool)
    flags['rev_metadata_contributed'] = flags['rev_found'].fillna(False).astype(bool) & flags['rev_metadata_same_primer'].fillna(False).astype(bool)
    flags['any_metadata'] = flags[['fwd_metadata_contributed','rev_metadata_contributed']].any(axis=1)
    flags['any_ambiguity'] = flags[['fwd_ambiguity_followup','rev_ambiguity_followup']].fillna(False).astype(bool).any(axis=1)
    selectors = [
        ('complete_informative_pair', flags['informative_pair'].fillna(False).astype(bool)),
        ('truncated_or_missing_reverse_site', flags['pair_status'].eq('no_rev')),
        ('terminal_only', flags['any_terminal'] & ~flags['any_metadata']),
        ('metadata_only', flags['any_metadata'] & ~flags['any_terminal']),
        ('terminal_and_metadata', flags['any_terminal'] & flags['any_metadata']),
        ('ambiguous_window_followup', flags['any_ambiguity']),
    ]
    rows=[]
    for label, mask in selectors:
        sub = flags[mask]
        if sub.empty:
            continue
        r = sub.iloc[0]
        seq = seqs.get((r.species, r.accession), '')
        p = pairs.loc[r.pair_name]
        rows.append({
            'artefact_class': label, 'species': r.species, 'accession': r.accession,
            'pair_name': r.pair_name, 'record_length': r.seq_len,
            'fwd_name': r.fwd_name, 'fwd_primer': p.fwd,
            'fwd_status': r.fwd_site_status, 'fwd_start': r.fwd_start, 'fwd_end': r.fwd_end,
            'fwd_window': r.fwd_window, 'fwd_context_plus_strand': context(seq,r.fwd_start,r.fwd_end),
            'rev_name': r.rev_name, 'rev_primer': p.rev,
            'rev_status': r.rev_site_status, 'rev_start': r.rev_start, 'rev_end': r.rev_end,
            'rev_window_primer_orientation': r.rev_window,
            'rev_context_plus_strand': context(seq,r.rev_start,r.rev_end),
            'record_start_context': seq[:60], 'record_end_context': seq[-60:],
            'fwd_metadata_match_fields': r.fwd_metadata_match_fields,
            'rev_metadata_match_fields': r.rev_metadata_match_fields,
        })
    outdir=Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    df=pd.DataFrame(rows)
    df.to_csv(outdir/'record_examples.csv', index=False)
    with open(outdir/'record_examples.txt','w',encoding='utf-8') as handle:
        for r in rows:
            handle.write(f"[{r['artefact_class']}] {r['species']} {r['accession']} {r['pair_name']}\n")
            handle.write(f"  Forward {r['fwd_name']} ({r['fwd_status']}): {r['fwd_context_plus_strand']}\n")
            handle.write(f"  Reverse {r['rev_name']} ({r['rev_status']}): {r['rev_context_plus_strand']}\n")
            handle.write(f"  Record start: {r['record_start_context']}\n")
            handle.write(f"  Record end:   {r['record_end_context']}\n\n")
    print(f'Wrote record examples to {outdir}')

if __name__ == '__main__':
    main()
