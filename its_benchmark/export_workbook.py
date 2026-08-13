import argparse
from pathlib import Path
import pandas as pd


def add_csv(writer, sheet_name, path):
    path = Path(path)
    if not path.exists():
        return False
    df = pd.read_csv(path)
    target = sheet_name[:31]
    df.to_excel(writer, sheet_name=target, index=False)
    ws = writer.sheets[target]
    for i, col in enumerate(df.columns):
        sample = [len(str(col))] + [len(str(x)) for x in df[col].head(300).fillna('').astype(str)]
        ws.set_column(i, i, min(max(10, min(max(sample), 35) + 2), 36))
    ws.freeze_panes(1, 0)
    ws.autofilter(0, 0, len(df), max(0, len(df.columns) - 1))
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results-root', required=True)
    parser.add_argument('--pair-catalog', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    results = Path(args.results_root)
    files = [
        ('Main_summary', results / 'main/species_pair_summary.csv'),
        ('Denominators', results / 'main/denominator_summary.csv'),
        ('Pair_calls', results / 'main/pair_calls.csv'),
        ('Site_flags', results / 'main/site_flags.csv'),
        ('Deduplicated_calls', results / 'main/calls_deduplicated.csv'),
        ('Flag_source_summary', results / 'audit/flag_source_summary.csv'),
        ('Flag_source_wide', results / 'audit/flag_source_summary_wide.csv'),
        ('Metadata_fields', results / 'audit/metadata_field_summary.csv'),
        ('Metadata_context', results / 'audit/metadata_context_field_summary.csv'),
        ('Metadata_evidence', results / 'audit/metadata_match_evidence.csv'),
        ('Metadata_availability', results / 'audit/metadata_availability_summary.csv'),
        ('Ambiguity_followup', results / 'audit/ambiguity_followup_summary.csv'),
        ('Dedup_sensitivity', results / 'audit/deduplication_sensitivity.csv'),
        ('Mismatch_sensitivity', results / 'audit/mismatch_sensitivity.csv'),
        ('Terminal_scenarios', results / 'audit/terminal_scenario_totals.csv'),
        ('Artefact_examples', results / 'audit/record_examples.csv'),
        ('Rarefaction_summary', results / 'audit/rarefaction_summary.csv'),
        ('Pair_catalog', Path(args.pair_catalog)),
    ]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out, engine='xlsxwriter') as writer:
        readme = pd.DataFrame({
            'Item': [
                'Description',
                'Main analysis',
                'Metadata exclusion rule',
                'Flag-source counting rule',
                'Quality note',
            ],
            'Value': [
                'Source tables and row-level outputs for informative-site ITS primer benchmarking.',
                'Primer-length terminal filtering with structured GenBank metadata flagging.',
                'Only trusted primer-specific metadata fields can exclude a detected site; context-only matches are retained for audit.',
                'Metadata contributes to a terminal/metadata category only when the corresponding primer site was detected in that record.',
                'Public FASTA records lack per-base quality scores. Binding windows with >5% ambiguous bases are flagged for follow-up but are not automatically excluded.',
            ],
        })
        readme.to_excel(writer, sheet_name='README', index=False)
        ws = writer.sheets['README']
        ws.set_column(0, 0, 28)
        ws.set_column(1, 1, 115)

        guide = pd.DataFrame([
            {'Sheet': sheet, 'Contents': Path(path).name}
            for sheet, path in files
            if Path(path).exists()
        ])
        guide.to_excel(writer, sheet_name='Sheet_guide', index=False)
        writer.sheets['Sheet_guide'].set_column(0, 0, 28)
        writer.sheets['Sheet_guide'].set_column(1, 1, 58)

        for sheet, path in files:
            add_csv(writer, sheet, path)

    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
