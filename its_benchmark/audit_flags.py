import argparse
from pathlib import Path
import pandas as pd

CATEGORIES = [
    'terminal_only',
    'metadata_only',
    'both_terminal_and_metadata',
    'neither',
]


def _as_bool(series: pd.Series) -> pd.Series:
    """Coerce common CSV boolean representations safely."""
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return (
        series.fillna(False)
        .astype(str)
        .str.strip()
        .str.lower()
        .map({'true': True, 'false': False, '1': True, '0': False})
        .fillna(False)
        .astype(bool)
    )


def classify(terminal: bool, metadata: bool) -> str:
    if terminal and metadata:
        return 'both_terminal_and_metadata'
    if terminal:
        return 'terminal_only'
    if metadata:
        return 'metadata_only'
    return 'neither'


def aggregate_categories(df: pd.DataFrame, unit_cols, species_col='species') -> pd.DataFrame:
    grouped = (
        df.groupby(unit_cols, as_index=False, dropna=False)
        .agg(
            any_terminal=('any_terminal', 'max'),
            any_metadata=('any_metadata', 'max'),
        )
    )
    grouped['flag_source_category'] = [
        classify(t, m) for t, m in zip(grouped['any_terminal'], grouped['any_metadata'])
    ]

    rows = []
    groups = [('All taxa', grouped)] + [(sp, sub) for sp, sub in grouped.groupby(species_col, dropna=False)]
    for species_label, sub in groups:
        counts = sub['flag_source_category'].value_counts()
        for category in CATEGORIES:
            rows.append({
                'species': species_label,
                'flag_source_category': category,
                'count': int(counts.get(category, 0)),
                'total_units': int(len(sub)),
            })
    return pd.DataFrame(rows)


def _split_fields(value):
    if pd.isna(value) or str(value).strip() in ('', 'None'):
        return []
    return [x for x in str(value).split(';') if x]


def _add_all_taxa_row(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    numeric = [c for c in df.columns if c != 'species']
    total = {'species': 'All taxa'}
    for col in numeric:
        total[col] = int(df[col].sum())
    return pd.concat([df, pd.DataFrame([total])], ignore_index=True)


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Summarize the actual contribution of terminal filtering and trusted '
            'GenBank primer metadata to non-informative-site classification.'
        )
    )
    parser.add_argument('--site-flags', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.site_flags)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    required = [
        'species', 'accession', 'pair_name',
        'fwd_name', 'rev_name',
        'fwd_found', 'rev_found',
        'fwd_terminal', 'rev_terminal',
        'fwd_metadata_same_primer', 'rev_metadata_same_primer',
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f'Site-flags table is missing required columns: {missing}')

    for col in [
        'fwd_found', 'rev_found',
        'fwd_terminal', 'rev_terminal',
        'fwd_metadata_same_primer', 'rev_metadata_same_primer',
    ]:
        df[col] = _as_bool(df[col])

    # A metadata match contributes to classification only when that primer site
    # was actually detected in the record. Metadata associated with an absent
    # site remains useful audit evidence but cannot itself exclude a detected
    # primer-template site.
    df['fwd_metadata_contributed'] = df['fwd_found'] & df['fwd_metadata_same_primer']
    df['rev_metadata_contributed'] = df['rev_found'] & df['rev_metadata_same_primer']
    df['any_terminal'] = df['fwd_terminal'] | df['rev_terminal']
    df['any_metadata'] = df['fwd_metadata_contributed'] | df['rev_metadata_contributed']

    pair_summary = aggregate_categories(df, ['species', 'accession', 'pair_name'])
    pair_summary.insert(0, 'counting_unit', 'accession_pair')

    accession_summary = aggregate_categories(df, ['species', 'accession'])
    accession_summary.insert(0, 'counting_unit', 'accession')

    site_frames = []
    for side in ['fwd', 'rev']:
        site = pd.DataFrame({
            'species': df['species'],
            'accession': df['accession'],
            'primer_name': df[f'{side}_name'],
            'direction': side,
            'start': df.get(f'{side}_start'),
            'end': df.get(f'{side}_end'),
            'found': df[f'{side}_found'],
            'any_terminal': df[f'{side}_terminal'],
            'any_metadata': df[f'{side}_metadata_contributed'],
        })
        # Primer-site counting is meaningful only for detected sites.
        site_frames.append(site[site['found']].copy())

    sites = pd.concat(site_frames, ignore_index=True)
    sites = sites.drop_duplicates(
        ['species', 'accession', 'primer_name', 'direction', 'start', 'end']
    )
    site_summary = aggregate_categories(
        sites,
        ['species', 'accession', 'primer_name', 'direction', 'start', 'end'],
    )
    site_summary.insert(0, 'counting_unit', 'accession_primer_site')

    flag_summary = pd.concat(
        [accession_summary, site_summary, pair_summary],
        ignore_index=True,
    )
    flag_summary.to_csv(outdir / 'flag_source_summary.csv', index=False)

    # Wide-format summary for downstream reporting.
    wide = (
        flag_summary.pivot_table(
            index=['counting_unit', 'species', 'total_units'],
            columns='flag_source_category',
            values='count',
            aggfunc='first',
            fill_value=0,
        )
        .reset_index()
    )
    for category in CATEGORIES:
        if category not in wide.columns:
            wide[category] = 0
    wide = wide[[
        'counting_unit', 'species', 'total_units',
        'terminal_only', 'metadata_only',
        'both_terminal_and_metadata', 'neither',
    ]]
    wide.to_csv(outdir / 'flag_source_summary_wide.csv', index=False)

    # Retain both trusted exclusion evidence and context-only evidence.
    trusted_rows = []
    context_rows = []
    for side in ['fwd', 'rev']:
        trusted_field_col = f'{side}_metadata_match_fields'
        trusted_term_col = f'{side}_metadata_match_terms'
        context_field_col = f'{side}_metadata_context_match_fields'
        context_term_col = f'{side}_metadata_context_match_terms'

        if trusted_field_col in df.columns:
            trusted_cols = [
                'species', 'accession', f'{side}_name',
                trusted_field_col, trusted_term_col,
            ]
            for rec in df[trusted_cols].drop_duplicates().itertuples(index=False):
                species, accession, primer, fields, terms = rec
                for field in _split_fields(fields):
                    trusted_rows.append({
                        'species': species,
                        'accession': accession,
                        'primer_name': primer,
                        'direction': side,
                        'metadata_field': field,
                        'match_terms': terms,
                        'evidence_class': 'trusted_primer_field',
                    })

        if context_field_col in df.columns:
            context_cols = [
                'species', 'accession', f'{side}_name',
                context_field_col, context_term_col,
            ]
            for rec in df[context_cols].drop_duplicates().itertuples(index=False):
                species, accession, primer, fields, terms = rec
                for field in _split_fields(fields):
                    context_rows.append({
                        'species': species,
                        'accession': accession,
                        'primer_name': primer,
                        'direction': side,
                        'metadata_field': field,
                        'match_terms': terms,
                        'evidence_class': 'context_only_not_excluded',
                    })

    evidence_columns = [
        'species', 'accession', 'primer_name', 'direction',
        'metadata_field', 'match_terms', 'evidence_class',
    ]
    trusted_df = pd.DataFrame(trusted_rows, columns=evidence_columns)
    context_df = pd.DataFrame(context_rows, columns=evidence_columns)
    evidence = pd.concat([trusted_df, context_df], ignore_index=True)
    evidence.to_csv(outdir / 'metadata_match_evidence.csv', index=False)

    def field_summary(source: pd.DataFrame, filename: str):
        if source.empty:
            summary = pd.DataFrame(
                columns=['species', 'primer_name', 'direction', 'metadata_field', 'count']
            )
        else:
            summary = (
                source.groupby(
                    ['species', 'primer_name', 'direction', 'metadata_field'],
                    as_index=False,
                )
                .size()
                .rename(columns={'size': 'count'})
            )
        summary.to_csv(outdir / filename, index=False)

    field_summary(trusted_df, 'metadata_field_summary.csv')
    field_summary(context_df, 'metadata_context_field_summary.csv')

    # Availability is reported at accession level, not repeated accession-pair rows.
    df['trusted_metadata_evidence'] = (
        df['fwd_metadata_same_primer'] | df['rev_metadata_same_primer']
    )
    availability = (
        df.groupby(['species', 'accession'], as_index=False)
        .agg(trusted_primer_metadata=('trusted_metadata_evidence', 'max'))
        .groupby('species', as_index=False)
        .agg(
            total_accessions=('accession', 'nunique'),
            accessions_with_trusted_primer_metadata=('trusted_primer_metadata', 'sum'),
        )
    )
    availability['accessions_without_trusted_primer_metadata'] = (
        availability['total_accessions']
        - availability['accessions_with_trusted_primer_metadata']
    )
    availability = _add_all_taxa_row(availability)
    availability.to_csv(outdir / 'metadata_availability_summary.csv', index=False)

    # Ambiguity follow-up audit (> threshold; does not alter informative status).
    ambiguity_rows = []
    for side in ['fwd', 'rev']:
        follow_col = f'{side}_ambiguity_followup'
        if follow_col not in df.columns:
            continue
        follow = _as_bool(df[follow_col])
        sub = df[follow]
        if not sub.empty:
            summary = (
                sub.groupby(['species', f'{side}_name'], as_index=False)
                .size()
                .rename(columns={f'{side}_name': 'primer_name', 'size': 'count'})
            )
            summary['direction'] = side
            ambiguity_rows.append(summary)

    ambiguity = (
        pd.concat(ambiguity_rows, ignore_index=True)
        if ambiguity_rows
        else pd.DataFrame(columns=['species', 'primer_name', 'count', 'direction'])
    )
    ambiguity.to_csv(outdir / 'ambiguity_followup_summary.csv', index=False)
    print(f'Wrote corrected audit summaries to {outdir}')


if __name__ == '__main__':
    main()
