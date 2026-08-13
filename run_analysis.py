from __future__ import annotations
import argparse
import json
import subprocess
import sys
from pathlib import Path
import pandas as pd


def run(command):
    command = [str(x) for x in command]
    print('\n>>> ' + ' '.join(command))
    subprocess.run(command, check=True)


def totals(summary_path: Path, scenario: str):
    summary = pd.read_csv(summary_path)
    return {
        'scenario': scenario,
        'accession_pair_rows': int(summary['N'].sum()),
        'informative_pair_rows': int(summary['informative'].sum()),
        'eligible_rows': int(summary['eligible'].sum()),
        'covered_rows': int(summary['covered'].sum()),
    }


def main():
    parser = argparse.ArgumentParser(description='Run informative-site ITS primer benchmarking and numerical audits.')
    parser.add_argument('--fasta-root', required=True)
    parser.add_argument('--pair-catalog', default='data/pair_catalog_auto.csv')
    parser.add_argument('--genbank-dir', default='')
    parser.add_argument('--outdir', default='results')
    parser.add_argument('--site-variant-targets', default='data/site_variant_targets.csv')
    parser.add_argument('--skip-tests', action='store_true')
    parser.add_argument('--skip-sensitivity', action='store_true')
    parser.add_argument('--skip-site-variants', action='store_true')
    parser.add_argument('--validate-study', action='store_true', help='Validate sequences and summary counts against bundled study manifests.')
    args = parser.parse_args()

    root = Path.cwd()
    python = sys.executable
    fasta_root = Path(args.fasta_root).resolve()
    pair_catalog = Path(args.pair_catalog).resolve()
    genbank_dir = Path(args.genbank_dir).resolve() if args.genbank_dir else None
    results = Path(args.outdir).resolve()
    main_dir = results / 'main'
    audit_dir = results / 'audit'
    main_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_tests:
        run([python, root / 'tests' / 'run_minimal_test.py', '--project-root', root])

    if args.validate_study:
        run([python, '-m', 'its_benchmark.validate_sequences', '--fasta-root', fasta_root,
             '--manifest', root / 'data' / 'accession_sequence_manifest.csv',
             '--out', audit_dir / 'sequence_manifest_validation.csv'])

    accessions = results / 'accessions.tsv'
    flags = main_dir / 'site_flags.csv'
    calls = main_dir / 'pair_calls.csv'
    run([python, '-m', 'its_benchmark.make_accession_list', '--fasta-root', fasta_root, '--out', accessions])

    flag_command = [python, '-m', 'its_benchmark.flag_informative_sites',
                    '--fasta-root', fasta_root, '--pair-catalog', pair_catalog,
                    '--terminal-mode', 'primer_length', '--out', flags]
    if genbank_dir:
        flag_command.extend(['--genbank-dir', genbank_dir])
    run(flag_command)
    run([python, '-m', 'its_benchmark.scan_pairs', '--fasta-root', fasta_root,
         '--pair-catalog', pair_catalog, '--site-flags', flags, '--out', calls])
    run([python, '-m', 'its_benchmark.summarize', '--calls', calls, '--outdir', main_dir])
    run([python, '-m', 'its_benchmark.audit_flags', '--site-flags', flags, '--outdir', audit_dir])
    run([python, '-m', 'its_benchmark.dedup_sensitivity', '--calls', calls,
         '--out', audit_dir / 'deduplication_sensitivity.csv'])
    run([python, '-m', 'its_benchmark.sensitivity', '--calls', calls,
         '--out', audit_dir / 'mismatch_sensitivity.csv'])
    run([python, '-m', 'its_benchmark.record_examples', '--fasta-root', fasta_root,
         '--site-flags', flags, '--pair-catalog', pair_catalog, '--outdir', audit_dir])

    scenario_rows = [totals(main_dir / 'species_pair_summary.csv', 'primer_length')]
    if not args.skip_sensitivity:
        for width in [20, 30, 50]:
            scenario_dir = results / f'terminal_{width}bp'
            scenario_dir.mkdir(parents=True, exist_ok=True)
            scenario_flags = scenario_dir / 'site_flags.csv'
            scenario_calls = scenario_dir / 'pair_calls.csv'
            command = [python, '-m', 'its_benchmark.flag_informative_sites',
                       '--fasta-root', fasta_root, '--pair-catalog', pair_catalog,
                       '--terminal-mode', 'fixed', '--terminal-bp', width, '--out', scenario_flags]
            if genbank_dir:
                command.extend(['--genbank-dir', genbank_dir])
            run(command)
            run([python, '-m', 'its_benchmark.scan_pairs', '--fasta-root', fasta_root,
                 '--pair-catalog', pair_catalog, '--site-flags', scenario_flags, '--out', scenario_calls])
            run([python, '-m', 'its_benchmark.summarize', '--calls', scenario_calls, '--outdir', scenario_dir])
            scenario_rows.append(totals(scenario_dir / 'species_pair_summary.csv', f'terminal_{width}bp'))
    pd.DataFrame(scenario_rows).to_csv(audit_dir / 'terminal_scenario_totals.csv', index=False)

    rarefaction_rows = []
    if not args.skip_site_variants:
        targets = pd.read_csv(args.site_variant_targets)
        rarefaction_dir = audit_dir / 'site_variants'
        rarefaction_dir.mkdir(parents=True, exist_ok=True)
        for target in targets.itertuples(index=False):
            run([python, '-m', 'its_benchmark.site_variants', '--fasta-root', fasta_root,
                 '--pair-catalog', pair_catalog, '--site-flags', flags,
                 '--species', target.species, '--primer-name', target.primer_name,
                 '--direction', target.direction, '--outdir', rarefaction_dir])
        for path in rarefaction_dir.glob('*_rarefaction_summary.csv'):
            rarefaction_rows.append(pd.read_csv(path))
        if rarefaction_rows:
            pd.concat(rarefaction_rows, ignore_index=True).to_csv(audit_dir / 'rarefaction_summary.csv', index=False)

    workbook = results / 'benchmark_results.xlsx'
    run([python, '-m', 'its_benchmark.export_workbook', '--results-root', results,
         '--pair-catalog', pair_catalog, '--out', workbook])

    settings = {
        'terminal_filter': 'primer_length',
        'metadata_exclusion': 'trusted_primer_fields',
        'site_search_max_mismatches_per_primer': 4,
        'strict_total_mismatches': 2,
        'strict_terminal_mismatch_allowed': False,
        'ambiguity_followup_threshold': 0.05,
        'rarefaction_replicates': 50,
        'rarefaction_seed': 17,
    }
    (results / 'analysis_settings.json').write_text(json.dumps(settings, indent=2), encoding='utf-8')

    if args.validate_study:
        run([python, '-m', 'its_benchmark.validate_results',
             '--summary', main_dir / 'species_pair_summary.csv',
             '--expected', root / 'data' / 'study_expected_metrics.json',
             '--out', audit_dir / 'study_result_validation.csv'])

    run([python, '-m', 'its_benchmark.record_environment', '--results-dir', results])
    print(f'\nAnalysis complete: {results}')


if __name__ == '__main__':
    main()
