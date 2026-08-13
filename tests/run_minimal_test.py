import argparse
import subprocess
import sys
import tempfile
from pathlib import Path
import pandas as pd

COMPARE = {
    'site_flags.csv': [
        'accession', 'pair_status', 'fwd_site_status', 'rev_site_status',
        'informative_pair', 'fwd_metadata_same_primer',
        'fwd_metadata_context_name_match', 'fwd_ambiguity_followup',
    ],
    'pair_calls.csv': [
        'accession', 'informative_pair', 'eligible', 'strict_hit', 'elig_reason',
    ],
    'species_pair_summary.csv': [
        'species', 'pair_name', 'N', 'informative', 'eligible', 'covered',
    ],
}


def run(command):
    subprocess.run([str(x) for x in command], check=True)


def compare(actual: Path, expected: Path, columns):
    left = pd.read_csv(actual)[columns].sort_values(columns[:2]).reset_index(drop=True)
    right = pd.read_csv(expected)[columns].sort_values(columns[:2]).reset_index(drop=True)
    pd.testing.assert_frame_equal(left, right, check_dtype=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project-root', required=True)
    args = parser.parse_args()
    root = Path(args.project_root).resolve()
    study_catalog = pd.read_csv(root / 'data' / 'pair_catalog_auto.csv')
    its2_kyo2 = set(study_catalog.loc[study_catalog['rev_name'].eq('ITS2_KYO2'), 'rev'].dropna())
    if its2_kyo2 != {'TTYRCTRCGTTCTTCATC'}:
        raise AssertionError(f'Unexpected ITS2_KYO2 sequence: {sorted(its2_kyo2)}')
    data = root / 'tests' / 'minimal_data'
    expected = root / 'tests' / 'expected'
    python = sys.executable

    with tempfile.TemporaryDirectory() as temporary:
        temporary = Path(temporary)
        flags = temporary / 'site_flags.csv'
        calls = temporary / 'pair_calls.csv'
        summary_dir = temporary / 'summary'
        run([python, '-m', 'its_benchmark.flag_informative_sites',
             '--fasta-root', data / 'raw_fasta', '--pair-catalog', data / 'pair_catalog.csv',
             '--genbank-dir', data / 'genbank', '--terminal-mode', 'primer_length', '--out', flags])
        run([python, '-m', 'its_benchmark.scan_pairs',
             '--fasta-root', data / 'raw_fasta', '--pair-catalog', data / 'pair_catalog.csv',
             '--site-flags', flags, '--out', calls])
        run([python, '-m', 'its_benchmark.summarize', '--calls', calls, '--outdir', summary_dir])

        compare(flags, expected / 'site_flags_expected.csv', COMPARE['site_flags.csv'])
        compare(calls, expected / 'pair_calls_expected.csv', COMPARE['pair_calls.csv'])
        compare(summary_dir / 'species_pair_summary.csv', expected / 'species_pair_summary_expected.csv', COMPARE['species_pair_summary.csv'])
    print('Minimal test dataset: PASS')


if __name__ == '__main__':
    main()
