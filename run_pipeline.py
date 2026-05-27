"""Run the ITS primer benchmark pipeline.

This script assumes FASTA files are organized by species folders and a pair catalog
is provided. It produces CSV outputs and an optional consolidated Excel workbook.
"""
import argparse
import subprocess
import sys
from pathlib import Path


def run(cmd):
    print('\n' + ' '.join(cmd))
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description='Run informative-site ITS primer benchmarking.')
    ap.add_argument('--fasta-root', required=True)
    ap.add_argument('--pair-catalog', default='data/pair_catalog_auto.csv')
    ap.add_argument('--genbank-dir', default='')
    ap.add_argument('--outdir', default='results')
    ap.add_argument('--terminal-mode', choices=['fixed', 'primer_length'], default='primer_length')
    ap.add_argument('--terminal-bp', type=int, default=30)
    ap.add_argument('--export-workbook', action='store_true')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    tables = outdir / 'tables'
    main_subdir = tables / 'ncbi_primerlength_metadata'
    tables.mkdir(parents=True, exist_ok=True)
    main_subdir.mkdir(parents=True, exist_ok=True)

    py = sys.executable
    accessions = tables / 'accessions.tsv'
    flags = tables / 'site_flags_ncbi_primerlength_metadata.csv'
    calls = tables / 'calls_ncbi_primerlength_metadata.csv'

    run([py, '-m', 'its_benchmark.make_accession_list', '--fasta-root', args.fasta_root, '--out', str(accessions)])
    flag_cmd = [py, '-m', 'its_benchmark.flag_informative_sites', '--fasta-root', args.fasta_root,
                '--pair-catalog', args.pair_catalog, '--terminal-mode', args.terminal_mode,
                '--terminal-bp', str(args.terminal_bp), '--out', str(flags)]
    if args.genbank_dir:
        flag_cmd.extend(['--genbank-dir', args.genbank_dir])
    run(flag_cmd)
    run([py, '-m', 'its_benchmark.scan_pairs', '--fasta-root', args.fasta_root, '--pair-catalog', args.pair_catalog,
         '--site-flags', str(flags), '--out', str(calls)])
    run([py, '-m', 'its_benchmark.summarize', '--calls', str(calls), '--outdir', str(main_subdir)])

    if args.export_workbook:
        run([py, '-m', 'its_benchmark.export_workbook', '--tables-dir', str(tables),
             '--main-subdir', 'ncbi_primerlength_metadata', '--out', str(outdir / 'Supplementary_Data_S1.xlsx')])


if __name__ == '__main__':
    main()
