import argparse
from pathlib import Path
import pandas as pd
from .common import iter_species_fastas, sequence_hash


def main():
    parser = argparse.ArgumentParser(description='Validate local FASTA records against an accession and SHA-256 manifest.')
    parser.add_argument('--fasta-root', required=True)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--out', default='')
    args = parser.parse_args()

    expected = pd.read_csv(args.manifest)
    observed_rows = []
    for species, accession, path, header, sequence in iter_species_fastas(Path(args.fasta_root)):
        observed_rows.append({
            'species': species,
            'accession': accession,
            'sequence_length_observed': len(sequence),
            'sequence_sha256_observed': sequence_hash(sequence),
        })
    observed = pd.DataFrame(observed_rows)
    merged = expected.merge(observed, on=['species', 'accession'], how='outer', indicator=True)
    merged['length_matches'] = merged['sequence_length'].eq(merged['sequence_length_observed'])
    merged['hash_matches'] = merged['sequence_sha256'].eq(merged['sequence_sha256_observed'])
    merged['passed'] = merged['_merge'].eq('both') & merged['length_matches'] & merged['hash_matches']

    if args.out:
        output = Path(args.out)
        output.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output, index=False)

    failures = merged[~merged['passed']]
    print(f'Expected accessions: {len(expected)}')
    print(f'Observed accessions: {len(observed)}')
    print(f'Validated accessions: {int(merged.passed.sum())}')
    if not failures.empty:
        print(f'Validation failures: {len(failures)}')
        raise SystemExit(1)
    print('Sequence manifest validation: PASS')


if __name__ == '__main__':
    main()
