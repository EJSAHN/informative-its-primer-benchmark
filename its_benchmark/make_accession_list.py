import argparse
from pathlib import Path
from .common import iter_species_fastas


def main():
    parser = argparse.ArgumentParser(description='Create an accession list from species-organized FASTA folders.')
    parser.add_argument('--fasta-root', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    fasta_root = Path(args.fasta_root).resolve()
    rows = []
    for species, accession, path, header, sequence in iter_species_fastas(fasta_root):
        rows.append((species, accession, str(Path(path).resolve().relative_to(fasta_root)), len(sequence)))

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', encoding='utf-8') as handle:
        handle.write('accession\tspecies\tsequence_length\tsource_fasta\n')
        for species, accession, path, length in sorted(rows):
            handle.write(f'{accession}\t{species}\t{length}\t{path}\n')
    print(f'Wrote {len(rows)} accession rows to {output}')


if __name__ == '__main__':
    main()
