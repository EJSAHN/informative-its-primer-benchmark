import argparse
from pathlib import Path
from .common import iter_species_fastas

def main():
    ap = argparse.ArgumentParser(description='Create an accession list from species-organized FASTA folders.')
    ap.add_argument('--fasta-root', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    rows = []
    for species, acc, fp, header, seq in iter_species_fastas(Path(args.fasta_root)):
        rows.append((species, acc, str(fp), len(seq)))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w', encoding='utf-8') as out:
        for species, acc, fp, n in sorted(rows):
            out.write(f'{acc}\t{species}\t{n}\t{fp}\n')
    print(f'Wrote {len(rows)} accession rows to {args.out}')

if __name__ == '__main__':
    main()
