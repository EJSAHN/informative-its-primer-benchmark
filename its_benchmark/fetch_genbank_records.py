import argparse
import csv
import time
from pathlib import Path
from Bio import Entrez, SeqIO

FOLDER_NAMES = {
    'M. perniciosa': 'Mperniciosa',
    'M. roreri': 'Mroreri',
    'P. palmivora': 'Ppalmivora',
    'P. megakarya': 'Pmegakarya',
}


def read_accessions(path: Path):
    text = path.read_text(encoding='utf-8-sig')
    first = text.splitlines()[0] if text.splitlines() else ''
    delimiter = ',' if ',' in first else '\t'
    rows = []
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    if reader.fieldnames and {'accession', 'species'}.issubset(reader.fieldnames):
        for row in reader:
            if row.get('accession') and row.get('species'):
                rows.append((row['accession'].strip(), row['species'].strip()))
        return rows
    # Backward-compatible headerless TSV: accession, species, length, source.
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) >= 2 and parts[0].lower() != 'accession':
            rows.append((parts[0].strip(), parts[1].strip()))
    return rows


def main():
    parser = argparse.ArgumentParser(description='Fetch FASTA and GenBank records for a versioned accession list.')
    parser.add_argument('--accessions', required=True, help='CSV or TSV with accession and species columns.')
    parser.add_argument('--fasta-root', required=True)
    parser.add_argument('--genbank-dir', required=True)
    parser.add_argument('--email', required=True, help='Contact email required by NCBI Entrez.')
    parser.add_argument('--api-key', default='', help='Optional NCBI API key.')
    parser.add_argument('--batch-size', type=int, default=100)
    parser.add_argument('--sleep', type=float, default=0.34)
    args = parser.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    rows = read_accessions(Path(args.accessions))
    fasta_root = Path(args.fasta_root)
    genbank_root = Path(args.genbank_dir)
    fasta_root.mkdir(parents=True, exist_ok=True)
    genbank_root.mkdir(parents=True, exist_ok=True)

    species_by_accession = {accession: species for accession, species in rows}
    for start in range(0, len(rows), args.batch_size):
        chunk = rows[start:start + args.batch_size]
        ids = [accession for accession, species in chunk]
        print(f'Fetching {start + 1}-{start + len(chunk)} of {len(rows)}')
        time.sleep(args.sleep)
        handle = Entrez.efetch(db='nuccore', id=','.join(ids), rettype='gb', retmode='text')
        records = list(SeqIO.parse(handle, 'genbank'))
        handle.close()
        for record in records:
            species = species_by_accession.get(record.id, species_by_accession.get(record.id.split('.')[0]))
            if species is None:
                print(f'WARN: species label not found for {record.id}; record skipped')
                continue
            folder = FOLDER_NAMES.get(species, species.replace('. ', '').replace(' ', ''))
            fasta_dir = fasta_root / folder
            gb_dir = genbank_root / folder
            fasta_dir.mkdir(parents=True, exist_ok=True)
            gb_dir.mkdir(parents=True, exist_ok=True)
            SeqIO.write(record, fasta_dir / f'{record.id}.fasta', 'fasta')
            SeqIO.write(record, gb_dir / f'{record.id}.gb', 'genbank')
    print(f'FASTA records: {fasta_root}')
    print(f'GenBank records: {genbank_root}')


if __name__ == '__main__':
    main()
