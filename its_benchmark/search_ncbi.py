"""Optional utility to download ITS/rDNA records from NCBI nucleotide for the four cacao taxa."""
import argparse, time
from pathlib import Path
from Bio import Entrez, SeqIO
from .common import NCBI_SPECIES

def main():
    ap = argparse.ArgumentParser(description='Search and download ITS/rDNA records from NCBI nucleotide.')
    ap.add_argument('--email', required=True, help='NCBI requires a real contact email.')
    ap.add_argument('--out-fasta-root', required=True)
    ap.add_argument('--out-genbank', required=True)
    ap.add_argument('--max-records-per-species', type=int, default=5000)
    ap.add_argument('--sleep', type=float, default=0.34)
    args = ap.parse_args()
    Entrez.email = args.email
    out_fasta = Path(args.out_fasta_root); out_fasta.mkdir(parents=True, exist_ok=True)
    out_gb = Path(args.out_genbank); out_gb.mkdir(parents=True, exist_ok=True)
    query_template = '({organism}[Organism]) AND ("internal transcribed spacer" OR ITS OR ribosomal RNA OR rDNA) NOT mitochondrion'
    for short, organism in NCBI_SPECIES.items():
        term = query_template.format(organism=organism)
        print(f'\nSearching {organism}: {term}')
        h = Entrez.esearch(db='nuccore', term=term, retmax=args.max_records_per_species)
        rec = Entrez.read(h); h.close()
        ids = rec.get('IdList', [])
        folder = short.replace('. ', '').replace(' ', '')
        fasta_dir = out_fasta / folder; fasta_dir.mkdir(parents=True, exist_ok=True)
        gb_dir = out_gb / folder; gb_dir.mkdir(parents=True, exist_ok=True)
        for i in range(0, len(ids), 100):
            chunk = ids[i:i+100]
            time.sleep(args.sleep)
            h = Entrez.efetch(db='nuccore', id=','.join(chunk), rettype='gb', retmode='text')
            records = list(SeqIO.parse(h, 'genbank')); h.close()
            for rec in records:
                SeqIO.write(rec, gb_dir / f'{rec.id}.gb', 'genbank')
                SeqIO.write(rec, fasta_dir / f'{rec.id}.fasta', 'fasta')
        print(f'Downloaded {len(ids)} records into {fasta_dir} and {gb_dir}')

if __name__ == '__main__':
    main()
