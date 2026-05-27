import argparse, time
from pathlib import Path
from Bio import Entrez, SeqIO

def main():
    ap = argparse.ArgumentParser(description='Fetch GenBank flatfiles for accessions listed by make_accession_list.')
    ap.add_argument('--accessions', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--email', required=True)
    ap.add_argument('--batch-size', type=int, default=100)
    ap.add_argument('--sleep', type=float, default=0.34)
    args = ap.parse_args()
    Entrez.email = args.email
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    with open(args.accessions, encoding='utf-8') as handle:
        for line in handle:
            if not line.strip():
                continue
            acc, species, n, fp = line.rstrip('\n').split('\t')[:4]
            rows.append((acc, species))
    for i in range(0, len(rows), args.batch_size):
        chunk = rows[i:i+args.batch_size]
        ids = [r[0] for r in chunk]
        print(f'Fetching {i+1}-{i+len(chunk)} of {len(rows)}')
        time.sleep(args.sleep)
        h = Entrez.efetch(db='nuccore', id=','.join(ids), rettype='gb', retmode='text')
        records = list(SeqIO.parse(h, 'genbank'))
        h.close()
        for rec in records:
            SeqIO.write(rec, out / f'{rec.id}.gb', 'genbank')
    print(f'Done. GenBank files in {out}')

if __name__ == '__main__':
    main()
