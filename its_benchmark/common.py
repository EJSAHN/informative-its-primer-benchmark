from __future__ import annotations
from pathlib import Path
from typing import List
import re

IUPAC = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'T'},
    'R': {'A','G'}, 'Y': {'C','T'}, 'S': {'G','C'}, 'W': {'A','T'},
    'K': {'G','T'}, 'M': {'A','C'}, 'B': {'C','G','T'}, 'D': {'A','G','T'},
    'H': {'A','C','T'}, 'V': {'A','C','G'}, 'N': {'A','C','G','T'}
}
REGEX_BASE = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'U':'T',
    'R':'[AG]', 'Y':'[CT]', 'S':'[GC]', 'W':'[AT]',
    'K':'[GT]', 'M':'[AC]', 'B':'[CGT]', 'D':'[AGT]',
    'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'
}
COMP = str.maketrans('ACGTRYSWKMBDHVNacgtryswkmbdhvnUu', 'TGCAYRSWMKVHDBNtgcayrswmkvhdbnAa')

SPECIES_DIRS = {
    'Mperniciosa': 'M. perniciosa',
    'Mroreri': 'M. roreri',
    'Ppalmivora': 'P. palmivora',
    'Pmegakarya': 'P. megakarya',
    'Moniliophthora_perniciosa': 'M. perniciosa',
    'Moniliophthora_roreri': 'M. roreri',
    'Phytophthora_palmivora': 'P. palmivora',
    'Phytophthora_megakarya': 'P. megakarya',
}
NCBI_SPECIES = {
    'M. perniciosa': 'Moniliophthora perniciosa',
    'M. roreri': 'Moniliophthora roreri',
    'P. palmivora': 'Phytophthora palmivora',
    'P. megakarya': 'Phytophthora megakarya',
}

def primer_to_regex(primer: str) -> str:
    return ''.join(REGEX_BASE.get(b.upper(), '[ACGT]') for b in primer.upper().replace('U','T'))

def rc(seq: str) -> str:
    return seq.upper().translate(COMP)[::-1].replace('U','T')

def fasta_records(path: Path):
    header = None
    chunks: List[str] = []
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks).upper().replace('U','T')
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, ''.join(chunks).upper().replace('U','T')

def accession_from_header(header: str, fallback: str = '') -> str:
    token = header.split()[0] if header else fallback
    if '|' in token:
        parts = [p for p in token.split('|') if p and p not in ('gb','emb','dbj','ref')]
        if parts:
            token = parts[0]
    return token.strip()

def infer_species_from_folder(name: str):
    if name in SPECIES_DIRS:
        return SPECIES_DIRS[name]
    low = name.lower()
    if 'perniciosa' in low:
        return 'M. perniciosa'
    if 'roreri' in low:
        return 'M. roreri'
    if 'palmivora' in low:
        return 'P. palmivora'
    if 'megakarya' in low:
        return 'P. megakarya'
    return None

def iter_species_fastas(fasta_root: Path):
    fasta_root = Path(fasta_root)
    for sub in sorted(fasta_root.iterdir()):
        if not sub.is_dir():
            continue
        species = infer_species_from_folder(sub.name)
        if species is None:
            continue
        files = []
        for pattern in ('*.fa', '*.fasta', '*.fas', '*.fna'):
            files.extend(sub.glob(pattern))
        for fp in sorted(files):
            for header, seq in fasta_records(fp):
                acc = accession_from_header(header, fp.stem)
                yield species, acc, fp, header, seq

def allowed(base: str):
    return IUPAC.get(base.upper(), {'A','C','G','T'})

def mismatch_count_limited(primer: str, window: str, max_mm: int) -> int:
    primer = primer.upper().replace('U','T')
    window = window.upper().replace('U','T')
    mm = 0
    for p, w in zip(primer, window):
        if w not in allowed(p):
            mm += 1
            if mm > max_mm:
                return mm
    return mm

def three_prime_mismatch(primer: str, window: str, k: int = 1) -> bool:
    primer = primer.upper().replace('U','T')
    window = window.upper().replace('U','T')
    if k <= 0:
        return False
    return any(w not in allowed(p) for p, w in zip(primer[-k:], window[-k:]))

def scan_best(seq: str, primer: str, max_mm: int = 4):
    seq = seq.upper().replace('U','T')
    primer = primer.upper().replace('U','T')
    L = len(primer)
    if len(seq) < L:
        return None
    pat = primer_to_regex(primer)
    m = re.search(pat, seq)
    if m:
        w = seq[m.start():m.start()+L]
        return {'pos': m.start(), 'mm': 0, 'window': w, 'len': L, 'any_3p': False}
    nblocks = max_mm + 1
    cuts = [round(i * L / nblocks) for i in range(nblocks + 1)]
    candidates = set()
    for b in range(nblocks):
        a, c = cuts[b], cuts[b+1]
        seed = primer[a:c]
        if not seed:
            continue
        for sm in re.finditer(primer_to_regex(seed), seq):
            start = sm.start() - a
            if 0 <= start <= len(seq) - L:
                candidates.add(start)
    best = None
    for i in sorted(candidates):
        w = seq[i:i+L]
        mm = mismatch_count_limited(primer, w, max_mm)
        if mm <= max_mm:
            rec = {'pos': i, 'mm': mm, 'window': w, 'len': L, 'any_3p': three_prime_mismatch(primer, w)}
            if best is None or (mm, i) < (best['mm'], best['pos']):
                best = rec
                if mm == 0:
                    break
    return best

def scan_pair(seq: str, fwd: str, rev: str, max_mm: int = 4):
    f = scan_best(seq, fwd, max_mm=max_mm)
    r_on_rc = scan_best(rc(seq), rev, max_mm=max_mm)
    if r_on_rc:
        r_start_plus = len(seq) - (r_on_rc['pos'] + len(rev))
        r_end_plus = len(seq) - r_on_rc['pos']
        r = dict(r_on_rc)
        r['start_plus'] = r_start_plus
        r['end_plus'] = r_end_plus
    else:
        r = None
    if f:
        f = dict(f)
        f['start_plus'] = f['pos']
        f['end_plus'] = f['pos'] + len(fwd)
    return f, r

def is_terminal_hit(start: int, end: int, seq_len: int, terminal_bp: int) -> bool:
    return start < terminal_bp or end > (seq_len - terminal_bp)

def primer_appears_in_metadata(text: str, primer_name: str, primer_seq: str) -> bool:
    t = (text or '').lower()
    if primer_name and primer_name.lower() in t:
        return True
    p = (primer_seq or '').lower().replace('u','t')
    return bool(p and p in t)
