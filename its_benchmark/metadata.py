from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import re
from Bio import SeqIO


def _flatten(value) -> str:
    if value is None:
        return ''
    if isinstance(value, (list, tuple, set)):
        return ' '.join(_flatten(v) for v in value)
    return str(value)


def extract_metadata_fields(record) -> Dict[str, str]:
    """Extract structured GenBank metadata fields from a SeqRecord."""
    fields: Dict[str, str] = {'description': record.description or ''}
    for key, value in (record.annotations or {}).items():
        fields[f'annotation:{key}'] = _flatten(value)
    for idx, feature in enumerate(record.features or []):
        feature_type = getattr(feature, 'type', 'feature')
        for key, value in (getattr(feature, 'qualifiers', {}) or {}).items():
            fields[f'feature:{idx}:{feature_type}:{key}'] = _flatten(value)
    return fields


def load_genbank_metadata(genbank_dir: str | Path) -> Dict[str, Dict[str, str]]:
    """Load GenBank metadata and index records by accession with and without version."""
    root = Path(genbank_dir)
    metadata: Dict[str, Dict[str, str]] = {}
    if not root.exists():
        return metadata
    for path in root.rglob('*.gb'):
        try:
            record = next(SeqIO.parse(str(path), 'genbank'))
            fields = extract_metadata_fields(record)
            keys = {
                record.id,
                record.id.split('.')[0],
                getattr(record, 'name', record.id),
                getattr(record, 'name', record.id).split('.')[0],
            }
            for key in keys:
                metadata[key] = fields
        except Exception as exc:
            print(f'WARN: could not parse {path}: {exc}')
    return metadata


def _canonical_name_token(token: str) -> str:
    return re.sub(r'[_-]+', '', token.upper())


def _name_match(text: str, primer_name: str) -> bool:
    target = _canonical_name_token(primer_name)
    if not target:
        return False
    tokens = re.findall(r'[A-Za-z0-9_-]+', text or '')
    return any(_canonical_name_token(token) == target for token in tokens)


def _sequence_match(text: str, primer_seq: str) -> bool:
    primer = primer_seq.upper().replace('U', 'T').replace(' ', '')
    if not primer:
        return False
    separator = r"[\s\-′’']*"
    pattern = r'(?<![A-Za-z])' + separator.join(re.escape(base) for base in primer) + r'(?![A-Za-z])'
    return re.search(pattern, text or '', flags=re.IGNORECASE) is not None


def _qualifier_name(field_name: str) -> str:
    return field_name.split(':')[-1] if ':' in field_name else field_name


def _is_trusted_primer_field(field_name: str) -> bool:
    """Identify fields explicitly intended to record primer use.

    Generic descriptions, locus names, references, notes, and comments are kept
    as context-only audit evidence and never exclude a site.
    """
    key = re.sub(r'[^a-z0-9]', '', _qualifier_name(field_name).lower())
    return key in {
        'pcrprimers', 'pcrprimer',
        'forwardprimer', 'reverseprimer',
        'fwdprimer', 'revprimer',
        'primerforward', 'primerreverse',
        'sequencingprimer', 'sequencingprimers',
        'primer', 'primers',
    }


def match_primer_metadata(fields: Dict[str, str], primer_name: str, primer_seq: str) -> Dict[str, object]:
    """Return trusted exclusion evidence and context-only metadata matches.

    A site is excluded by metadata only when the primer name or sequence appears
    in a primer-specific qualifier such as ``/PCR_primers``. Matches in other
    metadata fields are recorded for audit but do not alter site classification.
    """
    if not fields:
        return {
            'metadata_available': False,
            'matched': False,
            'name_match': False,
            'sequence_match': False,
            'match_fields': '',
            'match_terms': '',
            'context_name_match': False,
            'context_sequence_match': False,
            'context_match_fields': '',
            'context_match_terms': '',
        }

    trusted_fields: List[str] = []
    context_fields: List[str] = []
    trusted_name = trusted_sequence = False
    context_name = context_sequence = False

    for field_name, text in fields.items():
        name_found = _name_match(text, primer_name)
        sequence_found = _sequence_match(text, primer_seq)
        if _is_trusted_primer_field(field_name):
            if name_found or sequence_found:
                trusted_fields.append(field_name)
            trusted_name = trusted_name or name_found
            trusted_sequence = trusted_sequence or sequence_found
        else:
            if name_found or sequence_found:
                context_fields.append(field_name)
            context_name = context_name or name_found
            context_sequence = context_sequence or sequence_found

    trusted_terms = []
    if trusted_name:
        trusted_terms.append('primer_name')
    if trusted_sequence:
        trusted_terms.append('primer_sequence')

    context_terms = []
    if context_name:
        context_terms.append('primer_name')
    if context_sequence:
        context_terms.append('primer_sequence')

    return {
        'metadata_available': True,
        'matched': bool(trusted_name or trusted_sequence),
        'name_match': bool(trusted_name),
        'sequence_match': bool(trusted_sequence),
        'match_fields': ';'.join(sorted(set(trusted_fields))),
        'match_terms': ';'.join(trusted_terms),
        'context_name_match': bool(context_name),
        'context_sequence_match': bool(context_sequence),
        'context_match_fields': ';'.join(sorted(set(context_fields))),
        'context_match_terms': ';'.join(context_terms),
    }
