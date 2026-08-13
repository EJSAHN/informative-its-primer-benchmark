# Informative-site ITS primer benchmark

This repository provides a data-quality-aware workflow for ITS primer benchmarking. It separates records that are informative for direct primer-template inference from records affected by missing sites, terminal sequence, or trusted metadata indicating use of the same primer.

## Analysis steps

1. Read species-organized FASTA records and a primer-pair catalog.
2. Locate primer sites with IUPAC-aware matching.
3. Exclude detected sites that overlap the terminal primer-length segment.
4. Exclude detected sites supported by primer-specific GenBank qualifiers such as `/PCR_primers`.
5. Retain ordinary ITS locus-name matches in descriptions, notes, and references as context-only audit evidence.
6. Score informative pairs under predefined mismatch and amplicon-length rules.
7. Export row-level calls, denominator summaries, filtering audits, sensitivity analyses, and an optional Excel workbook.

The repository produces numerical outputs only; figure rendering is not included.

## Repository contents

```text
run_analysis.py                         Complete numerical analysis
its_benchmark/                          Python modules
data/pair_catalog_auto.csv              Primer-pair catalog
data/accession_species_list.csv         Versioned accession list used in the case study
data/accession_sequence_manifest.csv    Sequence lengths and SHA-256 hashes
data/study_expected_metrics.json        Optional worked-example validation targets
tests/minimal_data/                     Synthetic test inputs
tests/expected/                         Expected test outputs
requirements.txt                        Exact software versions used for the analysis
```

## Installation

Python 3.10 is recommended.

```bash
python -m venv .venv
```

Linux/macOS:

```bash
source .venv/bin/activate
```

Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Minimal test

```bash
python tests/run_minimal_test.py --project-root .
```

The test covers an informative pair, a terminal site, a missing reverse site, trusted `/PCR_primers` evidence, an ambiguous binding window, and context-only ITS locus labels.

## Obtain the case-study records

The repository includes versioned accessions but not the downloaded sequence files. Fetch them from NCBI:

```bash
python -m its_benchmark.fetch_genbank_records \
  --accessions data/accession_species_list.csv \
  --fasta-root data/raw_fasta \
  --genbank-dir data/genbank \
  --email you@example.org
```

Validate the downloaded sequences against the frozen hash manifest:

```bash
python -m its_benchmark.validate_sequences \
  --fasta-root data/raw_fasta \
  --manifest data/accession_sequence_manifest.csv \
  --out results/sequence_manifest_validation.csv
```

## Run the complete numerical analysis

```bash
python run_analysis.py \
  --fasta-root data/raw_fasta \
  --pair-catalog data/pair_catalog_auto.csv \
  --genbank-dir data/genbank \
  --validate-study
```

Omit `--validate-study` when applying the workflow to another dataset.

## Main scoring defaults

- Site search: no more than 4 mismatches per primer.
- Strict call: no more than 2 total mismatches across the pair.
- Strict call: no 3′-terminal mismatch on either primer.
- Strict call: product length within the pair-specific window.
- Terminal filter: one primer length from either record end.
- Metadata exclusion: trusted primer-specific fields only.
- Ambiguity follow-up: binding windows with more than 5% ambiguous bases are flagged but not automatically excluded.

## Main outputs

```text
results/main/site_flags.csv
results/main/pair_calls.csv
results/main/species_pair_summary.csv
results/main/denominator_summary.csv
results/audit/flag_source_summary.csv
results/audit/mismatch_sensitivity.csv
results/audit/deduplication_sensitivity.csv
results/audit/terminal_scenario_totals.csv
results/benchmark_results.xlsx
results/sha256_manifest.tsv
```

## Interpretation

`coverage_pct_eligible` is the strict-hit percentage among eligible informative rows. `coverage_pct_all` uses all deduplicated records for the corresponding species and primer pair and is retained only to quantify denominator distortion. Rarefaction outputs summarize observed internal binding-site words and are not evidence of biological haplotypes without independent validation.
