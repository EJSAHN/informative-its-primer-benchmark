# Informative-site ITS primer benchmark

This repository provides a data-quality-aware workflow for ITS primer benchmarking. It flags terminal, primer-derived, and otherwise non-informative primer sites before scoring primer-template mismatches, so coverage estimates are based on records that are informative for direct primer-template inference.

## What the workflow does

1. Reads species-organized FASTA files and a primer-pair catalog.
2. Flags primer sites that are terminal, metadata-derived, missing, or informative internal sequence.
3. Scores informative primer pairs under strict mismatch and amplicon-length rules.
4. Summarizes coverage over informative/eligible denominators and archive-wide denominators.
5. Exports row-level calls, deduplicated calls, summary tables, and an optional consolidated Excel workbook.

The workflow exports numerical tables only; figure rendering is not included.

## Repository contents

```text
run_pipeline.py                  Pipeline runner
its_benchmark/                   Python modules
data/pair_catalog_auto.csv       Primer-pair catalog used by the benchmark
data/accession_species_list.csv  GenBank accessions and target labels used in the worked example
data/species_counts.csv          Species-level accession counts for the worked example
requirements.txt                 Python dependencies
```

The worked-example FASTA and GenBank records are not bundled in this repository. They can be regenerated from the accession list or supplied as local FASTA folders with the same species labels.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scriptsctivate
pip install -r requirements.txt
```

## Input organization

Place FASTA files in species folders under a common root, for example:

```text
data/raw_fasta/
  Mperniciosa/
  Mroreri/
  Ppalmivora/
  Pmegakarya/
```

The default pair catalog is included at:

```text
data/pair_catalog_auto.csv
```

## Run the workflow

```bash
python run_pipeline.py   --fasta-root data/raw_fasta   --pair-catalog data/pair_catalog_auto.csv   --genbank-dir data/genbank   --terminal-mode primer_length   --export-workbook
```

If GenBank metadata are not available, omit `--genbank-dir`.

## Main output files

```text
results/tables/site_flags_ncbi_primerlength_metadata.csv
results/tables/calls_ncbi_primerlength_metadata.csv
results/tables/ncbi_primerlength_metadata/Table_all_pairs_by_species_informative.csv
results/tables/ncbi_primerlength_metadata/calls_deduplicated.csv
results/Supplementary_Data_S1.xlsx
```

## Strict scoring defaults

- Primer-site search tolerance: <=4 mismatches per primer.
- Strict hit: <=2 total mismatches across the primer pair.
- Strict hit: no 3-prime terminal mismatch on either primer.
- Strict hit: product length within the predefined pair-specific window.
- Main terminal filter: primer-length terminal window.

## Interpretation

Archive-wide coverage is reported only as a denominator-error diagnostic. Recommended performance estimates should be interpreted using informative and eligible denominators.
