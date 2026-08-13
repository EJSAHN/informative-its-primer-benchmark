# Data files

- `pair_catalog_auto.csv`: primer pairs and product-size windows used by the benchmark.
- `accession_species_list.csv`: versioned accessions and taxon labels used in the cacao case study.
- `species_counts.csv`: accession counts by taxon.
- `accession_sequence_manifest.csv`: expected sequence lengths and SHA-256 hashes for the downloaded records.
- `site_variant_targets.csv`: primer sites used for the observed-variant rarefaction summaries.
- `study_expected_metrics.json`: optional regression targets for the worked example.

FASTA and GenBank files are not bundled. Use `its_benchmark.fetch_genbank_records` to retrieve the versioned accessions, then run `its_benchmark.validate_sequences` before reproducing the case-study analysis.
