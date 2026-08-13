# Output dictionary

## Denominators

- `N`: deduplicated number of records for a species and primer pair.
- `informative`: records in which both primer sites are retained as informative internal sites.
- `eligible`: informative records with valid orientation and product-size window.
- `covered`: eligible records passing the strict mismatch rules.
- `coverage_pct_eligible`: `100 × covered / eligible`.
- `coverage_pct_all`: `100 × covered / N`.
- `delta_coverage_pp`: `coverage_pct_eligible - coverage_pct_all`.

## Metadata evidence

- `metadata_same_primer`: trusted primer-specific metadata indicates use of the same primer.
- `metadata_context_*`: a name or sequence match occurred in a non-primer-specific field and was retained for audit only.

## Flag-source categories

A metadata flag contributes to classification only when the corresponding primer site was detected.

- `terminal_only`
- `metadata_only`
- `both_terminal_and_metadata`
- `neither`
