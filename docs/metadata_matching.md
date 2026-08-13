# GenBank metadata matching

A detected primer site is excluded by metadata only when the same primer name or sequence occurs in a primer-specific qualifier, including `/PCR_primers`, `/forward_primer`, `/reverse_primer`, or `/sequencing_primer`.

Exact ITS labels found in descriptions, references, notes, or ordinary feature annotations are recorded as context-only evidence. They do not exclude a site because terms such as `ITS1` and `ITS2` can denote rDNA regions rather than primers used to generate a record.

A missing metadata match is not evidence that a primer was absent during sequence generation. Terminal filtering is applied independently as a complementary safeguard.
