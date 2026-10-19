# Release notes

## Data and model licensing update — 2026-10-19

The authors authorized CC BY 4.0 for their publicly released original aggregate tables and frozen model definitions. The exact 32-file scope is recorded in DATA_LICENSE_MANIFEST.tsv; LICENSE-DATA.md provides the notice, legal-code link and attribution information. The original software remains under MIT. README.md, DATA_ACCESS.md, LICENSE_SCOPE.md and the public checksum inventory were synchronized. No scientific table, model value or analysis code was changed, and no patient-level or raw sequencing data were added. The repository now contains 61 files, including the two new licensing records. Patient-data access arrangements remain pending confirmation.

## First public code and aggregate-data release — 2026-10-19

This release was assembled from the verified internal delivery. Core scoring and clinical calculations were retained. The two reconstruction Python scripts now take explicit input paths instead of a workstation-specific default. Aggregate table source paths were reduced to basenames; all other table cells and scientific estimates were preserved. Individual accessions in the compatibility summary were omitted.

Included result tables are research aggregates. Internal participant-level data, raw sequencing matrices, case consent/ethics records, internal logs, mixed-license archives and third-party models/force fields were excluded. The full internal-replay scripts are retained as source code but their complete input package is not publicly distributed.

This release does not claim new model training, fresh R installation, raw-data-to-all-figures reproduction, or calculation of the unavailable complete-set external random-null P values. Public file integrity, Python syntax and preservation of aggregate values were checked during release preparation.
