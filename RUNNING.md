# Running the supplied code

## Verify this public download

`python code/verify_public_release.py` checks the exact public file inventory and SHA256 values without importing scientific libraries. A Git checkout's `.git` metadata and generated Python caches are ignored.

## Original clinical analysis

`code/frozen_scoring_core.py` retains the original numerical calculations. `config/scoring.json` defines the expected expression inputs and reference results. `code/clinical_cohort_recheck.R` uses the input columns specified in `config/clinical_cohorts.tsv`; R versions are recorded in `environments/R_packages_expected.tsv`.

The original internal entry points and `code/verify_package.py` are retained for provenance. They refer to the complete internal package and cannot be used as a successful public-package replay without obtaining the missing original inputs. Use `verify_public_release.py` for this release. No patient inputs are downloaded automatically.

## SIDISH input reconstruction (requires separate original inputs)

Prepare a separate Python environment matching `modules/SIDISH_input_bridge_reconstruction/input_bridge_reconstruction/audit/requirements_verified_environment.txt`. The input directory must contain `processed_adata_sidish.h5ad`, `processed_bulk_sidish.csv`, `sc_genes.tsv`, `sc_cells.tsv`, `sc_counts_genes_by_cells.mtx`, `sc_metadata.csv`, `bulk_expression_samples_by_gene.csv` and `bulk_survival.csv`.

```sh
python modules/SIDISH_input_bridge_reconstruction/scripts/reconstruct_sidish_inputs.py --input-root /path/to/original_inputs --output-root /path/to/new_reconstruction
python modules/SIDISH_input_bridge_reconstruction/scripts/verify_reconstructed_bridge.py --input-root /path/to/original_inputs --output-root /path/to/new_reconstruction
```

Use a new output directory outside both the original inputs and repository. Generated reconstruction outputs contain cell/participant-level information and should not be committed. `verify_count_export_binding.R` accepts the source RDS, count-export directory and existing audit-output directory as three positional arguments.
