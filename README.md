# DLBCL spatial mapping and survival analysis

Analysis code, frozen model definitions and aggregate results for **Spatial mapping of survival-associated cell states and prognostic evaluation of a frozen gene pool in diffuse large B-cell lymphoma**.

Authors: Kangzhe Ruan, Bingyue Cao and Yongchang Wei. Kangzhe Ruan and Bingyue Cao contributed equally. Correspondence: weiyongchang8@163.com.

## Public release scope

This first public release contains original scoring/clinical-analysis code, SIDISH input-reconstruction code, frozen coefficient and gene-weight definitions, and cohort/gene/random-pool aggregate result tables. No participant-level clinical records, per-patient expression matrices, raw sequencing data or individual case documents are included. It is not a complete raw-data-to-figures reproduction package. See [data access](DATA_ACCESS.md), [licensing](LICENSE_SCOPE.md), and [release notes](RELEASE_NOTES.md).

## Contents

|Directory|Contents|
|---|---|
|`code/`|Frozen scoring functions, clinical R analysis and internal-replay entry points|
|`config/`|Frozen cohort and scoring settings; patient-level input paths are requirements, not bundled files|
|`models/`|Frozen coefficients, development scales, ordered genes and module weights|
|`modules/SIDISH_input_bridge_reconstruction/scripts/`|Input reconstruction and independent verification scripts|
|`source_data/`|Aggregate clinical, bridge, calibration, stability and random-pool coverage results|
|`provenance/`|Required original input hashes and list of inputs omitted from this public subset|
|`checksums/`|SHA256 inventory of the public files|

## Verification and use

Python 3.12 is used by the verification utility. From the repository root, run:

```sh
python code/verify_public_release.py
```

The calculation core in `code/frozen_scoring_core.py` accepts a patient-by-gene expression DataFrame together with the supplied gene universe and frozen model-state JSON. Outputs are continuous research scores, not validated clinical predictions. The original analysis entry points `reproduce_frozen_scores.py`, `reproduce_clinical.py` and `run_scoped_replay.py` require separately available original inputs; they cannot complete using this public subset alone. `provenance/unbundled_inputs.tsv` lists the omitted frozen inputs. Do not substitute aggregate tables for patient-level inputs or disable hash checks.

The clinical replay environment is recorded in `requirements-lock.txt` and `environments/R_packages_expected.tsv`. SIDISH reconstruction uses a separate environment recorded in its module; do not combine these locks. See [running the code](RUNNING.md).

## Scientific interpretation

The clinical score remains frozen. All 499 random candidate pools completed training; 330 met external nonzero-gene coverage requirements and 169 had coverage gaps. Complete-set external random-null performance and empirical P values were not estimated. Negative findings and gene-removal sensitivity remain in the result tables.

SIDISH preprocessing was reconstructed and matched saved inputs (29,668 cells, 2,154 genes and 470 bulk samples). This is not recovery of the original study-specific training entry point and no retraining is claimed. Earlier scoped internal replay does not imply that all historical workflows or every figure can be regenerated from this public subset.

## License

Original software: MIT, copyright Kangzhe Ruan and Bingyue Cao. Data/model reuse terms and third-party exceptions are described in [LICENSE_SCOPE.md](LICENSE_SCOPE.md).
