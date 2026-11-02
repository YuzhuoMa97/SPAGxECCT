# Validation records and interpretation

This maintenance candidate changes production implementation structure while
preserving the documented statistical policy. The public source baseline is
recorded separately from the not-yet-pushed local maintenance source.
Read [release scope](../RELEASE_STATUS.md) and
[maintenance notes](../MAINTENANCE_NOTES.md).

## Current candidate

- Fresh Windows build and numerical units passed.
- Seventeen regression suites passed. An initial private driver path-escaping
  defect affected two reference invocations; only those two were rerun after
  fixing the driver, with the original diagnostics retained.
- Previous/current full walkthroughs produced 85 byte-identical files and
  145,985 identical numeric fields.
- Focused scanner, sample mapping, file-resource, null-model specification,
  outlier and CGF checks supplement the full application regressions.
- Linux results are recorded explicitly in ACCEPTANCE.json when complete.

Developer scripts regenerate their summaries under cpp/results/. Reference
tests may use Python or R; analysis and bundled tutorials do not require them.
Matching old/current implementations is reproducibility evidence. Independent
scalar/R/high-precision checks and population-level calibration answer
different questions.

## Remaining limits

macOS remains unvalidated. The current maintenance source has not yet been
pushed or checked by hosted CI. Binary compatibility depends on OS, CPU and
runtime libraries; no finite test set guarantees every input, extreme tail,
historical default P value or biobank-scale speed claim.

Preserved R code keeps its original dependency chain. The original full-CGF
R baseline and the accelerated default need not have identical P values in
every case. Historical results are not relabeled as new tests.
