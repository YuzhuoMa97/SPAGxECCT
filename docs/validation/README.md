# Validation records and interpretation

This maintenance candidate changes production implementation structure while
preserving the documented statistical policy. The historical source baseline, local validation and subsequently published
maintenance commit are recorded separately.
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

Published develop commit d90fb4c45849757c816c5d1ba5d191911443e83c was
checked by [hosted CI](https://github.com/YuzhuoMa97/SPAGxECCT/actions/runs/41416750904):
Windows and Ubuntu passed, while the macOS build/unit job failed.
macOS remains unvalidated. These results apply to that exact commit and do
not certify later changes. Binary compatibility depends on OS, CPU and
runtime libraries; no finite test set guarantees every input, extreme tail,
historical default P value or biobank-scale speed claim.

Preserved R code keeps its original dependency chain. The original full-CGF
R baseline and the accelerated default need not have identical P values in
every case. Historical results are not relabeled as new tests.

## Portable reference calculations

Five checks previously relied on long double having extra precision or
range, or on a cancellation-prone linear probability. Their independent
oracles now use exact integer conditional-HWE enumeration with 80/120-digit
Decimal checks, or 50-digit multiprecision calculations. The generation
script and HWE fixtures are included in cpp/tests/. Original tolerances
and all 53 log-probability cases remain. No analysis formula was changed.
The revised complete unit target passed locally on Windows; actual macOS
CI must be rerun. Platform diagnostics print the native floating-point
precision so assumptions are visible in future test logs.
