# Validation records and interpretation

SPAGxE uses the documented paper-oriented statistical policy with guarded
partial-normal CGF acceleration. The checked public commit and historical
local implementation comparisons are recorded separately. Read
[release scope](../RELEASE_STATUS.md) and
[maintenance notes](../MAINTENANCE_NOTES.md).

## Latest checked public source

Develop commit `2a14885bd9ebc1bcf51910db852a38e212793225` passed the
[2 November 2026 hosted workflow](https://github.com/YuzhuoMa97/SPAGxECCT/actions/runs/41438512290).
All Windows, Ubuntu and macOS jobs completed their source build, numerical
unit, reader/regression and source/extracted-runtime tutorial checks.
The output comparison checks completed successfully. Check Actions for
later commits instead of assuming that this record covers them.

This verifies the tested source and CI runner routes. It does not publish a
macOS binary or guarantee compatibility with all OS/CPU/compiler versions.
At this review there was no C++ Release with public runtime attachments;
users can build the complete develop source. A later Release must identify
its own source commit, platform packages and checksums.

## Portable reference calculations

Earlier macOS failures exposed test references that assumed extra native
long-double precision/range or used cancellation-prone linear probabilities.
The corrected oracles use exact conditional-HWE fixtures, multiprecision
calculations and independent reference identities. The original accuracy
criteria remain. Production analysis formulas were not changed.

The final log-probability suite has 54 cases: the original 53 and one focused
normal-reference underflow-band regression. The CGF suite retains all 28
cases and the original identity grid. There are no newly added skips or
tolerance relaxations. Existing guards for 13 wider-long-double CGF
measurements remain; these measurements do not run on hosts without the
required native precision/range. The successful hosted run is not presented
as evidence that they ran on such hosts.

Native floating-point diagnostics remain visible in test output. Generation
scripts and HWE fixtures are included in cpp/tests/. A portable oracle
correction is distinct from a change to the statistical model.

## Earlier local implementation checks

- Fresh Windows build and numerical units passed.
- Seventeen local regression suites passed. Two reference invocations were
  rerun after correcting a private driver path-escaping error; their original
  diagnostics were retained.
- Previous/current full walkthroughs produced 85 byte-identical files and
  145,985 identical numeric fields.
- Scanner, sample-mapping, file-resource, null-model, outlier and CGF checks
  supplement those application regressions.
- Linux local build, numerical units, regression routes, complete walkthrough
  and installed-prefix execution passed, as recorded in ACCEPTANCE.json.

These are earlier local records, not counts newly extracted from the latest
hosted run. Earlier public CI results are retained in the machine-readable
history instead of replacing them with the newer successful result.

Developer scripts regenerate summaries under cpp/results/. Reference tests
may use Python or R; analysis and bundled tutorials do not require them.
Matching implementations is reproducibility evidence. Independent reference
checks and population-level calibration answer different questions.

## Scope and remaining limits

Binary compatibility depends on OS, CPU and runtime libraries. No finite test
set guarantees every input, extreme tail, historical default P value or
biobank-scale speed claim. The partial-normal diagnostic threshold is not a
uniform bound on P-value error.

Preserved R code keeps its original dependency chain. The original full-CGF
R baseline and accelerated default need not give identical P values in every
case. Historical results are not relabeled as newly executed tests.
