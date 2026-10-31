# Validation records and interpretation

The current candidate is a local repair of public development commit
`ce9c2473649e0be2a8b70dd2f6186a03b25d21aa`. Read
[release status](../RELEASE_STATUS.md) for the distinction between the failed
public Actions run and locally passing repaired source.

The [acceptance summary](ACCEPTANCE.json) records the current Windows/Linux
scope and retains the earlier acceptance object under a historical key.
Historical source counts, hashes and test dates describe those snapshots;
they are not the identities of a later repaired release asset.

## Current repaired candidate

- CRLF/LF/no-final-newline BIM metadata and cross-format integration are
  explicitly tested in cpp/tests/run_integration.py.
- The Windows full-static link mode is checked from the compiler-driver
  command and PE imports, with native minimal-PATH tutorials.
- Windows numerical units, 17 regression suites, six tutorials and the full
  33-group walkthrough have passing local records.
- Ubuntu 24.04.5 x86-64/WSL2, GCC 13.3: 22 stages cover units, 14 CI suites,
  integration, tutorials, bundle extraction and PREFIX installation.
- Cross-platform association/CCT values agree within the documented floating
  point tolerance. This does not replace actual-R, independent scalar or
  high-precision references, nor population-level type-I error studies.

Developer scripts regenerate their own JSON summaries under cpp/results/.
Analysis and example execution do not require Python or R; developer/reference
tests may require them. Source-only checks and actual binary execution are
different evidence.

## Remaining limits

The original public CI run failed. A repaired local source tree is not a new
remote commit, and its checks must run again after the maintainer pushes it.
macOS remains unverified. No finite test collection proves universal extreme-
tail calibration, every historical default P value, every possible input or
a biobank-scale speed advantage. Native Linux testing under WSL2 does not
measure bare-metal biobank throughput or support all older glibc releases.

Original R/full-CGF comparisons apply only to matched common domains. The
guarded partial default deliberately adds an approximation; full-CGF candidate
diagnostics are not a P-value error bound. The historical R source and interface
remain available, with their own versions and dependencies.
