# SPAGxE 1.2.0 release scope

SPAGxE is the author-designated official standalone C++ distribution of the
SPAGxE method family, maintained by Yuzhuo Ma (yuzhuoma@amss.ac.cn).
It runs without R or Python. Preserved R sources retain their historical
interface, version and dependency chain.

## Methods and inputs

| Method | Scope |
|---|---|
| SPAGxE | Shared AF, unrelated subjects |
| SPAGxE+ | Shared AF, sparse-GRM score-variance correction |
| SPAGxEmix | Individual-specific AF, unrelated multi-ancestry subjects |
| SPAGxEmix_local | Ancestry-specific G/H conditional on known H |

All four methods use one paper-oriented statistical policy with guarded
partial-normal CGF acceleration. Included are additive genotype/dosage
BED/PGEN/BGEN/VCF/BCF, native linear/logistic/Breslow Cox/ordinal null models,
explicit CCT weights, AF diagnostics, local text/count/BGEN input and native
hard-count packing. See [feature coverage](FEATURE_COVERAGE.md).

SPAGxEmixPlus is planned for a future release and is not included. Current
mix/local modes reject GRM input. The homogeneous GRM workflow corrects score
variance and is not a general GLMM null fitter. Local inputs require known H;
the program does not infer ancestry or implement arbitrary uncertain-H joint
posterior tests. Arbitrary R model objects, dominant/recessive recoding and
general counting-process/stratified survival formulas are not native inputs.

## Maintenance revision

The current source includes maintainability changes in file resource handling,
subject mapping, null-model result assembly, outlier partitioning, numerical
helper selection and shared SIMD coefficients. Public result columns,
statistical formulas, arithmetic/reduction order and convergence rules remain
unchanged. See [maintenance notes](MAINTENANCE_NOTES.md).

The canonical sparse-GRM table option is --grm-file. An obsolete pre-release
GRM option alias was removed; update old development commands to --grm-file.
Both sparse-triplet and PLINK GRM input layouts remain supported.

## Validation interpretation

### Published source CI

Develop commit `2a14885bd9ebc1bcf51910db852a38e212793225` passed the
[2 November 2026 hosted workflow](https://github.com/YuzhuoMa97/SPAGxECCT/actions/runs/41438512290).
The Windows, Ubuntu and macOS jobs completed successfully, including source
build, numerical units, reader/regression checks, and source and extracted
runtime tutorials with matching output checks. This is actual hosted macOS
execution, in addition to the earlier local Windows/Linux evidence.

The final portable-reference update changed two unit-test reference files
and their source manifest; production analysis code did not change. The
log-probability suite retains its original 53 cases and adds one focused
reference regression (54 total). The CGF suite retains 28 cases and the
original identity grid and accuracy criteria. Existing guards for 13 CGF
measurements that require wider native long double remain; success on a
binary64-long-double host does not mean those wider-precision measurements
ran there. No new skips or tolerance relaxation were added by the update.

These results identify one checked commit and the three CI runner routes.
Later source or documentation revisions need their own CI. They do not
certify every operating system, CPU, compiler, input or binary distribution.
There is no published macOS precompiled package. Default x86 builds require
Haswell/AVX2; dynamic Linux binaries also require compatible glibc/libstdc++.

### Earlier local implementation checks

Fresh Windows source compilation and numerical units passed. All 17 local
regression suites passed, covering CLI/scope, global/local readers, AF, CCT,
independent references and stability. The six source tutorials passed.
The previous/current full 33-group walkthrough produced 85 byte-identical
files, including 145,985 numeric fields with maximum delta 0. Linux local
build, tests, install-prefix execution and walkthrough checks also passed.
The acceptance record preserves those historical stages and dates; these
earlier figures are not relabeled as newly executed hosted measurements.

Earlier public commits c7c1159, d90fb4c and 2548530 passed Windows/Ubuntu
but did not complete macOS acceptance. Their platform-dependent unit
references were subsequently corrected. The successful run above supersedes
their failed current-platform status while retaining their historical record.

The documentation website was deployed from commit
`d389e486035d94c1aceb6acc23acf52a6a10f505`. A source push or website
deployment does not create a C++ Release or publish runtime attachments.

Original R/full-CGF comparisons apply only to matched common domains.
Partial-normal SPA adds a declared approximation; a full-CGF diagnostic
threshold is not a bound on P-value error. Finite fixtures are implementation
checks, not universal population/extreme-tail calibration proofs. No
representative biobank-scale speed advantage is claimed.

## Distribution boundaries

The C++ version is 1.2.0; DESCRIPTION retains the historical R version.
R installation does not install the standalone executable.
The candidate preserves applicable component copyright and license files.

At review, public Releases contained only historical R v1.1.0, with no
cpp-v1.2.0 C++ tag or runtime assets. Development sources can be downloaded
and compiled. A future tagged source archive will contain R plus cpp/.
Prepared local archives and CI artifacts are not already published Releases.
Match the final reviewed commit, tag, complete runtime archives and checksums.
An initial Windows/Linux pre-release should state its tested platform scope.
See [distribution instructions](GITHUB_DISTRIBUTION.md).
