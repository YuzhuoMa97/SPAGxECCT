# SPAGxE 1.2.0 release scope

This standalone application implements SPAGxE, SPAGxE+, SPAGxEmix and
SPAGxEmix_local using one documented statistical policy: paper-defined
weights and guarded partial-normal CGF acceleration. Original R remains a
historical interface with its own version/dependencies, not a C++ runtime.

Included: global BED/PGEN/BGEN/VCF/BCF; native linear/logistic/Breslow Cox/
ordinal null models; homogeneous sparse-GRM correction; individual AF;
known-H local text/count/BGEN; C++ hardcall packing; ID/allele/QC checks;
explicit CCT weights; annotated tutorials and build instructions.

Excluded: unpublished individual-AF + cross-individual-GRM GxE, ancestry
inference, arbitrary uncertain-H/joint posterior tests, or direct support for
every other program's local binary format.

SPAGxEmixPlus is planned for a future release and is not included in the current C++ distribution.

The homogeneous GRM method corrects the score variance. Native null fitting
provides the listed GLM/Cox/ordinal models; it is not a general SAIGE/GLMM null
fitter. A suitable externally fitted null residual can be supplied explicitly.

Development validation proceeds on public `develop`, followed by a reviewed
merge to `main` and a fixed release tag. See the
[development and release procedure](GITHUB_DISTRIBUTION.md).

## Validation interpretation

Acceptance checks include Git source transport/fresh build, numerical units,
CLI/scope guards, global/local readers, AF/CCT, independent formula/R references,
same-kernel numerical regression and native four-method tutorials.
Developer tests under cpp/tests regenerate their own summaries.

Windows x86-64 and Ubuntu 24.04.5 x86-64 under WSL2 have local execution
evidence for the repaired candidate described below. The default x86 binary
needs Haswell/AVX2. macOS has not passed validation. Check actual Actions for
the released commit; the original development run failed and a local repair
does not change that run. An initial pre-release remains appropriate.

Original R/full-CGF baseline agreement is limited to stated common domains,
not every historical default/defective edge case. The accelerated default can
change P values; independent partial-policy validation and full-reference delta
comparisons serve different purposes. R preservation was checked with compatible
installed dependencies, not all fresh-machine dependency routes.

Partial-normal SPA adds a declared approximation; full-CGF diagnostics at
1e-3 do not bound P-value error. Full-CGF SPA is itself an approximation. Passing examples does not prove universal
population/GWAS-tail calibration or a measured biobank-scale speed advantage.
See [method policy](METHOD_POLICY.md) and [equations](METHOD_EQUATIONS.md).

The tagged source archive contains preserved R plus cpp/. The C++ version is
1.2.0; root DESCRIPTION remains the historical R version. Binary assets require
a separate maintainer upload; see [distribution](GITHUB_DISTRIBUTION.md).


## Official release channel

This is the standalone C++ distribution designated and maintained as official
by Yuzhuo Ma, the method designer and original software author.
Contact: yuzhuoma@amss.ac.cn. The root README gives the official installation
and example entry points. Author designation does not change the measured
platform, comparison-domain or calibration limits above.
See [implementation choices](IMPLEMENTATION_CHOICE.md).


## Model and distribution boundaries

The four methods use additive allele dosage/counts. Native survival input is
one TIME:EVENT pair per subject, with Breslow ties. Dominant/recessive genotype
recoding, stratified or counting-process survival formulas and arbitrary R
model objects are not native interfaces. Supplying residuals does not establish
that every external model is valid for this test or add a prospective Wald leg.
See [feature coverage](FEATURE_COVERAGE.md) for the complete public scope.

CI supplies build/test jobs and development artifacts. It does not create a
GitHub Release or upload the complete prepared Windows bundle automatically.
The maintainer uploads the verified bundle and its matching checksum separately.


## Development verification, 31 October 2026

The public source was downloaded from commit
`ce9c2473649e0be2a8b70dd2f6186a03b25d21aa`. Its
[Actions run](https://github.com/YuzhuoMa97/SPAGxECCT/actions/runs/41034171345)
failed on Ubuntu, Windows and macOS. These results must not be relabeled as
successful after local work.

Native Ubuntu testing reproduced a CRLF BIM-reader defect: the final carriage
return entered REF and contaminated tabular output. The repaired reader
removes that terminator. The integration test now explicitly compares LF,
CRLF and no-final-newline BIM files across four trait models. A separate
Windows Makefile repair keeps library lookup static when a full `-static`
build is requested, instead of switching implicit runtime libraries back to
dynamic lookup. This explains a possible hosted Windows dependency failure;
the hosted rerun is still required. The statistical equations are unchanged.

The repaired local candidate passed:

- Windows numerical units, 17 regression suites, the six bundled tutorials
  and the 33-group complete walkthrough. Native runs used a minimal PATH.
  PE imports contain Windows system libraries, without compiler-runtime DLLs.
- Ubuntu 24.04.5/GCC 13.3 numerical units, 14 CI regression suites, 51
  integration checks over 52,416 rows, tutorials, distribution extraction and
  PREFIX installation: 22 recorded stages in total. This was a Linux ELF
  running under WSL2, not Windows Bash.
- Source and extracted bundles reproduced 45 tutorial result files. Across
  Windows and Linux, 43 association/CCT tables with 139,709 numeric fields
  agreed within about 1.42e-14 absolute difference. This is portability and
  reproducibility evidence, not an independent calibration proof.

Apply the repaired source to a development branch, rerun the actual hosted
checks, and identify the final accepted commit before making a Release.
Ubuntu-built dynamic binaries require compatible glibc/libstdc++ and do not
promise compatibility with every older Linux distribution. macOS remains
unverified; its original build failure has not been diagnosed from a complete
log. Preserved R code uses its own historical dependency chain.
