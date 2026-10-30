# Supported features and boundaries

This is the author-designated standalone C++ distribution maintained by
Yuzhuo Ma (yuzhuoma@amss.ac.cn). It provides one statistical policy with paper-defined weights and guarded
partial-normal CGF acceleration.
The historical R package remains a separate interface for reproducibility.

## Analysis methods

| Method | Supported model | Boundary |
|---|---|---|
| SPAGxE | Additive dosage/count GxE, shared AF, unrelated subjects | Branch B uses E*R0; report the chosen E coding |
| SPAGxE+ | Shared AF with sparse-GRM score correction | SPA-only; individual AF plus GRM is not provided |
| SPAGxEmix | PC-based individual AF, unrelated subjects | Rejects GRM; failed AF fits produce NA/status 9 |
| SPAGxEmix_local | Known ancestry copy counts H and ancestry-specific G | Unrelated subjects; no arbitrary uncertain-H integration |

The original-paper statistical definition is the baseline, with documented
input and numerical repairs and the declared [partial SPA](PARTIAL_SPA.md)
acceleration. This does not preserve every alternative software
approximation, selector or original R extension. Native tests use additive
allele dosage/counts; dominant/recessive recoding is not a native test mode.

Native trait models are linear, logistic, Breslow-ties Cox and ordered
cumulative logit. Cox input is a TIME:EVENT pair per subject; stratified,
counting-process or arbitrary R survival formulas are not native interfaces.
An external residual table must come from a suitable model; accepting a table
does not validate every possible external model and supplies no Wald refit.

## Inputs and workflow

| Feature | Available scope |
|---|---|
| Global genotype input | PLINK BED/BIM/FAM, PGEN/PVAR/PSAM, BGEN, VCF/VCF.gz and BCF |
| Local input | Paired G/H text, packed count triplets and known-H local BGEN |
| Local packing | C++ validation and packing of integer hardcall G/H; no rounding of uncertain genotypes |
| Sample alignment | Unique IID matching and trait-specific missingness masks |
| Association workflow | Multiple phenotype/environment columns, chunked input, marker threads, QC and selection |
| GRM input | Headered ID1/ID2/VALUE table or PLINK2 sparse GRM and companion IDs |
| Results | Score/SPA/available Wald/final P, branch, status and settings receipts |
| Combination | Global/local CCT with explicit group or per-component weights |
| Utilities | AF coefficient export and inverse-normal phenotype transform |

AF exports are not an association-cache input. Local input does not infer
ancestry, import every other local binary layout, or fit local plus GRM.
Ordinary BIM/FAM sidecars are reusable only after the documented identity,
ordering and allele conventions are satisfied; CRCs bind their actual bytes.
See [CLI](CLI_REFERENCE.md), [local counts](LOCAL_ANCESTRY_FORMAT.md) and
[local BGEN](LOCAL_BGEN_FORMAT.md).

SPAGxEmixPlus is planned for a future release and is not included in the current C++ distribution.
The GRM-adjusted score workflow does not imply a native general GLMM null fitter.
See [local input preparation](LOCAL_INPUT_PREPARATION.md) for native conversion
capabilities and the limits of ancestry-specific data generation.

## Distribution and evidence

Windows source builds and the prepared complete Windows bundle have local
execution evidence. The binary requires x86-64 Haswell/AVX2 and runs without
R, Python or a compiler. The source route needs the documented C/C++ toolchain,
Make and Bash; optional developer/reference checks may need Python or R.

Linux/macOS build/test jobs are supplied. Their actual results for the release
commit must be checked before claiming platform validation. CI artifacts are
development outputs and the workflow does not publish Releases automatically.
Keep the full cpp/ source tree; a source archive is not a precompiled binary.

Reference comparisons are limited to stated matched domains. Original R default
settings, revised statistical definitions, failure states and boundary repairs
can produce differences. The partial-CGF acceleration introduces a further
approximation; full-CGF diagnostics do not establish a P-value error bound. Existing tests
do not prove universal tail calibration or a measured speed advantage at biobank
scale. Required incorporated-source provenance and licenses remain in the
[notices](../THIRD_PARTY_NOTICES.md).

Start with the [four-method tutorials](../cpp/examples/tutorials/README.md),
[installation](INSTALLATION.md), [method choices](IMPLEMENTATION_CHOICE.md),
[release status](RELEASE_STATUS.md) and [validation](validation/README.md).
