# SPAGxE C++ 1.2.0

**Complete official teaching manual:** [download, installation, every current workflow, outputs and troubleshooting](../docs/USER_GUIDE.md).

Run the extended feature tutorial after the six introductory examples: `examples/tutorials/07_feature_walkthrough.sh` or `examples/tutorials/07_feature_walkthrough.cmd` from `cpp/`. See the guide for executable selection and fresh output directories.

One standalone C++17 implementation of the published SPAGxE method family,
with documented guarded partial-normal CGF acceleration.
Maintainer: Yuzhuo Ma, yuzhuoma@amss.ac.cn.

## Official C++ distribution and support

This repository is the sole official standalone C++ distribution of the
SPAGxE method family designated and maintained by **Yuzhuo Ma**, the method
designer and author of the original SPAGxE software.
For questions, bug reports, installation help or result interpretation,
contact **[yuzhuoma@amss.ac.cn](mailto:yuzhuoma@amss.ac.cn)**.

The historical R package remains available for reproducibility and compatibility.
Published-paper coauthorship and incorporated-source authorship and licenses
remain as recorded in CITATION_CPP.cff, OWNERSHIP.md and THIRD_PARTY_NOTICES.md.

SPAGxEmixPlus is planned for a future release and is not included in the current C++ distribution.

Prepare known-H ancestry-specific input with the
[local input guide](../docs/LOCAL_INPUT_PREPARATION.md). The native packer writes
integer count triplets; local BGEN reading and optional simulated BGEN generation
have distinct scopes explained there.

## Build and run

From the repository root, enter cpp/ in a supported build shell:

~~~sh
cd cpp
mkdir -p tmp results
export TMPDIR="$PWD/tmp"
export TMP="$TMPDIR" TEMP="$TMPDIR"
make -j4
make -j4 test
build/bin/SPAGxE --help
bash examples/tutorials/run_all.sh
~~~

Use a C++17 GCC/Clang compiler, C compiler, GNU Make and Bash/coreutils.
Windows source builds need MinGW/MSYS2/Rtools and produce build/bin/SPAGxE.exe.
The tested Windows binary requires Haswell/AVX2. Bundled dependency sources
are included; make does not install R or download libraries.
Optional developer/reference tests may use Python or R.

For the precompiled Windows bundle, extract the complete release asset and run
the root run_examples.cmd. The analysis and tutorials need no R/Python/compiler.
That root entry belongs to the binary bundle; the source repository uses the
commands above. See [installation](../docs/INSTALLATION.md).

## Four analysis commands

Run these inside cpp/. Create results/ first. Example-specific relaxed QC is
for small fixtures; choose real-data QC in advance.

~~~sh
# SPAGxE: shared AF, unrelated subjects. Include E in the null-model covariates.
build/bin/SPAGxE --method SPAGxE --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E --out results/base

# SPAGxE+: sparse GRM correction. Keep the companion .grm.id file.
build/bin/SPAGxE --method SPAGxE+ --pfile examples/1kg \
  --pheno examples/1kg.pheno --pheno-name Quantitative \
  --covar-name MALE,PC1,PC2,PC3,PC4 --envir-name MALE \
  --sp-grm-plink2 examples/1kg.grm.sp --out results/plus

# SPAGxEmix: BGEN and individual AF fitted from PCs.
build/bin/SPAGxE --method SPAGxEmix \
  --bgen examples/synthetic/synthetic.bgen ref-first \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E --pc-cols PC1,PC2 --out results/mix

# SPAGxEmix_local: ancestry-specific BGEN expected G conditional on known H.
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-bgen-manifest examples/local/bgen_demo/bits16_zlib.localbgen.tsv \
  --pheno examples/local/samples.pheno --pheno-name Quantitative \
  --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 --out results/local
~~~

[Annotated SH/CMD tutorials](examples/tutorials/README.md) cover each method,
four trait families, global/local CCT and native binary local packing.
Set SPAGXE_BIN (SH) or SPAGXE_EXE (CMD) to use another executable.
Use a fresh output directory when repeating packing: existing count triplets
are intentionally not overwritten.

## Inputs and interpretation

Use --pheno-name for native fitting or --resid-name for appropriately fitted
residuals. Residual-only input has no prospective Wald leg.
Native interaction tests use additive allele dosage/counts. Dominant/recessive
recoding and arbitrary survival formula extensions are not native interfaces.
See [feature coverage](../docs/FEATURE_COVERAGE.md).
Include every environment in --covar-name; include SPAGxEmix PCs there too.
Samples are matched by unique IID, with trait-specific missing-data masks.

Global formats: BED, PGEN, BGEN with explicit ref-first/ref-last/ref-unknown,
VCF/VCF.gz and BCF. The --grm-file option reads a headered ID1 ID2 VALUE
table; --sp-grm-plink2 reads .grm.sp with .grm.id.
Local formats: paired G/H text, packed counts, known-H local BGEN.
Unknown local H uses the documented complete-H mask; missing G is Hq;
minor-direction flipping is H-G. Uncertain-H inference is not implemented.

The published Branch B uses E*R0 after removing intercept and G from R.
Beyond the normal/SPA threshold (default 2), all four methods use guarded
partial-normal SPA. Fixed type-7 Tukey fences at 1.5 IQR select exact binomial
tail terms; central terms contribute their exact mean and variance. Candidate
roots are checked against the full CGF at diagnostic tolerance 1e-3. Failed
diagnostics, degenerate partitions and finite-support boundary cases use full
evaluation automatically. This tolerance is not a P-value error bound.
See [partial SPA](../docs/PARTIAL_SPA.md) for the policy and limitations.
SPAGxE+ uses score variance correction and remains SPA-only.
SPAGxEmix and local analysis do not accept GRM input.
See [method policy](../docs/METHOD_POLICY.md) for environment-coding sensitivity,
model assumptions and reference-comparison domains.

Output LOG10P columns store -log10(P); P=10^(-LOG10P). Separate score,
SPA, Wald, branch and status fields distinguish test components.
Score-derived BETA/SE are not the prospective interaction-regression estimates.
A selected individual-AF logistic fit must converge. Otherwise P values are
NA with status 9 (AF_FIT_FAILED); no pooled-model P is substituted.
AF exports preserve failure as status -2.
Use --help options and [CLI reference](../docs/CLI_REFERENCE.md).

## Utilities and project layout

--pack-local validates integer G/H and writes .localbed/.localbim/.localfam.
--cal-af-coef exports AF coefficients; association currently fits AF online
and does not consume that export as a cache. --int-pheno transforms traits.

Source modules:
- src/interaction: shared/mixed interaction scans and prospective Wald tests.
- src/allele_model: AF estimation and conditional-binomial CGF.
- src/local: local ancestry readers, kernels, packing and workflow.
- src/genotype, src/io, src/engine, src/util: input and numerical infrastructure.

make install PREFIX=/your/path installs the executable to PREFIX/bin.
make dist creates a complete runnable bundle in a new dist directory.
Clean-rebuild after changing compiler, alignment, SIMD or link flags.

See [validation](../docs/validation/README.md) and [release status](../docs/RELEASE_STATUS.md).
The executable uses the GPL license in LICENSE; incorporated components retain
their notices in [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).
