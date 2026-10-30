# Official SPAGxE C++ user guide

SPAGxE 1.2.0 provides native C++ workflows for SPAGxE, SPAGxE+, SPAGxEmix
and SPAGxEmix_local. This repository is the sole official standalone C++
distribution designated and maintained by **Yuzhuo Ma**, the method designer
and author of the original SPAGxE software. Questions, installation help and
bug reports: **[yuzhuoma@amss.ac.cn](mailto:yuzhuoma@amss.ac.cn)**.

SPAGxEmixPlus is planned for a future release and is not included in the
current C++ distribution. Incorporated components retain their authorship,
licenses and notices. The original R package remains available for historical
reproducibility and compatibility; installing it does not install this C++
executable.

This guide follows a complete analysis, from choosing a download to reading
results. Its commands use shipped examples. The examples are small validation
fixtures; their permissive QC settings are not recommended cohort thresholds.

## Contents

1. [Choose and obtain the software](#choose-and-obtain-the-software)
2. [Set the executable and start the walkthrough](#set-the-executable-and-start-the-walkthrough)
3. [Choose a method and prepare input](#choose-a-method-and-prepare-input)
4. [Run the four methods](#run-the-four-methods)
5. [Change input format or phenotype workflow](#change-input-format-or-phenotype-workflow)
6. [Select samples and variants and set QC](#select-samples-and-variants-and-set-qc)
7. [Pack and analyze local ancestry input](#pack-and-analyze-local-ancestry-input)
8. [Combine global and local results](#combine-global-and-local-results)
9. [Use AF export and phenotype transformation](#use-af-export-and-phenotype-transformation)
10. [Control output and resources](#control-output-and-resources)
11. [Interpret and retain results](#interpret-and-retain-results)
12. [Resolve common problems](#resolve-common-problems)
13. [Find every command-line option](#find-every-command-line-option)

## Choose and obtain the software

| Situation | Download | Next step |
|---|---|---|
| Windows user who wants to run analysis | Complete Windows ZIP plus SHA256SUMS.txt from [Releases](https://github.com/YuzhuoMa97/SPAGxECCT/releases) | Verify, extract everything, run the EXE |
| Source user on a supported C/C++ toolchain | Complete repository source or tagged source archive | Build inside `cpp/` |
| Tester of a development update | Exact `develop` commit archive or matching complete Actions candidate | Record the commit and validate before relying on it |
| Historical R user | Original R package and its compatible dependencies | Follow the historical R instructions |

The C++ analysis executable and supplied native tutorials do not call R or
Python. A C++ source build needs a C++17 compiler, a C compiler, GNU Make and
Bash/coreutils. Optional developer/reference checks may use R or Python.

### Windows binary download

After the maintainer publishes the Release, download these **two** assets:

- `SPAGxE-1.2.0-windows-x86_64.zip`
- `SHA256SUMS.txt`

Do not confuse the automatically generated **Source code ZIP** with the
precompiled Windows ZIP. The prepared binary requires Windows x86-64 with a
Haswell/AVX2-class CPU. Keep the executable, example files, documentation and
licenses together.

In PowerShell, in a directory containing the two downloaded files:

```powershell
$archive = 'SPAGxE-1.2.0-windows-x86_64.zip'
$expected = ((Get-Content -LiteralPath SHA256SUMS.txt | Where-Object { $_ -match [regex]::Escape($archive) }) -split '\s+')[0].ToLowerInvariant()
$actual = (Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actual -ne $expected) { throw 'ZIP checksum mismatch.' }
Expand-Archive -LiteralPath $archive -DestinationPath unpacked
Set-Location unpacked\SPAGxE-1.2.0-windows-x86_64
.\cpp\SPAGxE.exe --version
.\run_examples.cmd
```

Use a fresh extraction directory. The root launcher runs six tutorials and
writes `cpp/results/tutorials`. Native packing refuses existing output
triplets; do not repeatedly run it against an already populated prefix.

### Build the complete source

Download the selected source archive and extract it, or clone the repository.
Use the approved tag once it exists. During development, select `develop`
explicitly; before the C++ update is merged, `main` can still contain only the
historical R package. See the [development/release procedure](GITHUB_DISTRIBUTION.md).

For example, after the maintainer has pushed `develop`:

```sh
git clone --branch develop https://github.com/YuzhuoMa97/SPAGxECCT.git
cd SPAGxECCT
git rev-parse HEAD
```

An approved tag checkout uses `git checkout cpp-v1.2.0` once that tag exists.
Record the displayed commit when comparing builds or reporting a problem.

In a supported POSIX build shell, from the repository root:

```sh
cd cpp
mkdir -p tmp results
export TMPDIR="$PWD/tmp"
export TMP="$TMPDIR" TEMP="$TMPDIR"
make -j4
make -j4 test
build/bin/SPAGxE --version
bash examples/tutorials/run_all.sh
```

Windows source builds use the documented MSYS2/MinGW or compatible Rtools
environment and produce `build/bin/SPAGxE.exe`. A compiler supplied by Rtools
does not create an R runtime dependency. Keep the complete `cpp/` directory,
including `third_party/`, `Makefile`, `VERSION`, examples and tests. `src/`
alone cannot build the package. `make` compiles bundled dependencies locally.

`make install PREFIX=/chosen/path` optionally installs the executable in
`PREFIX/bin`; it does not install all example data. `make dist` creates a
complete bundle in a new distribution directory. Read
[installation](INSTALLATION.md) and the [C++ build manual](../cpp/README.md)
before changing compilers, architecture flags or static linkage. Actual
Linux/macOS status must be checked for the selected commit; a Windows Bash
build is not a Linux run.

## Set the executable and start the walkthrough

Every analysis block below runs **from `cpp/`**. The POSIX-shell recipes use
two variables, `SPAGXE` and `OUT`. Select the correct executable first:

```sh
# POSIX source build:
SPAGXE=build/bin/SPAGxE
# Windows source build in a supported Bash shell: SPAGXE=build/bin/SPAGxE.exe
# Extracted Windows binary bundle, from its cpp/: SPAGXE=./SPAGxE.exe
OUT=results/official_walkthrough
```

The full runnable walkthrough is also provided as
[SH](../cpp/examples/tutorials/07_feature_walkthrough.sh) and
[CMD](../cpp/examples/tutorials/07_feature_walkthrough.cmd). Windows binary
users can use the CMD script without Bash, Python or R. The SH recipes below
use standard shell utilities only to prepare small selection lists.

Choose a fresh `OUT` directory for a second walkthrough. Keep each example's
distinct output prefix. Input fixtures under `examples/` are read without
modification. Ordinary association writers can replace files at an existing
prefix, whereas native packing explicitly refuses an existing triplet. Keep
results you need and use a new prefix for another run.

<!-- example: UG01_help -->
```sh
mkdir -p "$OUT"
"$SPAGXE" --version
"$SPAGXE" --help
"$SPAGXE" --help options
"$SPAGXE" --help phenotype
"$SPAGXE" --help output
```

`--version` should report 1.2.0. Help is printed to the console; an unknown
topic or option exits with an error. Use `--help FLAG`, for example
`--help bgen`, to inspect an individual option.

## Choose a method and prepare input

| Method | Appropriate cohort | Required method-specific input |
|---|---|---|
| SPAGxE | Unrelated subjects with shared allele frequency | Global genotype input |
| SPAGxE+ | Shared-AF cohort with relatedness | Global genotypes and sparse GRM |
| SPAGxEmix | Unrelated multi-ancestry/admixed cohort | Global genotypes and suitable ancestry PCs for individual AF |
| SPAGxEmix_local | Unrelated cohort with known ancestry copy counts | Complete ancestry-specific G/H input |

SPAGxE+ implements the ordinary native null plus GRM score-correction path.
It does not claim to be a general SAIGE/GLMM null fitter; a suitable externally
fitted residual can be supplied explicitly. SPAGxEmix and local methods reject
GRM input. They do not enable an unpublished combined extension.

All four methods cover the paper's additive allele count/dosage workflow.
The software does not silently change to dominant/recessive genotype models.

### Phenotype and covariate tables

Use a headered whitespace-delimited table. TAB is a convenient delimiter:

```text
IID    Y       AGE    E       PC1
S1     1.7     54     0.4     0.02
S2     NA      61     1.1    -0.01
```

Use actual separators between fields. `FID IID ...` is also accepted, but
matching uses unique string IID, not a joint FID/IID key. Preserve leading
zeros and ensure that each table has unique IIDs. Missing numeric values may
use `NA`, `na`, `NaN`, `nan`, `.` or `-` in ordinary phenotype input. Local
genotype files have their own documented missing-value contract.

- `--pheno-name` selects raw trait columns; explicitly naming them avoids
  accidentally treating numeric metadata as outcomes.
- `--covar-name` selects null-model covariates. They come from `--pheno`
  unless a separate `--covar` table is supplied.
- Every tested environment must enter the native null covariates.
- SPAGxEmix's AF PCs must also enter the native null covariates. The default
  PC list is PC1,PC2,PC3,PC4; examples explicitly use the PCs actually present.
- Phenotype missingness is handled per trait. Different masks can require
  different AF/null fits; do not assume an AF export from another cohort is
  an interchangeable cache.
- Choose environmental coding before inspecting P values. E is used as
  supplied, and the paper-defined Branch-B approximation can depend on its
  zero point. The fixture's `Eshift=E+3` is not an independent covariate.

Native traits are quantitative, binary, survival (`TIME:EVENT`, positive time
and 0/1 event) and ordinal ordered numeric categories. Native Cox uses Breslow
ties. Default family selection is `auto`; explicit family selectors are shown
below. Include an intercept through the native fitter; do not add a constant
duplicate column to the covariates.

The current automatic rule recognizes two distinct nonmissing values as
binary, recoding the smaller to 0 and the larger to 1. Three to ten contiguous
integer levels are recognized as ordinal and shifted to start at 0. Other
nonconstant single columns are quantitative. Explicit `--regression-model`
is checked against this inference; it does not force an incompatible encoding
into another family. Check the model/recode messages, use intended ordered
integer codes and document which binary value represents the event.

### Global genotype files

Choose exactly one source in an association command:

| Option | Files to keep together |
|---|---|
| `--bfile PREFIX` | `PREFIX.bed`, `.bim`, `.fam`; SNP-major BED |
| `--pfile PREFIX` | `PREFIX.pgen`, `.pvar`, `.psam` |
| `--bgen FILE MODE` | BGEN plus its companion `.sample` if IDs are not embedded |
| `--vcf FILE` | VCF or supported VCF.gz |
| `--bcf FILE` | BCF2 |

For BGEN, `ref-first` means the first listed allele is REF; `ref-last` means
the last is REF. `ref-unknown` uses a provisional last REF and warns. Establish
the actual allele convention from the data producer; do not choose a mode to
obtain a desired sign or P value. The examples use known orientation.

## Run the four methods

### SPAGxE: quantitative trait and saved null residual

<!-- example: UG02_base -->
```sh
"$SPAGXE" --method SPAGxE --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E --regression-model linear \
  --save-resid --mac 0 --maf 0 --geno 1 --out "$OUT/base"
```

This fits a quantitative null in C++, tests GxE, and writes
`base.Quantitative.SPAGxE`, `base.analysis.settings.tsv` and
`base.null.resid`. Residual columns retain the fitted phenotype names.
The simulated permissive QC exposes edge cases; choose scientific thresholds
for real data before the scan.

### SPAGxE+: sparse relationship table

<!-- example: UG03_plus_table -->
```sh
"$SPAGXE" --method SPAGxE+ --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --grm-file examples/synthetic/synthetic.grm.tsv \
  --mac 0 --maf 0 --geno 1 --out "$OUT/plus_table"
```

The table header is `ID1 ID2 VALUE`. VALUE is relationship, or twice kinship;
self relationship is 1+F. Do not pass kinship values while describing them as
relationship. The same method is selected by supplying a GRM with SPAGxE.
Related-sample output uses the SPA/normal score leg without an IID Wald leg;
its file suffix is `.SPAGxE`, with the plus method identified by the command.

### SPAGxEmix: individual AF from ancestry PCs

<!-- example: UG04_mix -->
```sh
"$SPAGXE" --method SPAGxEmix \
  --bgen examples/synthetic/synthetic.bgen ref-first \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E --pc-cols PC1,PC2 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/mix"
```

PC1/PC2 enter both individual-AF fitting and the phenotype null. Their order
and units matter. The minor-oriented pooled/OLS/carrier-logistic cascade is
estimated online. Invalid or nonconverged selected AF fits produce explicit
NA/status 9, rather than a pooled replacement P. An `.afc` export is not an
association-cache argument.

### SPAGxEmix_local: ancestry-specific G/H

<!-- example: UG05_local_text -->
```sh
"$SPAGXE" --method SPAGxEmix_local \
  --local-manifest examples/local/cohort.local.tsv \
  --pheno examples/local/samples.pheno \
  --pheno-name Quantitative,Binary,Time:Event,Ordinal \
  --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 \
  --out "$OUT/local_text"
```

This writes each trait's `.SPAGxE_local` ancestry-specific table,
`.SPAGxE_local_CCT` across-ancestry table and a `.local.settings.tsv` receipt.
H must be known ancestry chromosome copies 0/1/2. Supply the complete ancestry
collection, with sum(H)=2 per complete person/variant. The program does not
infer ancestry. See [input preparation](LOCAL_INPUT_PREPARATION.md).

## Change input format or phenotype workflow

### BGEN with last-allele REF, VCF and BCF

These fixtures contain equivalent global genotypes with known allele metadata.

<!-- example: UG06_bgen_ref_last -->
```sh
"$SPAGXE" --bgen examples/synthetic/synthetic.reflast.bgen ref-last \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/bgen_ref_last"
```

A BGEN without embedded IDs automatically uses its companion `.sample` file.
Keep `synthetic.noids.sample` next to this fixture:

<!-- example: UG33_bgen_external_ids -->
```sh
"$SPAGXE" --bgen examples/synthetic/synthetic.noids.bgen ref-first \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/bgen_external_ids"
```

<!-- example: UG07_vcf -->
```sh
"$SPAGXE" --vcf examples/synthetic/synthetic.vcf \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/vcf"
```

<!-- example: UG08_bcf -->
```sh
"$SPAGXE" --bcf examples/synthetic/synthetic.bcf \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/bcf"
```

### PGEN and PLINK2 sparse GRM

The bundled `examples/1kg` files have their own IIDs, phenotypes and PCs.
Do not mix them with the synthetic cohort's phenotype or GRM files.

<!-- example: UG09_pgen -->
```sh
"$SPAGXE" --pfile examples/1kg --pheno examples/1kg.pheno \
  --pheno-name Quantitative --covar-name MALE,PC1,PC2,PC3,PC4 \
  --envir-name MALE --out "$OUT/pgen"
```

<!-- example: UG10_plus_plink2 -->
```sh
"$SPAGXE" --method SPAGxE+ --pfile examples/1kg \
  --pheno examples/1kg.pheno --pheno-name Quantitative \
  --covar-name MALE,PC1,PC2,PC3,PC4 --envir-name MALE \
  --sp-grm-plink2 examples/1kg.grm.sp --out "$OUT/plus_plink2"
```

Keep `examples/1kg.grm.id` beside `.grm.sp`; it maps numeric indices to IIDs.
Supply one GRM format, not both. Study-specific GRM construction is an upstream
step; verify its relationship scale and subject identity before analysis.

### Explicit binary, survival and ordinal models

<!-- example: UG11_trait_models -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Binary \
  --regression-model logistic --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/binary"
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Time:Event \
  --regression-model cox --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/survival"
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Ordinal \
  --regression-model ordinal --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/ordinal"
```

Use Breslow ties when matching an external Cox reference. An external R Cox
default using Efron need not reproduce the same residuals/P values.

### Several traits and environments in one scan

<!-- example: UG12_multiple -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno \
  --pheno-name Quantitative,Binary,Time:Event,Ordinal --regression-model auto \
  --covar-name X1,E,PC1,PC2 --envir-name E,X1 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/multiple"
```

Here E and X1 are two toy tested environments, and both enter the null.
One output file is produced per trait, with blocks for both environments.
Trait-specific missingness masks are supported. Do not use E and Eshift as
independent null covariates: Eshift is E+3 and would duplicate information
already represented by E and the intercept.

### Supplied residuals and reloading saved residuals

<!-- example: UG13_external_residuals -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno \
  --resid-name R_quant,R_binary,R_survival,R_ordinal --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/external_residuals"
```

<!-- example: UG14_reload_residual -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno "$OUT/base.null.resid" --resid-name Quantitative \
  --covar examples/synthetic/synthetic.pheno --covar-name E --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/reloaded_residual"
```

Residual input does not refit a null or add a prospective Branch-B Wald leg.
Even when the residual and score leg agree, its final P can therefore differ
from raw-trait Branch-B CCT. The user must supply suitable model residuals and
consistent subject/covariate/environment definitions. `--save-resid` belongs
to raw-trait fitting, not residual mode.

## Select samples and variants and set QC

Selection lists use the **data's IDs**, not row numbers. This shell preparation
keeps the first 400 synthetic IIDs, removes two of those, selects four variant
IDs, and then excludes one selected variant. It writes only inside OUT.

<!-- example: UG15_prepare_filters -->
```sh
awk 'NR<=401 {print $1}' examples/synthetic/synthetic.pheno > "$OUT/keep.txt"
printf 'IID\nS0001\nS0002\n' > "$OUT/remove.txt"
printf 'rsOracle01\nrsOracle02\nrsOracle03\nrsOracle04\n' > "$OUT/extract.txt"
printf 'rsOracle02\n' > "$OUT/exclude.txt"
```

<!-- example: UG16_filter_qc -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --keep "$OUT/keep.txt" --remove "$OUT/remove.txt" \
  --extract "$OUT/extract.txt" --exclude "$OUT/exclude.txt" --chr 1 \
  --geno 0.1 --maf 0.01 --mac 10 --hwe 0.000001 \
  --hard-call-threshold 0.1 --out "$OUT/filtered"
```

Do not interpret these thresholds as optimized for every study. Defaults are
geno .1, maf 1e-5, mac 10 and hwe 0 (disabled). `--hard-call-threshold` affects
HWE counts only, not replacement of analysis dosage by hardcalls. Selection
and QC are different: selected markers can fail QC and have unavailable tests.
Check the result/status information rather than expecting a P for every input.

Global QC MAC uses observed nonmissing minor dosage. Individual-AF model
selection has its distinct fitting-cohort MAC gate. Local MAC uses estimated
minor AF times known-H copies after masking, including copies whose missing
G will be mean-imputed. HWE is a global diploid hardcall diagnostic, not an
ancestry-specific exact variable-H test. Non-point posterior expectations
are not made into hardcalls merely because their mean is integer.

## Pack and analyze local ancestry input

### Understand and prepare the local schemas

In paired-text input, G and H use the same five metadata fields followed by
string IIDs. Each ancestry has one G stream and one H stream. For example:

```text
CHROM  POS  ID   REF  ALT  IID1  IID2
1      101  rs1  C    A    1     0
```

Displayed columns must be separated by actual TAB characters. The paired-text
manifest has four columns:

```text
ancestry   genotype         hapcount         counted_allele
Ancestry1  Ancestry1.G.tsv   Ancestry1.H.tsv   ALT
Ancestry2  Ancestry2.G.tsv   Ancestry2.H.tsv   ALT
```

Resolve paths relative to the manifest directory. Every stream must have the
same IID set and matching variant metadata/order; IID order is aligned, but
the program does not silently intersect variant catalogs. G can be fractional
when read directly, H must be integer 0/1/2. Known H=0 is structural G=0. Unknown
H in any ancestry excludes the subject from all ancestries at that marker;
missing G with known H is filled with H*q. Never replace unknown H with a
fabricated integer. For haplotype-resolved ALT indicators a1,a2 and labels L1,L2,
prepare `H_k=I(L1=k)+I(L2=k)` and `G_k=I(L1=k)*a1+I(L2=k)*a2` upstream.

The local BGEN manifest has three columns instead:

```text
ancestry   bgen             counted_allele
Ancestry1  Ancestry1.bgen   ALT
Ancestry2  Ancestry2.bgen   ALT
```

The packed triplet is a dedicated two-bit G/H container with zstd blocks and
CRC-bound sidecars, not a standard BED stream. Its headerless `.localbim`
columns are `CHROM ID CM POS ALT REF`; `.localfam` uses
`FID IID PAT MAT SEX PHENO`, with IID as the matching key. PHENO in that
sidecar is not the analysis outcome. Original BIM/FAM can be reused only after
the packer verifies orientation, order and matching IIDs. This preserves
compatible metadata without assuming arbitrary PLINK A1/A2 means ALT/REF.
Full byte layout and missing-value contracts are in
[LOCAL_ANCESTRY_FORMAT.md](LOCAL_ANCESTRY_FORMAT.md).

### Read existing packed counts and known-H BGEN

<!-- example: UG17_local_existing_counts -->
```sh
"$SPAGXE" --method SPAGxEmix_local \
  --local-prefix examples/local/binary_demo/local \
  --pheno examples/local/binary_demo/samples.pheno --resid-name R \
  --envir-name E --mac 0 --maf 0 --geno 1 --out "$OUT/local_existing"
```

<!-- example: UG18_local_bgen -->
```sh
"$SPAGXE" --method SPAGxEmix_local \
  --local-bgen-manifest examples/local/bgen_demo/bits16_zlib.localbgen.tsv \
  --pheno examples/local/bgen_demo/samples.pheno --resid-name R \
  --envir-name E --mac 0 --maf 0 --geno 1 --out "$OUT/local_bgen"
```

Local BGEN stores known H in ploidy and uses conditional posterior mean G.
It does not integrate arbitrary joint ancestry/genotype uncertainty into a
new exact SPA law. The supported profile is layout2, unphased biallelic,
embedded IID, H0/1/2, raw/zlib/zstd blocks and precision1..32. See the
[BGEN contract](LOCAL_BGEN_FORMAT.md).

### Generate hardcall triplets in C++

<!-- example: UG19_pack_local -->
```sh
"$SPAGXE" --pack-local \
  --local-manifest examples/local/binary_demo/cohort.local.tsv \
  --local-bim examples/local/binary_demo/local.localbim \
  --local-fam examples/local/binary_demo/local.localfam \
  --local-block-variants 2 --local-zstd-level 5 \
  --local-max-block-bytes 268435456 --out "$OUT/local_counts"
```

This creates `local_counts.localbed/.localbim/.localfam`. The optional supplied
sidecars are checked for complete IDs, variant order and ALT/REF direction;
their exact bytes are retained and CRC-bound to the binary. Ordinary BIM
A1/A2 is not universally ALT/REF. Do not rename a global BED as localbed or
edit sidecars after packing.

<!-- example: UG20_local_new_counts -->
```sh
"$SPAGXE" --method SPAGxEmix_local --local-prefix "$OUT/local_counts" \
  --pheno examples/local/binary_demo/samples.pheno --resid-name R \
  --envir-name E --mac 0 --maf 0 --geno 1 --out "$OUT/local_new"
```

<!-- example: UG21_repack_local -->
```sh
"$SPAGXE" --pack-local --local-prefix "$OUT/local_counts" \
  --local-block-variants 1 --local-zstd-level 0 \
  --local-max-block-bytes 268435456 --out "$OUT/local_repacked"
```

`--local-block-variants` is 1..65536. `--local-zstd-level` is -5..22, with
0 selecting zstd's default. The decoded-block budget defaults to 256 MiB and
is capped at 1 GiB; it is not a total-process memory limit. Reduce the requested
block size when it cannot fit. Packing has no association options or P values.

Packing accepts integer hardcalls only. It can also pack a known-H BGEN
collection **if every observed posterior is a point mass**; the general
posterior demo is not such a collection. Fractional G and uncertain posteriors
are rejected, without rounding. The program reads local BGEN but does not
provide a general native BGEN writer. The included optional Python fixture
generators simulate examples; they are not required to analyze the shipped data.

For real data, begin with resolved ancestry-specific G/H. Global G and H alone
usually do not determine each ancestry's G. The [preparation guide](LOCAL_INPUT_PREPARATION.md)
contains a complete TAB schema and explains upstream allele/ancestry assignment,
common unknown-H masks, H*q missing-G filling, H-G flips and structural H0.

## Combine global and local results

CCT combines **compatible result tables**, not arbitrary files from different
cohorts or marker catalogs. The following global VCF is matched to UG05's
local data, unlike the separate synthetic or 1kg examples.

<!-- example: UG22_cct_global_input -->
```sh
"$SPAGXE" --vcf examples/local/global.vcf \
  --pheno examples/local/samples.pheno --pheno-name Quantitative \
  --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 \
  --out "$OUT/cct_global"
```

<!-- example: UG23_cct_default -->
```sh
"$SPAGXE" --method SPAGxEcct \
  --cct-global "$OUT/cct_global.Quantitative.SPAGxE" \
  --cct-local "$OUT/local_text.Quantitative.SPAGxE_local_CCT" \
  --cct-phenotype Quantitative --envir-name E --out "$OUT/cct_default"
```

<!-- example: UG24_cct_weights -->
```sh
"$SPAGXE" --method SPAGxEcct \
  --cct-global "$OUT/cct_global.Quantitative.SPAGxE" \
  --cct-local "$OUT/local_text.Quantitative.SPAGxE_local_CCT" \
  --cct-phenotype Quantitative --envir-name E --cct-weights 0.6,0.4 \
  --out "$OUT/cct_weighted"
```

<!-- example: UG25_cct_components -->
```sh
"$SPAGXE" --method SPAGxEcct \
  --cct-global "$OUT/cct_global.Quantitative.SPAGxE" \
  --cct-local "$OUT/local_text.Quantitative.SPAGxE_local_CCT" \
  --cct-phenotype Quantitative --envir-name E --cct-per-ancestry-equal \
  --out "$OUT/cct_equal_components"
```

Default group weights are .5 global and .5 local. Explicit weights are finite,
positive and normalized. Per-ancestry equal weighting gives one nominal weight
to global and one to each valid ancestry for that row/environment; it is
mutually exclusive with `--cct-weights`. Weights change the scientific test;
choose them in advance. Variant, allele, phenotype and environment identities
must match. The global scan here is a matched toy illustration, not a substitute
for appropriate ancestry-PC adjustment in a real admixed cohort.

## Use AF export and phenotype transformation

### Export individual-AF model coefficients

<!-- example: UG26_af_export -->
```sh
"$SPAGXE" --cal-af-coef --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pc-cols PC1,PC2 \
  --mac 0 --maf 0 --geno 1 --hwe 0 --threads 2 \
  --compression gz --out "$OUT/af_models"
```

The output is `af_models.afc.gz`. Do not add `--method`, raw/residual phenotype
selection, environment, GRM or `--compression-level` to this utility. It uses
the codec's default level. The model rows have status plus intercept/PC betas;
they do not independently identify variants or samples. Retain selected marker
order, REF/ALT, sample set, PC labels and units. Status -2 is failed AF fit,
-1 is unfitted/QC-filtered, 0 pooled, 1 OLS, 2 logistic, 3 complemented logistic.
Skip every `#` comment line when reading the text. Association still estimates
AF online; it has no AFC-cache input. See [the export schema](AF_MODEL_FORMAT.md).

### Inverse-normal transform selected phenotype columns

<!-- example: UG27_int -->
```sh
"$SPAGXE" --int-pheno --pheno examples/synthetic/synthetic.pheno \
  --pheno-name Quantitative --out "$OUT/int_quantitative"
```

This writes `int_quantitative.txt`, preserving identifiers and missing values
and outputting the selected transformed columns. With no `--pheno-name`, all
non-ID columns are transformed; explicitly select outcomes to avoid changing
covariates, events or metadata. Do not add genotype, environment, model or
compression options. A transform is a scientific preprocessing choice, not a
mandatory step for every SPAGxE analysis.

The transformed table does not contain the original covariates. Supply them
as a separate covariate table, aligned by IID:

<!-- example: UG28_separate_covariates -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno "$OUT/int_quantitative.txt" --pheno-name Quantitative \
  --covar examples/synthetic/synthetic.pheno --covar-name X1,E,PC1,PC2 \
  --envir-name E --mac 0 --maf 0 --geno 1 --out "$OUT/int_analysis"
```

If `--covar` is given without `--covar-name`, all its non-ID columns enter the
null. Use an actual covariate-only file or explicitly select columns; do not
accidentally include outcomes in that default.

## Control output and resources

### Gzip and zstd output; compressed residual input

<!-- example: UG29_gzip -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E --save-resid \
  --compression gz --compression-level 6 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/gzip"
```

<!-- example: UG30_zstd -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --compression zst --compression-level 3 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/zstd"
```

<!-- example: UG31_compressed_residual -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno "$OUT/gzip.null.resid.gz" --resid-name Quantitative \
  --covar examples/synthetic/synthetic.pheno --covar-name E --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/compressed_residual"
```

Association/residual gzip levels are 1..9 and zstd levels 1..22. Omit
`--compression` for plain text. Readers support their documented plain/gzip/zstd
text paths. Use a codec-aware viewer; compressed output is not an empty table.

### Threads, chunks and statistical thresholds

<!-- example: UG32_resources -->
```sh
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --threads 4 --chunk-ksnp 0.5 --seed 20261029 \
  --spagxe-marginal-cutoff 0.001 --spa-z-threshold 2 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/resources"
```

`--threads` is positive, default 1. `--chunk-ksnp` is a positive multiple of .5;
each unit means 1024 markers, so .5 is 512. Larger chunks can use more memory;
a tiny fixture is not a speed benchmark. `--seed` is a compatibility setting:
association calculations are deterministic and do not consume it.

The marginal P cutoff defaults .001, with P<=cutoff entering Branch B.
The normal/SPA |Z| cutoff defaults 2, minimum .01. The command shows defaults
explicitly; do not tune thresholds after examining association significance.
All four methods use one guarded partial-normal policy, with fixed 1.5-IQR
partitioning and full-CGF candidate-root diagnostics. Automatic full evaluation
is protection within that policy. It is not a user-selectable second method.
The 1e-3 diagnostic threshold is not a P-value error bound. See
[partial SPA](PARTIAL_SPA.md) and [method equations](METHOD_EQUATIONS.md).

## Interpret and retain results

| Output | Meaning |
|---|---|
| `PREFIX.TRAIT.SPAGxE` | Global SPAGxE or SPAGxE+ result |
| `PREFIX.TRAIT.SPAGxEmix` | Global individual-AF result |
| `PREFIX.TRAIT.SPAGxE_local` | Ancestry-specific local rows |
| `PREFIX.TRAIT.SPAGxE_local_CCT` | Across-ancestry local combination |
| `PREFIX.TRAIT.SPAGxE_global_local_CCT` | Separate global/local combination utility output |
| `.analysis.settings.tsv` / `.local.settings.tsv` | Version, policy and local processing receipt |
| `.null.resid` with optional codec suffix | Requested fitted null residuals |
| `.afc` with optional codec suffix | AF coefficient export |
| `.txt` from INT | Selected inverse-normalized outcomes |

`LOG10P` means **-log10(P)**, not P and not log10(P). Recover P as
`10^(-LOG10P)` when representable; retain log values for extreme tails.
Each environment has normal/SPA, available Wald and final GxE fields, signed
statistics, branch and status. Working-score BETA/SE are not the prospective
Wald coefficient/SE. Branch1=A and branch2=B. Unavailable Wald is NA, and GRM
and residual-mode final P use the score leg.

Read `SPA_STATUS` together with P:

| Code | Interpretation |
|---|---|
| 0 | Successful SPA, including valid support/endpoint handling |
| 1 | Normal branch |
| 2 | Successful near-zero saddlepoint reduction |
| 3..6 | Normal fallback after a numerical SPA failure; inspect the case |
| 7..8 | Unavailable/postprocessing/no-test status |
| 9 | Invalid or nonconverged individual-AF fit; P is NA |

Full-CGF protection of a partial candidate is different from normal fallback
after numerical failure. A finite CCT P does not erase a failed component's
status. QC and missingness fields, local H masks and tested-allele flags must
also be interpreted. See [the CLI reference](CLI_REFERENCE.md) for exact columns.

Keep the source tag or commit, executable SHA256, command line, input/sample/
variant provenance, AF/GRM conventions, settings receipts and result checksums.
Changing the phenotype null, Cox ties, AF routing, E coding, CCT weights or
endpoint policy can change P beyond the partial approximation alone. Matching
an old R result requires matched definitions and samples, not merely the same
method name. Numerical acceptance is not proof of all-population tail calibration.

## Resolve common problems

| Symptom | Action |
|---|---|
| EXE or example file missing | Extract/download the complete correct package; check current directory and executable path |
| Only R files appear after download | Select the C++ branch/commit or approved C++ tag; historical R tags do not contain the new executable |
| Output cannot be opened | Create the parent directory and select a fresh prefix; retain existing output rather than forcing overwrite |
| Environment/PC missing from null | Add each tested environment and each mix AF PC to `--covar-name` |
| Rank-deficient design | Remove duplicated/constant covariates; E and Eshift are collinear in the supplied fixture |
| No valid P for a marker | Inspect QC, usable sample/copy counts, variance and SPA_STATUS; mono/all-missing data may have no test |
| AF_FIT_FAILED | Inspect genotype/PC design and fit diagnostics; do not interpret NA as P=1 or silently substitute pooled AF |
| Mix/local plus GRM rejected | Choose a supported method/cohort; no combined extension is available |
| BGEN allele error or sign mismatch | Verify the producer's REF order and companion IDs; local BGEN has its separate known-H profile |
| Local CRC/metadata error | Use all three matching sidecars and repack validated inputs; do not edit or rename an unrelated BED |
| Local hardcall packing rejects input | Fractional dosage or non-point posteriors cannot be silently rounded; use supported direct input |
| Local block budget exceeded | Reduce block variants or explicitly choose a supported budget; this is not a total-process-memory guarantee |
| INT result lacks covariates | Its output contains identifiers and selected traits; supply a separate covariate table |
| A combined P looks finite despite a warning | Inspect component SPA_STATUS and unavailable fields before interpreting it |

For a bug report, email the version/commit, platform, exact command, console
error and a minimal **nonprivate** reproducing example. Do not send confidential
cohort data or private study materials. Dataset-specific model validity remains
part of study design.

## Find every command-line option

The executable is the definitive interface: `--help options` lists canonical
flags and `--help TOPIC` explains methods, genotype, phenotype, sparse GRM,
packing and output. Legacy compatible aliases do not create new algorithms.

| Function group | Options | Worked examples |
|---|---|---|
| Identify/help | `--version`, `--help` | UG01 |
| Select method | `--method` | UG02..UG05, UG23..UG25 |
| Global readers | `--bfile`, `--pfile`, `--bgen`, `--vcf`, `--bcf` | UG02, UG04, UG06..UG10, UG33 |
| Raw phenotype/null | `--pheno`, `--pheno-name`, `--covar`, `--covar-name`, `--envir-name`, `--regression-model` | UG02, UG11, UG12, UG28 |
| Residual workflow | `--resid-name`, `--save-resid` | UG02, UG13, UG14, UG29, UG31 |
| Individual AF | `--pc-cols`, `--cal-af-coef` | UG04, UG26 |
| Sparse relatedness | `--grm-file`, `--sp-grm-plink2` | UG03, UG10 |
| Local readers | `--local-manifest`, `--local-prefix`, `--local-bgen-manifest`, `--local-max-block-bytes` | UG05, UG17..UG21 |
| Local conversion | `--pack-local`, `--local-bim`, `--local-fam`, `--local-block-variants`, `--local-zstd-level` | UG19, UG21 |
| Sample/marker selection | `--keep`, `--remove`, `--extract`, `--exclude`, `--chr` | UG15, UG16 |
| QC | `--geno`, `--maf`, `--mac`, `--hwe`, `--hard-call-threshold` | UG16 |
| Global/local CCT | `--cct-global`, `--cct-local`, `--cct-phenotype`, `--cct-weights`, `--cct-per-ancestry-equal` | UG22..UG25 |
| Transform outcome | `--int-pheno` | UG27, UG28 |
| Output codec | `--out`, `--compression`, `--compression-level` | All analyses; UG26, UG29..UG31 |
| Compute controls | `--threads`, `--chunk-ksnp`, `--seed` | UG32 |
| Statistical routing | `--spagxe-marginal-cutoff`, `--spa-z-threshold` | UG32 |

Full definitions and limits are in [CLI_REFERENCE.md](CLI_REFERENCE.md),
[FEATURE_COVERAGE.md](FEATURE_COVERAGE.md), [METHOD_POLICY.md](METHOD_POLICY.md)
and [the validation record](validation/README.md).
