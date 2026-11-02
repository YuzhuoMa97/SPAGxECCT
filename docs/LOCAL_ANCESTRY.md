# Native local ancestry analysis

The canonical method name is **SPAGxEmix_local**; `SPAGxElocal` remains an accepted alias. File decoding, phenotype null models, residual processing, SPA, interaction Wald fits and CCT run in the SPAGxE C++ executable. SPAGxE methods and release maintainer: **Yuzhuo Ma <yuzhuoma@amss.ac.cn>**.

G is the ancestry-specific allele count or dosage. H is the **known integer** ancestry chromosome-copy count, 0, 1 or 2. Supply every ancestry, with `0 <= G <= H` and complete `sum(H)=2`. Input options are paired G/H text, a packed hardcall count container, or one known-H BGEN file per ancestry. See the [G/H format specification](LOCAL_ANCESTRY_FORMAT.md) and [local BGEN profile](LOCAL_BGEN_FORMAT.md).

## Native phenotype and residual scans

Run these examples from `cpp/`:

```sh
mkdir -p results/local
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-manifest examples/local/cohort.local.tsv \
  --pheno examples/local/samples.pheno \
  --pheno-name Quantitative,Binary,Ordinal,Time:Event \
  --covar-name X,E --envir-name E --regression-model auto \
  --maf 0 --mac 0 --geno 1 --threads 4 --save-resid \
  --out results/local/native
```

These permissive QC thresholds expose small simulated examples, including degenerate branches; choose appropriate thresholds for a real cohort. On Windows use `build\bin\SPAGxE.exe`, create the output directory, and put each command on one line or use the continuation syntax for your shell.

`auto` fits linear, logistic, ordinal or Cox null models according to each phenotype specification. Include the tested environments in the genotype-independent null model together with the intended covariates. The marginal genetic screen selects Branch A or Branch B. Raw phenotype mode can fit the full-model interaction Wald test in Branch B. The local Wald design includes independent ancestry-H nuisance main-effect columns; it is not a joint model of every ancestry's genotype and interaction.

This release uses guarded partial-normal CGF acceleration with the paper-defined
normal/SPA threshold. For each ancestry, central H-count terms contribute their
exact mean and variance and tail terms retain Bin(H,q) contributions. Full-CGF
diagnostics can select automatic full evaluation. This is one statistical
policy; see [partial SPA](PARTIAL_SPA.md).

For supplied residuals and multiple environments:

```sh
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-manifest examples/local/cohort.local.tsv \
  --pheno examples/local/samples.pheno \
  --resid-name R_quant,R_binary,R_ordinal,R_survival,R_masked \
  --envir-name E,Eshift --maf 0 --mac 0 --geno 1 \
  --out results/local/residual
```

Residual mode performs no prospective Wald refit because raw responses are unavailable. Each residual/phenotype has its own sample mask. `R_masked` illustrates that behavior; `Eshift` is a numerical stress case and does not imply that local Branch B is invariant to arbitrary environmental shifts.

## Packed counts and genotype probabilities

The independently simulated packed-count example contains 450 people, 18 markers and two ancestries:

```sh
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-prefix examples/local/binary_demo/local \
  --pheno examples/local/binary_demo/samples.pheno --resid-name R \
  --envir-name E --maf 0 --mac 0 --geno 1 \
  --out results/local/binary
```

`--local-binary` is an alias of `--local-prefix`. The prefix names matching `.localbed/.localbim/.localfam` files. These are ancestry G/H count files, not ordinary PLINK BED files.

The independently simulated BGEN example uses three ancestries, 16-bit quantized probabilities and zlib compression:

```sh
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-bgen-manifest examples/local/bgen_demo/bits16_zlib.localbgen.tsv \
  --pheno examples/local/samples.pheno \
  --pheno-name Quantitative,Binary,Ordinal,Time:Event \
  --covar-name X,E --envir-name E --regression-model auto \
  --maf 0 --mac 0 --geno 1 --threads 4 \
  --out results/local/bgen
```

The three-column BGEN manifest is also recognized by `--local-manifest`. Choose one local input option per run. Supported BGEN inputs are Layout 2, unphased and biallelic, with embedded IIDs, known per-ancestry ploidy H=0/1/2, raw/zlib/zstd blocks and probability precision 1 through 32 bits. The analysis consumes the posterior mean G with known H; it does not propagate the complete posterior into the SPA CGF. Uncertain ancestry counts and ancestry inference are outside this model.

The two- and three-ancestry paired-text examples use the same phenotype IDs. Ancestry labels and counts come from the manifest/header; analysis does not assume two ancestries. IIDs are matched as strings, independently of file sample order. Marker metadata and order must match exactly across ancestry files.

## Missingness, direction and results

Unknown H is never imputed. A person with any missing ancestry H is excluded from every ancestry test at that marker. Missing G with H=0 becomes structural zero; other missing G is filled with H*q using the ancestry frequency estimated from observed G/H.

Canonical input counts ALT. A REF-counted text value is converted using **H-G**, not 2-G. BGEN categories intrinsically count its second listed allele; a manifest REF declaration alone does not swap its actual alleles or probability categories. A real allele-order reversal must exchange both allele metadata and categories G versus H-G. The test then selects each ancestry's minor allele. `ALT_FREQ` retains the canonical ALT frequency, while `TESTED_ALLELE` and `MINOR_FLIPPED` identify the direction of signed test statistics. At frequency exactly 0.5, the ALT tie convention applies.

Each phenotype produces:

- `PREFIX.PHENOTYPE.SPAGxE_local`: one row per marker and ancestry, with QC, marginal screen, SPA, available Wald/CCT results, branch/status, allele direction and sample counts. The output suffix is retained independently of the canonical CLI method spelling. `LOG10P_*` is **-log10(P)**.
- `PREFIX.PHENOTYPE.SPAGxE_local_CCT`: an equal-weight combination of available ancestry-specific final GxE P values for each environment. Failed-QC/NA tests are omitted; inspect the contributing ancestry count. This is cross-ancestry CCT. Use the separate `SPAGxEcct` workflow when explicitly combining matched global and local results.
- `PREFIX.local.settings.tsv`: settings and completion receipt. `--save-resid` additionally exports fitted null residuals for reuse.

Other local options include `--keep`, `--remove`, `--extract`, `--exclude`, `--chr`, `--compression gz|zst`, `--threads`, `--spagxe-marginal-cutoff`, `--spa-z-threshold` and `--local-max-block-bytes`. Global individual-frequency PC and sparse-GRM options are separate interfaces; this local workflow does not claim a validated local+GRM covariance model.

## Numerical checks and comparison with the original R

Agreement with the preserved original R implementation is assessed only on a valid common domain: compatible additive coding, complete finite inputs, nondegenerate designs, matched null-model residuals/covariates/environments, and settings supported by both implementations. Direct original-source comparisons restrict ancestry frequency to its valid minor-allele domain. They do not establish equality for every original R extension or defective input branch.

The native implementation uses H-G for local allele reversal and H*q for missing local genotypes. Those corrections can differ from original R branches using 2-G or 2*q. Missing-H filtering, singular designs and dosage assumptions also require explicit interpretation. Corrected edge cases are checked against independent formulas rather than forced to reproduce defective outputs. Numerical agreement on simulated cases does not establish population-wide calibration for every ancestry inference or imputation mechanism.

Portable verification from `cpp/`:

```sh
python -B tests/verify_local.py build/bin/SPAGxE results/verify_local
python -B tests/verify_local_bgen.py build/bin/SPAGxE results/verify_local_bgen
```

The example count containers and BGEN probabilities are generated by this project's simulation scripts. Independent numerical references use `tests/local_oracle.R` and `tests/verify_local_oracle.py`; R is needed only for those development comparisons, not for delivered analyses.