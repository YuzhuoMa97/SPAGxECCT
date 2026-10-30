# Individual allele-frequency model exports

The `--cal-af-coef` utility exports the same individual-frequency policy used by native SPAGxEmix. SPAGxE methods and release maintainer: **Yuzhuo Ma <yuzhuoma@amss.ac.cn>**. Export is provided for inspection and reuse through the documented C++ reader. The association CLI estimates AF online; it does **not** import an AFC file as an association cache.

## Command and fitting cohort

From `cpp/`:

```sh
mkdir -p results/af
build/bin/SPAGxE --cal-af-coef \
  --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pc-cols PC1,PC2 \
  --maf 0 --mac 0 --geno 1 --hwe 0 --threads 4 \
  --out results/af/example
```

This writes `results/af/example.afc`; add `--compression gz` or `zst` for the matching native text codec. Windows uses `build\bin\SPAGxE.exe`; create the output directory and use your shell's continuation syntax or one command line. The permissive QC values illustrate the small simulated data, not recommended cohort thresholds.

The utility supports the global genotype formats accepted by the executable, sample keep/remove lists, marker extract/exclude/chromosome selection, QC and threading. Its PCs are matched by string IID, in the specified PC column order and original units. Missing-PC removal and finite/full-rank design checks precede fitting. Covariate-file mean filling follows the existing subject-table rules.

For each marker, QC uses observed nonmissing genotype statistics. G missing values are filled with `2*ALT_FREQ` only after QC passes. The AF model fits in the cohort's minor-allele direction and returns canonical ALT frequencies. Its minor-MAC gate is the minor count in that mean-imputed fitting cohort, mathematically `2*N*min(ALT_FREQ,1-ALT_FREQ)`. A model-gate minor MAC at most 20 selects pooled AF. The implementation sums the actual minor-oriented G vector to avoid cancellation in 1-ALT_FREQ moving an exact hardcall count of 20 above that boundary.

OLS is accepted when fewer than 10% of fitted genotypes lie outside [0,2]; accepted values are clipped to AF [0,1]. Otherwise, significant PC slopes select the minor-carrier logistic fit. There is no second rule requiring more than 20 carriers. Constant carrier outcomes retain the documented pooled routing. A selected AF fit that is invalid, cannot be solved, or does not converge is recorded as failed; it is not substituted by a pooled-model P value.

## Text schema and allele direction

Text files begin with:

```text
#STATUS BETA0 BETA1 ...
#AFMODEL_SCHEMA 2 -2=AF_fit_failed -1=unfitted_or_QC_filtered 0=pooled_ALT 1=ALT_OLS 2=logistic 3=complement_logistic
```

Actual fields are TAB separated. Skip **every** comment line beginning with `#`, then read one status and exactly `1+nPC` coefficients per data row. Coefficients are written with 17 significant digits. The first beta is the intercept; remaining betas follow `--pc-cols` order. There is one data row per marker selected by the genotype catalog, including a status -1 row when selected-marker QC fails. These files do not contain marker IDs or allele metadata; retain the matching selected marker list and its order.

Let `X=[1|PCs]`, `sigmoid(t)=1/(1+exp(-t))`, and `f(t)=1-sqrt(1-sigmoid(t))`.

The carrier-to-AF conversion is evaluated from the linear predictor using
`f(t)=-expm1(-softplus(t)/2)`, where
`softplus(t)=max(t,0)+log1p(exp(-abs(t)))`. This is the same mathematical
model. It avoids first rounding the carrier probability to one and losing
the corresponding allele-frequency complement. The complement used by
status 3 is evaluated as `exp(-softplus(t)/2)`. This numerical repair does
not alter the minor-MAC threshold, OLS eligibility, selected PCs, convergence
requirement or failed-fit status.

| STATUS | Meaning | Reconstructed ALT AF |
| --- | --- | --- |
| -2 | Selected AF fit failed or was invalid | No usable coefficients; reconstruction is rejected |
| -1 | Not fitted, filtered by QC, or an unwritten raw binary slot | No fitted AF; the decoder rejects it |
| 0 | Completed pooled-frequency fallback | Caller-supplied pooled `ALT_FREQ` |
| 1 | ALT-scale OLS coefficients | `clip(X*beta/2,0,1)` |
| 2 | Logistic where cohort ALT is the fitted minor direction | `f(X*beta)` |
| 3 | Logistic fitted in minor direction when cohort ALT is major | `1-f(X*beta_minor)` |

For a high-ALT OLS fit, stored coefficients are `beta0_ALT=2-beta0_minor` and `beta_j_ALT=-beta_j_minor`. For a high-ALT logistic fit, simply negating beta does not produce the required ALT complement. Status 3 keeps the minor-carrier coefficients and explicitly complements the nonlinear AF transform. GxE scores and global effects remain ALT-oriented.

Status 0 contains zero betas and requires the matching pooled ALT frequency for reconstruction. Status -1 also contains zero betas but means that no eligible model was fitted; status -2 contains zero placeholders for a failed fit. Neither negative status may be interpreted as status 0. Unknown statuses, invalid dimensions and nonfinite coefficients fail explicitly.

## Binary record layout

The C++ `IndivAFWriter` additionally accepts a filename ending in `.bin`. The current CLI writes the text `.afc` family.

A binary export has no header. Each raw genotype marker index has a fixed record:

```text
int8 STATUS
double BETA0
double BETA1
...
double BETA_nPC
```

On the supported builds, stride is `1+8*(1+nPC)` bytes and record offset is `raw_marker_index*stride`. Raw marker indices are zero based. Each unwritten slot is initialized to status -1 with zero coefficients, including markers omitted by selection. States 2 and 3 use the same record size. The representation follows the build's native double/byte order; no endian conversion, nPC header, marker metadata or automatic provenance check is included. Prefer text for interchange and retain the exact genotype/PC metadata.

A binary reader needs nPC and raw marker indices. A text reader needs nPC and the exact selected-marker order. Neither file independently identifies sample IIDs, PC labels/scaling or REF/ALT alleles. Regenerate older exports for this release: previous zero placeholders cannot retrospectively be distinguished from completed pooled fits, and older readers may not understand states -2, -1 and 3.

## QC MAC, fitting MAC and phenotype masks

QC and model selection have different MAC roles:

- `--mac` filters the **observed** nonmissing minor count from the genotype backend, before imputation.
- The AF-cascade threshold of 20 uses the minor count of the completed, mean-imputed **fitting cohort**. This is a model-selection rule, not an extra user QC threshold.
- `--geno`, `--maf` and `--hwe` are applied through the existing backend statistics and threshold definitions. HWE uses the native global hardcall/dosage-to-hardcall diagnostic. All-missing markers are unfitted.

For model/export equivalence, use the same IIDs, PC columns/order/values, genotype alleles, sample masks and missing-data rules. Online SPAGxEmix refits AF for each phenotype's analysis mask. An export made on all PC-complete people need not equal an online fit after a phenotype, environment or GRM exclusion. Reproduce that analysis mask with explicit selected IIDs when making a matched export.

The model is diploid and global. Known-H local ancestry input has a different G/H probability law; these coefficients do not replace its local ancestry counts.

## Verification

```sh
python -B tests/verify_af_utility.py build/bin/SPAGxE results/verify_af_utility
```

The standard-library test generates its own PLINK files, checks three explicit AF boundary cases and high-ALT complements against an independent binomial Newton fit, and exercises QC, thread equality, extract and gzip. The C++ `af_model_test` additionally checks model-to-vector equivalence, fractional dosage, the MAC=20 boundary, binary/text/gzip storage and rejection of unknown/unfitted statuses. These checks validate the measured numerical and I/O cases; they do not establish universal population calibration.
See [release scope](RELEASE_STATUS.md). The test commands above regenerate their own summaries.


## Failed AF fits

The text export uses AFMODEL_SCHEMA 2. Status -2 is AF_fit_failed; status -1
remains unfitted/QC_filtered, 0 pooled ALT, 1 ALT OLS, 2 logistic and 3
complement-logistic. Binary record size is unchanged (one status byte and
1+nPC doubles). A failed record stores zero placeholders, not usable coefficients.
Readers preserve that failed status and refuse to reconstruct AF from it.

IRLS uses tolerance 1e-6 and at most 100 iterations, with factorization and
finite-value checks. Finite coefficients alone do not establish convergence.
Failure does not silently switch ancestry adjustment to a pooled model.
In association output, that variant has NA marginal/interaction P values and
SPA_STATUS 9 (AF_FIT_FAILED). A valid monomorphic model retains NA_NO_TEST 8.
