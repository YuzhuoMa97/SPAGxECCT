# SPAGxE 1.2.0 command-line reference

SPAGxE methods and release maintainer: **Yuzhuo Ma, yuzhuoma@amss.ac.cn**.
Commands run from cpp/. Windows source builds use build/bin/SPAGxE.exe.
See [installation](INSTALLATION.md) and [tutorials](../cpp/examples/tutorials/README.md).

## Statistical definition

Branch B uses the published E*Rtilde weight. All four methods use guarded
partial-normal CGF acceleration beyond the normal/SPA Z threshold (default 2).
Fixed type-7 Tukey fences use 1.5 IQR; full-CGF candidate-root diagnostics use
tolerance 1e-3, with automatic full evaluation when needed. This tolerance is
not a bound on P-value error. Homogeneous GRM correction uses score adjustment. QC/thresholds remain configurable; there is no second
algorithm profile. Obsolete algorithm-profile, spa-mode, spa-exact,
grm-spa-adjust and outlier-iqr-multiplier selectors are rejected.
E is used as supplied, without silent centering. Choose its coding before
analysis; see [method policy](METHOD_POLICY.md).

## Methods and input

| Option | Meaning |
|---|---|
| --method | SPAGxE (default), SPAGxE+, SPAGxEmix, SPAGxEmix_local or SPAGxEcct |
| --bfile PREFIX | BED/BIM/FAM |
| --pfile PREFIX | PGEN/PVAR/PSAM |
| --bgen FILE MODE | BGEN; supply ref-first, ref-last or ref-unknown |
| --vcf FILE / --bcf FILE | VCF, VCF.gz or BCF |
| --pheno FILE | Headered raw phenotype or residual table |
| --pheno-name COLS | Comma-separated raw phenotypes; survival uses TIME:EVENT |
| --resid-name COLS | Supplied residuals; mutually exclusive with raw traits; no prospective Wald refit |
| --covar FILE / --covar-name COLS | Covariates, taken from pheno if no separate table is supplied |
| --envir-name COLS | Required tested environment column(s), comma separated |
| --regression-model | auto, linear, logistic, cox or ordinal |
| --save-resid | Export native residuals in raw-trait mode |
| --pc-cols COLS | Individual-AF PCs; default PC1,PC2,PC3,PC4 |
| --grm-file FILE | Sparse matrix with header ID1 ID2 VALUE; relationship=2*kinship |
| --sp-grm-plink2 FILE | Sparse .grm.sp and companion .grm.id; mutually exclusive with --grm-file |
| --local-manifest FILE | Paired G/H manifest or automatically recognized local-BGEN manifest |
| --local-prefix PREFIX | Packed .localbed/.localbim/.localfam count input |
| --local-bgen-manifest FILE | ancestry, bgen and counted_allele columns |
| --local-max-block-bytes | 268435456 by default; permitted 4096..1073741824 decoded bytes |

Choose exactly one global or one local genotype source. SPAGxE+ requires GRM;
supplying GRM to SPAGxE selects the same method. SPAGxEmix and local require
unrelated subjects and reject GRM. Local also rejects global PC-AF input.
SPAGxElocal and local-binary remain compatible spellings.

Unique string IIDs align samples; FID may also be present. Include all tested
environments in null-model covariates; mix PCs must also adjust that model.
Linear, logistic, Breslow-ties Cox and ordered cumulative-logit models run in C++.

## QC, selection and resources

| Option | Default/range |
|---|---|
| --spagxe-marginal-cutoff P | .001, [0,1]; P at or below the threshold enters Branch B |
| --spa-z-threshold Z | 2, at least .01 |
| --geno P | .1, [0,1], maximum missing fraction |
| --maf P | 1e-5, [0,.5] |
| --mac COUNT | 10, nonnegative |
| --hwe P | 0 disables; otherwise minimum diagnostic P in [0,1] |
| --hard-call-threshold X | .1, [0,.5), HWE hardcall distance, not an analysis-dosage replacement |
| --threads INT | 1, positive |
| --chunk-ksnp FLOAT | 8, positive multiple of .5; each unit is 1024 markers |
| --extract FILE / --exclude FILE | Include/exclude variant IDs |
| --chr SPEC | Chromosome set/range, for example 1-4,6,22,X |
| --keep FILE / --remove FILE | Include/exclude IIDs |
| --out PREFIX | Required; create the parent directory first |
| --compression gz or zst | Omit for plain text |
| --compression-level INT | gzip default 6/range 1..9; zstd default 3/range 1..22 |
| --seed UINT64 | Compatibility setting; association calculations are deterministic |

Global QC MAC uses observed nonmissing minor dosage. The individual-AF cascade
uses the completed fitting cohort for its distinct MAC=20 model-selection gate.
Local MAC uses estimated minor AF times known H after masking. Unknown H is not
imputed; any missing ancestry H removes that person from every ancestry branch
at the marker. Missing G is H*q and H=0 is structural zero.
HWE is a summed global diploid hardcall diagnostic, not an exact variable-H
ancestry test. Non-point posteriors cannot be treated as hardcalls merely
because their expected count is integer. See [local formats](LOCAL_ANCESTRY_FORMAT.md).

## Combination and utilities

--method SPAGxEcct uses --cct-global, --cct-local, --cct-phenotype and --out.
Default global/local group weights are .5,.5; --cct-weights supplies two positive
weights. --cct-per-ancestry-equal gives every available component equal weight.
Variant, allele, phenotype and environment identities must match.

--cal-af-coef exports the online minor-oriented AF policy. Unfitted/QC-filtered
rows have status -1; status 3 is complemented minor-carrier logistic AF.
Exports are not association-cache inputs. See [AF exports](AF_MODEL_FORMAT.md).

--pack-local creates validated hardcall G/H count triplets without fitting
association or inferring ancestry; see [packing](LOCAL_PACKING.md).

--int-pheno transforms specified non-ID columns to inverse-normal values,
preserving missing values. Without selected columns, all non-ID columns are
transformed. It rejects association options and output compression.

## Output interpretation

LOG10P is -log10(P). Outputs separate marginal statistics, normal score, SPA,
available Wald, final GxE P, working-score BETA/SE, branch and SPA_STATUS.
Working-score BETA/SE are not the prospective Wald coefficient/SE.
Branch 1=A and 2=B; unavailable Wald fields are NA. Related-sample final P is SPA-only.

SPA_STATUS 0 is success including support/endpoint handling; 1 is the normal
branch; 2 is successful reduction near zero saddlepoint. Codes 3..6 are guarded
normal fallbacks; codes at least 7 mean no test or postprocessing failure.
Interpret status with P. Numerical failure is not replaced by P=0.

Global scans write PREFIX.analysis.settings.tsv and local scans write
PREFIX.local.settings.tsv with version/policy receipts. The CGF policy is
`guarded_partial_normal`. Automatic full-CGF protection is part of this one
policy and is distinct from a normal fallback after a numerical SPA failure.
Both partial and full-CGF SPA remain approximations, not finite-sample enumeration.
See [partial SPA](PARTIAL_SPA.md). See
[equations](METHOD_EQUATIONS.md) and [method policy](METHOD_POLICY.md).


AF failure status: 9 AF_FIT_FAILED means the selected individual-AF
reference model did not converge, could not be solved or was invalid. Its P values are NA.
This is distinct from 8 NA_NO_TEST (a valid model with no test statistic).
--help plus, --help mix and --help SPAGxE+ are accepted aliases;
an unknown help topic exits with an error.
