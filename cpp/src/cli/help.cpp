// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE standalone help. Upstream C++ infrastructure attribution is retained
// in THIRD_PARTY_NOTICES; SPAGxE algorithms: Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
#include "util/version.hpp"
#include "cli/cli.hpp"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>
#include <stdexcept>

namespace cli {
namespace {
struct Option { const char *name; const char *argument; const char *description; };
const Option options[] = {
    {"method", "SPAGxE|SPAGxE+|SPAGxEmix|SPAGxEmix_local|SPAGxEcct", "Analysis method (default: SPAGxE); case-insensitive."},
    {"bfile", "PREFIX", "PLINK 1 .bed/.bim/.fam genotype prefix."},
    {"pfile", "PREFIX", "PLINK 2 .pgen/.pvar/.psam genotype prefix."},
    {"bgen", "FILE MODE", "BGEN; MODE is ref-first, ref-last, or ref-unknown (required)."},
    {"vcf", "FILE", "VCF or .vcf.gz genotype input; dosages preferred over GT."},
    {"bcf", "FILE", "BCF2 genotype input (htslib backend)."},
    {"pheno", "FILE", "Headered whitespace-delimited phenotype/residual table; IID or FID IID."},
    {"pheno-name", "COLS", "Comma-separated phenotype columns; Cox traits use TIME:EVENT."},
    {"resid-name", "COLS", "Precomputed residual columns; mutually exclusive with --pheno-name."},
    {"covar", "FILE", "Optional covariate table. Without --covar-name, all non-ID columns enter the null model."},
    {"covar-name", "COLS", "Null-model covariates from --covar, or from --pheno when --covar is absent."},
    {"envir-name", "COLS", "Environment columns (required); in fit mode also list each in --covar-name."},
    {"regression-model", "MODEL", "Null family: auto (default), linear, logistic, cox, ordinal."},
    {"save-resid", "", "Save fitted residuals in PREFIX.null.resid[.gz|.zst]."},
    {"pc-cols", "COLS", "SPAGxEmix/AF-model PCs (default: PC1,PC2,PC3,PC4); mix fit also requires them in --covar-name."},
    {"grm-file", "FILE", "Optional sparse GRM: header ID1 ID2 VALUE; engages SPAGxE+ with --method SPAGxE."},
    {"sp-grm-plink2", "FILE", "Optional sparse GRM in PLINK 2 .grm.sp format, with companion .grm.id."},
    {"pack-local", "", "Pack integer local G/H into .localbed/.localbim/.localfam; no analysis or R runtime."},
    {"local-bim", "FILE", "Optional six-column BIM metadata for --pack-local; allele/order validation required."},
    {"local-fam", "FILE", "Optional six-column FAM for --pack-local; preserves fields and matches samples by unique IID."},
    {"local-block-variants", "UINT", "Maximum packed variants per block (default64); 1..65536, bounded by decoded memory."},
    {"local-zstd-level", "INT", "Packed-count zstd level -5..22 (default3; 0 uses zstd default)."},
    {"local-manifest", "FILE", "SPAGxElocal paired ancestry G/H manifest; mutually exclusive with --local-prefix."},
    {"local-bgen-manifest", "FILE", "Local known-H BGEN layout-2 manifest: ancestry,bgen,counted_allele (TAB)."},
    {"local-prefix", "PREFIX", "SPAGxElocal SQRLOC01 binary prefix from the compatible hardcall layout; alias: --local-binary."},
    {"local-max-block-bytes", "UINT64", "Local binary decompression memory guard (default: 268435456; minimum: 4096)."},
    {"cct-global", "FILE", "SPAGxEcct global association table (SPAGxE/SPAGxEmix output)."},
    {"cct-local", "FILE", "SPAGxEcct across-ancestry SPAGxE_local_CCT table."},
    {"cct-phenotype", "NAME", "Required phenotype identity for SPAGxEcct; checked against both input filenames."},
    {"cct-weights", "GLOBAL,LOCAL", "SPAGxEcct positive group weights (default: 0.5,0.5; normalized internally)."},
    {"cct-per-ancestry-equal", "", "SPAGxEcct weights global:local as 1:K_valid using each row/environment's ancestry count."},
    {"spagxe-marginal-cutoff", "P", "Branch routing threshold in [0,1] (default: 0.001)."},
    {"spa-z-threshold", "Z", "Use SPA above |Z| threshold (default: 2; minimum: 0.01)."},
    {"threads", "INT", "Worker threads (default: 1; minimum: 1)."},
    {"chunk-ksnp", "FLOAT", "Markers per chunk in units of 1024 (default: 8; positive multiples of 0.5)."},
    {"geno", "P", "Maximum missing-genotype rate in [0,1] (default: 0.1)."},
    {"maf", "P", "Minimum minor-allele frequency in [0,0.5] (default: 1e-5)."},
    {"mac", "COUNT", "Minimum minor-allele count (default: 10)."},
    {"hwe", "P", "Minimum exact Hardy-Weinberg P in [0,1] (default: 0; disabled)."},
    {"hard-call-threshold", "X", "Maximum distance to 0/1/2 for HWE counts only (default: 0.1; [0,0.5))."},
    {"extract", "FILE", "Variant-ID include list."},
    {"exclude", "FILE", "Variant-ID exclude list."},
    {"chr", "SPEC", "Chromosome labels or numeric ranges, e.g. 1-4,6,22,X."},
    {"keep", "FILE", "Subject-IID include list (shared input parser also accepts FID IID)."},
    {"remove", "FILE", "Subject-IID exclude list."},
    {"out", "PREFIX", "Required output prefix; the parent directory must exist."},
    {"compression", "gz|zst", "Compress association/residual/AF output (default: plain text)."},
    {"compression-level", "INT", "Association/residual compression: gz 1..9 (default 6), zst 1..22 (default 3)."},
    {"seed", "UINT64", "Compatibility option; SPAGxE calculations are deterministic and do not consume a seed."},
    {"cal-af-coef", "", "Utility: estimate and export per-variant individual-AF model coefficients."},
    {"int-pheno", "", "Utility: inverse-normal transform phenotype columns into PREFIX.txt."},
    {"help", "[TOPIC]", "Show help: spagxe, spagxemix, spagxelocal, genotype, phenotype, sp-grm, output, options, or flag name."},
    {"version", "", "Print the SPAGxE version and exit."}
};
std::string lowercase(std::string value) {
    for (auto &ch : value) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return value;
}
void printOption(const Option &option) {
    std::string label = std::string("--") + option.name;
    if (*option.argument) label += std::string(" ") + option.argument;
    std::fprintf(stderr, "  %-42s %s\n", label.c_str(), option.description);
}
void shortHelp() {
    std::fprintf(stderr, "SPAGxE %s | Yuzhuo Ma | yuzhuoma@amss.ac.cn\n", SPAGXE_VERSION);
    std::fputs(R"(
Standalone C++ gene-environment interaction analysis.
SPAGxE algorithms: Yuzhuo Ma. Incorporated components: THIRD_PARTY_NOTICES.md.

Usage:
  SPAGxE --bfile GENO --pheno TABLE --pheno-name Y \
         --covar-name AGE,SEX,E --envir-name E --out results/run
  SPAGxE --method SPAGxEmix --bgen GENO.bgen ref-first \
         --pheno TABLE --pheno-name Y --pc-cols PC1,PC2 \
         --covar-name AGE,SEX,E,PC1,PC2 --envir-name E --out results/mix

  SPAGxE --method SPAGxEmix_local --local-manifest local/manifest.tsv \
         --pheno TABLE --pheno-name Y --covar-name AGE,E \
         --envir-name E --out results/local

Methods: SPAGxE (default), SPAGxE+, SPAGxEmix, SPAGxEmix_local; SPAGxEcct combines result tables.
Supply --grm-file or --sp-grm-plink2 with SPAGxE to select SPAGxE+. SPAGxEmix does not accept GRM input.
Utilities: --cal-af-coef; --int-pheno; --pack-local.
Global genotypes: exactly one of --bfile, --pfile, --bgen, --vcf, --bcf.
Local genotypes: --local-manifest, --local-prefix or --local-bgen-manifest.
Algorithm: one guarded partial-normal CGF policy for all four methods.
Fixed 1.5-IQR tails retain binomial terms; central terms match mean/variance.
Full-CGF root diagnostics (1e-3) select automatic full evaluation when needed.
This diagnostic tolerance is not a p-value error bound or a user profile.
The normal/SPA threshold defaults to 2; set it with --spa-z-threshold.
Tables: subjects are matched by IID; each environment enters the fitted null model.

Run 'SPAGxE --help options' for all flags, or '--help TOPIC' for details.
Examples and reproducible validation are provided in the repository.
)", stderr);
}
void phenotypeHelp() {
    std::fputs(R"(Phenotype/covariate input

Whitespace-delimited with a mandatory header. Accepted identifier layouts:
  IID Y AGE E             or      #FID IID Y AGE E
  S1  1.7 54 0.4                  F1   S1  1.7 54 0.4
FID is ignored; IID must identify a subject uniquely. Missing numeric cells:
NA, na, NaN, nan, . or -. Input supports plain text, gzip and zstd.

Fit mode: --pheno-name Y1,Y2 --regression-model auto|linear|logistic|cox|ordinal.
The null model is trait ~ intercept + X + E, fitted once per trait; include
all --envir-name columns in --covar-name. With --method SPAGxEmix, every
--pc-cols column must also enter --covar-name. For a survival trait specify
--pheno-name TIME:EVENT, with strictly positive times and 0/1 event indicators.
Ordinal traits use ordered numeric categories; check the validation report
for the scope of agreement with the original R implementation.

--covar-name selects columns from --covar; when --covar is absent it selects
from --pheno. If --covar is given without --covar-name, all non-ID columns
are used. Environment and PC columns may be supplied in either input table.

Residual mode: --resid-name R1,R2 uses precomputed null residuals and does not
refit a null model or add a Branch-B Wald leg. Supply environment columns
(and PCs for SPAGxEmix) in --pheno or --covar. --save-resid is fit-mode only.

Without --pheno-name and --resid-name, non-ID phenotype columns are selected
after excluding specified covariates, environments, and SPAGxEmix PCs.
Explicit --pheno-name is recommended for tables containing extra metadata.
Multiple traits are analyzed in one invocation, with trait-specific missing
values handled by the engine. Multiple environments are comma-separated.
)", stderr);
}
void genotypeHelp() {
    std::fputs(R"(Genotypes: choose exactly one input
  --bfile PREFIX       PLINK 1 .bed/.bim/.fam (SNP-major .bed).
  --pfile PREFIX       PLINK 2 .pgen/.pvar/.psam.
  --vcf FILE           VCF, .vcf.gz (htslib); use --bcf for a BCF file.
  --bcf FILE           BCF2 binary (htslib).
  --bgen FILE MODE     BGEN with explicit allele orientation:
                        ref-first: first allele is REF (UK Biobank convention).
                        ref-last: last allele is REF (PLINK-export convention).
                        ref-unknown: last allele is provisional REF; emits a warning.

BGEN sample identifiers are read from the embedded sample block or the
companion .sample file. Preserve the sample file beside its BGEN input when
required. Format-specific compression/dependency support is described in
README.md. The tested format matrix is in the validation report.

ALT is the tested allele. Missing genotypes are mean-imputed after QC.
--geno / --maf / --mac / --hwe apply variant QC; --hard-call-threshold affects
HWE hard calls only. --extract / --exclude / --chr subset variants; --keep /
--remove subset subjects. List files and allele ordering must match your data.
)", stderr);
}
void grmHelp() {
    std::fputs(R"(SPAGxE+ sparse GRM (optional)
  --grm-file FILE     Header: ID1 ID2 VALUE, one related pair per row.
  --sp-grm-plink2 FILE   Sparse .grm.sp and companion .grm.id (FID IID).
Supply one format, not both. VALUE is the genomic relationship (2*kinship),
with self-relationship 1+F. Subjects are aligned to genotype IIDs. Preserve
the companion .grm.id file to avoid ambiguity in numeric subject indices.
The GRM corrects retrospective score variances through quadratic forms;
no pairwise IBD input is needed. GRM input is supported only for SPAGxE+.
SPAGxEmix and local-ancestry methods require unrelated samples.

The centered score is rescaled by the square root of the independent/GRM
variance ratio before solving the guarded partial-normal SPA tails. This release
uses one score-based GRM adjustment convention.
)", stderr);
}
void outputHelp() {
    std::fputs(R"(Association output: PREFIX.<TRAIT>.SPAGxE[.gz|.zst]
               or PREFIX.<TRAIT>.SPAGxEmix[.gz|.zst]

Variant columns: CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC LOG10P_HWE
Marginal block: LOG10P_G Z_G BETA_G SE_G SPA_STATUS_G
Each environment appends:
  LOG10P_Gx<E>       Final interaction result: SPA/normal, or Branch-B CCT.
  LOG10P_Wald_Gx<E>  Branch-B Wald leg; NA when that leg is unavailable.
  Z_Gx<E>           Signed z corresponding to the final p-value.
  Z_Norm_Gx<E>      Uncalibrated score / sqrt(variance).
  BETA_Gx<E> SE_Gx<E>  Score-derived effect and SE (not the Wald coefficient).
  SPA_STATUS_Gx<E>  Status of the SPA/normal score-test leg.
  LOG10P_SPA_Gx<E>  Separate SPA/normal score-test p-value before CCT.
  BRANCH_Gx<E>      1=Branch A (marginal P>cutoff); 2=Branch B (P<=cutoff).

All LOG10P columns store -log10(P). Recover P=10^(-LOG10P); log-domain
results preserve extremely small p-values without underflow. Branch B
combines the score and Wald p-values by Cauchy combination only where a
supported raw-trait Wald fit exists. GRM and residual modes have no Wald leg.
A finite combined p-value does not erase a failed SPA leg: inspect its status.

SPA_STATUS codes:
  0 SPA_OK; 1 NORMAL; 2 SPA_W_SINGULAR;
  3 FALLBACK_MAXITER; 4 FALLBACK_GUARD_TEMP; 5 FALLBACK_GUARD_CURV;
  6 FALLBACK_NONFINITE; 7 NA_POST_FAIL; 8 NA_NO_TEST; 9 AF_FIT_FAILED.
Codes 3-6 indicate normal-approximation fallback after numerical SPA failure;
codes 7-9 indicate unavailable tests and NA. AF_FIT_FAILED means the individual
AF fit did not converge or was invalid; no pooled substitute or p-value is used.
Check status columns in analyses.

Every subject contributes to the working mean and variance. Central weights
use a Gaussian CGF and tail weights retain binomial terms. Automatic full-CGF
evaluation protects unsuitable candidates; it is distinct from a numerical
normal fallback. Both paths remain saddlepoint approximations, not exact
discrete p-value tests. See docs/PARTIAL_SPA.md.
)", stderr);
}
void localHelp() {
    std::fputs(R"(SPAGxEmix_local (SPAGxElocal alias): native local-ancestry GxE and across-ancestry CCT

Select exactly one paired ancestry genotype/haplotype-count input:
  --local-manifest FILE  TAB header: ancestry genotype hapcount counted_allele.
                         One row per ancestry; counted_allele is ALT or REF.
                         G/H tables have CHROM POS ID REF ALT then sample IIDs.
  --local-bgen-manifest FILE known-H BGEN layout-2 ancestry set.
  --local-prefix PREFIX SQRLOC01 binary prefix (alias --local-binary).
                         Preserve its companion variant/sample metadata.

Each ancestry uses G_i | H_i ~ Binomial(H_i,q), with integer H in {0,1,2}
and 0<=G<=H. REF-counted input is explicitly recoded to canonical ALT.
AF excludes missing G and its H denominator; missing G is filled with H*q.
If q>0.5 the kernel tests H-G at frequency 1-q. Output TESTED_ALLELE and
MINOR_FLIPPED identify this choice; effect/z fields refer to that tested
minor ancestry allele, rather than necessarily the global ALT. At q=0.5
ALT is retained. Variant-level unknown H removes that subject from every
ancestry test for the same variant; N_TESTED and N_H_MISSING report the mask.

Branch A uses lambda=sum(H*E*R^2)/sum(H*R^2). Branch B removes [1,G] from R
and uses E*R0. Guarded partial SPA uses H-weighted central Gaussian moments
and exact hapcount-binomial tail terms, with automatic full-CGF protection.
Raw-trait Branch B Wald adjusts for independent local H main-effect columns
(at most K-1, in addition to ordinary covariates); residual input has no Wald.
This follows the corrected original-R local test; H:E is not an additional
nuisance term. It does not claim a general conditional ancestry-by-environment
model with all ancestry interaction effects fitted simultaneously.

Each ancestry is reported separately in PREFIX.TRAIT.SPAGxE_local; the
SPAGxE_local_CCT file combines finite, QC-passing ancestry p-values with equal
weights. N_ANCESTRIES records how many entered. This is an across-ancestry
CCT, not a separate combination with a global SPAGxEmix test.

Use --pheno/--pheno-name or --resid-name, --covar-name, --envir-name and --out.
Subject/variant filters, QC, threading, chunking, compression and residual
export are supported. The HWE diagnostic uses summed global genotypes; it
does not impose diploid HWE on a variable-H ancestry component.
SPAGxElocal rejects ordinary global genotype inputs, sparse GRM, --pc-cols,
and --hard-call-threshold. Null-model PC adjustment uses --covar-name.
See docs/LOCAL_ANCESTRY.md for the complete data contract and examples.
)", stderr);
}
void cctHelp() {
    std::fputs(R"(SPAGxEcct: native Cauchy combination of global and local GxE results

  SPAGxE --method SPAGxEcct --cct-global global.Y.SPAGxEmix \
    --cct-local local.Y.SPAGxE_local_CCT --cct-phenotype Y \
    --envir-name E --out results/combined

Inputs may be plain text, .gz or .zst. Both must contain the same variants,
in the same order, with identical CHROM POS ID REF ALT. Duplicate keys,
different variant counts, missing environments and malformed P values fail
before the result file is opened. Filename suffixes must identify the declared
phenotype; this checks declared identity, not whether the source model was
scientifically specified for that phenotype.

Global input columns: LOG10P_Gx<ENV>.
Local input columns: LOG10P_ANCESTRY_CCT_Gx<ENV>, N_ANCESTRIES_Gx<ENV>.
Global filenames end in .PHENOTYPE.SPAGxE/.SPAGxEmix;
local filenames end in .PHENOTYPE.SPAGxE_local_CCT, before .gz/.zst.

Default --cct-weights 0.5,0.5 gives equal weight to the global and local groups.
Use --cct-weights 1,K for a fixed local-group weight K; this represents K
individual ancestry tests only if their valid count is constant. Alternatively
--cct-per-ancestry-equal uses 1:K_valid for every row/environment and thus
gives the global result and each valid ancestry result equal nominal weight.
These two weighting options are mutually exclusive.

If either group result is missing, the combined result is NA with a named
CCT_STATUS; missing legs are not silently dropped. Inputs store -log10(P);
the weighted calculation remains in the log domain for extreme tails. As in
the original-R CCT helper, P=1 is clamped to 0.999 before combination.
Output: PREFIX.PHENOTYPE.SPAGxE_global_local_CCT[.gz|.zst], containing both
input legs, normalized weights, valid ancestry counts and combination status.
No genotype/phenotype model fitting occurs; analysis-only flags are rejected.
)", stderr);
}
} // namespace
void printVersion() { std::fprintf(stdout, "SPAGxE %s\n", SPAGXE_VERSION); }
void printHelp(const std::string &rawTopic) {
    std::string topic = lowercase(rawTopic);
    if (topic.rfind("--", 0) == 0) topic.erase(0, 2);
    if (topic.empty() || topic == "__short__") { shortHelp(); return; }
    if (topic == "options") {
        std::fputs("SPAGxE command-line options:\n\n", stderr);
        for (const auto &option : options) printOption(option);
        return;
    }
    if (topic == "phenotype") { phenotypeHelp(); return; }
    if (topic == "genotype") { genotypeHelp(); return; }
    if (topic == "sp-grm" || topic == "grm" || topic == "plus" || topic == "spagxe+") { grmHelp(); return; }
    if (topic == "output") { outputHelp(); return; }
    if ((topic == "spagxelocal" || topic == "spagxemix_local") || topic == "local" || topic == "local-binary") { localHelp(); return; }
    if (topic == "spagxecct" || topic == "cct") { cctHelp(); return; }
    if (topic == "spagxe" || topic == "spagxemix" || topic == "mix") {
        shortHelp(); phenotypeHelp();
        grmHelp();
        outputHelp(); return;
    }
    if (topic == "cal-af-coef") {
        std::fputs(R"(SPAGxE --cal-af-coef --bfile GENO --covar PCS --pc-cols PC1,PC2 --out PREFIX
Estimate individual-AF model coefficients; writes PREFIX.afc[.gz|.zst].
AF status -2 records a failed/non-convergent fit; -1 is unfitted/QC filtered.
Failed coefficients cannot be reconstructed or used for association.
Supports all genotype formats, subject/variant filters, QC and threading.
No phenotype fitting or GRM options. Compression uses the codec default.
SPAGxEmix currently estimates individual AF on the fly; it does not read this
exported coefficient file. The utility is provided for inspection/reuse.
)", stderr); return;
    }
    if (topic == "pack-local") {
        std::fputs(R"(SPAGxE --pack-local --local-manifest INPUT.tsv --out PREFIX

Writes PREFIX.localbed, PREFIX.localbim and PREFIX.localfam using the documented
SQRLOC01 count format. Choose one input: --local-manifest, --local-prefix or
--local-bgen-manifest. Only integer G and known H (or explicit missing values)
are packed; fractional dosages/non-point posteriors are rejected, not rounded.
Known H=0 has structural G=0. Output direction is canonical ALT.

Optional --local-bim ORIGINAL.bim preserves its six-column variant metadata,
but CHROM/ID/POS and column5 ALT / column6 REF must match the G/H stream in
order. --local-fam ORIGINAL.fam preserves six columns and aligns by unique IID;
the same sample set is required. Ordinary BIM uses A1/A2: confirm orientation.
The packed footer binds actual sidecar bytes by CRC; do not swap them later.

--local-block-variants 64 and --local-zstd-level 3 are defaults. The memory budget
--local-max-block-bytes bounds decoded block capacity. Packing is refused if the
requested capacity exceeds the budget; reduce --local-block-variants or raise the budget.
Existing outputs are refused. No R, Python, null-model fit or P calculation.
See docs/LOCAL_PACKING.md for complete examples and input contracts.
)", stderr); return;
    }
    if (topic == "int-pheno") {
        std::fputs(R"(SPAGxE --int-pheno --pheno TABLE --pheno-name Y1,Y2 --out PREFIX
Writes inverse-normal transformed traits to PREFIX.txt. With no column list,
all non-ID columns are transformed. Missing values are retained as NA.
This mode accepts phenotype input and output options, not association inputs.
)", stderr); return;
    }
    for (const auto &option : options) if (topic == option.name) { printOption(option); return; }
    std::fprintf(stderr, "Unknown help topic: '%s'\n\n", rawTopic.c_str());
    throw std::invalid_argument("Unknown help topic: " + rawTopic);
}
} // namespace cli
