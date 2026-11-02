#!/usr/bin/env bash
# Complete feature walkthrough. Run in a fresh output directory.
# No R or Python runtime is used; list preparation uses awk and printf.
set -euo pipefail
CPP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$CPP_DIR"
SPAGXE="${SPAGXE_EXE:-${SPAGXE_BIN:-}}"
if [[ -z "$SPAGXE" ]]; then
    for candidate in "$CPP_DIR/SPAGxE" "$CPP_DIR/SPAGxE.exe" "$CPP_DIR/build/bin/SPAGxE" "$CPP_DIR/build/bin/SPAGxE.exe"; do
        if [[ -f "$candidate" ]]; then SPAGXE="$candidate"; break; fi
    done
fi
if command -v cygpath >/dev/null 2>&1; then SPAGXE="$(cygpath -u "$SPAGXE")"; fi
if [[ ! -f "$SPAGXE" ]]; then
    printf 'Set SPAGXE_EXE to a built SPAGxE executable.\n' >&2
    exit 1
fi
OUT="${1:-results/official_walkthrough}"
if [[ -e "$OUT" && ( ! -d "$OUT" || -n "$(ls -A "$OUT")" ) ]]; then
    printf 'Use a new or empty output directory: %s\n' "$OUT" >&2
    exit 1
fi
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"


# UG01_help: Executable version, general help, all options, phenotype and output help
printf '\n%s\n' 'UG01_help'
mkdir -p "$OUT"
"$SPAGXE" --version
"$SPAGXE" --help
"$SPAGXE" --help options
"$SPAGXE" --help phenotype
"$SPAGXE" --help output

# UG02_base: SPAGxE native quantitative null, PLINK BED, residual saving
printf '\n%s\n' 'UG02_base'
"$SPAGXE" --method SPAGxE --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E --regression-model linear \
  --save-resid --mac 0 --maf 0 --geno 1 --out "$OUT/base"

# UG03_plus_table: SPAGxE+ shared-AF GRM correction using ID1 ID2 VALUE table
printf '\n%s\n' 'UG03_plus_table'
"$SPAGXE" --method SPAGxE+ --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --grm-file examples/synthetic/synthetic.grm.tsv \
  --mac 0 --maf 0 --geno 1 --out "$OUT/plus_table"

# UG04_mix: SPAGxEmix individual-specific AF, ancestry PCs, embedded-ID BGEN ref-first
printf '\n%s\n' 'UG04_mix'
"$SPAGXE" --method SPAGxEmix \
  --bgen examples/synthetic/synthetic.bgen ref-first \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E --pc-cols PC1,PC2 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/mix"

# UG05_local_text: SPAGxEmix_local paired G/H, four native families, across-ancestry CCT
printf '\n%s\n' 'UG05_local_text'
"$SPAGXE" --method SPAGxEmix_local \
  --local-manifest examples/local/cohort.local.tsv \
  --pheno examples/local/samples.pheno \
  --pheno-name Quantitative,Binary,Time:Event,Ordinal \
  --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 \
  --out "$OUT/local_text"

# UG06_bgen_ref_last: Global BGEN last-listed REF orientation
printf '\n%s\n' 'UG06_bgen_ref_last'
"$SPAGXE" --bgen examples/synthetic/synthetic.reflast.bgen ref-last \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/bgen_ref_last"

# UG33_bgen_external_ids: Global BGEN without embedded IIDs, companion sample resolution
printf '\n%s\n' 'UG33_bgen_external_ids'
"$SPAGXE" --bgen examples/synthetic/synthetic.noids.bgen ref-first \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/bgen_external_ids"

# UG07_vcf: VCF global reader
printf '\n%s\n' 'UG07_vcf'
"$SPAGXE" --vcf examples/synthetic/synthetic.vcf \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/vcf"

# UG08_bcf: BCF2 global reader
printf '\n%s\n' 'UG08_bcf'
"$SPAGXE" --bcf examples/synthetic/synthetic.bcf \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/bcf"

# UG09_pgen: PGEN/PVAR/PSAM global reader with matched 1kg phenotype fixture
printf '\n%s\n' 'UG09_pgen'
"$SPAGXE" --pfile examples/1kg --pheno examples/1kg.pheno \
  --pheno-name Quantitative --covar-name MALE,PC1,PC2,PC3,PC4 \
  --envir-name MALE --out "$OUT/pgen"

# UG10_plus_plink2: SPAGxE+ PLINK2 sparse GRM plus grm.id
printf '\n%s\n' 'UG10_plus_plink2'
"$SPAGXE" --method SPAGxE+ --pfile examples/1kg \
  --pheno examples/1kg.pheno --pheno-name Quantitative \
  --covar-name MALE,PC1,PC2,PC3,PC4 --envir-name MALE \
  --sp-grm-plink2 examples/1kg.grm.sp --out "$OUT/plus_plink2"

# UG11_trait_models: Explicit native logistic, Breslow Cox and ordinal fitting
printf '\n%s\n' 'UG11_trait_models'
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

# UG12_multiple: Multiple raw traits and two environments with automatic families
printf '\n%s\n' 'UG12_multiple'
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno \
  --pheno-name Quantitative,Binary,Time:Event,Ordinal --regression-model auto \
  --covar-name X1,E,PC1,PC2 --envir-name E,X1 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/multiple"

# UG13_external_residuals: Several precomputed residual columns
printf '\n%s\n' 'UG13_external_residuals'
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno \
  --resid-name R_quant,R_binary,R_survival,R_ordinal --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/external_residuals"

# UG14_reload_residual: Native saved-residual reload with separate environment table
printf '\n%s\n' 'UG14_reload_residual'
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno "$OUT/base.null.resid" --resid-name Quantitative \
  --covar examples/synthetic/synthetic.pheno --covar-name E --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/reloaded_residual"

# UG15_prepare_filters: Create IID and marker include/exclude lists
printf '\n%s\n' 'UG15_prepare_filters'
awk 'NR<=401 {print $1}' examples/synthetic/synthetic.pheno > "$OUT/keep.txt"
printf 'IID\nS0001\nS0002\n' > "$OUT/remove.txt"
printf 'rsOracle01\nrsOracle02\nrsOracle03\nrsOracle04\n' > "$OUT/extract.txt"
printf 'rsOracle02\n' > "$OUT/exclude.txt"

# UG16_filter_qc: Combined keep/remove/extract/exclude/chromosome selection and all QC flags
printf '\n%s\n' 'UG16_filter_qc'
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --keep "$OUT/keep.txt" --remove "$OUT/remove.txt" \
  --extract "$OUT/extract.txt" --exclude "$OUT/exclude.txt" --chr 1 \
  --geno 0.1 --maf 0.01 --mac 10 --hwe 0.000001 \
  --hard-call-threshold 0.1 --out "$OUT/filtered"

# UG17_local_existing_counts: Existing native packed G/H with known H and supplied residual
printf '\n%s\n' 'UG17_local_existing_counts'
"$SPAGXE" --method SPAGxEmix_local \
  --local-prefix examples/local/binary_demo/local \
  --pheno examples/local/binary_demo/samples.pheno --resid-name R \
  --envir-name E --mac 0 --maf 0 --geno 1 --out "$OUT/local_existing"

# UG18_local_bgen: Local known-H BGEN layout2 zlib/16-bit profile
printf '\n%s\n' 'UG18_local_bgen'
"$SPAGXE" --method SPAGxEmix_local \
  --local-bgen-manifest examples/local/bgen_demo/bits16_zlib.localbgen.tsv \
  --pheno examples/local/bgen_demo/samples.pheno --resid-name R \
  --envir-name E --mac 0 --maf 0 --geno 1 --out "$OUT/local_bgen"

# UG19_pack_local: Native C++ integer G/H packer, original BIM/FAM reuse, block/level/budget controls
printf '\n%s\n' 'UG19_pack_local'
"$SPAGXE" --pack-local \
  --local-manifest examples/local/binary_demo/cohort.local.tsv \
  --local-bim examples/local/binary_demo/local.localbim \
  --local-fam examples/local/binary_demo/local.localfam \
  --local-block-variants 2 --local-zstd-level 5 \
  --local-max-block-bytes 268435456 --out "$OUT/local_counts"

# UG20_local_new_counts: Analyze newly generated native G/H triplet
printf '\n%s\n' 'UG20_local_new_counts'
"$SPAGXE" --method SPAGxEmix_local --local-prefix "$OUT/local_counts" \
  --pheno examples/local/binary_demo/samples.pheno --resid-name R \
  --envir-name E --mac 0 --maf 0 --geno 1 --out "$OUT/local_new"

# UG21_repack_local: Repack existing native G/H triplet with new block size and default zstd level
printf '\n%s\n' 'UG21_repack_local'
"$SPAGXE" --pack-local --local-prefix "$OUT/local_counts" \
  --local-block-variants 1 --local-zstd-level 0 \
  --local-max-block-bytes 268435456 --out "$OUT/local_repacked"

# UG22_cct_global_input: Global scan matching local fixture variant catalog for CCT
printf '\n%s\n' 'UG22_cct_global_input'
"$SPAGXE" --vcf examples/local/global.vcf \
  --pheno examples/local/samples.pheno --pheno-name Quantitative \
  --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 \
  --out "$OUT/cct_global"

# UG23_cct_default: Global/local CCT with default equal group weights
printf '\n%s\n' 'UG23_cct_default'
"$SPAGXE" --method SPAGxEcct \
  --cct-global "$OUT/cct_global.Quantitative.SPAGxE" \
  --cct-local "$OUT/local_text.Quantitative.SPAGxE_local_CCT" \
  --cct-phenotype Quantitative --envir-name E --out "$OUT/cct_default"

# UG24_cct_weights: Global/local CCT custom positive normalized group weights
printf '\n%s\n' 'UG24_cct_weights'
"$SPAGXE" --method SPAGxEcct \
  --cct-global "$OUT/cct_global.Quantitative.SPAGxE" \
  --cct-local "$OUT/local_text.Quantitative.SPAGxE_local_CCT" \
  --cct-phenotype Quantitative --envir-name E --cct-weights 0.6,0.4 \
  --out "$OUT/cct_weighted"

# UG25_cct_components: Global/local CCT equal nominal weight per valid component
printf '\n%s\n' 'UG25_cct_components'
"$SPAGXE" --method SPAGxEcct \
  --cct-global "$OUT/cct_global.Quantitative.SPAGxE" \
  --cct-local "$OUT/local_text.Quantitative.SPAGxE_local_CCT" \
  --cct-phenotype Quantitative --envir-name E --cct-per-ancestry-equal \
  --out "$OUT/cct_equal_components"

# UG26_af_export: Individual-AF export and gzip codec with default compression level
printf '\n%s\n' 'UG26_af_export'
"$SPAGXE" --cal-af-coef --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pc-cols PC1,PC2 \
  --mac 0 --maf 0 --geno 1 --hwe 0 --threads 2 \
  --compression gz --out "$OUT/af_models"

# UG27_int: Native inverse-normal transform of one selected quantitative phenotype
printf '\n%s\n' 'UG27_int'
"$SPAGXE" --int-pheno --pheno examples/synthetic/synthetic.pheno \
  --pheno-name Quantitative --out "$OUT/int_quantitative"

# UG28_separate_covariates: Analyze transformed phenotype with separately supplied original covariates
printf '\n%s\n' 'UG28_separate_covariates'
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno "$OUT/int_quantitative.txt" --pheno-name Quantitative \
  --covar examples/synthetic/synthetic.pheno --covar-name X1,E,PC1,PC2 \
  --envir-name E --mac 0 --maf 0 --geno 1 --out "$OUT/int_analysis"

# UG29_gzip: Gzip association output plus saved residual output
printf '\n%s\n' 'UG29_gzip'
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E --save-resid \
  --compression gz --compression-level 6 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/gzip"

# UG30_zstd: Zstd association output and explicit level
printf '\n%s\n' 'UG30_zstd'
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --compression zst --compression-level 3 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/zstd"

# UG31_compressed_residual: Read gzip saved residual directly and separate environment table
printf '\n%s\n' 'UG31_compressed_residual'
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno "$OUT/gzip.null.resid.gz" --resid-name Quantitative \
  --covar examples/synthetic/synthetic.pheno --covar-name E --envir-name E \
  --mac 0 --maf 0 --geno 1 --out "$OUT/compressed_residual"

# UG32_resources: Threads/chunks/seed compatibility plus explicit branch and SPA thresholds
printf '\n%s\n' 'UG32_resources'
"$SPAGXE" --bfile examples/synthetic/synthetic \
  --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative \
  --covar-name X1,E,PC1,PC2 --envir-name E \
  --threads 4 --chunk-ksnp 0.5 --seed 20261029 \
  --spagxe-marginal-cutoff 0.001 --spa-z-threshold 2 \
  --mac 0 --maf 0 --geno 1 --out "$OUT/resources"
printf '\nFeature walkthrough completed: %s\n' "$OUT"
