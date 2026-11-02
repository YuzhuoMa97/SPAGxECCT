#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/../.."
exe="${SPAGXE_BIN:-build/bin/SPAGxE}"
[ ! -f "$exe.exe" ] || exe="$exe.exe"
mkdir -p results/tutorials
# SPAGxE+: sparse-GRM covariance correction; SPA-only final P.
# Use bundled public PGEN/phenotype/GRM inputs. Keep the companion .grm.id.
"$exe" --method SPAGxE+ --pfile examples/1kg --pheno examples/1kg.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name MALE,PC1,PC2,PC3,PC4 --envir-name MALE --sp-grm-plink2 examples/1kg.grm.sp --threads 4 --out results/tutorials/plus
