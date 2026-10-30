#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/../.."
exe="${SPAGXE_BIN:-build/bin/SPAGxE}"
[ ! -f "$exe.exe" ] || exe="$exe.exe"
mkdir -p results/tutorials
# SPAGxE: unrelated subjects; paper weights with guarded partial-normal SPA.
# Fit four traits from PLINK BED. Include E in the null-model covariates.
"$exe" --method SPAGxE --bfile examples/synthetic/synthetic --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name X1,E,PC1,PC2 --envir-name E --out results/tutorials/spagxe
