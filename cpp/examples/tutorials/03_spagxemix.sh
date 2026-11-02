#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/../.."
exe="${SPAGXE_BIN:-build/bin/SPAGxE}"
[ ! -f "$exe.exe" ] || exe="$exe.exe"
mkdir -p results/tutorials
# SPAGxEmix: individual AF from genotype PCs.
# The same PCs enter AF and null models. BGEN allele orientation is explicit.
"$exe" --method SPAGxEmix --bgen examples/synthetic/synthetic.bgen ref-first --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name X1,E,PC1,PC2 --envir-name E --pc-cols PC1,PC2 --threads 4 --out results/tutorials/mix
