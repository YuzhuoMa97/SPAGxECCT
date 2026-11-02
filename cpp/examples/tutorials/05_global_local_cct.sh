#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/../.."
exe="${SPAGXE_BIN:-build/bin/SPAGxE}"
[ ! -f "$exe.exe" ] || exe="$exe.exe"
mkdir -p results/tutorials
# Optional CCT; run tutorial 04 first. Match metadata, phenotype and environment.
# Matched total G is a toy sum of local dosages. Real admixed analyses require suitable SNP PCs.
"$exe" --vcf examples/local/global.vcf --pheno examples/local/samples.pheno --pheno-name Quantitative --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 --out results/tutorials/global
# Dynamic weights give one nominal weight to each valid ancestry and one to global.
"$exe" --method SPAGxEcct --cct-global results/tutorials/global.Quantitative.SPAGxE --cct-local results/tutorials/local_text.Quantitative.SPAGxE_local_CCT --cct-phenotype Quantitative --envir-name E --cct-per-ancestry-equal --out results/tutorials/global_local
