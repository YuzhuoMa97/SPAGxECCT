#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/../.."
exe="${SPAGXE_BIN:-build/bin/SPAGxE}"
[ ! -f "$exe.exe" ] || exe="$exe.exe"
mkdir -p results/tutorials
# Convert integer G/H in C++; no ancestry inference or P-value calculation.
# The three original-format sidecars are validated; exact bytes are retained.
# Existing outputs are refused, so run the tutorial in a fresh results folder.
"$exe" --pack-local --local-manifest examples/local/binary_demo/cohort.local.tsv --local-bim examples/local/binary_demo/local.localbim --local-fam examples/local/binary_demo/local.localfam --out results/tutorials/packed_demo
# Read the newly created counts using the same published local association kernel.
"$exe" --method SPAGxEmix_local --local-prefix results/tutorials/packed_demo --pheno examples/local/binary_demo/samples.pheno --resid-name R --envir-name E --mac 0 --maf 0 --geno 1 --out results/tutorials/local_packed
