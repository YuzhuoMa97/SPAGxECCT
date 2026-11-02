#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/../.."
exe="${SPAGXE_BIN:-build/bin/SPAGxE}"
[ ! -f "$exe.exe" ] || exe="$exe.exe"
mkdir -p results/tutorials
# SPAGxEmix_local: all local examples are independently generated SPAGxE simulations.
# Paired G/H tables; native four-family null models and Wald/CCT. Relaxed QC exposes tiny-example edge cases.
"$exe" --method SPAGxEmix_local --local-manifest examples/local/cohort.local.tsv --pheno examples/local/samples.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 --threads 4 --out results/tutorials/local_text
# Packed hardcalls, H=0/1/2. External residuals have no prospective Wald leg.
"$exe" --method SPAGxEmix_local --local-prefix examples/local/binary_demo/local --pheno examples/local/binary_demo/samples.pheno --resid-name R --envir-name E --mac 0 --maf 0 --geno 1 --out results/tutorials/local_binary
# Known-H BGEN layout 2 per ancestry; use expected G from stored posteriors. No R runtime.
"$exe" --method SPAGxEmix_local --local-bgen-manifest examples/local/bgen_demo/bits16_zlib.localbgen.tsv --pheno examples/local/samples.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 --threads 4 --out results/tutorials/local_bgen
