# Local genotype posteriors in native BGEN input

SPAGxE reads one known-haplotype-count BGEN file per ancestry and performs the association scan entirely in C++. SPAGxE methods and release maintainer: **Yuzhuo Ma <yuzhuoma@amss.ac.cn>**. The reader, format checks and simulation encoder are independently implemented in this distribution. Existing dependency licenses remain applicable.

See [preparing local G/H input](LOCAL_INPUT_PREPARATION.md) for the upstream
data requirements and the native hardcall conversion route. The current C++
application reads this BGEN profile but does not provide a general BGEN export
command. Its optional Python fixture encoder simulates validation examples;
it does not reconstruct a posterior distribution from a real-data dosage.

## Data meaning

For each person, variant and ancestry, H is a **known integer chromosome-copy count** in 0, 1 or 2. The BGEN ploidy byte stores H. Conditional genotype probabilities describe G on those H copies. Every complete ancestry collection must satisfy sum(H)=2 for each person and marker. Supply every ancestry represented in the diploid autosomal cohort; ancestry count is not fixed at two.

The stored posterior has support G=0..H. The count is always for the **second allele listed in BGEN**, labeled ALT by this profile. Its expected count is P1+2P2. For H=1, P2 is structurally zero. For H=0, G=0 is structurally known and there are no stored probability integers. A genuinely missing posterior is marked in BGEN while its H remains known; G is unavailable except for the structural H=0 case.

The association engine receives the expected G and the integer H. Its underlying conditional-binomial working law is evaluated with the
[guarded partial-normal CGF policy](PARTIAL_SPA.md). It does not integrate the full posterior uncertainty into a new exact SPA distribution. A dosage alone does not identify its underlying posterior; the synthetic examples simulate probabilities directly and export expected values from their stored quantized probabilities.

H is a known ancestry chromosome-copy count, not a continuous ancestry probability. Analysis uses the conditional posterior mean G stored in BGEN. This interface neither infers ancestry nor exactly integrates arbitrary joint genotype/ancestry uncertainty.

## Manifest and identity matching

`--local-bgen-manifest FILE` accepts exactly these three TAB-separated columns:

```text
ancestry    bgen                    counted_allele
Ancestry1   ancestry1.bgen          ALT
Ancestry2   ancestry2.bgen          ALT
Ancestry3   ancestry3.bgen          ALT
```

The spaces shown above denote columns; use actual TAB delimiters. Ancestry labels must be unique, nonempty and contain letters, digits, underscores or hyphens; `global` is reserved. Relative file paths resolve against the manifest directory. The same three-column header is also recognized automatically under `--local-manifest`; the four-column G/H manifest remains a separate accepted schema.

Every ancestry file must contain embedded string IIDs and the same complete IID set. Sample order may differ between files: the decoder aligns them by IID, preserving leading zeros. Variant count, order and CHROM/POS/ID/REF/ALT must agree exactly. Duplicate IIDs or variant IDs, mismatched alleles, incomplete ancestry collections and mixed counted directions are rejected. The decoder does not silently intersect or reorder variants.

`counted_allele` records ALT/REF count interchange direction; all rows must agree. BGEN's probability categories intrinsically count its second listed allele regardless of this declaration. The direct decoder preserves that canonical ALT expectation, avoiding an unnecessary floating-point H-(H-G) round trip when REF interchange is declared. In paired text input, where G values themselves may count REF, canonicalization uses H-G. Changing only this manifest declaration does not swap the BGEN allele metadata or posterior categories.

The local statistical test subsequently selects the ancestry-specific minor allele, as reported in `TESTED_ALLELE` and `MINOR_FLIPPED`. Interpret signed statistics using those columns, while `ALT_FREQ` remains the original canonical ALT frequency.

## Supported BGEN profile

- BGEN Layout 2, unphased biallelic variants, embedded sample IDs and known H=0/1/2.
- Raw, zlib or zstd probability blocks; precision from 1 through 32 bits.
- Per-sample missing posterior flags with known H, including H=0 structural zeros.
- Streaming association without BGI indexing, an Oxford sample sidecar, a count-container sidecar, or another analysis package at runtime.

External-sample-only BGEN, phased probabilities, more than two alleles, uncertain H, joint ancestry/allele posteriors and non-diploid total ancestry counts are rejected by this profile. The global-genotype BGEN interface has its own scope and is not used to decode haploid/zero-copy ancestry probability blocks.

For a sample with H copies, Layout 2 stores H probability integers for categories 0..H-1. The remaining category H is recovered by subtracting their sum from 2^bits-1. Integer probabilities are packed least-significant bit first. Missing samples store zero probability integers; absent/reserved flags and support are checked. The exact decoded size is `10+N+ceil(bits*sum(H)/8)` bytes.

The reader validates header and sample offsets, counts, reserved flags, IID block size, allele metadata, ploidy extrema, precision, probability sums, packed padding, exact decoded lengths and complete compressed-stream consumption. Truncated files and trailing blocks/bytes fail explicitly. BGEN does not add a SPAGxE CRC footer; zlib/zstd integrity checks and independent provenance hashes are different forms of verification.

`--local-max-block-bytes` bounds a decoded probability block before allocation. This is separate from total process memory: identifiers, alignment maps, per-ancestry G/H vectors, compressed buffers and duplicate-variant indexes also consume memory. Genotypes are streamed; metadata/index memory can still grow with sample, ancestry and variant counts.

An uncertain posterior may have an exactly integer expected count. Such a marker is still flagged as dosage-derived, so ordinary hardcall HWE is unavailable. Default `--hwe 0` permits association; positive HWE filtering reports `FAIL_HWE_UNAVAILABLE` for these rows rather than inventing hardcall genotypes from expectations.

## Runnable SPAGxE examples

From `cpp/`, create `results/local_bgen` and run:

```sh
mkdir -p results/local_bgen
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-bgen-manifest examples/local/bgen_demo/bits16_zlib.localbgen.tsv \
  --pheno examples/local/samples.pheno \
  --pheno-name Quantitative,Binary,Ordinal,Time:Event \
  --covar-name X,E --envir-name E --regression-model auto \
  --mac 0 --maf 0 --geno 1 --threads 4 --out results/local_bgen/native
```

The example QC thresholds deliberately expose small-fixture edge cases. `SPAGxElocal` remains a compatible method spelling. Windows uses `build\bin\SPAGxE.exe`; put a command on one line or use the appropriate shell's continuation syntax.

For a residual-only scan:

```sh
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-bgen-manifest examples/local/bgen_demo/bits32_zstd.localbgen.tsv \
  --pheno examples/local/bgen_demo/samples.pheno --resid-name R \
  --envir-name E --mac 0 --maf 0 --geno 1 --out results/local_bgen/residual
```

The `bgen_demo` collection is generated entirely from this project's simulation code with seed 20261025. It includes three ancestries, sample-order permutations, H=0/1/2, missing probabilities, point calls, fractional expectations, integer expectations with uncertainty and ancestry AF above 0.5. The `binary_demo` count-container examples are independently simulated by the same project and do not require access to another research project.

Optional development verification uses only Python's standard library:

```sh
python -B tests/make_local_bgen_fixtures.py
python -B tests/verify_local_bgen.py build/bin/SPAGxE results/verify_local_bgen
```

The tests compare native association output against independently encoded posterior expectations exported as paired G/H text, exercise four raw trait families and residual input, and reject malformed files. These are numerical and I/O checks for the listed cases, not population-scale calibration or resource benchmarks.
