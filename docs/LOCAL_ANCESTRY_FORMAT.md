# Native local ancestry G/H formats

SPAGxE local analysis runs entirely in C++. SPAGxE methods and release maintainer: **Yuzhuo Ma <yuzhuoma@amss.ac.cn>**. This distribution independently implements its readers and generates its public examples from simulation. Existing dependency licenses remain applicable.

## G, H and allele direction

For a person, marker and ancestry, **H** is a known integer chromosome-copy count: 0, 1 or 2. **G** counts an allele on those H chromosomes, with `0 <= G <= H`. Supply every ancestry, so each complete person/marker has `sum(H)=2`. H is not a continuous ancestry probability. A global genotype is `sum(G)` when all counts and directions are available.

Hardcall G is integer. Paired text also supports finite fractional G in [0,H]. Local BGEN supplies the expected G from genotype probabilities conditional on known H. The association CGF uses the working `Bin(H,q)` model; it does not integrate arbitrary genotype/ancestry posterior uncertainty. Ancestry inference, uncertain H and non-diploid total ancestry counts are outside this input model.

For REF-counted text or packed counts, canonicalization to ALT is **H-G**. All ancestry streams must declare the same source counted direction. For BGEN, categories always count the second listed allele, irrespective of an ALT/REF interchange declaration in its manifest. The direct BGEN reader retains that canonical ALT expectation rather than applying a numerically lossy double complement; changing only the manifest declaration does not change stored allele order or posterior meaning.

A real allele-order reversal must exchange both REF/ALT metadata and posterior categories: new G corresponds to old H-G. Every ancestry file must use the same resulting metadata; mixed allele orders are rejected. The reader does not silently swap or strand-correct one ancestry to match another.

The kernel subsequently selects the ancestry-specific minor allele using H-G when canonical ALT frequency exceeds 0.5. Outputs preserve original `ALT_FREQ` and identify the tested direction through `TESTED_ALLELE` (`ALT`/`REF`) and `MINOR_FLIPPED`. Signed statistics refer to the tested allele. Exactly 0.5 uses the ALT tie convention.

## Paired text G/H streams

Use TAB delimiters, one marker per row, and these columns in both G and H:

```text
CHROM  POS  ID   REF  ALT  IID1  IID2  IID3 ...
1      101  rs1  C    A    0     1     0    ...
```

Displayed spaces represent columns; actual delimiters must be TAB. `CHR` is an alias of `CHROM`. The five metadata columns are required in that order. POS must be a positive integer no larger than 2^32-1; REF and ALT must be distinct biallelic allele labels. String IIDs retain leading zeros. Duplicate IIDs, ancestry labels and marker IDs are rejected.

The four-column SPAGxE manifest is:

```text
ancestry   genotype         hapcount         counted_allele
Ancestry1  ancestry1.G.tsv   ancestry1.H.tsv   ALT
Ancestry2  ancestry2.G.tsv   ancestry2.H.tsv   ALT
Ancestry3  ancestry3.G.tsv   ancestry3.H.tsv   ALT
```

Pass it through `--local-manifest FILE`. Paths resolve relative to the manifest directory unless absolute. Native `.gz` and `.zst` streams are supported. Ancestry labels must be unique and nonempty, contain letters/digits/underscores/hyphens, and not equal the reserved label `global`.

Every G/H file must have the same IID set, but each file may order samples independently: the reader aligns by IID. Marker count, order and CHROM/POS/ID/REF/ALT must agree exactly. There is no silent intersection or marker reordering. Global G and ancestry H alone do not determine ancestry-specific G without resolved ancestry alleles.

Missing tokens are an empty field, `NA`, `.`, `-9`, `NaN`, `nan` or `NAN`. Known H=0 restores missing G to structural zero. Otherwise, q is estimated from observed G/H and missing G is filled with H*q for that ancestry. H is never imputed. If any ancestry H is unknown, that person is excluded from **all ancestry branches at that marker**. Outputs report `N_H_MISSING` and `N_TESTED`. Phenotype/residual masks and keep/remove filters precede this marker-specific mask. Every fully known ancestry-H collection must sum to two.

## Native BGEN probabilities

`--local-bgen-manifest FILE` accepts a separate three-column schema:

```text
ancestry   bgen              counted_allele
Ancestry1  ancestry1.bgen    ALT
Ancestry2  ancestry2.bgen    ALT
Ancestry3  ancestry3.bgen    ALT
```

The same header is auto-detected under `--local-manifest`. Choose one local input option per run. The supported profile is BGEN Layout 2, unphased biallelic probabilities, embedded IIDs, known per-ancestry H=0/1/2, raw/zlib/zstd blocks and 1..32-bit precision. All ancestry files must agree on marker metadata/order and contain the same IID set. Sample order may differ.

External-ID-only BGEN, phased/multiallelic probabilities and uncertain H are rejected. The reader streams expected G to the native kernel; it does not turn the full posterior into a new exact test. See [LOCAL_BGEN_FORMAT.md](LOCAL_BGEN_FORMAT.md) for the probability encoding, checks, missingness and allele-interchange semantics.

## Packed hardcall count container

`--local-prefix PREFIX` (alias `--local-binary PREFIX`) reads `.localbed/.localbim/.localfam` directly. This is a G/H count format, not ordinary PLINK BED. Do not rename a different binary format to these extensions. The `binary_demo` example is independently generated by this project. The native `--pack-local` utility creates these files from integer G/H; see [packing and original BIM/FAM reuse](LOCAL_PACKING.md).

All integers are little endian:

| Offset | Bytes | Meaning |
| --- | --- | --- |
| 0 | 8 | ASCII magic `SQRLOC01` |
| 8 | 8 | uint64 sample count N |
| 16 | 8 | uint64 marker count M |
| 24 | 4 | uint32 ancestry count K |
| 28 | 4 | uint32 maximum rows per block B |
| 32 | 8 | uint64 absolute footer offset |
| 40 | 1 | source counted allele: 0 ALT, 1 REF |
| 41 | variable | K ancestry labels, each uint16 byte length followed by label bytes |

Each block stores uint32 row count, uint64 decoded byte count, uint64 compressed byte count, uint32 decoded CRC32 and an independent zstd frame. Each marker contains planes `G1,H1,G2,H2,...,GK,HK`. A plane uses `ceil(N/4)` bytes. Two-bit codes 0/1/2 are counts; code 3 is missing. The first sample occupies the least significant pair; unused padding bits are zero. Decoded length is exactly `rows * 2*K*ceil(N/4)`. Fractional G cannot be represented.

The footer contains uint64 block count, one uint64 absolute offset per block, uint32 CRC32 of the complete `.localbim`, and uint32 CRC32 of the complete `.localfam`. The footer ends exactly at EOF. Sidecars are headerless:

- `.localbim`: `CHROM ID CM POS ALT REF`, one marker per row.
- `.localfam`: `FID IID PAT MAT SEX PHENO`, one sample per row. Unique IID is the matching key; PHENO here is not the analysis phenotype.

Checks include magic, dimensions/overflow, file/index extents, complete compressed-frame consumption, decoded length, CRCs, padding, labels, duplicate IDs, G<=H and complete H sums. CRC32 detects accidental corruption and mismatched sidecars; it does not authenticate provenance. Keep independent hashes and rebuild matching files after any sidecar edit. Original PLINK A1/A2 is not automatically ALT/REF; the packer validates orientation and can preserve supplied sidecar bytes.

## Resource limits and QC

`--local-max-block-bytes` limits a decoded binary or BGEN block before allocation: default 268435456 bytes, permitted 4096..1073741824. Repack oversized blocks or explicitly raise the budget. This is not a total process-memory limit. Compressed buffers, identifiers, alignment maps, marker indexes, phenotype arrays and worker buffers also consume memory; metadata/index memory can grow with sample/marker count. Genotypes are streamed in batches.

Ancestry-specific MAC/MAF/G-missingness drive `--mac`, `--maf` and `--geno`. Failed rows retain their QC status and NA test values. Missingness is evaluated after the common complete-H mask; `--geno` is not a threshold on unknown-H frequency. MAC is the ancestry minor frequency times all known H copies after masking, consistent with mean-imputed G.

`--hwe` is a **global diploid hardcall diagnostic** based on observed sums of local G, not ordinary diploid HWE within a variable-H ancestry branch. Fractional text G or non-point BGEN posteriors make this diagnostic unavailable; this includes an uncertain posterior whose expected G happens to be integer. Default `--hwe 0` permits the marker; positive HWE filtering reports `FAIL_HWE_UNAVAILABLE`.

This local interface does not infer ancestry, ingest arbitrary joint ancestry/allele posteriors or provide a validated local+GRM covariance model. Global SPAGxE+ provides the published homogeneous related-sample method.

## Reproducible examples and numerical scope

`cpp/tests/make_local_fixtures.py` creates the paired-text examples and the `cpp/examples/local/binary_demo/` packed-count example. `cpp/tests/make_local_bgen_fixtures.py` independently simulates and encodes `cpp/examples/local/bgen_demo/`; expected text G is calculated from stored quantized probabilities, rather than exported by the C++ decoder.

`verify_local.py`, `verify_local_absolute.py` and `verify_local_bgen.py` check the listed formats, sample alignment, real allele-order reversals, masks, four trait families, residual scans, threads, budgets and malformed inputs. Independent R/formula references are development checks; runtime analysis requires no R.

Original R agreement is claimed only in the valid shared additive domain with complete finite inputs, supported frequency direction, nondegenerate designs and matched residual/null-model settings. Native H-G/H*q corrections and explicit missing-H policies are validated separately. Passing format and numerical examples is not a claim of universal calibration across all ancestry/imputation methods or equality with every original R branch.