# Native local G/H packing and reuse of BIM/FAM

SPAGxE 1.2.0 adds a pure-C++ data utility. It does not fit a model, compute P values or infer local ancestry; all four association methods use the single documented policy with paper-defined weights and guarded partial-normal SPA.

Start with [preparing local G/H input](LOCAL_INPUT_PREPARATION.md) for upstream
data requirements, a minimal text schema and the complete text-to-packed
analysis route. The native utility generates hardcall count triplets; it is not
a general local-BGEN exporter.

From cpp after building:

```sh
mkdir -p results/packed
build/bin/SPAGxE --pack-local \
  --local-manifest examples/local/binary_demo/cohort.local.tsv \
  --out results/packed/cohort

build/bin/SPAGxE --method SPAGxEmix_local \
  --local-prefix results/packed/cohort \
  --pheno examples/local/binary_demo/samples.pheno \
  --resid-name R --envir-name E --mac 0 --maf 0 --geno 1 \
  --out results/packed/association
```

This tiny fixture relaxes QC to expose edge cases; choose real-data QC in advance. See the runnable [packing tutorial](../cpp/examples/tutorials/06_local_packing.sh).

## Outputs and supported inputs

The three outputs are PREFIX.localbed, PREFIX.localbim and PREFIX.localfam. The wire format remains SQRLOC01. Use exactly one of --local-manifest, --local-prefix, --local-bgen-manifest. Text must contain validated G/H; existing counts may be repacked to different blocks. BGEN input is accepted only when all observed posteriors are point masses and expected G is integer. Fractional text dosage or non-point posterior is rejected, even if its expectation is an integer; no silent rounding occurs.

Output stores canonical ALT counts. REF-counted text/count input is recoded as H-G. Known H=0 has structural G=0. Other explicit missing G/H remains missing. The ancestry-inference step is upstream of this utility.

Existing output files are refused. Owned temporary files are created beside the target, validated and committed without replacing existing files. A failed conversion does not overwrite a pre-existing sidecar or binary.

## Can original .bim and .fam be reused?

**Their six-column layouts can be reused, subject to validation. Their mere extensions do not establish compatible contents.** Use the native utility:

```sh
build/bin/SPAGxE --pack-local --local-manifest cohort.local.tsv \
  --local-bim original.bim --local-fam original.fam \
  --out cohort_local
```

- BIM must have the same number/order of variants as G/H and matching CHROM, ID, POS. Column5 must equal the G/H canonical ALT and column6 its REF. CM is preserved.
- Ordinary PLINK1 BIM describes A1/A2; A1 is not universally VCF ALT. Verify the upstream counted allele and metadata. The tool refuses an unverified mismatch; it does not guess strand/orientation.
- FAM must contain the same complete sample set with unique string IID in column2. FID/PAT/MAT/SEX/PHENO are preserved. A different FAM sample order is supported: binary data is aligned to that order. It does not silently subset samples.
- When provided sidecars pass validation, their original bytes and newline style are copied into the outputs. The new footer CRCs bind those exact bytes.
- Do not edit/reorder/replace sidecars after packing. Repack with the validated inputs. The FAM PHENO field is not used as the association phenotype; supply --pheno separately.

Without those options, metadata is generated from the G/H stream; FAM uses FID 0, parents 0, sex 0 and PHENO -9, while retaining IID. K/ancestry labels are stored in the binary header, not in extra BIM/FAM columns.

## Association frequency and MAC

The association workflow estimates q from observed G and its corresponding H,
then fills missing G with H*q. Its reported/filtering local MAC is
`min(q,1-q) * sum(H)` across the common complete-H sample mask, including copies
with missing G. This mean-imputed count can differ from the observed-only minor
count. Packing preserves missing G; it does not pre-impute it. See the
[frequency and missingness explanation](LOCAL_INPUT_PREPARATION.md#local-frequency-mac-and-missing-g).

## Blocks and memory

--local-block-variants defaults64, allowed1..65536; --local-zstd-level defaults3, allowed-5..22 (0 selects zstd's default). The requested decoded block capacity must fit --local-max-block-bytes; otherwise packing is refused. Reduce --local-block-variants or raise the budget explicitly. The current supported decoded-block budget defaults256MiB and is capped1GiB; it is not a total-memory bound. Repack overly large source blocks upstream before importing them. See [wire specification](LOCAL_ANCESTRY_FORMAT.md).

## Paths and platform scope

Manifests use UTF8 paths. Relative entries resolve against the manifest parent; Windows local/packing CLI preserves wide Unicode arguments. Plain/gzip/zstd local text, local count sidecars, local BGEN, phenotype headers and filter files use the Unicode-aware layer. This does not make a universal Unicode promise for every legacy global-genotype third-party reader.

The public tests independently parse planes/header/footer/CRCs, verify optional sidecar bytes and rejection/preservation behavior, and compare native text/binary association results. Runtime packing requires only the C++ executable. Python/R are development validation dependencies only.
