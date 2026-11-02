# Preparing local-ancestry G/H input

SPAGxEmix_local analyzes ancestry-specific allele counts or dosages conditional
on known ancestry chromosome-copy counts. Analysis and hardcall packing run in
the standalone C++ executable. The prepared examples are included in the source
and Windows bundle; running them requires no R or Python.

## What must be known upstream

For person i, variant j and ancestry k, H is the number of the person's two
chromosomes assigned to ancestry k. G counts the ALT alleles on those H
chromosomes. For resolved haplotypes with ALT indicators a1,a2 and ancestry
labels L1,L2:

```text
H_k = I(L1 = k) + I(L2 = k)
G_k = I(L1 = k) * a1 + I(L2 = k) * a2
```

Each allele indicator is 0 or 1, and each chromosome's ancestry label must
be resolved and aligned to that allele at the same variant. Every complete
person/variant has sum(H_k)=2 and 0<=G_k<=H_k. H=0 implies structural G=0.
Supply all ancestries represented in the cohort, not only the ancestry being
tested. String IIDs and allele metadata must be consistent across files.

Prepare these counts from your validated ancestry-resolved genotype output.
The program does not infer ancestry or supply a universal parser for every
ancestry-inference tool. Global diploid G and ancestry copy counts H alone do
not generally identify each ancestry-specific G: a heterozygote with one copy
of each ancestry still needs the ALT allele assigned to the correct chromosome.
Continuous ancestry probabilities are not known H and must not be rounded into
this model. See the [local format contract](LOCAL_ANCESTRY_FORMAT.md).

## Choose the input representation

| Available data | Supported route | Native conversion/generation |
|---|---|---|
| Integer G/H hardcalls | Paired text, or packed count triplet | `--pack-local` creates `.localbed/.localbim/.localfam` |
| Fractional G dosage with known integer H | Paired text directly | Hardcall packing refuses fractional G |
| Genotype probabilities conditional on known H | Supported local BGEN directly | Reader is available; no general C++ BGEN export command is provided |
| Existing packed hardcalls | Packed count triplet directly | `--pack-local --local-prefix` can repack blocks |
| Point-mass BGEN with known H | BGEN directly or packed counts | `--pack-local --local-bgen-manifest` can pack point calls |
| Uncertain ancestry H or unresolved ancestry-specific alleles | Outside this interface | No ancestry inference or hard-calling of uncertain H |

Packed counts and local BGEN preserve different information. A posterior may
have an integer expectation without being a point mass. The packer rejects
that uncertain posterior rather than silently discarding uncertainty. A dosage
alone does not specify a unique probability distribution, so no posterior is
invented from a dosage for BGEN export.

## Minimal paired-text schema

Use actual TAB characters. The following two-person illustration describes one
variant; use the supplied cohort example below for runnable analysis.

`Ancestry1.G.tsv`:

```text
CHROM	POS	ID	REF	ALT	001	002
1	101	rs1	C	A	1	1
```

`Ancestry1.H.tsv`:

```text
CHROM	POS	ID	REF	ALT	001	002
1	101	rs1	C	A	2	1
```

`Ancestry2.G.tsv`:

```text
CHROM	POS	ID	REF	ALT	001	002
1	101	rs1	C	A	0	0
```

`Ancestry2.H.tsv`:

```text
CHROM	POS	ID	REF	ALT	001	002
1	101	rs1	C	A	0	1
```

The four-column manifest `cohort.local.tsv` is:

```text
ancestry	genotype	hapcount	counted_allele
Ancestry1	Ancestry1.G.tsv	Ancestry1.H.tsv	ALT
Ancestry2	Ancestry2.G.tsv	Ancestry2.H.tsv	ALT
```

Paths are relative to the manifest directory. Every stream must contain the
same complete IID set and matching variant metadata/order; IID order may
differ because it is aligned explicitly. REF-counted text is supported with
a consistent `REF` declaration and is canonicalized by H-G. For new files,
keeping canonical ALT direction throughout preparation is the clearest route.

G may be missing (`NA`, `.`, `-9` or another documented token). H is never
imputed. If any ancestry H is missing, that person is excluded from every
ancestry branch at that variant. With known H, missing G is filled by H*q;
known H=0 remains structural zero. Do not replace uncertain H with a fabricated
integer or fill missing G with zero. Full rules are in the
[format specification](LOCAL_ANCESTRY_FORMAT.md).

## Runnable C++ packing and analysis

Run from `cpp/` after compiling the source. On Windows use
`build\bin\SPAGxE.exe`; in the complete Windows Release bundle the executable
is `SPAGxE.exe` inside `cpp/`. The supplied
[SH](../cpp/examples/tutorials/06_local_packing.sh) and
[CMD](../cpp/examples/tutorials/06_local_packing.cmd) tutorials use the correct
shell syntax and can also be launched by `run_all`.

The existing `examples/local/binary_demo/cohort.local.tsv` points to complete
integer G/H streams for a simulated cohort. First create a new output prefix:

```sh
mkdir -p results/local_preparation
build/bin/SPAGxE --pack-local \
  --local-manifest examples/local/binary_demo/cohort.local.tsv \
  --out results/local_preparation/cohort
```

This creates:

```text
results/local_preparation/cohort.localbed
results/local_preparation/cohort.localbim
results/local_preparation/cohort.localfam
```

Read the newly generated triplet in the native local association workflow:

```sh
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-prefix results/local_preparation/cohort \
  --pheno examples/local/binary_demo/samples.pheno \
  --resid-name R --envir-name E --mac 0 --maf 0 --geno 1 \
  --out results/local_preparation/packed_analysis
```

The same input can be analyzed directly as text:

```sh
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-manifest examples/local/binary_demo/cohort.local.tsv \
  --pheno examples/local/binary_demo/samples.pheno \
  --resid-name R --envir-name E --mac 0 --maf 0 --geno 1 \
  --out results/local_preparation/text_analysis
```

The relaxed QC thresholds expose small-example cases; specify appropriate
real-data QC before analysis. Here R is a supplied null residual column and
E is the environmental variable. For raw-trait native fitting, follow the
[four-method tutorials](../cpp/examples/tutorials/README.md). Inspect the local
result, ancestry-combination result and `.local.settings.tsv` receipt. The
receipt identifies the guarded partial-normal policy and missingness rules.

Existing packing outputs are refused. Use a fresh prefix for a repeat, or retain
the previous files and choose a new destination. Never edit or swap sidecars
after packing: CRCs bind their exact bytes to the container.

## Reusing ordinary BIM/FAM metadata

`--local-bim original.bim --local-fam original.fam` can preserve compatible
six-column sidecars while packing. The packer validates the complete sample set,
variant order/coordinates/IDs and allele direction. Ordinary BIM columns A1/A2
are not universally ALT/REF. Do not merely rename ordinary BED or sidecars to
`.localbed/.localbim/.localfam`. Follow the
[packing and sidecar instructions](LOCAL_PACKING.md).

## Reading known-H local BGEN

One BGEN file per ancestry contains Layout-2 unphased biallelic probabilities
and embedded IIDs. Per-sample ploidy stores known H=0/1/2. Genotype probability
categories count BGEN's second listed allele; do not change that meaning by
editing only a manifest declaration. The reader checks the required collection
and receives expected G from its stored probabilities.

From `cpp/`, a supplied runnable collection is:

```sh
build/bin/SPAGxE --method SPAGxEmix_local \
  --local-bgen-manifest examples/local/bgen_demo/bits16_zlib.localbgen.tsv \
  --pheno examples/local/bgen_demo/samples.pheno --resid-name R \
  --envir-name E --mac 0 --maf 0 --geno 1 \
  --out results/local_preparation/bgen_analysis
```

This uses conditional posterior mean G with the Bin(H,q) working law and
guarded partial-normal SPA. It does not integrate arbitrary joint
genotype/ancestry uncertainty into an exact test. Read the
[BGEN profile](LOCAL_BGEN_FORMAT.md) before preparing real BGEN collections.

## Local frequency, MAC and missing G

After phenotype/sample filters and the common complete-H mask, q is the
observed ALT sum divided by the H sum for people with observed G. Define
`minor_q=min(q,1-q)`. The reported/filtering local MAC is
`minor_q * sum(H)` across all known H in that mask, including copies whose G
will be mean-imputed. It is the mean-imputed count convention, not the
observed-only minor count `minor_q * sum(H for observed G)`.

Consequently, missing G can make those two MAC definitions differ. H=0 adds no
copies to either definition. `--geno` measures missing G after the common H
mask; it is not a maximum fraction of unknown H. Outputs report `N_H_MISSING`
separately. Ancestry-specific testing can flip to REF when ALT q>0.5, while
`ALT_FREQ` retains the canonical ALT frequency. Use `TESTED_ALLELE` and
`MINOR_FLIPPED` for signed-statistic interpretation.

## Optional developer regeneration

All required examples are already shipped. A developer with Python can
optionally regenerate independent simulated fixtures from the source checkout:

```sh
python -B tests/make_local_fixtures.py
python -B tests/make_local_bgen_fixtures.py
```

These scripts write their documented example fixtures, so use a disposable
checkout or preserve local edits first. The BGEN fixture script simulates
posteriors and derives reference G from quantized stored probabilities. It is
not a general real-data BGEN exporter. Python is not needed for the native
packing, association commands or shipped tutorials above.
