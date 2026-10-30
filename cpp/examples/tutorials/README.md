# Annotated standalone C++ examples

The six focused tutorials and the complete feature walkthrough use the same
SPAGxE executable and the same numerical policy. They do not invoke R or Python.
Start with the [official user guide](../../../docs/USER_GUIDE.md).

| Scripts (SH and CMD) | Purpose |
|---|---|
| [01_spagxe.sh](01_spagxe.sh) / [CMD](01_spagxe.cmd) | SPAGxE: PLINK BED, native quantitative/binary/survival/ordinal null models |
| [02_spagxe_plus.sh](02_spagxe_plus.sh) / [CMD](02_spagxe_plus.cmd) | SPAGxE+: PGEN and sparse GRM, four trait families |
| [03_spagxemix.sh](03_spagxemix.sh) / [CMD](03_spagxemix.cmd) | SPAGxEmix: BGEN and PCs, individual allele frequencies |
| [04_spagxemix_local.sh](04_spagxemix_local.sh) / [CMD](04_spagxemix_local.cmd) | SPAGxEmix_local: text G/H, count binary and known-H BGEN |
| [05_global_local_cct.sh](05_global_local_cct.sh) / [CMD](05_global_local_cct.cmd) | Global/local CCT and explicit ancestry grouping |
| [06_local_packing.sh](06_local_packing.sh) / [CMD](06_local_packing.cmd) | Native C++ count packing and packed-input analysis |
| [07_feature_walkthrough.sh](07_feature_walkthrough.sh) / [CMD](07_feature_walkthrough.cmd) | All 33 runnable guide groups and all 53 canonical CLI parameters |

## Working directory and executable

For a source checkout, first build inside `cpp/`; the executable is
`build/bin/SPAGxE` on POSIX or `build/bin/SPAGxE.exe` on Windows. For the Windows
binary Release, use `cpp/SPAGxE.exe` without a compiler. The six focused scripts default to the source-build executable. The binary
Release's root launcher sets their executable path automatically. The optional
07 walkthrough locates either layout. `SPAGXE_BIN` (SH) or `SPAGXE_EXE` (CMD)
can override the selected path.
Do not copy the executable away from the accompanying examples and licenses.

## Six focused tutorials

From a built source checkout, inside `cpp/` on POSIX:

```sh
bash examples/tutorials/run_all.sh
```

From a built source checkout, inside `cpp/` on Windows:

```powershell
.\examples\tutorials\run_all.cmd
```

For a Windows binary Release, start at the extracted bundle root instead:

```powershell
.\run_examples.cmd
```

Both routes write under `cpp/results/tutorials`. The root launcher belongs
to the binary Release, and does not belong to the source checkout. Use a fresh output location when repeating
packing: `--pack-local` refuses existing output triplets.

## Complete feature walkthrough

This optional, separate tutorial executes the same commands as the marked
UG01 through UG33 groups in the guide, in document order. It is not part of
`run_all`, so a quick first run remains small. Default output is
`cpp/results/official_walkthrough`. Use a new or empty directory for each run.

```sh
# From cpp/. The optional first argument selects a fresh output directory.
bash examples/tutorials/07_feature_walkthrough.sh results/my_walkthrough
```

```powershell
# From cpp/. Windows uses built-in PowerShell for file-list preparation.
.\examples\tutorials\07_feature_walkthrough.cmd -OutputDirectory results/my_walkthrough
```

The Windows CMD launcher calls the adjacent
[PowerShell script](07_feature_walkthrough.ps1). The POSIX script uses standard
`awk` and `printf` for four small selection lists; all analysis, AF export, INT,
CCT and local count packing are native C++ operations.

The walkthrough includes PLINK BED/PGEN, BGEN with and without embedded IDs,
VCF/BCF, all four methods, all native trait families, input/QC filters, compressed
output and residual reuse, local conversion, CCT weights, and resource controls.
It stops on a failed command. Check result status codes before interpreting P.

## Numerical interpretation

No statistical profile selection is required. All four methods use the paper's
weights and guarded partial-normal SPA. Central terms use matched Gaussian
moments; inaccurate partial candidates fall back to full-CGF evaluation.
Environment coding is retained as supplied. Branch-B SPA need not be invariant
to shifting E's zero point. Environmental main effects enter the null model.
Tiny examples relax QC to illustrate edge-case handling; choose justified QC
thresholds for real data. Read the guide's output/status and numerical sections.
