#!/usr/bin/env python3
"""SPAGxE CLI contract checks, standard library only; no analysis dependencies.
Run from cpp/: python tests/cli_validation.py --binary build/bin/SPAGxE[.exe]
"""
import argparse
import json
import pathlib
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error("Executable does not exist: " + str(binary))
    cases = []

    def check(name, command, code, expected, stdout_only=False):
        result = subprocess.run([str(binary), *command], capture_output=True,
                                text=True, encoding="utf-8", errors="replace", timeout=20)
        output = result.stdout if stdout_only else result.stdout + result.stderr
        passed = result.returncode == code and expected in output
        cases.append({"name": name, "passed": passed, "expected_exit": code,
                      "actual_exit": result.returncode, "expected_text": expected,
                      "command_arguments": command,
                      "diagnostic": "" if passed else output[-4000:]})
        print(("PASS " if passed else "FAIL ") + name)

    check("version", ["--version"], 0, "SPAGxE 1.2.0", True)
    check("short-version", ["-V"], 0, "SPAGxE 1.2.0", True)
    check("help-plus-alias", ["--help", "plus"], 0, "SPAGxE+ sparse GRM")
    check("help-mix-alias", ["--help", "mix"], 0, "SPAGxEmix")
    check("unknown-help-fails", ["--help", "unknown-topic"], 1, "Unknown help topic")
    check("short-help-ownership", ["--help"], 0, "yuzhuoma@amss.ac.cn")
    check("help-guarded-partial-cgf", ["--help"], 0, "guarded partial-normal")
    check("help-output-provenance", ["--help", "output"], 0, "LOG10P_SPA_Gx")
    check("help-admixture", ["--help", "SPAGxEmix"], 0, "SPAGxEmix")
    check("help-grm", ["--help", "sp-grm"], 0, "SPAGxE+ sparse GRM")
    check("help-local", ["--help", "SPAGxElocal"], 0, "Binomial(H_i,q)")
    check("help-cct", ["--help", "SPAGxEcct"], 0, "--cct-per-ancestry-equal")
    check("help-packing", ["--help", "pack-local"], 0, "PREFIX.localbed")

    with tempfile.TemporaryDirectory(prefix="cli_validation_", dir=".") as directory:
        base = ["--bfile", "unused_genotype", "--pheno", "unused_phenotype",
                "--pheno-name", "Y", "--envir-name", "E", "--covar-name", "X,E",
                "--out", str(pathlib.Path(directory) / "test")]
        invalid = [
            ("unknown-option", ["--unknown-option"], "Unknown option"),
            ("unsupported-method", ["--method", "SPACox"], "--method must be SPAGxE, SPAGxE+, SPAGxEmix or SPAGxEmix_local"),
            ("threads-zero", ["--threads", "0"], "--threads must be at least 1"),
            ("threads-negative", ["--threads", "-2"], "--threads must be at least 1"),
            ("threads-fraction", ["--threads", "1.5"], "requires an integer"),
            ("threads-overflow", ["--threads", "9999999999999"], "out of range"),
            ("nan-threshold", ["--maf", "NaN"], "finite numeric"),
            ("inf-threshold", ["--maf", "Inf"], "finite numeric"),
            ("trailing-numeric-text", ["--maf", "0.1abc"], "finite numeric"),
            ("maf-negative", ["--maf", "-0.01"], "--maf must be in [0,0.5]"),
            ("maf-large", ["--maf", "0.51"], "--maf must be in [0,0.5]"),
            ("missing-rate-large", ["--geno", "1.1"], "--geno must be in [0,1]"),
            ("mac-negative", ["--mac", "-1"], "--mac must be non-negative"),
            ("hwe-large", ["--hwe", "1.1"], "--hwe must be in [0,1]"),
            ("hard-call-boundary", ["--hard-call-threshold", "0.5"], "--hard-call-threshold must be in [0,0.5)"),
            ("marginal-cutoff-large", ["--spagxe-marginal-cutoff", "1.1"], "--spagxe-marginal-cutoff must be in [0,1]"),
            ("spa-threshold-singularity", ["--spa-z-threshold", "0"], "--spa-z-threshold must be at least 0.01"),
            ("chunk-zero", ["--chunk-ksnp", "0"], "positive multiple of 0.5"),
            ("chunk-nonmultiple", ["--chunk-ksnp", "0.7"], "positive multiple of 0.5"),
            ("chunk-overflow", ["--chunk-ksnp", "999999999999"], "supported integer range"),
            ("duplicate-flag", ["--threads", "1", "--threads", "2"], "specified more than once"),
            ("genotype-mutual-exclusion", ["--pfile", "unused_pgen"], "Exactly one of --bfile"),
            ("duplicate-environment-list", [], "--envir-name contains duplicate entry"),
            ("empty-environment-list", [], "--envir-name contains an empty list entry"),
            ("trailing-environment-comma", [], "--envir-name contains an empty list entry"),
            ("residual-phenotype-mutual-exclusion", ["--resid-name", "R"], "--pheno-name and --resid-name are mutually exclusive"),
            ("base-pc-option-rejected", ["--pc-cols", "PC1"], "--pc-cols is for SPAGxEmix"),
            ("grm-mutual-exclusion", ["--grm-file", "a", "--sp-grm-plink2", "b"], "are mutually exclusive"),
            ("local-input-requires-method", ["--local-manifest", "a"], "require --method SPAGxElocal"),
            ("local-memory-requires-method", ["--local-max-block-bytes", "4096"], "require --method SPAGxElocal"),
            ("cct-input-requires-method", ["--cct-global", "a"], "require --method SPAGxEcct"),
            ("bad-compression", ["--compression", "zip"], "--compression must be gz or zst"),
            ("compression-level-without-codec", ["--compression-level", "3"], "requires --compression"),
            ("gzip-level-range", ["--compression", "gz", "--compression-level", "10"], "1..9 for gz"),
            ("negative-seed", ["--seed", "-1"], "non-negative integer"),
            ("chromosome-trailing-text", ["--chr", "1x-3"], "Invalid --chr numeric range"),
            ("chromosome-reversed-range", ["--chr", "3-1"], "Invalid or excessive --chr range"),
            ("utility-with-method", ["--cal-af-coef", "--method", "SPAGxE"], "cannot be combined with --method"),
            ("utility-mutual-exclusion", ["--cal-af-coef", "--int-pheno"], "are mutually exclusive"),
        ]
        for name, extra, message in invalid:
            command = base + extra
            if name in ("duplicate-environment-list", "empty-environment-list", "trailing-environment-comma"):
                command = base.copy()
                command[command.index("--envir-name") + 1] = {
                    "duplicate-environment-list": "E,E", "empty-environment-list": "E,,F",
                    "trailing-environment-comma": "E,"}[name]
            check(name, command, 1, message)
        check("missing-value", ["--threads"], 1, "Missing value for --threads")
        check("flag-as-missing-value", ["--pheno", "--out", "unused"], 1, "Missing value for --pheno")
        check("missing-bgen-allele-mode", ["--bgen", "cohort.bgen"], 1, "ref-first, ref-last, or ref-unknown")
        check("bad-bgen-allele-mode", ["--bgen", "cohort.bgen", "unknown"], 1, "ref-first, ref-last, or ref-unknown")
        local = ["--method", "SPAGxElocal", "--local-manifest", "unused_local",
                 "--pheno", "unused_phenotype", "--pheno-name", "Y", "--envir-name", "E",
                 "--covar-name", "X,E", "--out", str(pathlib.Path(directory) / "local")]
        for name, extra, message in [
            ("local-global-genotype-rejected", ["--bfile", "a"], "does not accept --bfile"),
            ("local-grm-rejected", ["--grm-file", "a"], "GRM input is not supported"),
            ("local-pc-af-rejected", ["--pc-cols", "PC1"], "does not fit individual AF"),
            ("local-hardcall-rejected", ["--hard-call-threshold", "0.1"], "does not use --hard-call-threshold"),
            ("local-input-mutual-exclusion", ["--local-prefix", "b"], "requires one of"),
            ("local-memory-small", ["--local-max-block-bytes", "100"], "must be at least 4096"),
        ]:
            check(name, local + extra, 1, message)
        pack = ["--pack-local", "--local-prefix", "unused_counts", "--out", str(pathlib.Path(directory)/"pack")]
        for name, flags, fragment in [
            ("pack-block-zero", ["--local-block-variants","0"], "1..65536"),
            ("pack-block-large", ["--local-block-variants","65537"], "1..65536"),
            ("pack-zstd-low", ["--local-zstd-level","-6"], "-5..22"),
            ("pack-zstd-high", ["--local-zstd-level","23"], "-5..22"),
            ("pack-analysis-option", ["--pheno","unused"], "does not accept analysis option"),
            ("pack-mutual-utility", ["--int-pheno"], "are mutually exclusive"),
        ]:
            check(name,pack+flags,1,fragment)
        check("packing-sidecar-needs-utility", base+["--local-bim","unused.bim"],1,"require --pack-local")
        check("local-input-required", ["--method", "SPAGxElocal"], 1, "requires one of")
        check("local-alias-conflict", ["--method", "SPAGxElocal", "--local-prefix", "a", "--local-binary", "b"],
              1, "select the same input")

    report = {"binary": str(binary), "total": len(cases),
              "passed": sum(item["passed"] for item in cases), "cases": cases}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("{passed}/{total} CLI checks passed".format(**report))
    return 0 if report["passed"] == report["total"] else 1


if __name__ == "__main__":
    sys.exit(main())

