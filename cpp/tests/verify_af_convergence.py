"""Check failed AF fits produce explicit NA rather than plausible p-values.

Run from cpp/: python -B tests/verify_af_convergence.py BINARY NEW_OUTPUT_DIR
Optionally add --baseline PREVIOUS_BINARY to check converged AF/model invariants
and report CGF-dependent differences separately.
Only the Python standard library is required.
"""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import subprocess
from reference_comparison import comparison, depends_on_spa_cgf, valid_changed_p


def write_fixture(prefix):
    n = 400
    pc = [-2 + 4 * i / (n - 1) for i in range(n)]
    separated = [int(x > 0) for x in pc]
    ordinary = [int((i * .618033988749895) % 1 < 1 / (1 + math.exp(-3 * x)))
                for i, x in enumerate(pc)]
    columns = {
        "Separated": separated,
        "Separated_ALT_complement": [2 - g for g in separated],
        "Nonseparated": ordinary,
        "Nonseparated_ALT_complement": [2 - g for g in ordinary],
    }
    ids = ["s%04d" % i for i in range(n)]
    prefix.with_suffix(".fam").write_text(
        "".join("%s %s 0 0 0 -9\n" % (i, i) for i in ids), encoding="utf-8")
    prefix.with_suffix(".bim").write_text(
        "".join("1 %s 0 %d A G\n" % (name, j + 1)
                for j, name in enumerate(columns)), encoding="utf-8")
    prefix.with_suffix(".pheno").write_text(
        "IID\tPC\tE\tY\n" + "".join(
            "%s\t%.17g\t%.17g\t%.17g\n" %
            (iid, x, math.cos(i * .13) + math.sin(i * .071),
             .2 * x + .3 * (math.cos(i * .13) + math.sin(i * .071)) +
             math.sin(i * 1.371) + .4 * math.cos(i * .837))
            for i, (iid, x) in enumerate(zip(ids, pc))), encoding="utf-8")
    raw = bytearray((0x6c, 0x1b, 1))
    for values in columns.values():
        for start in range(0, n, 4):
            byte = 0
            for offset, value in enumerate(values[start:start + 4]):
                byte |= (3, 2, 0)[value] << (2 * offset)
            raw.append(byte)
    prefix.with_suffix(".bed").write_bytes(raw)
    return list(columns)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--baseline", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = output / "cohort"
    names = write_fixture(prefix)
    checks = []
    commands = []

    def check(name, passed, detail=""):
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    def run(executable, arguments, label):
        result = subprocess.run([str(executable), *map(str, arguments)],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                timeout=120)
        (output / (label + ".console.log")).write_bytes(result.stdout)
        commands.append({"label": label, "exit": result.returncode})
        if result.returncode:
            raise RuntimeError(label + " failed: " +
                               result.stdout.decode("utf-8", errors="replace")[-2000:])

    common = ["--bfile", prefix, "--pheno", prefix.with_suffix(".pheno"),
              "--pc-cols", "PC", "--mac", "0", "--maf", "0", "--geno", "1", "--hwe", "0"]
    association = ["--method", "SPAGxEmix", "--pheno-name", "Y",
                   "--covar-name", "PC,E", "--envir-name", "E"]

    def scan(executable, label, threads):
        destination = output / label
        run(executable, common + association + ["--threads", str(threads), "--out", destination], label)
        files = [p for p in output.glob(label + ".*")
                 if "SPAGxEmix" in p.name and p.is_file() and not p.name.endswith(".log")]
        if len(files) != 1:
            raise RuntimeError("Expected one association table for " + label)
        with files[0].open(encoding="utf-8") as stream:
            return {row["ID"]: row for row in csv.DictReader(stream, delimiter="\t")}

    rows = scan(binary, "current", 1)
    parallel = scan(binary, "parallel", 2)
    check("all selected variants remain visible", set(rows) == set(names))
    for name in names[:2]:
        row = rows[name]
        check(name + " marginal AF failure status", row["SPA_STATUS_G"] == "9")
        check(name + " interaction AF failure status", row["SPA_STATUS_GxE"] == "9")
        for column in ("LOG10P_G", "LOG10P_GxE", "LOG10P_SPA_GxE", "LOG10P_Wald_GxE"):
            check(name + " " + column + " is NA", row[column] == "NA", row[column])
    for name in names[2:]:
        row = rows[name]
        check(name + " valid score P", math.isfinite(float(row["LOG10P_SPA_GxE"])))
        check(name + " valid final P", math.isfinite(float(row["LOG10P_GxE"])))
        check(name + " no AF failure", row["SPA_STATUS_GxE"] != "9")
    check("thread count leaves complete rows unchanged", rows == parallel)

    current_models = None
    for compression, suffix in (("", ".afc"), ("gz", ".afc.gz")):
        destination = output / ("models" + ("_gz" if compression else ""))
        extra = ["--compression", compression] if compression else []
        run(binary, ["--cal-af-coef", *common, *extra, "--out", destination],
            "export_" + (compression or "plain"))
        path = Path(str(destination) + suffix)
        opener = gzip.open if compression else open
        with opener(path, "rt", encoding="utf-8") as stream:
            lines = list(stream)
        check("export explicitly documents failed-fit status " + suffix,
              any("-2=AF_fit_failed" in line for line in lines))
        models = [line.strip().split("\t") for line in lines if line.strip() and not line.startswith("#")]
        if not compression:
            current_models = models
        check("export preserves failed and converged states " + suffix,
              [int(row[0]) for row in models] == [-2, -2, 2, 3])
        check("failed export has no usable coefficients " + suffix,
              all(float(value) == 0 for row in models[:2] for value in row[1:]))

    cgf_differences = []
    if args.baseline:
        baseline = scan(args.baseline.resolve(), "baseline", 1)
        count = 0
        maximum = 0.0
        for name in names[2:]:
            for column, value in rows[name].items():
                old = baseline[name][column]
                passed, delta = comparison(value, old, 1e-10)
                if depends_on_spa_cgf(column):
                    cgf_differences.append({"ID": name, "field": column,
                                            "current": value, "full_reference": old,
                                            "absolute_difference": delta})
                    check(name + " valid changed P " + column,
                          valid_changed_p(column, value, old), value + " vs " + old)
                else:
                    if delta is not None:
                        maximum = max(maximum, delta)
                        count += 1
                    check(name + " invariant baseline " + column, passed, value + " vs " + old)
        check("baseline comparison checked numerical fields", count > 0)
        destination = output / "baseline_models"
        run(args.baseline.resolve(), ["--cal-af-coef", *common, "--out", destination], "baseline_export")
        with Path(str(destination)+".afc").open(encoding="utf-8") as stream:
            old_models = [line.strip().split("\t") for line in stream
                          if line.strip() and not line.startswith("#")]
        check("baseline AF export row count", len(old_models) == len(current_models))
        for index in range(2, len(current_models)):
            current_model, old_model = current_models[index], old_models[index]
            check(names[index] + " AF export schema", len(current_model) == len(old_model))
            for column, (value, old) in enumerate(zip(current_model, old_model)):
                passed, delta = comparison(value, old, 1e-10)
                check(names[index] + " AF coefficient " + str(column), passed, value + " vs " + old)
    else:
        maximum = None
    run(binary, ["--help", "output"], "output_help")
    help_text = (output / "output_help.console.log").read_text(encoding="utf-8", errors="replace")
    check("output help describes status9", "9 AF_FIT_FAILED" in help_text)
    result = {"passed": all(x["passed"] for x in checks), "commands": commands,
              "checks": len(checks), "max_common_domain_difference": maximum,
              "baseline_comparison_scope": "AF/model invariants; SPA-dependent differences are descriptive",
              "cgf_dependent_differences": cgf_differences,
              "failures": [x for x in checks if not x["passed"]]}
    (output / "AF_CONVERGENCE.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return int(not result["passed"])


if __name__ == "__main__":
    raise SystemExit(main())
