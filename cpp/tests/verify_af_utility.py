"""Independent AF-export/QC integration checks; Python standard library only.
Run from cpp: python -B tests/verify_af_utility.py BINARY NEW_OUTPUT_DIRECTORY
"""
import gzip
import json
import math
from pathlib import Path
import subprocess
import sys


def write_fixture(prefix, pcs, columns):
    n = len(pcs)
    ids = ["%06d" % (i + 1) for i in range(n)]
    prefix.with_suffix(".pheno").write_text(
        "IID\tPC\n" + "".join("%s\t%.17g\n" % pair for pair in zip(ids, pcs)), encoding="utf-8")
    prefix.with_suffix(".fam").write_text(
        "".join("%s %s 0 0 0 -9\n" % (iid, iid) for iid in ids), encoding="utf-8")
    prefix.with_suffix(".bim").write_text(
        "".join("1 %s 0 %d A G\n" % (name, j + 1) for j, name in enumerate(columns)), encoding="utf-8")
    raw = bytearray((0x6c, 0x1b, 1))
    for values in columns.values():
        for start in range(0, n, 4):
            byte = 0
            for offset, value in enumerate(values[start:start + 4]):
                code = 1 if value is None else (3, 2, 0)[value]
                byte |= code << (2 * offset)
            raw.append(byte)
    prefix.with_suffix(".bed").write_bytes(raw)


def read_models(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [(int(fields[0]), [float(v) for v in fields[1:]])
                for fields in (line.split() for line in handle)
                if fields and not fields[0].startswith("#")]


def logistic_mle(pcs, carriers):
    """Independent two-parameter score/Hessian Newton solution of binomial MLE."""
    beta = [0.0, 0.0]
    for _ in range(200):
        s0 = s1 = h00 = h01 = h11 = 0.0
        for x, y in zip(pcs, carriers):
            eta = beta[0] + beta[1] * x
            ex = math.exp(-abs(eta))
            mu = 1 / (1 + ex) if eta >= 0 else ex / (1 + ex)
            w = mu * (1 - mu)
            s0 += y - mu
            s1 += (y - mu) * x
            h00 += w
            h01 += w * x
            h11 += w * x * x
        det = h00 * h11 - h01 * h01
        delta = [(h11 * s0 - h01 * s1) / det, (h00 * s1 - h01 * s0) / det]
        beta = [b + d for b, d in zip(beta, delta)]
        if max(map(abs, delta)) < 1e-12:
            return beta
    raise RuntimeError("independent logistic MLE did not converge")


def reconstructed(beta, pcs, complement=False):
    values = []
    for x in pcs:
        eta = beta[0] + beta[1] * x
        ex = math.exp(-abs(eta))
        mu = 1 / (1 + ex) if eta >= 0 else ex / (1 + ex)
        minor = 1 - math.sqrt(1 - mu)
        values.append(1 - minor if complement else minor)
    return values


def main():
    binary = Path(sys.argv[1]).resolve()
    output = Path(sys.argv[2]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    checks = []

    def check(name, ok, detail=""):
        checks.append({"name": name, "pass": bool(ok), "detail": detail})

    def run(prefix, tag, extra=()):
        destination = output / tag
        command = [str(binary), "--cal-af-coef", "--bfile", str(prefix),
                   "--pheno", str(prefix.with_suffix(".pheno")), "--pc-cols", "PC",
                   "--maf", "0", "--mac", "0", "--geno", "1", "--hwe", "0",
                   "--out", str(destination)] + list(extra)
        # Each option is emitted once: override defaults before invocation.
        for flag in ("--maf", "--mac", "--geno", "--hwe"):
            if flag in extra:
                first = command.index(flag)
                del command[first:first + 2]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output / (tag + ".console.txt")).write_bytes(result.stdout)
        check(tag + " completed", result.returncode == 0, str(result.returncode))
        if result.returncode:
            return []
        suffix = ".afc.gz" if "--compression" in extra and extra[extra.index("--compression") + 1] == "gz" else ".afc"
        return read_models(Path(str(destination) + suffix))

    n = 400
    boundary_cases = [
        ([3.] * 60 + [1.] * 40 + [0.] * 300,
         [2] * 58 + [0] * 2 + [2] * 38 + [0] * 2 + [1] * 10 + [0] * 290),
        ([-1.] * 40 + [0.] * 320 + [1.] * 40,
         [1] * 3 + [0] * 37 + [1] * 10 + [0] * 310 + [2] * 35 + [0] * 5),
        ([-1 + 2 * i / (n - 1) for i in range(n)], [0] * n)]
    boundary_cases[2][1][0] = 2
    for i in range(372, 399, 2):
        boundary_cases[2][1][i] = 2
    for j, (pcs, g) in enumerate(boundary_cases):
        mean = sum(pcs) / n
        sd = math.sqrt(sum((x - mean) ** 2 for x in pcs) / (n - 1))
        pcs = [(x - mean) / sd for x in pcs]
        prefix = output / ("boundary%d" % j)
        write_fixture(prefix, pcs, {"low": g, "high": [2 - value for value in g]})
        models = run(prefix, "boundary%d_export" % j)
        check("boundary%d statuses 2/3" % j, [row[0] for row in models] == [2, 3])
        expected = logistic_mle(pcs, [float(v > .5) for v in g])
        if len(models) == 2:
            for allele, (status, beta) in enumerate(models):
                error = max(abs(a - b) for a, b in zip(
                    reconstructed(beta, pcs, allele == 1), reconstructed(expected, pcs, allele == 1)))
                check("boundary%d allele%d independent AF" % (j, allele), error < 2e-7, error)
        parallel = run(prefix, "boundary%d_threads" % j, ("--threads", "4"))
        check("boundary%d threading equality" % j, parallel == models)

    pcs = [math.sin(i * .37) for i in range(n)]
    common = [(0, 1, 1, 2)[i % 4] for i in range(n)]
    columns = {"common": common, "missing": [None] * 80 + common[80:],
               "rare": [1] * 4 + [0] * (n - 4), "lowmac": [1] * 8 + [0] * (n - 8),
               "hwe": [1] * n, "mono": [0] * n, "allmissing": [None] * n}
    prefix = output / "qc"
    write_fixture(prefix, pcs, columns)
    base = run(prefix, "qc_base")
    check("allmissing explicitly unfitted", len(base) == 7 and base[-1][0] == -1)
    check("monomorphic fitted pooled", len(base) == 7 and base[-2][0] == 0)
    for tag, flag, threshold, failed in (
            ("missing", "--geno", ".1", {"missing", "allmissing"}),
            ("maf", "--maf", ".02", {"rare", "lowmac", "mono", "allmissing"}),
            ("mac", "--mac", "10", {"rare", "lowmac", "mono", "allmissing"}),
            ("hwe", "--hwe", "1e-6", {"hwe", "allmissing"})):
        rows = run(prefix, "qc_" + tag, (flag, threshold))
        actual = {name for name, row in zip(columns, rows) if row[0] == -1}
        check("QC " + tag + " exact filtered set", len(rows) == len(columns) and actual == failed, sorted(actual))
    include = output / "extract.txt"
    include.write_text("common\n", encoding="ascii")
    selected = run(prefix, "qc_extract", ("--extract", str(include)))
    check("extract retains only selected marker", len(selected) == 1 and selected[0] == base[0])
    compressed = run(prefix, "qc_gzip", ("--compression", "gz"))
    check("gzip export matches plain", compressed == base)
    summary = {"checks": len(checks), "passed": sum(row["pass"] for row in checks),
               "failed": sum(not row["pass"] for row in checks), "results": checks}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in ("checks", "passed", "failed")}, ensure_ascii=True))
    return int(summary["failed"] != 0)


if __name__ == "__main__":
    sys.exit(main())
