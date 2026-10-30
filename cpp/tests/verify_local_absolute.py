"""Windows native-path and opaque Unicode/string IID regression.

Uses absolute paths including the containing project's Chinese directories.
All generated files are under the explicitly selected output directory.
"""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2]).resolve()
FIX = ROOT / "examples/local"
OUT.mkdir(parents=True, exist_ok=True)
checks = []


def run(tag, binary=False, fixture=FIX):
    prefix = OUT / tag
    args = [str(BIN), "--method", "SPAGxElocal", "--pheno", str(fixture / "samples.pheno"), "--resid-name", "R" if binary else "R_quant", "--envir-name", "E", "--maf", "0", "--mac", "0", "--geno", "1", "--out", str(prefix)]
    args += ["--local-binary", str(fixture / "local")] if binary else ["--local-manifest", str(fixture / "cohort.local.tsv")]
    process = subprocess.run(args, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", timeout=120)
    (OUT / (tag + ".console.txt")).write_text(process.stdout, encoding="utf-8")
    checks.append({"case": tag, "pass": process.returncode == 0, "exit": process.returncode})
    return Path(str(prefix) + (".R" if binary else ".R_quant") + ".SPAGxE_local")


baseline = run("absolute_text")
run("absolute_binary", binary=True, fixture=FIX / "binary_demo")
unicode_fixture = OUT / "\u7956\u6e90\u683c\u5f0f_\u524d\u5bfc\u96f6_ID"
shutil.copytree(FIX, unicode_fixture, ignore=shutil.ignore_patterns("binary_demo", "bgen_demo", "three_ancestry"), dirs_exist_ok=True)
replace = {"s001": "0001", "s002": "0002", "s003": "\u6837\u672c003"}
for file in unicode_fixture.glob("*.tsv"):
    if file.name == "cohort.local.tsv":
        continue
    with open(file, encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    rows[0] = [replace.get(token, token) for token in rows[0]]
    with open(file, "w", encoding="utf-8", newline="") as stream:
        csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)
pheno = unicode_fixture / "samples.pheno"
with open(pheno, encoding="utf-8", newline="") as stream:
    rows = list(csv.reader(stream, delimiter="\t"))
for row in rows[1:]:
    row[0] = replace.get(row[0], row[0])
with open(pheno, "w", encoding="utf-8", newline="") as stream:
    csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)
relabel = run("absolute_Unicode_IID", fixture=unicode_fixture)
if baseline.exists() and relabel.exists():
    checks.append({"case": "leading_zero_Unicode_IID_numeric_identity", "pass": baseline.read_bytes() == relabel.read_bytes()})
summary = {"checks": len(checks), "passed": sum(item["pass"] for item in checks), "details": checks, "scope": "Absolute Chinese Windows directory paths and opaque string IID labels. Manifest filenames themselves use ASCII G/H basenames."}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2, ensure_ascii=False))
raise SystemExit(any(not item["pass"] for item in checks))
