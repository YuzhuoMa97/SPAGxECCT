"""Native local GxE / Wald / CCT comparison with independent R Bin(H,q) oracle.
Run from cpp: python tests/verify_local_oracle.py [binary] [validation-dir].
Generate the reference first with Rscript tests/local_oracle.R.
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "build/bin" / ("SPAGxE.exe" if os.name == "nt" else "SPAGxE")
OUT = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else ROOT.parent.parent / "validation/local_oracle"
FIX = ROOT / "examples/local"
OUT.mkdir(parents=True, exist_ok=True)


def tab(path):
    with open(path, encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def number(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return math.nan


def relative(path):
    return os.path.relpath(path, ROOT)


reference = tab(FIX / "local_reference.tsv")
combined = tab(FIX / "local_combined_reference.tsv")
traits = {"quantitative": ("Quantitative", "R_quant"), "binary": ("Binary", "R_binary"), "survival": ("Time_Event", "R_survival"), "ordinal": ("Ordinal", "R_ordinal"), "masked": (None, "R_masked")}
failures, comparisons, summaries, native_residual_checks = [], [], [], []
outputs = {}


def check(tag, detail, expected, actual, tolerance=2e-5):
    if (math.isnan(expected) and math.isnan(actual)) or expected == actual:
        ok, error = True, 0.0
    else:
        error = abs(actual - expected)
        ok = math.isfinite(error) and error <= tolerance
    if not ok:
        failures.append(f"{tag}:{detail}: {actual} vs {expected}, error={error}")
    comparisons.append({"run": tag, "detail": detail, "expected": expected, "actual": actual, "abs_error": error, "pass": ok})
    return error


for dataset in ("two", "three"):
    dataset_traits = traits
    for mode in ("native", "residual"):
        for cgf in ("guarded_partial_normal",):
            for threads in (1, 4):
                tag = f"{dataset}_{mode}_{cgf}_t{threads}"
                prefix = OUT / tag
                fixture = FIX if dataset == "two" else FIX / "three_ancestry"
                pheno = FIX / "samples.pheno"
                command = [str(BIN), "--method", "SPAGxElocal", "--local-manifest", relative(fixture / "cohort.local.tsv"), "--pheno", relative(pheno), "--mac", "0", "--threads", str(threads), "--out", relative(prefix)]
                if mode == "native":
                    command += ["--pheno-name", "Quantitative,Binary,Time:Event,Ordinal", "--regression-model", "auto", "--covar-name", "X,E", "--envir-name", "E", "--save-resid"]
                else:
                    command += ["--resid-name", "R_quant,R_binary,R_survival,R_ordinal,R_masked", "--envir-name", "E,Eshift"]
                process = subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", timeout=60)
                (OUT / (tag + ".console.txt")).write_text(process.stdout, encoding="utf-8")
                if process.returncode:
                    failures.append(f"{tag}: process exit {process.returncode}: {process.stdout[-2000:]}")
                    continue
                start, bad_start = len(comparisons), len(failures)
                for trait, (native_col, residual_col) in dataset_traits.items():
                    if mode == "native" and native_col is None:
                        continue
                    column = native_col if mode == "native" else residual_col
                    local_path = Path(str(prefix) + "." + column + ".SPAGxE_local")
                    combined_path = Path(str(prefix) + "." + column + ".SPAGxE_local_CCT")
                    rows = tab(local_path)
                    result = {(r["ID"], r["ANCESTRY"]): r for r in rows}
                    actual_combined = {r["ID"]: r for r in tab(combined_path)}
                    expected_rows = [r for r in reference if r["DATASET"] == dataset and r["MODE"] == mode and r["CGF"] == cgf and r["TRAIT"] == trait]
                    expected_keys = {(r["ID"], r["ANCESTRY"]) for r in expected_rows}
                    if len(rows) != len(result) or set(result) != expected_keys:
                        failures.append(f"{tag}:{trait}: missing, extra or duplicate ancestry/variant rows")
                    for oracle in expected_rows:
                        ident, ancestry, en = oracle["ID"], oracle["ANCESTRY"], oracle["ENV"]
                        row = result.get((ident, ancestry))
                        if row is None:
                            continue
                        fields = {"LOG10P_G": "LOG10P_G", "Z_G": "Z_G", "Z_Norm": "Z_Norm_Gx" + en, "LOG10P_SPA": "LOG10P_SPA_Gx" + en, "LOG10P_WALD": "LOG10P_Wald_Gx" + en, "LOG10P_FINAL": "LOG10P_Gx" + en, "N": "N_TESTED", "N_H_MISSING": "N_H_MISSING", "AF": "ALT_FREQ", "MISS_RATE": "MISS_RATE", "MAC": "MAC", "FLIPPED": "MINOR_FLIPPED"}
                        for expected_field, actual_field in fields.items():
                            check(tag, f"{trait}:{ident}:{ancestry}:{en}:{expected_field}", number(oracle[expected_field]), number(row[actual_field]))
                        if math.isfinite(number(oracle["VAR"])) and number(oracle["VAR"]) > 0:
                            check(tag, f"{trait}:{ident}:{ancestry}:{en}:SE", 1 / math.sqrt(number(oracle["VAR"])), number(row["SE_Gx" + en]))
                            check(tag, f"{trait}:{ident}:{ancestry}:{en}:branch", 2 if oracle["BRANCH"] == "B" else 1, number(row["BRANCH_Gx" + en]), 0)
                    expected_combined = [r for r in combined if r["DATASET"] == dataset and r["MODE"] == mode and r["CGF"] == cgf and r["TRAIT"] == trait]
                    for oracle in expected_combined:
                        row = actual_combined.get(oracle["ID"])
                        if row is None:
                            failures.append(f"{tag}:{trait}:{oracle['ID']}: combined row absent")
                            continue
                        en = oracle["ENV"]
                        check(tag, f"{trait}:{oracle['ID']}:{en}:ANCESTRY_CCT", number(oracle["LOG10P_ANCESTRY_CCT"]), number(row["LOG10P_ANCESTRY_CCT_Gx" + en]))
                        check(tag, f"{trait}:{oracle['ID']}:{en}:N_ANCESTRIES", number(oracle["N_ANCESTRIES"]), number(row["N_ANCESTRIES_Gx" + en]), 0)
                    digest = (hashlib.sha256(local_path.read_bytes()).hexdigest(), hashlib.sha256(combined_path.read_bytes()).hexdigest())
                    key = (dataset, mode, cgf, trait)
                    if threads == 1:
                        outputs[key] = digest
                    elif outputs.get(key) != digest:
                        failures.append(f"{tag}:{trait}: results differ between 1 and 4 threads")
                if mode == "native":
                    source_residual = {r["IID"]: r for r in tab(FIX / "native_R_residuals.tsv")}
                    fitted = {r["IID"]: r for r in tab(str(prefix) + ".null.resid")}
                    for name, reference_col in (("Quantitative", "R_quant_native"), ("Binary", "R_binary_native"), ("Time_Event", "R_survival_native"), ("Ordinal", "R_ordinal_native")):
                        error = max(abs(number(fitted[iid][name]) - number(ref[reference_col])) for iid, ref in source_residual.items())
                        native_residual_checks.append({"run": tag, "trait": name, "subjects": len(fitted), "max_abs_error": error})
                        if not math.isfinite(error) or error > 1e-5:
                            failures.append(f"{tag}:{name}: native null residual error {error}")
                summaries.append({"run": tag, "comparisons": len(comparisons) - start, "failures": len(failures) - bad_start})

with open(OUT / "local_cpp_comparison.tsv", "w", encoding="utf-8", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=["run", "detail", "expected", "actual", "abs_error", "pass"], delimiter="\t")
    writer.writeheader()
    writer.writerows(comparisons)
summary = {"cgf_policy": "guarded_partial_normal", "passed": not failures, "binary": str(BIN), "binary_sha256": hashlib.sha256(BIN.read_bytes()).hexdigest(), "reference_rows": len(reference), "run_count": len(summaries), "comparison_count": len(comparisons), "max_abs_error": max((r["abs_error"] for r in comparisons), default=0), "native_residual_checks": native_residual_checks, "summaries": summaries, "failure_count": len(failures), "failures": failures}
(OUT / "local_cpp_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps({k: summary[k] for k in ("passed", "reference_rows", "run_count", "comparison_count", "max_abs_error", "failure_count")}, indent=2))
if failures:
    print("\n".join(failures[:30]))
sys.exit(0 if not failures else 1)
