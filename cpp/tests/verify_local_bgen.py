"""Native known-H BGEN parsing, numerical equivalence and rejection checks.

Run: python -B tests/verify_local_bgen.py [binary] [output-directory]
Only SPAGxE-owned deterministic simulations and Python's stdlib are used.
"""
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys

from make_local_bgen_fixtures import probabilities, simulated_records, table, write_bgen

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "build/bin" / ("SPAGxE.exe" if os.name == "nt" else "SPAGxE")
OUT = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else ROOT.parent.parent / "validation/local_bgen"
FIX = ROOT / "examples/local/bgen_demo"
OUT.mkdir(parents=True, exist_ok=True)
simulation = json.loads((FIX / "SIMULATION.json").read_text())
recipes = {row["stem"]: row for row in simulation["recipes"]}
checks, failures = [], []


def check(name, passed, detail=""):
    checks.append({"name": name, "pass": bool(passed), "detail": str(detail)})
    if not passed:
        failures.append(name + ": " + str(detail))


def rows(path):
    with open(path, encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run(tag, manifest, text=False, automatic=False, extra=(), success=True, native=False, pheno=None, error=None, absolute=False):
    prefix = OUT / tag
    relative = lambda p: str(p) if absolute else os.path.relpath(p, ROOT)
    cmd = [str(BIN), "--method", "SPAGxElocal", "--local-manifest" if text or automatic else "--local-bgen-manifest", relative(manifest),
           "--pheno", relative(pheno or (FIX.parent / "samples.pheno" if native else FIX / "samples.pheno")),
           "--envir-name", "E", "--maf", "0", "--mac", "0", "--geno", "1", "--out", relative(prefix)]
    if native:
        cmd += ["--pheno-name", "Quantitative,Binary,Ordinal,Time:Event", "--covar-name", "X,E", "--regression-model", "auto", "--spagxe-marginal-cutoff", "1"]
    else:
        cmd += ["--resid-name", "R"]
    cmd += list(extra)
    process = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", timeout=120)
    (OUT / (tag + ".console.txt")).write_text(process.stdout, encoding="utf-8")
    check(tag + "_exit", (process.returncode == 0) == success, process.stdout[-1200:] if (process.returncode == 0) != success else process.returncode)
    if error:
        check(tag + "_diagnostic", error in process.stdout, process.stdout[-1200:] if error not in process.stdout else error)
    return prefix, process


def result(prefix, phenotype="R", combined=False):
    return Path(str(prefix) + "." + phenotype + (".SPAGxE_local_CCT" if combined else ".SPAGxE_local"))


def equivalent(name, left, right, phenotype="R", exact=False):
    worst, okay = 0.0, True
    for combined in (False, True):
        a, b = rows(result(left, phenotype, combined)), rows(result(right, phenotype, combined))
        okay &= len(a) == len(b)
        for x, y in zip(a, b):
            okay &= x.keys() == y.keys()
            for key in x:
                # A posterior can remain uncertain even when its expectation
                # is exactly an integer. Text has lost that uncertainty, so
                # only the BGEN source reliably marks HWE unavailable there.
                if key == "LOG10P_HWE_GLOBAL":
                    continue
                if x[key] == y[key]:
                    continue
                try:
                    u, v = float(x[key]), float(y[key])
                except ValueError:
                    okay = False
                    continue
                if math.isnan(u) and math.isnan(v):
                    continue
                error = abs(u - v)
                worst = max(worst, error)
                okay &= math.isfinite(error) and error <= (0 if exact else 2e-9) * max(1, abs(u), abs(v))
    check(name, okay, "maximum absolute difference %.17g" % worst)


baseline = {}
for stem, recipe in recipes.items():
    native, process = run(stem, FIX / (stem + ".localbgen.tsv"))
    text, text_process = run(stem + "_text", FIX / (stem + ".local.tsv"), text=True)
    if process.returncode == 0 and text_process.returncode == 0:
        equivalent(stem + "_posterior_expectation_equivalence", native, text, exact=True)
        output = rows(result(native))
        check(stem + "_all_ancestry_rows", len(output) == 30)
        check(stem + "_zero_H_structural_missing", all(row["MISS_RATE"] == "0" for row in output if row["ID"] == "H_zero" and row["ANCESTRY"] != "Ancestry3"))
    baseline[stem] = native

chosen = "bits16_zlib"
automatic, _ = run("automatic_manifest_detection", FIX / (chosen + ".localbgen.tsv"), automatic=True)
equivalent("automatic_manifest_numerical_identity", baseline[chosen], automatic, exact=True)
threaded, _ = run("threads4_repeat", FIX / (chosen + ".localbgen.tsv"), extra=["--threads", "4"])
check("threads4_byte_identity", result(baseline[chosen]).read_bytes() == result(threaded).read_bytes() and result(baseline[chosen], combined=True).read_bytes() == result(threaded, combined=True).read_bytes())
for suffix, text in (("native", False), ("native_text", True)):
    prefix, process = run(suffix, FIX / (chosen + (".local.tsv" if text else ".localbgen.tsv")), text=text, native=True)
    if not text:
        native_prefix = prefix
    elif process.returncode == 0:
        for phenotype in ("Quantitative", "Binary", "Ordinal", "Time_Event"):
            equivalent("native_" + phenotype + "_equivalence", native_prefix, prefix, phenotype)

ref_manifest = OUT / "REF.localbgen.tsv"
ref = rows(FIX / (chosen + ".localbgen.tsv"))
table(ref_manifest, [["ancestry", "bgen", "counted_allele"], *[[r["ancestry"], os.path.relpath(FIX / r["bgen"], OUT), "REF"] for r in ref]])
ref_output, _ = run("REF_direction", ref_manifest)
equivalent("REF_H_minus_G_canonicalization", baseline[chosen], ref_output)
integer_uncertainty = [r for r in rows(result(baseline[chosen])) if r["ID"] == "integer_mean"]
check("integer_expectation_retains_posterior_HWE_unavailability", all(r["LOG10P_HWE_GLOBAL"] == "NA" for r in integer_uncertainty))

unicode_dir = OUT / "Unicode_IIDs_\u7956\u6e90"
unicode_dir.mkdir(exist_ok=True)
unicode_ids, unicode_records = simulated_records(160, 3)
renamed_ids = [{"s001": "0001", "s002": "0002", "s003": "\u6837\u672c003"}.get(iid, iid) for iid in unicode_ids]
unicode_manifest = [["ancestry", "bgen", "counted_allele"]]
for a, branch in enumerate(unicode_records):
    label = "Ancestry%d" % (a+1)
    write_bgen(unicode_dir / (label + ".bgen"), renamed_ids, branch, 16, "zlib")
    unicode_manifest.append([label, label + ".bgen", "ALT"])
table(unicode_dir / "cohort.localbgen.tsv", unicode_manifest)
table(unicode_dir / "samples.pheno", [["IID", "R", "E"], *[[new, source["R"], source["E"]] for new, source in zip(renamed_ids, rows(FIX / "samples.pheno"))]])
unicode_out, _ = run("absolute_BGEN_Unicode_string_IIDs", unicode_dir / "cohort.localbgen.tsv", pheno=unicode_dir / "samples.pheno", absolute=True)
equivalent("Unicode_leading_zero_IID_numeric_identity", baseline[chosen], unicode_out, exact=True)

two_bgen, _ = run("two_ancestry_bgen", FIX / "two_ancestry/cohort.localbgen.tsv")
two_text, _ = run("two_ancestry_text", FIX / "two_ancestry/cohort.local.tsv", text=True)
equivalent("two_ancestry_equivalence", two_bgen, two_text, exact=True)
check("two_ancestry_rows", len(rows(result(two_bgen))) == 20)

# Test an actual physical allele permutation, not only a manifest declaration.
# A/G becomes G/A and P_new(g)=P_old(H-g), on the stored probability grid.
alleles = OUT / "actual_allele_permutation"
alleles.mkdir(exist_ok=True)
permutation_ids, permutation_records = simulated_records(160, 3, seed=20261027)
original_manifest = [["ancestry", "bgen", "counted_allele"]]
swapped_manifest = [["ancestry", "bgen", "counted_allele"]]
for a, branch in enumerate(permutation_records):
    original, swapped = [], []
    for j, source in enumerate(branch[:3]):
        posterior, reversed_posterior = [], []
        for i, copies in enumerate(source["H"]):
            q = (.19, .31, .73)[a] + .015 * math.sin(i+j)
            p = [math.comb(copies, g) * q**g * (1-q)**(copies-g) for g in range(copies+1)]
            stored = probabilities(p, 16)
            posterior.append([v/65535 for v in stored])
            reversed_posterior.append([v/65535 for v in reversed(stored)])
        metadata = ["1", str(9001+j), "allele_%d" % j, "A", "G"]
        original.append({"metadata": metadata, "H": source["H"], "P": posterior})
        swapped.append({"metadata": metadata[:3] + ["G", "A"], "H": source["H"], "P": reversed_posterior})
    label = "Ancestry%d" % (a+1)
    original_file, swapped_file = label + ".AG.bgen", label + ".GA.bgen"
    write_bgen(alleles / original_file, permutation_ids, original, 16, "zlib")
    write_bgen(alleles / swapped_file, permutation_ids, swapped, 16, "zlib")
    original_manifest.append([label, original_file, "ALT"])
    swapped_manifest.append([label, swapped_file, "REF"])
table(alleles / "AG.localbgen.tsv", original_manifest)
table(alleles / "GA.localbgen.tsv", swapped_manifest)
ag, _ = run("actual_AG", alleles / "AG.localbgen.tsv", extra=["--spagxe-marginal-cutoff", "1"])
ga, _ = run("actual_GA_REF", alleles / "GA.localbgen.tsv", extra=["--spagxe-marginal-cutoff", "1"])
left, right = rows(result(ag)), rows(result(ga))
metadata_ok, frequency_ok, direction_ok, statistics_ok = True, True, True, True
largest = 0.0
for x, y in zip(left, right):
    metadata_ok &= (x["REF"], x["ALT"], y["REF"], y["ALT"]) == ("A", "G", "G", "A")
    frequency_ok &= abs(float(x["ALT_FREQ"]) + float(y["ALT_FREQ"]) - 1) <= 2e-14
    direction_ok &= x["TESTED_ALLELE"] != y["TESTED_ALLELE"]
    for key in x:
        if not (key.startswith(("LOG10P_", "Z_", "BETA_", "SE_", "SPA_STATUS_", "BRANCH_")) or key == "MAC"):
            continue
        if x[key] == y[key]:
            continue
        u, v = float(x[key]), float(y[key])
        if math.isnan(u) and math.isnan(v):
            continue
        error = abs(u-v)
        largest = max(largest, error)
        statistics_ok &= math.isfinite(error) and error <= 2e-9 * max(1, abs(u), abs(v))
check("actual_allele_permutation_metadata", metadata_ok)
check("actual_allele_permutation_ALT_frequency", frequency_ok)
check("actual_allele_permutation_physical_minor_direction", direction_ok)
check("actual_allele_permutation_minor_statistics", statistics_ok, "maximum absolute difference %.17g" % largest)
mixed = copy.deepcopy(original_manifest)
mixed[2] = [swapped_manifest[2][0], swapped_manifest[2][1], "ALT"]
table(alleles / "mixed.localbgen.tsv", mixed)
run("reject_actual_mixed_allele_order", alleles / "mixed.localbgen.tsv", success=False, error="metadata/order")


def altered(tag, change, stem="bits16_raw", error=None):
    directory = OUT / ("input_" + tag)
    directory.mkdir(exist_ok=True)
    entries = rows(FIX / (stem + ".localbgen.tsv"))
    for entry in entries:
        shutil.copyfile(FIX / entry["bgen"], directory / entry["bgen"])
    path = directory / "cohort.localbgen.tsv"
    table(path, [["ancestry", "bgen", "counted_allele"], *[[r["ancestry"], r["bgen"], r["counted_allele"]] for r in entries]])
    change(directory, entries, recipes[stem])
    return run("reject_" + tag, path, success=False, error=error)


def binary_change(directory, entries, fn, ancestry=0):
    path = directory / entries[ancestry]["bgen"]
    data = bytearray(path.read_bytes())
    fn(data)
    path.write_bytes(data)


def set_u32(data, offset, value):
    data[offset:offset + 4] = struct.pack("<I", value)


def raw_offset(recipe, marker=0, ancestry="Ancestry1"):
    return recipe["offsets"][ancestry][marker]["raw_probability"]


def raw_change(fn):
    return lambda d, e, r: binary_change(d, e, lambda data: fn(data, raw_offset(r)))


for tag, offset, value, error in (
    ("bad_magic", 16, ord("X"), "magic"),
    ("reserved_flags", 20, 72, "reserved flags"),
    ("wrong_layout", 20, 4, "Layout 2"),
    ("unknown_codec", 20, 11, "Layout 2"),
):
    altered(tag, lambda d, e, r, off=offset, val=value: binary_change(d, e, lambda b: b.__setitem__(off, val)), error=error)
altered("external_sample_IDs", lambda d, e, r: binary_change(d, e, lambda b: b.__setitem__(23, b[23] & 127)), error="embedded sample IDs")
altered("empty_variants", lambda d, e, r: binary_change(d, e, lambda b: set_u32(b, 8, 0)), error="empty dimensions")
altered("sample_header_mismatch", lambda d, e, r: binary_change(d, e, lambda b: set_u32(b, 12, 159)), error="sample block")
def dishonest_samples(directory, entries, recipe):
    def change(data):
        set_u32(data, 12, 100000000)
        set_u32(data, 28, 100000000)
    binary_change(directory, entries, change)


altered("sample_count_vs_embedded_bytes", dishonest_samples, error="sample block")
altered("variant_count_mismatch", lambda d, e, r: binary_change(d, e, lambda b: set_u32(b, 8, 11), ancestry=1), error="variant counts differ")
altered("truncated", lambda d, e, r: binary_change(d, e, lambda b: b.__delitem__(slice(len(b)-3, None))), error="bounds")
altered("trailing_bytes", lambda d, e, r: binary_change(d, e, lambda b: b.extend(b"trailing")), error="trailing bytes")
altered("probability_N", raw_change(lambda b, o: set_u32(b, o, 159)), error="dimensions mismatch")
altered("probability_alleles", raw_change(lambda b, o: b.__setitem__(slice(o+4, o+6), struct.pack("<H", 3))), error="dimensions mismatch")
altered("ploidy_three", raw_change(lambda b, o: b.__setitem__(o+7, 3)), error="ploidy H")
altered("ploidy_reserved", raw_change(lambda b, o: b.__setitem__(o+8, b[o+8] | 64)), error="ploidy/reserved")
altered("ploidy_bounds", raw_change(lambda b, o: b.__setitem__(o+7, 1)), error="ploidy/reserved")
altered("phased", raw_change(lambda b, o: b.__setitem__(o+8+160, 1)), error="phased")
altered("bits_zero", raw_change(lambda b, o: b.__setitem__(o+9+160, 0)), error="1..32")
altered("bits_above32", raw_change(lambda b, o: b.__setitem__(o+9+160, 33)), error="1..32")


def probability_failure(data, start, missing=False):
    h = [v & 63 for v in data[start+8:start+8+160]]
    i = next(i for i, copies in enumerate(h) if copies == (1 if missing else 2))
    if missing:
        data[start+8+i] |= 128
    location = start + 10 + 160 + 2 * sum(h[:i])
    data[location:location+2] = b"\xff\xff"
    if not missing:
        data[location+2:location+4] = b"\xff\xff"


altered("probability_sum", raw_change(lambda b, o: probability_failure(b, o)), error="sum above one")
altered("missing_nonzero_probabilities", raw_change(lambda b, o: probability_failure(b, o, True)), error="missing posterior")


def padding_failure(directory, entries, recipe):
    def change(data):
        for record in recipe["offsets"]["Ancestry1"]:
            start = record["raw_probability"]
            bits = sum(v & 63 for v in data[start+8:start+8+160])
            if bits % 8:
                data[start+record["decoded_bytes"]-1] |= 128
                return
        raise AssertionError("No padding fixture")
    binary_change(directory, entries, change)


altered("nonzero_padding", padding_failure, stem="bits1_raw", error="padding")
altered("metadata_order", lambda d, e, r: binary_change(d, e, lambda b: set_u32(b, r["offsets"]["Ancestry1"][0]["position"], 9)), error="metadata/order")


def iid_change(directory, entries, recipe, duplicate):
    def change(data):
        first_length = int.from_bytes(data[32:34], "little")
        if duplicate:
            second = 34 + first_length
            length = int.from_bytes(data[second:second+2], "little")
            assert length == first_length
            data[second+2:second+2+length] = data[34:34+first_length]
        else:
            data[34] = ord("x")
    binary_change(directory, entries, change)


altered("duplicate_IID", lambda d, e, r: iid_change(d, e, r, True), error="duplicate sample IID")
altered("mismatched_IID_set", lambda d, e, r: iid_change(d, e, r, False), error="omits IID")


def manifest_failure(directory, entries, recipe, kind):
    if kind == "incomplete":
        entries = entries[:1]
    if kind == "direction":
        entries[1]["counted_allele"] = "REF"
    if kind == "label":
        entries[1]["ancestry"] = entries[0]["ancestry"]
    table(directory / "cohort.localbgen.tsv", [["ancestry", "bgen", "counted_allele"], *[[r["ancestry"], r["bgen"], r["counted_allele"]] for r in entries]])


altered("incomplete_ancestry_H", lambda d, e, r: manifest_failure(d, e, r, "incomplete"), error="sum to two")
altered("mixed_directions", lambda d, e, r: manifest_failure(d, e, r, "direction"), error="same counted direction")
altered("duplicate_ancestry", lambda d, e, r: manifest_failure(d, e, r, "label"), error="duplicate ancestry")


def compressed_failure(directory, entries, recipe, trailing):
    def change(data):
        info = recipe["offsets"]["Ancestry1"][-1]
        if trailing:
            data.extend(b"junk")
            set_u32(data, info["block"], int.from_bytes(data[info["block"]:info["block"]+4], "little") + 4)
        else:
            data[-1] ^= 255
    binary_change(directory, entries, change)


altered("zlib_integrity", lambda d, e, r: compressed_failure(d, e, r, False), stem="bits16_zlib", error="zlib")
altered("zlib_trailing_compressed_bytes", lambda d, e, r: compressed_failure(d, e, r, True), stem="bits16_zlib", error="zlib")
altered("zstd_trailing_compressed_bytes", lambda d, e, r: compressed_failure(d, e, r, True), stem="bits16_zstd", error="zstd frame")

large = OUT / "large_block"
large.mkdir(exist_ok=True)
ids = ["large%04d" % i for i in range(1000)]
record = {"metadata": ["1", "1", "large", "C", "A"], "H": [2] * len(ids), "P": [[.25, .5, .25]] * len(ids)}
write_bgen(large / "large.bgen", ids, [record], 32, "zlib")
table(large / "cohort.localbgen.tsv", [["ancestry", "bgen", "counted_allele"], ["Only", "large.bgen", "ALT"]])
table(large / "samples.pheno", [["IID", "R", "E"], *[[iid, math.sin(i*.7), math.cos(i*.3)] for i, iid in enumerate(ids)]])
run("large_block_default", large / "cohort.localbgen.tsv", pheno=large / "samples.pheno")
run("reject_decoded_block_budget", large / "cohort.localbgen.tsv", pheno=large / "samples.pheno", extra=["--local-max-block-bytes", "4096"], success=False, error="decoded block exceeds")

summary = {"checks": len(checks), "passed": sum(r["pass"] for r in checks), "failures": failures,
           "binary_sha256": hashlib.sha256(BIN.read_bytes()).hexdigest(),
           "scope": "Known-H BGEN 1..32-bit raw/zlib/zstd: independent encoded probabilities, native four-trait/text equivalence, masks/direction, integrity and malformed inputs. No calibration claim."}
(OUT / "checks.json").write_text(json.dumps(checks, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
(OUT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2, ensure_ascii=True))
raise SystemExit(bool(failures))
