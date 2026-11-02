"""End-to-end native local ancestry I/O, mask, QC and format verification.

No R or third-party Python packages are required. Run from any directory:
  python -B tests/verify_local.py [SPAGxE-binary] [new-output-directory]
Uses independently simulated count-container fixtures and paired TSV.
Numerical R oracle comparisons are a separate verification stage.
"""
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import zlib

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "build/bin" / ("SPAGxE.exe" if os.name == "nt" else "SPAGxE")
OUT = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else ROOT.parent.parent / "validation/local_io"
FIX = ROOT / "examples/local"
OUT.mkdir(parents=True, exist_ok=True)
checks, failures = [], []


def check(name, ok, detail=""):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    if not ok:
        failures.append(name + ": " + str(detail))


def read_table(path):
    with open(path, encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def command(tag, manifest=None, binary=None, native=False, extra=(), pheno=None, success=True):
    prefix = OUT / tag
    pheno = pheno or (FIX / "samples.pheno")
    relative = lambda path: os.path.relpath(path, ROOT)
    cmd = [str(BIN), "--method", "SPAGxElocal", "--pheno", relative(pheno), "--envir-name", "E", "--maf", "0", "--mac", "0", "--geno", "1", "--out", relative(prefix)]
    if native:
        cmd += ["--pheno-name", "Quantitative,Binary,Ordinal,Time:Event", "--covar-name", "X,E", "--regression-model", "auto"]
    else:
        cmd += ["--resid-name", "R"] if pheno.name == "samples.pheno" and pheno.parent.name == "binary_demo" else ["--resid-name", "R_quant,R_binary,R_ordinal,R_survival,R_masked"]
    if binary:
        cmd += ["--local-binary", relative(binary)]
    else:
        cmd += ["--local-manifest", relative(manifest or FIX / "cohort.local.tsv")]
    cmd += [relative(value) if Path(value).exists() and not str(value).startswith("--") else value for value in extra]
    result = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", timeout=120)
    (OUT / (tag + ".console.txt")).write_text(result.stdout, encoding="utf-8")
    check(tag + "_exit", (result.returncode == 0) == success, result.stdout[-1200:] if (result.returncode == 0) != success else result.returncode)
    return prefix, result


def local(prefix, trait="R_quant"):
    return Path(str(prefix) + "." + trait + ".SPAGxE_local")


def compare(tag, a, b, traits, tolerance=None, sort_rows=False):
    def same_rows(left, right):
        if len(left) != len(right):
            return False
        if sort_rows:
            key = lambda row: (row["CHROM"], int(row["POS"]), row["ID"], row.get("ANCESTRY", ""))
            left, right = sorted(left, key=key), sorted(right, key=key)
        for x, y in zip(left, right):
            if x.keys() != y.keys():
                return False
            for key in x:
                if x[key] == y[key]:
                    continue
                try:
                    u, v = float(x[key]), float(y[key])
                except ValueError:
                    return False
                if not ((math.isnan(u) and math.isnan(v)) or abs(u-v) <= (tolerance or 0) * max(1, abs(u), abs(v))):
                    return False
        return True
    for trait in traits:
        check(tag + "_" + trait, same_rows(read_table(local(a, trait)), read_table(local(b, trait))) if tolerance is not None or sort_rows else local(a, trait).read_bytes() == local(b, trait).read_bytes())
        ac = Path(str(a) + "." + trait + ".SPAGxE_local_CCT")
        bc = Path(str(b) + "." + trait + ".SPAGxE_local_CCT")
        check(tag + "_CCT_" + trait, same_rows(read_table(ac), read_table(bc)) if tolerance is not None or sort_rows else ac.read_bytes() == bc.read_bytes())


def mutate_case(tag, mutator):
    directory = OUT / ("input_" + tag)
    shutil.copytree(FIX, directory, ignore=shutil.ignore_patterns("binary_demo", "bgen_demo", "three_ancestry"), dirs_exist_ok=True)
    mutator(directory)
    return directory / "cohort.local.tsv"


def transform_table(path, fn):
    with open(path, encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    fn(rows)
    with open(path, "w", encoding="utf-8", newline="") as stream:
        csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)


def single_raw_block_fixture(source, target):
    """Write legal SQRLOC01 with all supplied rows in one zstd RAW frame.

    This independent stdlib encoder deliberately avoids a compression-library
    dependency. Descriptor 0xA0 means single segment, four-byte content size,
    no dictionary and no frame checksum. The last RAW block header is the
    little-endian uint24 (raw_length << 3) | 1.
    """
    target.mkdir(parents=True, exist_ok=True)
    entries = read_table(source / "cohort.local.tsv")
    def counts(path):
        return [{("CHROM" if key == "CHR" else key): value for key, value in row.items()}
                for row in read_table(path)]
    branches = [(counts(source / entry["genotype"]),
                 counts(source / entry["hapcount"])) for entry in entries]
    keys = ["CHROM", "POS", "ID", "REF", "ALT"]
    ids = list(branches[0][0][0])[5:]
    n, m, k = len(ids), len(branches[0][0]), len(entries)
    plane = (n + 3) // 4
    raw = bytearray()
    bim_rows = []
    for j in range(m):
        metadata = [branches[0][0][j][key] for key in keys]
        bim_rows.append("\t".join([metadata[0], metadata[2], "0", metadata[1], metadata[4], metadata[3]]) + "\n")
        for genotype, hapcount in branches:
            for rows in (genotype, hapcount):
                assert len(rows) == m and [rows[j][key] for key in keys] == metadata
                assert set(rows[j]) == set(keys + ids)
                packed = bytearray(plane)
                for i, iid in enumerate(ids):
                    token = rows[j][iid]
                    if token in ("", "NA", ".", "-9", "NaN", "nan", "NAN"):
                        value = 3
                    else:
                        number = float(token)
                        assert number in (0, 1, 2), "RAW fixture is hardcall-only"
                        value = int(number)
                    packed[i // 4] |= value << (2 * (i % 4))
                raw.extend(packed)
    assert len(raw) == m * 2 * k * plane and len(raw) < 128 * 1024
    frame = (b"\x28\xb5\x2f\xfd\xa0" + struct.pack("<I", len(raw))
             + ((len(raw) << 3) | 1).to_bytes(3, "little") + raw)
    labels = b"".join(struct.pack("<H", len(entry["ancestry"].encode("utf-8")))
                      + entry["ancestry"].encode("utf-8") for entry in entries)
    directions = {entry["counted_allele"] for entry in entries}
    assert len(directions) == 1 and next(iter(directions)) in ("ALT", "REF")
    direction = int(next(iter(directions)) == "REF")
    offset = 41 + len(labels)
    footer = offset + 24 + len(frame)
    header = b"SQRLOC01" + struct.pack("<QQIIQB", n, m, k, m, footer, direction) + labels
    bim = "".join(bim_rows).encode("utf-8")
    fam = "".join("0\t%s\t0\t0\t0\t-9\n" % iid for iid in ids).encode("utf-8")
    block = struct.pack("<IQQI", m, len(raw), len(frame), zlib.crc32(raw)) + frame
    index = struct.pack("<QQII", 1, offset, zlib.crc32(bim), zlib.crc32(fam))
    prefix = target / "single_block"
    prefix.with_suffix(".localbed").write_bytes(header + block + index)
    prefix.with_suffix(".localbim").write_bytes(bim)
    prefix.with_suffix(".localfam").write_bytes(fam)
    return prefix, len(raw)


resid_traits = ["R_quant", "R_binary", "R_ordinal", "R_survival", "R_masked"]
native_traits = ["Quantitative", "Binary", "Ordinal", "Time_Event"]
a, process = command("text_threads1", extra=["--threads", "1"])
if process.returncode:
    print(json.dumps({"checks": len(checks), "failures": failures}, indent=2))
    raise SystemExit(1)
b, _ = command("text_threads3", extra=["--threads", "3", "--chunk-ksnp", "0.5"])
compare("threads_and_chunk", a, b, resid_traits)
native, native_process = command("native_four_traits", native=True, extra=["--threads", "2", "--save-resid"])
if native_process.returncode == 0:
    for trait in native_traits:
        rows = read_table(local(native, trait))
        check("native_" + trait + "_rows", len(rows) == 16)
        check("native_" + trait + "_has_p", any(r["QC_STATUS"] == "PASS" and r["LOG10P_GxE"] != "NA" for r in rows))
three, _ = command("three_ancestry", manifest=FIX / "three_ancestry/cohort.local.tsv")
check("three_ancestry_rows", len(read_table(local(three))) == 24)
rows = read_table(local(a))
by = {(r["ID"], r["ANCESTRY"]): r for r in rows}
check("high_AF_flips", all(by["high_AF", anc]["TESTED_ALLELE"] == "REF" and float(by["high_AF", anc]["ALT_FREQ"]) > 0.5 for anc in ("EUR", "EAS")))
check("missing_H_joint_mask", all(int(by["H_missing", anc]["N_H_MISSING"]) == 2 and int(by["H_missing", anc]["N_TESTED"]) == 158 for anc in ("EUR", "EAS")))
check("structural_zero_G_restore", float(by["structural_zero", "EUR"]["MISS_RATE"]) == 0 and by["structural_zero", "EUR"]["QC_STATUS"] == "FAIL_QC")
check("dosage_HWE_unavailable", all(by["dosage", anc]["LOG10P_HWE_GLOBAL"] == "NA" for anc in ("EUR", "EAS")))
check("dosage_default_passes", all(by["dosage", anc]["QC_STATUS"] == "PASS" for anc in ("EUR", "EAS")))
masked = read_table(local(a, "R_masked"))
check("per_phenotype_mask", int(masked[0]["N_TESTED"]) == 142)
hw, _ = command("dosage_HWE_filter", extra=["--hwe", "0.0001"])
check("dosage_HWE_filter_status", all(r["QC_STATUS"] == "FAIL_HWE_UNAVAILABLE" for r in read_table(local(hw)) if r["ID"] == "dosage"))
filtered, _ = command("marker_filters", extra=["--extract", str(FIX / "extract.txt"), "--exclude", str(FIX / "exclude.txt"), "--chr", "1-2"])
check("extract_exclude", {r["ID"] for r in read_table(local(filtered))} == {"high_AF", "G_missing"})
kept, _ = command("sample_filters", extra=["--keep", str(FIX / "keep.txt"), "--remove", str(FIX / "remove.txt")])
check("keep_remove", int(read_table(local(kept))[0]["N_TESTED"]) == 104)

def reordered(d):
    for name in ("EUR.G.tsv", "EUR.H.tsv", "EAS.G.tsv", "EAS.H.tsv"):
        def reverse(rows):
            for row in rows:
                row[5:] = list(reversed(row[5:]))
        transform_table(d / name, reverse)

perm, _ = command("permuted_IIDs", manifest=mutate_case("permuted", reordered))
compare("IID_alignment", a, perm, resid_traits, tolerance=2e-12)

def compressed_inputs(d):
    for name in ("EUR.G.tsv", "EUR.H.tsv", "EAS.G.tsv", "EAS.H.tsv"):
        with gzip.open(d / (name + ".gz"), "wb") as stream:
            stream.write((d / name).read_bytes())
    transform_table(d / "cohort.local.tsv", lambda rows: [row.__setitem__(j, row[j] + ".gz") for row in rows[1:] for j in (1, 2)])

gz, _ = command("gzip_inputs", manifest=mutate_case("gzip", compressed_inputs))
compare("gzip_input_equivalence", a, gz, resid_traits)
gzout, _ = command("gzip_output", extra=["--compression", "gz"])
for trait in resid_traits:
    check("gzip_output_" + trait, gzip.decompress(Path(str(local(gzout, trait)) + ".gz").read_bytes()) == local(a, trait).read_bytes())

def complement_ref(d):
    for anc in ("EUR", "EAS"):
        with open(d / (anc + ".H.tsv"), encoding="utf-8", newline="") as stream:
            h = list(csv.reader(stream, delimiter="\t"))
        def complement(rows):
            for j, row in enumerate(rows[1:], 1):
                for i in range(5, len(row)):
                    if row[i] != "NA" and h[j][i] != "NA":
                        row[i] = "%.17g" % (float(h[j][i]) - float(row[i]))
        transform_table(d / (anc + ".G.tsv"), complement)
    transform_table(d / "cohort.local.tsv", lambda rows: [row.__setitem__(3, "REF") for row in rows[1:]])

ref, _ = command("REF_inputs", manifest=mutate_case("ref", complement_ref))
# Integer fields are exact. Fractional REF->ALT incurs the inevitable
# subtraction round trip; scientific values must agree to tight tolerance.
for trait in resid_traits:
    left, right = read_table(local(a, trait)), read_table(local(ref, trait))
    exact, close = True, True
    for x, y in zip(left, right):
        for key in x:
            if x["ID"] != "dosage":
                exact &= x[key] == y[key]
            else:
                try:
                    u, v = float(x[key]), float(y[key])
                    close &= (math.isnan(u) and math.isnan(v)) or abs(u-v) <= 2e-9 * max(1, abs(u), abs(v))
                except ValueError:
                    close &= x[key] == y[key]
    check("REF_HminusG_integer_" + trait, exact)
    check("REF_HminusG_dosage_" + trait, close)

def greater_than_h(d):
    with open(d / "EUR.H.tsv", encoding="utf-8", newline="") as stream:
        h = list(csv.reader(stream, delimiter="\t"))
    index = next(i for i in range(5, len(h[1])) if h[1][i] in ("0", "1"))
    transform_table(d / "EUR.G.tsv", lambda rows: rows[1].__setitem__(index, str(int(h[1][index]) + 1)))


bad_cases = {
    "G_greater_H": greater_than_h,
    "G_out_of_range": lambda d: transform_table(d / "EUR.G.tsv", lambda rows: rows[1].__setitem__(5, "3")),
    "negative_G": lambda d: transform_table(d / "EUR.G.tsv", lambda rows: rows[1].__setitem__(5, "-1")),
    "fractional_H": lambda d: transform_table(d / "EUR.H.tsv", lambda rows: rows[1].__setitem__(5, "0.5")),
    "duplicate_IID": lambda d: transform_table(d / "EUR.G.tsv", lambda rows: rows[0].__setitem__(6, rows[0][5])),
    "absent_IID": lambda d: transform_table(d / "EUR.G.tsv", lambda rows: rows[0].__setitem__(6, "not_in_other_files")),
    "metadata_mismatch": lambda d: transform_table(d / "EAS.H.tsv", lambda rows: rows[1].__setitem__(2, "different_variant")),
    "duplicate_variant": lambda d: transform_table(d / "EUR.G.tsv", lambda rows: rows[2].__setitem__(2, rows[1][2])),
    "mixed_counted_alleles": lambda d: transform_table(d / "cohort.local.tsv", lambda rows: rows[1].__setitem__(3, "REF")),
    "H_sum_not_two": lambda d: transform_table(d / "EUR.H.tsv", lambda rows: rows[1].__setitem__(5, "0" if rows[1][5] != "0" else "1")),
    "row_shape": lambda d: transform_table(d / "EUR.G.tsv", lambda rows: rows[1].pop()),
    "allele_mismatch": lambda d: transform_table(d / "EUR.H.tsv", lambda rows: rows[1].__setitem__(4, "T")),
    "nonfinite_G": lambda d: transform_table(d / "EUR.G.tsv", lambda rows: rows[1].__setitem__(5, "Inf")),
    "invalid_header": lambda d: transform_table(d / "EUR.G.tsv", lambda rows: rows[0].__setitem__(1, "POSITION")),
    "duplicate_ancestry": lambda d: transform_table(d / "cohort.local.tsv", lambda rows: rows[2].__setitem__(0, rows[1][0])),
    "reserved_ancestry": lambda d: transform_table(d / "cohort.local.tsv", lambda rows: rows[2].__setitem__(0, "global")),
}
for tag, fn in bad_cases.items():
    command("reject_" + tag, manifest=mutate_case(tag, fn), success=False)

count_demo = FIX / "binary_demo"
count_text, _ = command("count_demo_text", manifest=count_demo / "cohort.local.tsv", pheno=count_demo / "samples.pheno")
count_binary, _ = command("count_demo_binary", binary=count_demo / "local", pheno=count_demo / "samples.pheno")
compare("simulated_count_binary_text", count_text, count_binary, ["R"], tolerance=0, sort_rows=True)
count_rows = len(read_table(local(count_binary, "R")))
check("simulated_count_binary_rows", count_rows == 36)
single_block, raw_bytes = single_raw_block_fixture(count_demo, OUT / "binary_single_raw_block")
check("single_raw_block_payload_size", raw_bytes == 8136 and raw_bytes > 4096, raw_bytes)
raw_output, _ = command("count_demo_single_RAW_block", binary=single_block, pheno=count_demo / "samples.pheno")
compare("single_RAW_block_text", count_text, raw_output, ["R"], tolerance=0, sort_rows=True)
_, budget_error = command("reject_binary_block_budget", binary=single_block, pheno=count_demo / "samples.pheno", extra=["--local-max-block-bytes", "4096"], success=False)
check("decoded_block_budget_error", budget_error.returncode != 0 and "decoded block exceeds" in budget_error.stdout,
      budget_error.stdout[-1200:] if "decoded block exceeds" not in budget_error.stdout else "decoded block exceeds")

def binary_bad(tag, change):
    target = OUT / ("binary_" + tag)
    target.mkdir(exist_ok=True)
    for ext in (".localbed", ".localbim", ".localfam"):
        shutil.copyfile(count_demo / ("local" + ext), target / ("local" + ext))
    change(target)
    command("reject_binary_" + tag, binary=target / "local", pheno=count_demo / "samples.pheno", success=False)


def bytes_change(directory, fn):
    p = directory / "local.localbed"
    data = bytearray(p.read_bytes())
    fn(data)
    p.write_bytes(data)


binary_bad("sample_sidecar_CRC", lambda d: (d / "local.localfam").write_bytes((d / "local.localfam").read_bytes() + b"\n"))
binary_bad("variant_sidecar_CRC", lambda d: (d / "local.localbim").write_bytes((d / "local.localbim").read_bytes() + b"\n"))
binary_bad("truncated", lambda d: bytes_change(d, lambda b: b.__delitem__(slice(len(b)-10, None))))
binary_bad("trailing_bytes", lambda d: bytes_change(d, lambda b: b.extend(b"bad")))
binary_bad("bad_magic", lambda d: bytes_change(d, lambda b: b.__setitem__(0, ord("X"))))
binary_bad("huge_N", lambda d: bytes_change(d, lambda b: b.__setitem__(slice(8, 16), struct.pack("<Q", 2**40))))
binary_bad("bad_direction", lambda d: bytes_change(d, lambda b: b.__setitem__(40, 2)))
binary_bad("corrupt_frame", lambda d: bytes_change(d, lambda b: b.__setitem__(100, b[100] ^ 0xA5)))
binary_bad("footer_location", lambda d: bytes_change(d, lambda b: b.__setitem__(slice(32, 40), struct.pack("<Q", 10))))

summary = {"checks": len(checks), "passed": sum(c["pass"] for c in checks), "failures": failures, "binary_sha256": hashlib.sha256(BIN.read_bytes()).hexdigest(), "count_demo_binary_sha256": hashlib.sha256((count_demo / "local.localbed").read_bytes()).hexdigest(), "scope": "Native end-to-end format, masks, deterministic threading and QC. Scientific R numerical oracle is separate."}
(OUT / "checks.json").write_text(json.dumps(checks, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
(OUT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2, ensure_ascii=True))
raise SystemExit(bool(failures))
