"""Native local-path regression: UTF-8 manifests and Chinese/emoji folders.

Usage: python -B tests/verify_local_unicode.py BINARY OUTPUT_DIRECTORY
All fixtures are independently generated SPAGxE examples, copied under OUTPUT.
No genotype inference or numerical-method changes are tested here.
"""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
from make_local_bgen_fixtures import raw_zstd

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "examples/local"
BINARY = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "tmp").mkdir(exist_ok=True)
checks = []
commands = 0


def check(case, passed, detail=""):
    checks.append({"case": case, "passed": bool(passed), "detail": detail})


def text_table(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.reader(stream, delimiter="\t"))


def save_table(path, rows, codec="plain"):
    # Fixture fields have no TAB/newline, and LF is part of the text profile.
    raw = ("\n".join("\t".join(row) for row in rows) + "\n").encode("utf-8")
    payload = gzip.compress(raw, mtime=0) if codec == "gz" else raw_zstd(raw) if codec == "zst" else raw
    path.write_bytes(payload)


def argument(path, style):
    path = Path(path)
    return str(path.resolve()) if style == "absolute" else os.path.relpath(path, ROOT)


def run(case, options, prefix, style="relative"):
    global commands
    commands += 1
    env = os.environ.copy()
    for name in ("TMP", "TEMP", "TMPDIR"):
        env[name] = str(OUT / "tmp")
    result = subprocess.run([str(BINARY), *options, "--out", argument(prefix, style)],
                            cwd=ROOT, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=120)
    (OUT / (case + ".console.txt")).write_bytes(result.stdout)
    check(case + ":exit", result.returncode == 0,
          result.stdout.decode("utf-8", errors="replace")[-2000:] if result.returncode else "")
    check(case + ":probe_cleanup", not list(prefix.parent.glob(prefix.name + "*.spagxe-probe-*")))
    return result.returncode == 0


def analyze(case, source, path, pheno, trait, style="relative", compression=None, unicode_output=True, extra=()):
    output_dir = OUT / ("\u7ed3\u679c_🧪" if unicode_output else "baseline")
    output_dir.mkdir(exist_ok=True)
    prefix = output_dir / case
    options = ["--method", "SPAGxEmix_local", source, argument(path, style),
               "--pheno", argument(pheno, style), "--resid-name", trait,
               "--envir-name", "E", "--geno", "1", "--maf", "0", "--mac", "0"]
    options += list(extra)
    if compression:
        options += ["--compression", compression]
    success = run(case, options, prefix, style)
    suffix = "." + compression if compression else ""
    paths = [Path(str(prefix) + "." + trait + "." + kind + suffix)
             for kind in ("SPAGxE_local", "SPAGxE_local_CCT")]
    for path in paths:
        check(case + ":exists:" + path.suffix, success and path.is_file())
    return paths


def equal(case, actual, expected):
    for index, (got, want) in enumerate(zip(actual, expected)):
        check(case + ":numeric_and_metadata_identity:" + str(index),
              got.is_file() and want.is_file() and got.read_bytes() == want.read_bytes())


baseline = {}
baseline["text"] = analyze("baseline_text", "--local-manifest", FIX / "cohort.local.tsv",
                           FIX / "samples.pheno", "R_quant", unicode_output=False)
baseline["binary"] = analyze("baseline_binary", "--local-prefix", FIX / "binary_demo/local",
                             FIX / "binary_demo/samples.pheno", "R", unicode_output=False)
baseline["bgen"] = analyze("baseline_bgen", "--local-bgen-manifest",
                           FIX / "bgen_demo/bits16_zlib.localbgen.tsv",
                           FIX / "samples.pheno", "R_quant", unicode_output=False)

# Chinese plus supplementary-plane characters in parent, manifest, table,
# phenotype and output names. Manifest contents are always UTF-8.
for codec in ("plain", "gz", "zst"):
    extension = "" if codec == "plain" else "." + codec
    folder = OUT / ("\u7956\u6e90\u6570\u636e_🧬_" + codec)
    folder.mkdir(exist_ok=True)
    pheno = folder / ("\u8868\u578b_🧪.tsv" + extension)
    save_table(pheno, text_table(FIX / "samples.pheno"), codec)
    original = text_table(FIX / "cohort.local.tsv")
    local_rows = [original[0]]
    for index, row in enumerate(original[1:]):
        entries = []
        for column, role in ((1, "\u57fa\u56e0\u578b"), (2, "\u7956\u6e90\u6570")):
            target = folder / (role + "_🧬_" + str(index) + ".tsv" + extension)
            save_table(target, text_table(FIX / row[column]), codec)
            entries.append(target)
        local_rows.append([row[0], *entries, row[3]])
    for entries_style in ("relative", "absolute"):
        rows = [local_rows[0]]
        for label, geno, hap, allele in local_rows[1:]:
            names = [path.as_posix() if entries_style == "absolute" else path.name for path in (geno, hap)]
            rows.append([label, *names, allele])
        manifest = folder / ("\u6e05\u5355_📁_" + entries_style + ".tsv" + extension)
        save_table(manifest, rows, codec)
        for cli_style in ("relative", "absolute"):
            case = "text_" + codec + "_" + entries_style + "_" + cli_style
            actual = analyze(case, "--local-manifest", manifest, pheno, "R_quant", cli_style)
            equal(case, actual, baseline["text"])

folder = OUT / "\u4e8c\u8fdb\u5236_🧬"
folder.mkdir(exist_ok=True)
prefix = folder / "\u5c40\u90e8\u7956\u6e90_📁"
for suffix in (".localbed", ".localbim", ".localfam"):
    shutil.copyfile(Path(str(FIX / "binary_demo/local") + suffix), Path(str(prefix) + suffix))
pheno = folder / "\u8868\u578b_🧪.tsv"
shutil.copyfile(FIX / "binary_demo/samples.pheno", pheno)
for style in ("relative", "absolute"):
    case = "binary_" + style
    equal(case, analyze(case, "--local-prefix", prefix, pheno, "R", style), baseline["binary"])

folder = OUT / "\u540e\u9a8cBGEN_🧬"
folder.mkdir(exist_ok=True)
pheno = folder / "\u8868\u578b_🧪.tsv"
shutil.copyfile(FIX / "samples.pheno", pheno)
original = text_table(FIX / "bgen_demo/bits16_zlib.localbgen.tsv")
files = []
for index, row in enumerate(original[1:]):
    target = folder / ("\u7956\u6e90\u540e\u9a8c_🧬_" + str(index) + ".bgen")
    shutil.copyfile(FIX / "bgen_demo" / row[1], target)
    files.append((row[0], target, row[2]))
for entries_style in ("relative", "absolute"):
    rows = [original[0]]
    for label, path, allele in files:
        rows.append([label, path.as_posix() if entries_style == "absolute" else path.name, allele])
    manifest = folder / ("\u540e\u9a8c\u6e05\u5355_📁_" + entries_style + ".tsv")
    save_table(manifest, rows)
    for cli_style in ("relative", "absolute"):
        case = "bgen_" + entries_style + "_" + cli_style
        equal(case, analyze(case, "--local-bgen-manifest", manifest, pheno, "R_quant", cli_style), baseline["bgen"])

# Exercise gzip/zstd writers at a Unicode output path and feed the bytes back
# through native TextReader. The CCT transport uses ASCII CLI basenames so
# this test does not enlarge the Unicode-CLI promise to non-local methods.
if baseline["text"][1].is_file():
    rows = text_table(baseline["text"][1])
    meta = ["CHROM", "POS", "ID", "REF", "ALT"]
    indices = [rows[0].index(name) for name in meta]
    global_file = OUT / "transport.R_quant.SPAGxE"
    save_table(global_file, [meta + ["LOG10P_GxE"]] +
               [[row[i] for i in indices] + ["0.3"] for row in rows[1:]])
    def combine(case, local):
        prefix = OUT / case
        options = ["--method", "SPAGxEcct", "--cct-global", argument(global_file, "relative"),
                   "--cct-local", argument(local, "relative"), "--cct-phenotype", "R_quant",
                   "--envir-name", "E"]
        ok = run(case, options, prefix)
        result = Path(str(prefix) + ".R_quant.SPAGxE_global_local_CCT")
        return result if ok else result
    expected = combine("transport_plain", baseline["text"][1])
    for codec in ("gz", "zst"):
        produced = analyze("unicode_output_" + codec, "--local-manifest", FIX / "cohort.local.tsv",
                           FIX / "samples.pheno", "R_quant", "absolute", codec)
        local = OUT / ("transport.R_quant.SPAGxE_local_CCT." + codec)
        if produced[1].is_file():
            shutil.copyfile(produced[1], local)
        else:
            continue
        actual = combine("transport_" + codec, local)
        check("unicode_output_" + codec + ":native_roundtrip_identity",
              actual.is_file() and expected.is_file() and actual.read_bytes() == expected.read_bytes())

# Header inference must use the same compressed Unicode-aware reader as
# explicit phenotype mode. Test phenotype and covariate auto-selection.
folder = OUT / "\u81ea\u52a8\u8868\u578b_🧪"
folder.mkdir(exist_ok=True)
original = text_table(FIX / "samples.pheno")
positions = [original[0].index(name) for name in ("Quantitative", "X", "E")]
combined = [[original[0][0], "Y", "X", "E"]] + [[row[0], *[row[i] for i in positions]] for row in original[1:]]
plain = OUT / "auto_reference.pheno"
save_table(plain, combined)

def fit_analysis(case, pheno, covar=None, auto_pheno=False, auto_covar=False):
    prefix = OUT / "\u7ed3\u679c_🧪" / case
    options = ["--method", "SPAGxEmix_local", "--local-manifest", argument(FIX / "cohort.local.tsv", "relative"),
               "--pheno", argument(pheno, "relative"), "--envir-name", "E",
               "--regression-model", "linear", "--geno", "1", "--maf", "0", "--mac", "0"]
    if not auto_pheno:
        options += ["--pheno-name", "Y"]
    if covar is not None:
        options += ["--covar", argument(covar, "relative")]
    if not auto_covar:
        options += ["--covar-name", "X,E"]
    ok = run(case, options, prefix)
    result = [Path(str(prefix) + ".Y." + kind) for kind in ("SPAGxE_local", "SPAGxE_local_CCT")]
    for path in result:
        check(case + ":exists:" + path.suffix, ok and path.is_file())
    return result

fit_baseline = fit_analysis("auto_header_baseline", plain)
for codec in ("gz", "zst"):
    pheno = folder / ("\u81ea\u52a8\u9009\u62e9_🧬.tsv." + codec)
    save_table(pheno, combined, codec)
    case = "auto_pheno_" + codec
    equal(case, fit_analysis(case, pheno, auto_pheno=True), fit_baseline)
pheno = folder / "\u7ed3\u5c40_🧬.tsv.gz"
covar = folder / "\u534f\u53d8\u91cf_📁.tsv.zst"
save_table(pheno, [[original[0][0], "Y"]] + [[row[0], row[positions[0]]] for row in original[1:]], "gz")
save_table(covar, [[original[0][0], "X", "E"]] + [[row[0], row[positions[1]], row[positions[2]]] for row in original[1:]], "zst")
equal("auto_covar_header", fit_analysis("auto_covar_header", pheno, covar, auto_covar=True), fit_baseline)

# Keep/remove and extract/exclude each exercise a Unicode file independently,
# then all four absolute Unicode paths are combined. Identical ASCII lists
# provide the numerical baseline for precisely the same selected cohort/SNPs.
folders = {"ascii": OUT / "filter_ascii", "unicode": OUT / "\u7b5b\u9009_🧬"}
for folder in folders.values():
    folder.mkdir(exist_ok=True)
ids = [row[0] for row in original[1:]]
manifest = text_table(FIX / "cohort.local.tsv")
geno = text_table(FIX / manifest[1][1])
id_column = geno[0].index("ID")
variants = [row[id_column] for row in geno[1:]]
contents = {"keep": ids[:-2], "remove": [ids[1]],
            "extract": variants[:-1] or variants,
            "exclude": [variants[1] if len(variants) > 2 else variants[-1]]}
files = {}
for role, values in contents.items():
    files[role] = {}
    for kind, folder in folders.items():
        target = folder / ((role if kind == "ascii" else "\u5217\u8868_📁_" + role) + ".txt")
        target.write_text("\n".join(values) + "\n", encoding="utf-8")
        files[role][kind] = target

def filter_flags(unicode_roles, style="relative"):
    options = []
    for role in ("keep", "remove", "extract", "exclude"):
        options += ["--" + role, argument(files[role]["unicode" if role in unicode_roles else "ascii"], style)]
    return options

filtered = analyze("filter_baseline", "--local-manifest", FIX / "cohort.local.tsv",
                   FIX / "samples.pheno", "R_quant", unicode_output=False, extra=filter_flags(set()))
for role in ("keep", "remove", "extract", "exclude"):
    case = "unicode_filter_" + role
    equal(case, analyze(case, "--local-manifest", FIX / "cohort.local.tsv",
                        FIX / "samples.pheno", "R_quant", extra=filter_flags({role})), filtered)
case = "unicode_filters_absolute"
equal(case, analyze(case, "--local-manifest", FIX / "cohort.local.tsv",
                    FIX / "samples.pheno", "R_quant", "absolute",
                    extra=filter_flags(set(contents), "absolute")), filtered)

# Preserve the existing global native-ACP contract for ambiguous byte paths.
# Use a relative CLI path so the project's Chinese absolute parent cannot
# accidentally make UTF-8 auto-detection fall back to ACP and hide the defect.
acp_collision = False
if os.name == "nt":
    try:
        encoded = "\u5784".encode("mbcs")
        acp_collision = encoded.decode("utf-8") != "\u5784"
    except (UnicodeError, LookupError):
        pass
if acp_collision:
    folder = OUT / "\u5784"
    folder.mkdir(exist_ok=True)
    pheno = folder / "pheno.tsv"
    shutil.copyfile(ROOT / "examples/synthetic/synthetic.pheno", pheno)
    global_tables = []
    for case, data, prefix in (
            ("global_native_ascii", ROOT / "examples/synthetic/synthetic.pheno", OUT / "global_native_ascii"),
            ("global_native_acp_collision", pheno, folder / "global_native_acp")):
        options = ["--method", "SPAGxE", "--bfile", "examples/synthetic/synthetic",
                   "--pheno", argument(data, "relative"), "--resid-name", "R_quant",
                   "--envir-name", "E", "--maf", "0", "--mac", "0"]
        ok = run(case, options, prefix)
        table = Path(str(prefix) + ".R_quant.SPAGxE")
        check(case + ":table", ok and table.is_file())
        global_tables.append(table)
    check("global_native_acp_collision:identity", all(x.is_file() for x in global_tables) and
          global_tables[0].read_bytes() == global_tables[1].read_bytes())

report = {"passed": all(row["passed"] for row in checks), "commands": commands,
          "checks": len(checks), "failure_count": sum(not row["passed"] for row in checks),
          "binary_sha256": hashlib.sha256(BINARY.read_bytes()).hexdigest(),
          "global_native_ACP_collision_case_exercised": acp_collision,
          "scope": "Local CLI/text/binary/known-H BGEN UTF-8 paths, Chinese/emoji names, absolute/relative manifest entries and CLI paths, plain/gzip/zstd text, automatic phenotype/covariate headers, keep/remove/extract/exclude, compressed-output native roundtrip; numerical identity, not calibration",
          "details": checks}
(OUT / "unicode_summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({key: report[key] for key in ("passed", "commands", "checks", "failure_count", "binary_sha256")}, indent=2))
if not report["passed"]:
    for row in checks:
        if not row["passed"]:
            print(json.dumps(row, ensure_ascii=True))
raise SystemExit(0 if report["passed"] else 1)
