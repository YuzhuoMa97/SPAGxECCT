"""Deterministic native local ancestry examples; stdlib only, all writes explicit.

Run `python -B tests/make_local_fixtures.py [output]`.
All examples are generated from this program's fixed simulation seeds.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import struct
import sys
import zlib

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "examples/local"
OUT.mkdir(parents=True, exist_ok=True)
rng = random.Random(20261024)
n = 160
ids = ["s%03d" % (i + 1) for i in range(n)]
env = [rng.gauss(0, 1) for _ in ids]
x = [rng.gauss(0, 1) for _ in ids]


def write(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as stream:
        csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)


def binomial(h, q):
    return sum(rng.random() < q for _ in range(h))


def data(k):
    labels = ["EUR", "EAS", "AFR"][:k]
    g, h = [], []
    names = ["common_null", "main_effect", "high_AF", "rare", "G_missing", "H_missing", "structural_zero", "dosage"]
    for j, name in enumerate(names):
        hj = [[0] * n for _ in labels]
        gj = [[0] * n for _ in labels]
        for i in range(n):
            for _ in range(2):
                a = rng.randrange(k)
                hj[a][i] += 1
            if name == "structural_zero":
                hj[0][i] = 0
                for a in range(1, k):
                    hj[a][i] = 0
                hj[-1][i] = 2
            for a in range(k):
                q = 0.8 if name == "high_AF" else 0.025 if name == "rare" else 0.22 + 0.12 * a
                gj[a][i] = binomial(hj[a][i], q)
                if name == "dosage" and hj[a][i] > 0:
                    gj[a][i] = hj[a][i] * (0.2 + 0.1 * a + 0.03 * math.sin(i))
        if name == "G_missing":
            for i in (1, 7, 22, 59):
                gj[0][i] = "NA"
        if name == "H_missing":
            for i in (6, 13):
                hj[0][i] = "NA"
        if name == "structural_zero":
            gj[0][3] = "NA"
        g.append(gj)
        h.append(hj)
    return labels, names, g, h


def save(directory, k):
    directory.mkdir(parents=True, exist_ok=True)
    labels, names, g, h = data(k)
    manifest = [["ancestry", "genotype", "hapcount", "counted_allele"]]
    for a, label in enumerate(labels):
        for array, tag in ((g, "G"), (h, "H")):
            rows = [["CHROM", "POS", "ID", "REF", "ALT", *ids]]
            rows += [[str(1 + j % 2), str(100 + j), name, "C", "A", *array[j][a]] for j, name in enumerate(names)]
            write(directory / (label + "." + tag + ".tsv"), rows)
        manifest.append([label, label + ".G.tsv", label + ".H.tsv", "ALT"])
    write(directory / "cohort.local.tsv", manifest)
    return g, h


g, h = save(OUT, 2)
save(OUT / "three_ancestry", 3)
quantitative, binary, ordinal, time, event = [], [], [], [], []
for i in range(n):
    main = g[1][0][i]
    quantitative.append(4 * main + 0.4 * x[i] + 0.5 * env[i] + rng.gauss(0, 0.8))
    prob = 1 / (1 + math.exp(-(-2 + 2.8 * main + 0.4 * x[i] + 0.4 * env[i])))
    binary.append(int(rng.random() < prob))
    latent = 2.2 * main + 0.4 * x[i] + 0.3 * env[i] + math.log((u := rng.random()) / (1 - u))
    ordinal.append(int(latent > 0.2) + int(latent > 2.2))
    time.append(-math.log(rng.random()) * math.exp(-1.0 * main - 0.2 * env[i] - 0.1 * x[i]))
    event.append(int(rng.random() < 0.8))
rows = [["IID", "Quantitative", "Binary", "Ordinal", "Time", "Event", "X", "E", "Eshift", "R_quant", "R_binary", "R_ordinal", "R_survival", "R_masked"]]
for i, iid in enumerate(ids):
    r = quantitative[i] - 0.4 * x[i] - 0.5 * env[i]
    rb = binary[i] - 1 / (1 + math.exp(-(-1 + 0.4 * x[i] + 0.4 * env[i])))
    ro = ordinal[i] - 0.3 * x[i] - 0.2 * env[i] - 0.7
    rs = event[i] - time[i] * math.exp(0.2 * env[i] + 0.1 * x[i])
    rows.append([iid, quantitative[i], binary[i], ordinal[i], time[i], event[i], x[i], env[i], env[i] + 8, r, rb, ro, rs, "NA" if i % 9 == 0 else r])
write(OUT / "samples.pheno", rows)
write(OUT / "keep.txt", [[iid] for i, iid in enumerate(ids) if i % 3 != 0])
write(OUT / "remove.txt", [[iid] for iid in ids[:4]])
write(OUT / "extract.txt", [["main_effect"], ["high_AF"], ["G_missing"]])
write(OUT / "exclude.txt", [["main_effect"]])

def zstd_raw(payload):
    assert len(payload) < 128 * 1024
    return (b"\x28\xb5\x2f\xfd\xa0" + struct.pack("<I", len(payload))
            + ((len(payload) << 3) | 1).to_bytes(3, "little") + payload)


def binary_demo():
    """Independent count-container encoder, with no external source files."""
    target = OUT / "binary_demo"
    target.mkdir(exist_ok=True)
    simulation = random.Random(20261025)
    sample_ids = ["demo%04d" % (i + 1) for i in range(450)]
    labels, m, block = ["Ancestry1", "Ancestry2"], 18, 4
    tables = [[["CHROM", "POS", "ID", "REF", "ALT", *sample_ids]] for _ in range(4)]
    payloads, bim_rows = [], []
    plane = (len(sample_ids) + 3) // 4
    for j in range(m):
        metadata = [str(j % 2 + 1), str(j + 1001), "count_%02d" % (j + 1), "G", "T"]
        values = [[] for _ in range(4)]
        for i in range(len(sample_ids)):
            h1 = sum(simulation.random() < .45 for _ in range(2))
            h2 = 2 - h1
            g1 = sum(simulation.random() < .12 + .03 * (j % 7) for _ in range(h1))
            g2 = sum(simulation.random() < .24 + .02 * (j % 5) for _ in range(h2))
            for row, value in zip(values, (g1, h1, g2, h2)):
                row.append(value)
        if j == 3:
            values[0][7] = "NA"
        packed = bytearray()
        for table, row in zip(tables, values):
            table.append(metadata + row)
            codes = bytearray(plane)
            for i, value in enumerate(row):
                codes[i // 4] |= (3 if value == "NA" else value) << (2 * (i % 4))
            packed.extend(codes)
        payloads.append(packed)
        bim_rows.append("\t".join([metadata[0], metadata[2], "0", metadata[1], metadata[4], metadata[3]]) + "\n")
    manifest = [["ancestry", "genotype", "hapcount", "counted_allele"]]
    for a, label in enumerate(labels):
        write(target / (label + ".G.tsv"), tables[2 * a])
        write(target / (label + ".H.tsv"), tables[2 * a + 1])
        manifest.append([label, label + ".G.tsv", label + ".H.tsv", "ALT"])
    write(target / "cohort.local.tsv", manifest)
    fam = "".join("0\t%s\t0\t0\t0\t-9\n" % iid for iid in sample_ids).encode("utf-8")
    bim = "".join(bim_rows).encode("utf-8")
    labels_raw = b"".join(struct.pack("<H", len(label)) + label.encode("ascii") for label in labels)
    offsets, blocks = [], []
    cursor = 41 + len(labels_raw)
    for start in range(0, m, block):
        payload = b"".join(payloads[start:start + block])
        frame = zstd_raw(payload)
        encoded = struct.pack("<IQQI", min(block, m - start), len(payload), len(frame), zlib.crc32(payload)) + frame
        offsets.append(cursor)
        blocks.append(encoded)
        cursor += len(encoded)
    header = b"SQRLOC01" + struct.pack("<QQIIQB", len(sample_ids), m, len(labels), block, cursor, 0) + labels_raw
    index = struct.pack("<Q", len(offsets)) + b"".join(struct.pack("<Q", offset) for offset in offsets)
    index += struct.pack("<II", zlib.crc32(bim), zlib.crc32(fam))
    (target / "local.localbed").write_bytes(header + b"".join(blocks) + index)
    (target / "local.localbim").write_bytes(bim)
    (target / "local.localfam").write_bytes(fam)
    write(target / "samples.pheno", [["IID", "R", "E"], *[[iid, math.sin(i * .7) + .2 * math.cos(i * .2), math.cos(i * .31) + .01 * i] for i, iid in enumerate(sample_ids)]])
    files = ["local.localbed", "local.localbim", "local.localfam"]
    (target / "SIMULATION.json").write_text(json.dumps({"generator": "SPAGxE tests/make_local_fixtures.py", "seed": 20261025, "samples": len(sample_ids), "variants": m, "ancestries": labels, "sha256": {name: hashlib.sha256((target / name).read_bytes()).hexdigest() for name in files}}, indent=2) + "\n", encoding="utf-8")


binary_demo()
print(OUT)
