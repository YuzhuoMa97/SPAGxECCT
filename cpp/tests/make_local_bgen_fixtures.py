"""Create SPAGxE-owned synthetic known-H BGEN collections using stdlib only.

This is a fixture encoder, not an imputation tool. Probability vectors are
simulated directly. Expected text values are calculated from the quantized
stored probabilities rather than reconstructed from pre-existing dosages.
"""
import csv
import json
import math
from pathlib import Path
import random
import struct
import sys
import zlib

ROOT = Path(__file__).resolve().parents[1]


def table(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as stream:
        csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)


def raw_zstd(payload):
    assert len(payload) < 128 * 1024
    return (b"\x28\xb5\x2f\xfd\xa0" + struct.pack("<I", len(payload))
            + ((len(payload) << 3) | 1).to_bytes(3, "little") + payload)


def probabilities(p, bits):
    denominator = (1 << bits) - 1
    scaled = [x * denominator / sum(p) for x in p]
    result = [math.floor(x) for x in scaled]
    for i in sorted(range(len(p)), key=lambda i: (scaled[i] - result[i], -i), reverse=True)[:denominator - sum(result)]:
        result[i] += 1
    if denominator >= sum(x > 0 for x in p):
        for i in range(len(p)):
            if p[i] > 0 and result[i] == 0:
                donor = max(range(len(p)), key=lambda j: result[j])
                assert result[donor] > 1
                result[donor] -= 1
                result[i] += 1
    assert sum(result) == denominator
    return result


def string(value, width):
    encoded = value.encode("utf-8")
    return len(encoded).to_bytes(width, "little") + encoded


def write_bgen(path, ids, records, bits, codec="zlib", embedded=True):
    sample_block = struct.pack("<II", 8 + sum(2 + len(i.encode("utf-8")) for i in ids), len(ids))
    sample_block += b"".join(string(i, 2) for i in ids)
    if not embedded:
        sample_block = b""
    compression = {"raw": 0, "zlib": 1, "zstd": 2}[codec]
    flags = compression | (2 << 2) | (int(embedded) << 31)
    encoded = bytearray(struct.pack("<IIII", 20 + len(sample_block), 20, len(records), len(ids)) + b"bgen" + struct.pack("<I", flags) + sample_block)
    offsets, expected = [], []
    for record in records:
        chrom, pos, ident, ref, alt = record["metadata"]
        begin = len(encoded)
        encoded += string(ident + "_snp", 2) + string(ident, 2)
        chromosome_offset = len(encoded) + 2
        encoded += string(chrom, 2)
        position_offset = len(encoded)
        encoded += struct.pack("<IH", int(pos), 2) + string(ref, 4) + string(alt, 4)
        h, posterior = record["H"], record["P"]
        ploidy = bytes(int(x) | (128 if p is None else 0) for x, p in zip(h, posterior))
        raw = bytearray(struct.pack("<IHBB", len(ids), 2, min(h), max(h)) + ploidy + bytes([0, bits]))
        values, mean = [], []
        for copies, p in zip(h, posterior):
            if p is None:
                stored = [0] * copies
                mean.append("NA" if copies else 0)
            else:
                quantized = probabilities(p, bits)
                assert len(quantized) == copies + 1
                stored = quantized[:-1]
                mean.append(sum(g * q for g, q in enumerate(quantized)) / ((1 << bits) - 1))
            values.extend(stored)
        packed, accumulator, available = bytearray(), 0, 0
        for value in values:
            accumulator |= value << available
            available += bits
            while available >= 8:
                packed.append(accumulator & 255)
                accumulator >>= 8
                available -= 8
        if available:
            packed.append(accumulator)
        raw += packed
        block_offset = len(encoded)
        if codec == "raw":
            encoded += struct.pack("<I", len(raw))
            probability_offset = len(encoded)
            encoded += raw
        else:
            compressed = zlib.compress(raw, 3) if codec == "zlib" else raw_zstd(raw)
            encoded += struct.pack("<II", len(compressed) + 4, len(raw))
            probability_offset = None
            encoded += compressed
        offsets.append({"identification": begin, "chromosome": chromosome_offset, "position": position_offset,
                        "block": block_offset, "raw_probability": probability_offset, "decoded_bytes": len(raw)})
        expected.append(mean)
    path.write_bytes(encoded)
    return expected, offsets


def simulated_records(n, k, seed=20261025):
    rng = random.Random(seed)
    ids = ["s%03d" % (i + 1) for i in range(n)]
    names = ["uncertain", "point_calls", "high_AF", "G_missing", "H_zero", "integer_mean", "zero_alleles", "rare", "uncertain_2", "point_calls_2"]
    branches = [[] for _ in range(k)]
    for j, name in enumerate(names):
        h = [[0] * n for _ in range(k)]
        p = [[] for _ in range(k)]
        for i in range(n):
            for _ in range(2):
                h[rng.randrange(k)][i] += 1
            if name == "H_zero":
                for a in range(k):
                    h[a][i] = 2 if a == k - 1 else 0
            if name == "integer_mean":
                for a in range(k):
                    h[a][i] = 2 if a == 0 else 0
            for a in range(k):
                copies = h[a][i]
                q = .8 if name == "high_AF" else .008 if name == "rare" else .15 + .12 * a + .025 * math.sin(i + j)
                if name == "zero_alleles":
                    prob = [1] + [0] * copies
                elif name == "integer_mean" and copies == 2:
                    prob = [.25, .5, .25]
                elif name.startswith("point_calls"):
                    g = sum(rng.random() < q for _ in range(copies))
                    prob = [int(t == g) for t in range(copies + 1)]
                else:
                    prob = [math.comb(copies, g) * q ** g * (1-q) ** (copies-g) for g in range(copies + 1)]
                if name == "G_missing" and i % 31 == 0:
                    prob = None
                p[a].append(prob)
        for a in range(k):
            branches[a].append({"metadata": [str(j % 2 + 1), str(20001+j), name, "C", "A"], "H": h[a], "P": p[a]})
    return ids, branches


def generate(output):
    output.mkdir(parents=True, exist_ok=True)
    ids, branches = simulated_records(160, 3)
    labels = ["Ancestry1", "Ancestry2", "Ancestry3"]
    table(output / "samples.pheno", [["IID", "R", "E"], *[[iid, math.sin(i * .73) + .3 * math.cos(i * .31), math.cos(i * .17)] for i, iid in enumerate(ids)]])
    recipes = []
    for bits in (1, 8, 16, 24, 32):
        for codec in ("raw", "zlib", "zstd"):
            stem = "bits%d_%s" % (bits, codec)
            manifest = [["ancestry", "bgen", "counted_allele"]]
            text_manifest = [["ancestry", "genotype", "hapcount", "counted_allele"]]
            info = {}
            for a, label in enumerate(labels):
                order = list(range(len(ids))) if a == 0 else list(reversed(range(len(ids)))) if a == 1 else list(range(17, len(ids))) + list(range(17))
                ordered = [{"metadata": row["metadata"], "H": [row["H"][i] for i in order], "P": [row["P"][i] for i in order]} for row in branches[a]]
                binary = stem + "." + label + ".bgen"
                expected, offsets = write_bgen(output / binary, [ids[i] for i in order], ordered, bits, codec)
                reverse = {original: i for i, original in enumerate(order)}
                gt, ht = [["CHROM", "POS", "ID", "REF", "ALT", *ids]], [["CHROM", "POS", "ID", "REF", "ALT", *ids]]
                for j, row in enumerate(branches[a]):
                    gt.append(row["metadata"] + [expected[j][reverse[i]] for i in range(len(ids))])
                    ht.append(row["metadata"] + row["H"])
                gfile, hfile = stem + "." + label + ".G.tsv", stem + "." + label + ".H.tsv"
                table(output / gfile, gt)
                table(output / hfile, ht)
                manifest.append([label, binary, "ALT"])
                text_manifest.append([label, gfile, hfile, "ALT"])
                info[label] = offsets
            table(output / (stem + ".localbgen.tsv"), manifest)
            table(output / (stem + ".local.tsv"), text_manifest)
            recipes.append({"stem": stem, "bits": bits, "codec": codec, "offsets": info})
    (output / "SIMULATION.json").write_text(json.dumps({"generator": "SPAGxE tests/make_local_bgen_fixtures.py", "seed": 20261025, "samples": len(ids), "variants": len(branches[0]), "ancestries": labels, "recipes": recipes}, indent=2) + "\n", encoding="utf-8")
    two = output / "two_ancestry"
    two.mkdir(exist_ok=True)
    two_ids, two_records = simulated_records(160, 2, seed=20261026)
    bgen_manifest = [["ancestry", "bgen", "counted_allele"]]
    pair_manifest = [["ancestry", "genotype", "hapcount", "counted_allele"]]
    for a, records in enumerate(two_records):
        label = "Ancestry%d" % (a+1)
        expected, _ = write_bgen(two / (label + ".bgen"), two_ids, records, 16, "zlib")
        gt, ht = [["CHROM", "POS", "ID", "REF", "ALT", *two_ids]], [["CHROM", "POS", "ID", "REF", "ALT", *two_ids]]
        for j, record in enumerate(records):
            gt.append(record["metadata"] + expected[j])
            ht.append(record["metadata"] + record["H"])
        table(two / (label + ".G.tsv"), gt)
        table(two / (label + ".H.tsv"), ht)
        bgen_manifest.append([label, label + ".bgen", "ALT"])
        pair_manifest.append([label, label + ".G.tsv", label + ".H.tsv", "ALT"])
    table(two / "cohort.localbgen.tsv", bgen_manifest)
    table(two / "cohort.local.tsv", pair_manifest)
    return recipes


if __name__ == "__main__":
    generate(Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "examples/local/bgen_demo")
    print("Generated SPAGxE synthetic known-H BGEN fixtures")
