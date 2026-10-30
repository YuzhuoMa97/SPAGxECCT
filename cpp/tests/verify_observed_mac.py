#!/usr/bin/env python3
"""Observed-dosage MAC/QC across VCF DS, BGEN GP and phenotype masks.

Run: python tests/verify_observed_mac.py BINARY OUTPUT_DIRECTORY
Uses only the standard library; generates independent tiny VCF/BGEN inputs.
It checks the actual SPAGxE product path and AF coefficient utility. It does
not claim to execute the unused generic fused-GEMM MethodBase adapter.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import zlib


def string16(value):
    raw = value.encode()
    return struct.pack('<H', len(raw)) + raw


def string32(value):
    raw = value.encode()
    return struct.pack('<I', len(raw)) + raw


def fixtures(folder):
    ids = ['mac%03d' % i for i in range(12)]
    # Quantization-exact GP: p(0),p(1),p(2) are multiples of 1/255;
    # DS uses the decoded GP expectation, not the pre-encoding target.
    source = {
        'low_missing': [.2] * 8 + [None] * 4,
        'edge_missing': [.2, .4, .6, .2, .4, .2, .2, .2] + [None] * 4,
        'high_missing': [.8, .8, .8, .4, .4, .4, .2, .2] + [None] * 4,
        'no_missing': [.2] * 12,
        'opposite_missing': [None] * 4 + [.2, .4, .6, .2, .4, .2, .2, .2],
    }
    probabilities = {key: [(0, 0) if x is None else
                          (255 - round(255*x), round(255*x)) for x in values]
                     for key, values in source.items()}
    dosage = {key: [None if x is None else probabilities[key][i][1]/255.
                    for i, x in enumerate(values)] for key, values in source.items()}
    vcf = folder/'observed.vcf'
    with vcf.open('w', newline='\n') as out:
        out.write('##fileformat=VCFv4.2\n##contig=<ID=1>\n')
        out.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ids) + '\n')
        for j, (key, values) in enumerate(dosage.items(), 1):
            out.write('\t'.join(['1', str(j), key, 'A', 'G', '.', 'PASS', '.', 'DS'] +
                                ['.' if x is None else format(x, '.17g') for x in values]) + '\n')
    sample_data = struct.pack('<I', len(ids)) + b''.join(string16(i) for i in ids)
    sample_data = struct.pack('<I', 4+len(sample_data)) + sample_data
    bgen = folder/'observed.bgen'
    with bgen.open('wb') as out:
        out.write(struct.pack('<IIII4sI', 20+len(sample_data), 20, len(dosage), len(ids),
                              b'bgen', (1 << 31) | 8 | 1))
        out.write(sample_data)
        for j, (key, values) in enumerate(dosage.items(), 1):
            out.write(string16(key) + string16(key) + string16('1') +
                      struct.pack('<IH', j, 2) + string32('A') + string32('G'))
            raw = struct.pack('<IHBB', len(ids), 2, 2, 2)
            raw += bytes(130 if x is None else 2 for x in values) + bytes((0, 8))
            raw += bytes(x for pair in probabilities[key] for x in pair)
            compressed = zlib.compress(raw)
            out.write(struct.pack('<II', len(compressed)+4, len(raw)) + compressed)
    masks = {'Rall': list(range(12)), 'Rhead': list(range(8)), 'Rtail': list(range(4, 12))}
    pheno = folder/'phenotypes.tsv'
    with pheno.open('w', newline='\n') as out:
        out.write('IID\tE\tPC1\tRall\tRhead\tRtail\n')
        for i, iid in enumerate(ids):
            e = math.sin(.53*i) + .13*i
            r = math.cos(.37*i) + .23*math.sin(.91*i)
            out.write('\t'.join([iid, format(e, '.17g'), format((i-5.5)/12, '.17g')] +
                                [format(r, '.17g') if i in masks[t] else 'NA' for t in masks]) + '\n')
    for trait, mask in masks.items():
        (folder/(trait+'.keep')).write_text(''.join(ids[i]+'\n' for i in mask))
    return pheno, dosage, masks, {'vcf': ['--vcf', str(vcf)],
                                'bgen': ['--bgen', str(bgen), 'ref-first']}


def main():
    binary = Path(sys.argv[1]).resolve()
    output = Path(sys.argv[2]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    pheno, dosage, masks, formats = fixtures(output)
    comparisons = 0
    commands = 0
    failures = []
    records = []

    def check(condition, label):
        nonlocal comparisons
        comparisons += 1
        if not condition:
            failures.append(label)

    def run(tag, flags):
        nonlocal commands
        prefix = output/tag
        command = [str(binary), *flags, '--out', str(prefix)]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output/(tag+'.console.log')).write_bytes(result.stdout)
        commands += 1
        if result.returncode:
            raise RuntimeError('%s exited %d; see console log' % (tag, result.returncode))
        return prefix

    def table(path):
        with path.open() as stream:
            return {row['ID']: row for row in csv.DictReader(stream, delimiter='\t')}

    def same_row(left, right):
        if left.keys() != right.keys():
            return False
        for key in left:
            try:
                if abs(float(left[key])-float(right[key])) >= 2e-6:
                    return False
            except ValueError:
                if left[key] != right[key]:
                    return False
        return True

    qc = ['--maf', '0', '--geno', '1', '--hwe', '0']
    shared = ['--pheno', str(pheno), '--envir-name', 'E', '--spagxe-marginal-cutoff', '0',
              '--spa-z-threshold', '100', '--chunk-ksnp', '.5', *qc]
    outputs = {}
    for fmt, genotype in formats.items():
        for cutoff in (1.5, 2.5, 3.5):
            tag = '%s_mac%g' % (fmt, cutoff)
            multi = run(tag+'_multi', [*genotype, *shared, '--resid-name', ','.join(masks),
                                      '--mac', str(cutoff), '--threads', '2'])
            for trait, mask in masks.items():
                single = run(tag+'_'+trait, [*genotype, *shared, '--resid-name', trait,
                                            '--mac', str(cutoff), '--threads', '1'])
                a = table(Path(str(multi)+'.'+trait+'.SPAGxE'))
                b = table(Path(str(single)+'.'+trait+'.SPAGxE'))
                check(set(a) == set(dosage) == set(b), tag+'/'+trait+'/variant identity')
                # AF export on the same subject mask has to reject the same
                # variants; -1 means QC-filtered, not a pooled AF fit.
                af = run(tag+'_'+trait+'_af', [*genotype, '--cal-af-coef', '--covar', str(pheno),
                         '--pc-cols', 'PC1', '--keep', str(output/(trait+'.keep')), '--mac', str(cutoff), *qc])
                with Path(str(af)+'.afc').open() as stream:
                    lines = (line for line in stream if not line.startswith('#') or line.startswith('#STATUS\t'))
                    afrows = list(csv.DictReader(lines, delimiter='\t'))
                check(len(afrows) == len(dosage), tag+'/'+trait+'/AF record count')
                for (key, values), afrow in zip(dosage.items(), afrows):
                    known = [values[i] for i in mask if values[i] is not None]
                    total = math.fsum(known)
                    mac = min(total, 2*len(known)-total)
                    frequency = total/(2*len(known))
                    missing = 1-len(known)/len(mask)
                    passes = mac >= cutoff
                    row = a[key]
                    check(same_row(row, b[key]), tag+'/'+trait+'/'+key+'/single vs multiple masks')
                    check(abs(float(row['MAC'])-mac) < 2e-6, tag+'/'+trait+'/'+key+'/observed MAC')
                    check(abs(float(row['ALT_FREQ'])-frequency) < 2e-7, tag+'/'+trait+'/'+key+'/observed AF')
                    check(abs(float(row['MISS_RATE'])-missing) < 2e-7, tag+'/'+trait+'/'+key+'/missing rate')
                    tested = row['LOG10P_GxE'] != 'NA'
                    check(tested == passes, tag+'/'+trait+'/'+key+'/actual test availability')
                    if passes:
                        check(math.isfinite(float(row['LOG10P_GxE'])) and row['SPA_STATUS_GxE'] == '1',
                              tag+'/'+trait+'/'+key+'/finite normal result/status')
                    else:
                        check(row['LOG10P_GxE'] == 'NA', tag+'/'+trait+'/'+key+'/MAC filtered result NA')
                    check((int(afrow['#STATUS']) != -1) == passes, tag+'/'+trait+'/'+key+'/AF QC -1 agreement')
                    records.append({'format': fmt, 'mac_cutoff': cutoff, 'trait': trait, 'variant': key,
                                    'observed_subjects': len(known), 'observed_mac': mac,
                                    'mean_imputed_mac_would_be': mac*len(mask)/len(known),
                                    'expected_test_available': passes, 'actual_test_available': tested,
                                    'af_status': int(afrow['#STATUS'])})
                outputs[(fmt, cutoff, trait)] = a
    for cutoff in (1.5, 2.5, 3.5):
        for trait in masks:
            a = outputs[('vcf', cutoff, trait)]
            b = outputs[('bgen', cutoff, trait)]
            for key in dosage:
                for name in a[key]:
                    x, y = a[key][name], b[key][name]
                    try:
                        check(abs(float(x)-float(y)) < 2e-6, 'VCF DS vs BGEN GP/'+trait+'/'+key+'/'+name)
                    except ValueError:
                        check(x == y, 'VCF DS vs BGEN GP/'+trait+'/'+key+'/'+name)
    summary = {'passed': not failures, 'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
               'commands': commands, 'checks': comparisons, 'failures': failures, 'records': records,
               'scope': 'Actual product genotype/QC path and AF export; generic fused MethodBase path is not exercised.'}
    (output/'observed_mac_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k: summary[k] for k in ('passed', 'commands', 'checks', 'failures')}))
    raise SystemExit(bool(failures))


if __name__ == '__main__':
    main()
