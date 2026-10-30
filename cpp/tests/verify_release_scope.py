#!/usr/bin/env python3
"""Check the four public analysis methods and their relatedness boundaries.

Run: python tests/verify_release_scope.py BINARY OUTPUT_DIRECTORY
Uses only the Python standard library and the bundled synthetic examples.
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys


def main():
    binary = Path(sys.argv[1]).resolve()
    output = Path(sys.argv[2]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    cpp = Path(__file__).resolve().parents[1]
    example = cpp / 'examples' / 'synthetic'
    local_example = cpp / 'examples' / 'local' / 'bgen_demo'
    failures, records = [], []
    commands = checks = 0

    def check(condition, label):
        nonlocal checks
        checks += 1
        if not condition:
            failures.append(label)

    def run(tag, flags, expected_message=None):
        nonlocal commands
        commands += 1
        prefix = output / tag
        command = [str(binary), *flags, '--out', str(prefix)]
        env = dict(os.environ)
        for name in ('TEMP', 'TMP', 'TMPDIR'):
            env[name] = str(output)
        result = subprocess.run(command, cwd=cpp, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output / (tag + '.console.log')).write_bytes(result.stdout)
        console = result.stdout.decode('utf-8', errors='replace')
        if expected_message is None:
            check(result.returncode == 0, tag + ': successful analysis')
        else:
            check(result.returncode != 0, tag + ': rejected')
            check(expected_message.lower() in console.lower(), tag + ': explicit reason')
            check(not list(output.glob(tag + '.*.SPAGxE*')), tag + ': no association table')
        records.append({'case': tag, 'exit_code': result.returncode,
                        'expected': 'rejected' if expected_message else 'success'})
        return prefix, console

    def table(prefix, suffix, value_column='LOG10P_GxE'):
        filename = Path(str(prefix) + suffix)
        check(filename.is_file(), prefix.name + ': output table exists')
        if not filename.is_file():
            return []
        with filename.open(newline='', encoding='utf-8') as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        check(bool(rows), prefix.name + ': nonempty output')
        check(bool(rows) and value_column in rows[0], prefix.name + ': interaction column')
        finite = [row for row in rows if row.get(value_column) not in (None, 'NA', 'NaN', 'nan')]
        check(any(math.isfinite(float(row[value_column])) for row in finite),
              prefix.name + ': at least one finite interaction result')
        return rows

    # Encode the same sparse GRM in the two accepted input layouts.
    with (example / 'synthetic.pheno').open(newline='') as stream:
        ids = [row['IID'] for row in csv.DictReader(stream, delimiter='\t')]
    indices = {iid: i for i, iid in enumerate(ids)}
    plink_grm = output / 'scope.grm.sp'
    (output / 'scope.grm.id').write_text(''.join(iid + '\t' + iid + '\n' for iid in ids))
    with (example / 'synthetic.grm.tsv').open(newline='') as source, plink_grm.open('w') as target:
        for row in csv.DictReader(source, delimiter='\t'):
            target.write('%d\t%d\t%s\n' %
                         (indices[row['ID1']], indices[row['ID2']], row['Value']))
    triplet_flags = ['--grm-file', str(example / 'synthetic.grm.tsv')]
    plink_flags = ['--sp-grm-plink2', str(plink_grm)]
    common = ['--bfile', str(example / 'synthetic'), '--pheno', str(example / 'synthetic.pheno'),
              '--resid-name', 'R_quant', '--envir-name', 'E', '--maf', '0', '--mac', '0', '--geno', '1']
    local = ['--method', 'SPAGxEmix_local', '--local-bgen-manifest',
             str(local_example / 'bits16_zlib.localbgen.tsv'),
             '--pheno', str(local_example / 'samples.pheno'), '--resid-name', 'R',
             '--envir-name', 'E', '--maf', '0', '--mac', '0', '--geno', '1']

    _, help_text = run('help', ['--help'])
    for method in ('SPAGxE', 'SPAGxE+', 'SPAGxEmix', 'SPAGxEmix_local'):
        check(method in help_text, 'help documents ' + method)
    prefix, _ = run('base', ['--method', 'SPAGxE', *common])
    table(prefix, '.R_quant.SPAGxE')
    prefix, _ = run('related_triplet', ['--method', 'SPAGxE+', *common, *triplet_flags])
    related_triplet = table(prefix, '.R_quant.SPAGxE')
    prefix, _ = run('related_plink', ['--method', 'SPAGxE+', *common, *plink_flags])
    related_plink = table(prefix, '.R_quant.SPAGxE')
    check(len(related_triplet) == len(related_plink), 'GRM layouts preserve row count')
    for left, right in zip(related_triplet, related_plink):
        check(left == right, 'GRM layouts identical result: ' + left.get('ID', '?'))
    unrelated = ['--method', 'SPAGxEmix', *common, '--pc-cols', 'PC1,PC2']
    prefix, _ = run('unrelated', unrelated)
    table(prefix, '.R_quant.SPAGxEmix')
    prefix, _ = run('local', local)
    table(prefix, '.R.SPAGxE_local')
    for layout, grm_flags in (('grab', triplet_flags), ('plink', plink_flags)):
        run('unrelated_reject_' + layout + '_grm', [*unrelated, *grm_flags],
            'supports unrelated samples; GRM input is not supported')
        run('local_reject_' + layout + '_grm', [*local, *grm_flags],
            'GRM input is not supported')
    run('related_requires_grm', ['--method', 'SPAGxE+', *common],
        'SPAGxE+ requires a sparse GRM')

    summary = {'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
               'command_count': commands, 'check_count': checks,
               'failure_count': len(failures), 'failures': failures, 'cases': records}
    (output / 'release_scope_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: summary[k] for k in ('command_count', 'check_count', 'failure_count')}))
    return bool(failures)


if __name__ == '__main__':
    sys.exit(main())
