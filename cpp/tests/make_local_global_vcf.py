"""Rebuild the matched global VCF dosage example by summing local G over ancestries.
This deterministic tutorial conversion uses only Python's standard library.
"""
import csv,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'examples/local'
def table(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
manifest=table(ROOT/'cohort.local.tsv');gs=[table(ROOT/r['genotype']) for r in manifest];hs=[table(ROOT/r['hapcount']) for r in manifest]
ids=list(gs[0][0])[5:];keys=['CHROM','POS','ID','REF','ALT']
def num(s):
 try:return float(s)
 except ValueError:return math.nan
with (ROOT/'global.vcf').open('w',newline='\n',encoding='utf-8') as f:
 f.write('##fileformat=VCFv4.2\n##FORMAT=<ID=DS,Number=1,Type=Float,Description="Total ALT dosage summed from ancestry-specific G; missing if any G/H is unknown">\n')
 for chrom in dict.fromkeys(r['CHROM'] for r in gs[0]):f.write('##contig=<ID='+chrom+'>\n')
 f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(ids)+'\n')
 for j in range(len(gs[0])):
  meta=[gs[0][j][k] for k in keys];ds=[]
  for iid in ids:
   g=[num(x[j][iid]) for x in gs];h=[num(x[j][iid]) for x in hs]
   for a in range(len(g)):
    if h[a]==0 and math.isnan(g[a]):g[a]=0
    if manifest[a]['counted_allele']=='REF':g[a]=h[a]-g[a]
   ds.append(format(sum(g),'.17g') if all(math.isfinite(x) for x in g+h) else '.')
  f.write('\t'.join(meta+['.','PASS','.','DS']+ds)+'\n')
print('Wrote matched examples/local/global.vcf')
