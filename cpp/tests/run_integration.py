"""Reproducible integration checks. Run from cpp; standard-library Python only."""
from pathlib import Path
import argparse,csv,subprocess,json,time,gzip,hashlib,math,os
p=argparse.ArgumentParser();p.add_argument('--binary',default='build/bin/SPAGxE.exe' if os.name=='nt' else 'build/bin/SPAGxE');p.add_argument('--out-dir',default='examples_output/integration');a=p.parse_args()
bin=Path(a.binary).resolve();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True);records=[]
base=['--pheno','examples/synthetic/synthetic.pheno','--pheno-name','Quantitative,Binary,Time:Event,Ordinal','--covar-name','X1,E,PC1,PC2','--envir-name','E','--mac','1','--maf','0','--geno','0.1','--chunk-ksnp','0.5']

def run(label,opts,common=base):
 cmd=[str(bin)]+common+opts+['--out',str(out/label)];t=time.perf_counter()
 cp=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=120)
 (out/(label+'.log')).write_bytes(cp.stdout)
 rec={'case':label,'seconds':time.perf_counter()-t,'exit':cp.returncode};records.append(rec)
 if cp.returncode: raise RuntimeError(label+' failed: '+cp.stdout.decode(errors='replace')[-1500:])
 files=sorted(out.glob(label+'.*'));return [f for f in files if f.suffix not in ('.log','.resid') and ('SPAGxE' in f.name)]

def tables(files):
 result={}
 for f in files:
  name=f.name.split('.',1)[1].replace('.gz','').replace('.zst','')
  raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
  result[name]=raw
 return result

def same(label,got,expected):
 assert got.keys()==expected.keys(),(label,got.keys(),expected.keys())
 for k in got: assert got[k]==expected[k],label+' differs at '+k
 records.append({'case':label,'check':'byte-identical','tables':len(got)})

bed=tables(run('bed',['--bfile','examples/synthetic/synthetic','--save-resid']))
# BIM allele fields must be independent of the producing system's line endings.
# Compare all four trait outputs with LF, CRLF and a final line without a newline.
source_prefix=Path('examples/synthetic/synthetic')
bim_rows=source_prefix.with_suffix('.bim').read_bytes().splitlines()
assert bim_rows and all(b'\r' not in row and b'\n' not in row for row in bim_rows)
for label,separator,terminate in [('lf',b'\n',True),('crlf',b'\r\n',True),('crlf_no_final_newline',b'\r\n',False)]:
 prefix=out/('bim_'+label)
 prefix.with_suffix('.bed').write_bytes(source_prefix.with_suffix('.bed').read_bytes())
 prefix.with_suffix('.fam').write_bytes(source_prefix.with_suffix('.fam').read_bytes())
 prefix.with_suffix('.bim').write_bytes(separator.join(bim_rows)+(separator if terminate else b''))
 got=tables(run('bim_'+label,['--bfile',str(prefix)]))
 same('BIM_line_endings_'+label,got,bed)
 assert all(b'\r' not in raw for raw in got.values()), 'BIM line terminator leaked into association output'
for label,opts in [('vcf',['--vcf','examples/synthetic/synthetic.vcf']),('bgen',['--bgen','examples/synthetic/synthetic.bgen','ref-first']),('bgen_raw',['--bgen','examples/synthetic/synthetic.raw.bgen','ref-first']),('bgen16',['--bgen','examples/synthetic/synthetic.16bit.bgen','ref-first']),('bgen_raw16',['--bgen','examples/synthetic/synthetic.raw.16bit.bgen','ref-first'])]:
 same('cross_format_'+label,tables(run(label,opts)),bed)
for label,opts in [('bcf',['--bcf','examples/synthetic/synthetic.bcf']),('bgen_noids',['--bgen','examples/synthetic/synthetic.noids.bgen','ref-first']),('bgen_reflast',['--bgen','examples/synthetic/synthetic.reflast.bgen','ref-last']),('bgen_phased',['--bgen','examples/synthetic/synthetic.phased.bgen','ref-first']),('bgen_v11',['--bgen','examples/synthetic/synthetic.v11.bgen','ref-first']),('bgen_v11raw',['--bgen','examples/synthetic/synthetic.v11raw.bgen','ref-first'])]:
 same('cross_format_'+label,tables(run(label,opts)),bed)
same('gzip_output',tables(run('gzip',['--bfile','examples/synthetic/synthetic','--compression','gz'])),bed)
# Repeat a panel to force several independent chunks and real worker parallelism.
src=Path('examples/synthetic/synthetic'); tiled=out/'tiled'; n=800; per=(n+3)//4; raw=src.with_suffix('.bed').read_bytes()[3:]; bim=src.with_suffix('.bim').read_text().splitlines()
tiled.with_suffix('.fam').write_bytes(src.with_suffix('.fam').read_bytes());tiled.with_suffix('.bed').write_bytes(bytes((108,27,1))+raw*32)
tiled.with_suffix('.bim').write_text('\n'.join('\t'.join([r.split()[0],r.split()[1]+'_tile'+str(t)]+r.split()[2:]) for t in range(32) for r in bim)+'\n')
t1=tables(run('thread1',['--bfile',str(tiled),'--threads','1']));t4=tables(run('thread4',['--bfile',str(tiled),'--threads','4']));same('thread_reproducibility',t4,t1)
# Multiple environments, residual path, exact CGF, GRM and admixed-AF paths.
run('multi_env',['--pheno','examples/synthetic/synthetic.pheno','--resid-name','R_quant,R_binary,R_survival','--bfile','examples/synthetic/synthetic','--envir-name','E,Eshift','--mac','1'],[])
run('mix',['--method','SPAGxEmix','--bfile','examples/synthetic/synthetic','--pc-cols','PC1,PC2'])
run('plus',['--bfile','examples/synthetic/synthetic','--grm-file','examples/synthetic/synthetic.grm.tsv'])

run('residual',['--pheno','examples/synthetic/synthetic.pheno','--resid-name','R_quant,R_binary,R_survival,R_offset','--envir-name','E','--bfile','examples/synthetic/synthetic','--mac','1'],[])
real=['--pheno','examples/1kg.pheno','--pheno-name','Quantitative,Binary,Time:Event,Ordinal','--covar-name','MALE,PC1,PC2,PC3,PC4','--envir-name','MALE','--pfile','examples/1kg','--threads','4','--chunk-ksnp','0.5']
run('1kg_base',[],real);run('1kg_mix',['--method','SPAGxEmix'],real);run('1kg_plus',['--sp-grm-plink2','examples/1kg.grm.sp'],real)
# Independent phenotype masks: multi-trait output must match separate fits.
with Path('examples/synthetic/synthetic.pheno').open() as f:
 panel=list(csv.DictReader(f,delimiter='\t')); fields=list(panel[0])
for j,r in enumerate(panel):
 if j%23==0:r['Quantitative']='NA'
 if j%29==0:r['Binary']='NA'
 if j%31==0:r['Time']='NA'
missing=out/'missing.pheno'
with missing.open('w',newline='') as f:
 w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(panel)
changed=base[:];changed[1]=str(missing)
multi=tables(run('missing_multi',['--bfile','examples/synthetic/synthetic'],changed))
for tr in ('Quantitative','Binary','Time:Event','Ordinal'):
 one=changed[:];one[3]=tr
 single=tables(run('missing_'+tr.replace(':','_'),['--bfile','examples/synthetic/synthetic'],one))
 for key,value in single.items():assert multi[key]==value,('phenotype mask',key)
records.append({'case':'separate_vs_multi_missing_masks','check':'byte-identical','tables':4})
# Matching subject/variant filters across the required PLINK and BGEN readers.
(out/'keep.txt').write_text('IID\n'+'\n'.join(r['IID'] for r in panel[:400])+'\n')
(out/'remove.txt').write_text('IID\n'+'\n'.join(r['IID'] for r in panel[:5])+'\n')
(out/'extract.txt').write_text('\n'.join('rsOracle%02d'%i for i in range(1,21))+'\n')
(out/'exclude.txt').write_text('rsOracle02\nrsOracle04\n')
filters=['--keep',str(out/'keep.txt'),'--remove',str(out/'remove.txt'),'--extract',str(out/'extract.txt'),'--exclude',str(out/'exclude.txt'),'--chr','1']
fbed=tables(run('filter_bed',['--bfile','examples/synthetic/synthetic']+filters))
fbgen=tables(run('filter_bgen',['--bgen','examples/synthetic/synthetic.bgen','ref-first']+filters))
same('filters_bed_vs_bgen',fbgen,fbed)
for table in fbed.values():assert len(table.splitlines())==19, 'expected 18 filtered variants'
# All numerical output must be finite or explicitly NA, with statuses recorded.
rows=0
for f in out.glob('*SPAGxE*'):
 if f.suffix in ('.gz','.zst'): continue
 with f.open() as h:
  reader=csv.DictReader(h,delimiter='\t')
  for row in reader:
   rows+=1
   assert all(v not in ('inf','-inf','Inf','-Inf','nan') for v in row.values()),(str(f),row['ID'])
records.append({'case':'finite_output','rows_checked':rows})
(out/'integration_summary.json').write_text(json.dumps({'passed':True,'records':records},indent=2))
print(json.dumps({'passed':True,'checks':len(records),'rows_checked':rows,'summary':str(out/'integration_summary.json')}))

