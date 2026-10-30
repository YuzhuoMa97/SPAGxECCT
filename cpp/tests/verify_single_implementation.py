"""One guarded partial-normal policy: CLI migration, receipts and alias checks."""
from pathlib import Path
import csv,hashlib,json,os,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];BIN=Path(sys.argv[1]).resolve();OUT=Path(sys.argv[2]).resolve();OUT.mkdir(parents=True,exist_ok=True)
checks=[]
def check(name,ok,detail=''):checks.append(dict(name=name,passed=bool(ok),detail=detail))
def run(name,options,base=None,code=0,fragment=None):
 pref=OUT/name
 if base is None:base=['--bfile','examples/synthetic/synthetic','--pheno','examples/synthetic/synthetic.pheno','--resid-name','R_quant','--envir-name','E','--mac','0','--maf','0']
 cp=subprocess.run([str(BIN),*base,*options,'--out',str(pref)],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
 (OUT/(name+'.log')).write_bytes(cp.stdout)
 text=cp.stdout.decode('utf-8',errors='replace');check(name+':exit',cp.returncode==code,text[-1200:]if cp.returncode!=code else'')
 if fragment:check(name+':message',fragment in text)
 return pref
def receipt(pref,suffix='.analysis.settings.tsv'):
 with Path(str(pref)+suffix).open(encoding='utf-8')as h:return dict(list(csv.reader(h,delimiter='\t'))[1:])
p=run('default',[]);r=receipt(p)
check('single_settings',r['version']==(ROOT/'VERSION').read_text().strip()and r['cgf']=='guarded_partial_normal'and r['grm_adjustment']=='score'and'algorithm_profile'not in r)
for flag,value in [('--algorithm-profile','paper'),('--algorithm-profile','unsupported'),('--spa-mode','full'),('--spa-mode','hybrid'),('--spa-exact',None),('--grm-spa-adjust','score'),('--grm-spa-adjust','tail'),('--outlier-iqr-multiplier','1.5')]:
 opts=[flag]+([]if value is None else[value]);run('removed_'+flag[2:]+'_'+str(value),opts,code=1,fragment='was removed')
grm=['--grm-file','examples/synthetic/synthetic.grm.tsv']
p=run('plus_explicit',['--method','SPAGxE+']+grm);q=run('plus_implicit',grm)
check('plus_alias_equal',Path(str(p)+'.R_quant.SPAGxE').read_bytes()==Path(str(q)+'.R_quant.SPAGxE').read_bytes())
run('plus_requires_grm',['--method','SPAGxE+'],code=1,fragment='requires a sparse GRM')
local=['--method','SPAGxEmix_local','--local-bgen-manifest','examples/local/bgen_demo/bits16_zlib.localbgen.tsv','--pheno','examples/local/bgen_demo/samples.pheno','--resid-name','R','--envir-name','E','--mac','0','--maf','0','--geno','1']
p=run('local_canonical',[],base=local);r=receipt(p,'.local.settings.tsv')
check('local_single_settings',r['version']==(ROOT/'VERSION').read_text().strip()and r['cgf']=='guarded_partial_normal'and'algorithm_profile'not in r)
legacy=local.copy();legacy[1]='SPAGxElocal';q=run('local_old_spelling',[],base=legacy)
check('local_alias_same_algorithm',Path(str(p)+'.R.SPAGxE_local').read_bytes()==Path(str(q)+'.R.SPAGxE_local').read_bytes())
run('local_sources_exclusive',['--local-prefix','unused'],base=local,code=1,fragment='requires one of')
cp=subprocess.run([str(BIN),'--help'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
check('help_one_implementation',b'one' in cp.stdout and b'guarded partial-normal' in cp.stdout and b'1.5' in cp.stdout and b'1e-3' in cp.stdout)
report=dict(passed=all(x['passed']for x in checks),total=len(checks),failure_count=sum(not x['passed']for x in checks),binary_sha256=hashlib.sha256(BIN.read_bytes()).hexdigest(),checks=checks)
(OUT/'single_implementation_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:report[k]for k in ['passed','total','failure_count']},indent=2))
sys.exit(0 if report['passed']else 1)
