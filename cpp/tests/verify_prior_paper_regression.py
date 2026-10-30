"""Audit the current guarded partial policy against a frozen full-CGF binary.
Model/AF/normal-score/Wald/metadata fields are strict acceptance invariants.
SPA/final P, P-derived Z and SPA status differences form a delta study.
"""
from pathlib import Path
import argparse,csv,hashlib,json,math,os,subprocess
from reference_comparison import comparison,depends_on_spa_cgf,valid_changed_p
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',required=True,type=Path);p.add_argument('--baseline',required=True,type=Path);p.add_argument('--out-dir',required=True,type=Path);a=p.parse_args()
 root=Path(__file__).resolve().parents[1];out=a.out_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 current=a.binary.resolve();old=a.baseline.resolve();records=[];fail=[];count=0;maximum=0
 deltas=[];strict_count=0;strict_maximum=0
 def rows(path):
  with path.open(encoding='utf-8')as h:return list(csv.DictReader(h,delimiter='\t'))
 def compare(x,y,tag):
  nonlocal count,maximum,strict_count,strict_maximum
  rx,ry=rows(x),rows(y)
  if len(rx)!=len(ry):fail.append(tag+':row count');return
  for row1,row2 in zip(rx,ry):
   if row1.keys()!=row2.keys():fail.append(tag+':columns');return
   for name in row1:
    v,w=row1[name],row2[name]
    ok,delta=comparison(v,w,1e-8)
    count+=1
    if delta is not None:maximum=max(maximum,delta)
    ident=':'.join(filter(None,[row1.get('ID'),row1.get('ANCESTRY')]))
    if depends_on_spa_cgf(name):
     deltas.append(dict(case=tag,ID=ident,field=name,current=v,full_reference=w,absolute_difference=delta,equal_within_1e_8=ok))
     if not valid_changed_p(name,v,w):fail.append(tag+':'+ident+':invalid or lost P:'+name)
    else:
     strict_count+=1
     if delta is not None:strict_maximum=max(strict_maximum,delta)
     if not ok:fail.append(tag+':'+ident+':'+name+':'+v+' vs '+w)
 cases=[]
 for method in ['base','plus','mix']:
  for mode in ['native','residual']:
   for env in ['E','Eshift']:
    cmd=['--bfile','examples/synthetic/synthetic','--pheno','examples/synthetic/synthetic.pheno','--envir-name',env,'--mac','0']
    cmd+=['--pheno-name','Quantitative,Binary,Time:Event,Ordinal','--covar-name','X1,'+env+',PC1,PC2']if mode=='native'else['--resid-name','R_quant,R_binary,R_survival,R_ordinal']
    if method in['mix']:cmd+=['--method','SPAGxEmix','--pc-cols','PC1,PC2']
    if method in['plus']:cmd+=['--grm-file','examples/synthetic/synthetic.grm.tsv']
    cases.append((method+'_'+mode+'_'+env,cmd))
 for fmt,source in [('text',['--local-manifest','examples/local/cohort.local.tsv']),('bgen',['--local-bgen-manifest','examples/local/bgen_demo/bits16_zlib.localbgen.tsv'])]:
  for mode in ['native','residual']:
   cmd=['--method','SPAGxEmix_local',*source,'--pheno','examples/local/samples.pheno','--envir-name','E','--mac','0','--maf','0','--geno','1']
   cmd+=['--pheno-name','Quantitative,Binary,Time:Event,Ordinal','--covar-name','X,E']if mode=='native'else['--resid-name','R_quant,R_binary,R_survival,R_ordinal']
   cases.append(('local_'+fmt+'_'+mode,cmd))
 for tag,cmd in cases:
  prefixes=[]
  for label,binary in [('current',current),('previous_paper',old)]:
   pref=out/(tag+'_'+label)
   actual_cmd=[('--sp-grm-grab'if x=='--grm-file'and label=='previous_paper'else x)for x in cmd]
   cp=subprocess.run([str(binary),*actual_cmd,'--out',str(pref)],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=120)
   (out/(pref.name+'.log')).write_bytes(cp.stdout)
   if cp.returncode:raise RuntimeError(cp.stdout.decode(errors='replace')[-1000:])
   prefixes.append(pref)
  files=sorted(x for x in out.glob(prefixes[0].name+'.*')if 'SPAGxE'in x.name and not x.name.endswith('.log'))
  for f in files:
   suffix=f.name[len(prefixes[0].name):];compare(f,out/(prefixes[1].name+suffix),tag+suffix)
  records.append(dict(case=tag,tables=len(files)))
 with(out/'full_reference_differences.tsv').open('w',encoding='utf-8',newline='')as h:
  writer=csv.DictWriter(h,fieldnames=['case','ID','field','current','full_reference','absolute_difference','equal_within_1e_8'],delimiter='\t');writer.writeheader();writer.writerows(deltas)
 finite=sorted(r['absolute_difference']for r in deltas if r['absolute_difference']is not None)
 quantile=lambda p:finite[round((len(finite)-1)*p)]if finite else None
 result=dict(passed=not fail,comparison_role='strict invariants plus full-reference delta study',cgf_policy='guarded_partial_normal',binary_sha256=hashlib.sha256(current.read_bytes()).hexdigest(),baseline_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),case_count=len(cases),commands=2*len(cases),field_comparisons=count,max_abs_numeric_difference=maximum,invariant_field_comparisons=strict_count,max_invariant_numeric_difference=strict_maximum,cgf_dependent_fields=len(deltas),cgf_fields_different_above_1e_8=sum(not r['equal_within_1e_8']for r in deltas),cgf_absolute_difference_quantiles={str(p):quantile(p)for p in [.5,.9,.99,1]},failures=fail[:100],failure_count=len(fail),records=records)
 (out/'prior_paper_regression.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
 print(json.dumps({k:result[k]for k in ['passed','case_count','commands','field_comparisons','max_abs_numeric_difference']},indent=2))
 return int(bool(fail))
if __name__=='__main__':raise SystemExit(main())
