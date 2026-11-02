#!/usr/bin/env python3
"""Native global/local CCT numerical, alignment and failure-contract checks."""
import argparse,csv,gzip,json,math,pathlib,subprocess,tempfile

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',required=True,type=pathlib.Path)
    p.add_argument('--out-dir',required=True,type=pathlib.Path)
    p.add_argument('--rscript',type=pathlib.Path)
    a=p.parse_args();binary=a.binary.resolve();root=a.out_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    cases=[]
    def check(name,condition,detail=''):
        cases.append(dict(name=name,passed=bool(condition),detail=detail));print(('PASS 'if condition else'FAIL ')+name)
    def write(path,header,rows):
        text='\t'.join(header)+'\n'+''.join('\t'.join(map(str,row))+'\n' for row in rows)
        if path.suffix=='.gz':
            with gzip.open(path,'wt',encoding='utf-8',newline='') as h:h.write(text)
        else:path.write_text(text,encoding='utf-8')
    pairs=[(.01,.2,2),(.9995,.9999,3),(.6,.7,1),(1e-10,.999999,2),(.1,.1,3),(.5,.5,1),(.99999,.99999,2)]
    logs=[(-math.log10(x),-math.log10(y),n)for x,y,n in pairs]
    logs.extend([(1000.,2.,2),(400.,400.,3),(None,2.,2),(2.,None,0)])
    keys=[['1',str(i+100),'v'+str(i),'A','G']for i in range(len(logs))]
    gh=['CHROM','POS','ID','REF','ALT','LOG10P_GxE']
    lh=['CHROM','POS','ID','REF','ALT','LOG10P_ANCESTRY_CCT_GxE','N_ANCESTRIES_GxE']
    gr=[key+['NA'if g is None else format(g,'.17g')]for key,(g,l,n)in zip(keys,logs)]
    lr=[key+['NA'if l is None else format(l,'.17g'),n]for key,(g,l,n)in zip(keys,logs)]
    global_file=root/'global.Y.SPAGxEmix';local_file=root/'local.Y.SPAGxE_local_CCT'
    write(global_file,gh,gr);write(local_file,lh,lr)
    def run(name,extra=(),g=global_file,l=local_file,success=True,expect=''):
        out=root/name
        cmd=[str(binary),'--method','SPAGxEcct','--cct-global',str(g),'--cct-local',str(l),
             '--cct-phenotype','Y','--envir-name','E','--out',str(out),*map(str,extra)]
        r=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=30)
        check(name,(r.returncode==0)==success and (not expect or expect in r.stderr),r.stderr[-1500:]if r.returncode else'')
        return out,r
    def read(out,suffix=''):
        file=pathlib.Path(str(out)+'.Y.SPAGxE_global_local_CCT'+suffix)
        handle=gzip.open(file,'rt',encoding='utf-8')if suffix=='.gz'else file.open(encoding='utf-8')
        with handle:return list(csv.DictReader(handle,delimiter='\t'))
    def expected(x,y,wg,wl):
        if x is None or y is None:return None
        if x>=300 and y>=300:return max(x,y)+math.log10((wg*10**(x-max(x,y))+wl*10**(y-max(x,y)))/(wg+wl))
        if max(x,y)>=300:
            z,w=(x,wg)if x>y else(y,wl)
            return z+math.log10(w/(wg+wl))
        p1,p2=10**-x,10**-y
        if p1==1:p1=.999
        if p2==1:p2=.999
        t=(wg/math.tan(math.pi*p1)+wl/math.tan(math.pi*p2))/(wg+wl)
        probability=math.atan(1/t)/math.pi if t>0 else .5-math.atan(t)/math.pi
        return -math.log10(probability)
    numeric_rows=[]
    for name,extra,weights in [('default',(),(.5,.5)),('weighted',('--cct-weights','1,4'),(1,4)),
                               ('ancestry_equal',('--cct-per-ancestry-equal',),None)]:
        out,r=run(name,extra)
        if r.returncode:continue
        values=read(out);check(name+'-row-count',len(values)==len(logs))
        for i,(row,(x,y,n))in enumerate(zip(values,logs)):
            wg,wl=weights if weights else(1,n)
            want=expected(x,y,wg,wl);got=row['LOG10P_GLOBAL_LOCAL_CCT_GxE']
            if want is None:check(name+'-missing-'+str(i),got=='NA'and row['CCT_STATUS_GxE'].startswith('MISSING_'))
            else:
                got=float(got);check(name+'-numeric-'+str(i),abs(got-want)<=2e-10*(1+abs(want)),str((got,want)))
                if i<len(pairs):numeric_rows.append([x,y,wg,wl,got])
            check(name+'-weights-'+str(i),abs(float(row['WEIGHT_GLOBAL_GxE'])-wg/(wg+wl))<1e-14)
    gg=root/'global.Y.SPAGxEmix.gz';ll=root/'local.Y.SPAGxE_local_CCT.gz';write(gg,gh,gr);write(ll,lh,lr)
    out,r=run('compressed',('--compression','gz'),g=gg,l=ll)
    if r.returncode==0:check('gzip-input-output-identical',read(out,'.gz')==read(root/'default'))
    malformed=[('mismatch',[row[:]for row in lr],'metadata/order mismatch'),('short',lr[:-1],'different numbers'),
               ('duplicate',[lr[0],lr[0],*lr[2:]],'metadata/order mismatch'),
               ('badvalue',[row[:]for row in lr],'invalid -log10(P)'),('badcount',[row[:]for row in lr],'zero contributing')]
    malformed[0][1][0][4]='C';malformed[3][1][0][5]='-2';malformed[4][1][0][6]=0
    for name,rows,expect in malformed:
        path=root/(name+'.Y.SPAGxE_local_CCT');write(path,lh,rows)
        output=root/(name+'.Y.SPAGxE_global_local_CCT');output.write_text('SENTINEL\n',encoding='utf-8')
        run(name,l=path,success=False,expect=expect)
        check(name+'-preserves-existing-output',output.read_text(encoding='utf-8')=='SENTINEL\n')
    wrong=root/'wrong.Other.SPAGxEmix';write(wrong,gh,gr);run('trait-mismatch',g=wrong,success=False,expect='global filename')
    env_missing=root/'env.Y.SPAGxEmix';write(env_missing,gh[:-1]+['LOG10P_GxOther'],gr);run('env-mismatch',g=env_missing,success=False,expect='missing column')
    run('mutually-exclusive-weights',('--cct-weights','1,2','--cct-per-ancestry-equal'),success=False,expect='mutually exclusive')
    run('unsupported-threads',('--threads','2'),success=False,expect='does not accept analysis option')
    run('negative-weight',('--cct-weights','-1,2'),success=False,expect='finite positive')
    run('equal-explicit-weights',('--cct-weights','1,1'))
    if a.rscript:
        oracle_in=root/'r_input.tsv';oracle_out=root/'r_expected.tsv'
        write(oracle_in,['GLOBAL','LOCAL','WG','WL','CPP'],numeric_rows)
        # ASCII relative arguments avoid legacy Windows Rscript argument-codepage
        # conversion of the Chinese workspace path; all files stay in out-dir.
        (root/'reference_cct.R').write_bytes(pathlib.Path(__file__).with_name('reference_cct.R').read_bytes())
        r=subprocess.run([str(a.rscript.resolve()),'reference_cct.R','r_input.tsv','r_expected.tsv'],
                          cwd=root,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=60)
        check('independent-R-oracle-runs',r.returncode==0,r.stderr[-1000:])
        if r.returncode==0:
            with oracle_out.open()as h:ref=list(csv.DictReader(h,delimiter='\t'))
            for i,(cpp,rr)in enumerate(zip(numeric_rows,ref)):
                want=float(rr['EXPECTED']);check('independent-R-'+str(i),abs(cpp[-1]-want)<=2e-10*(1+abs(want)),str((cpp[-1],want)))
    summary=dict(total=len(cases),passed=sum(c['passed']for c in cases),cases=cases,binary=str(binary))
    (root/'cct_validation.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(f"{summary['passed']}/{summary['total']} CCT checks passed")
    return 0 if summary['passed']==summary['total']else 1
if __name__=='__main__':raise SystemExit(main())
