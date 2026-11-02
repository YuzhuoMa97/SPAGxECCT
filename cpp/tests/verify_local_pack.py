"""Black-box local count packer checks. Only independent simulated inputs.

Usage: python -B tests/verify_local_pack.py BINARY OUTPUT_DIRECTORY
The CRC checks derive expected count planes from input values independently;
round trips use association only as a reader check, not a calibration proof.
"""
from pathlib import Path
import csv, hashlib, json, math, os, struct, subprocess, sys, zlib

CPP=Path(__file__).resolve().parents[1]
BINARY=Path(sys.argv[1]).resolve()
OUT=Path(sys.argv[2]).resolve()
OUT.mkdir(parents=True,exist_ok=True)
checks=[];commands=0

def check(name,ok,detail=''):
    checks.append(dict(name=name,passed=bool(ok),detail=detail))

def run(name,flags,prefix=None,success=True,fragment=None):
    global commands
    commands+=1
    if prefix is None:prefix=OUT/name
    cp=subprocess.run([str(BINARY),*flags,'--out',str(prefix)],cwd=CPP,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
    (OUT/(name+'.log')).write_bytes(cp.stdout)
    text=cp.stdout.decode('utf-8',errors='replace')
    check(name+':exit',(cp.returncode==0)==success,text[-1200:] if ((cp.returncode==0)!=success) else '')
    if fragment:check(name+':message',fragment in text)
    if not success:
        check(name+':no_new_triplet',not any(Path(str(prefix)+ext).exists() for ext in ('.localbed','.localbim','.localfam')))
        check(name+':temporary_cleanup',not list(prefix.parent.glob(prefix.name+'.pack-tmp-*')))
    return prefix

ids=['0001','0002','0003','0004','0005','0006','sample7']
meta=[('1','101','own1','C','A'),('1','102','own2','G','T'),('2','201','own3','A','C')]
H=[[[0,1,2,1,2,0,1],[1,1,0,2,1,0,2],[2,0,1,1,0,2,1]],
   [[2,1,0,1,0,2,1],[1,1,2,0,1,2,0],[0,2,1,1,2,0,1]]]
G=[[[math.nan,0,1,1,2,0,math.nan],[0,1,0,2,0,0,1],[1,0,1,0,0,2,0]],
   [[1,1,0,0,0,1,1],[1,0,1,0,1,1,0],[0,1,0,1,2,0,1]]]

def write_input(folder,ref=False,fractional=False):
    folder.mkdir(exist_ok=True)
    files=[]
    for a,label in enumerate(['A','B']):
        paths=[]
        for kind in ['G','H']:
            file=folder/(label+'.'+kind+'.tsv');paths.append(file)
            # Different source orders exercise G/H alignment before packing.
            order=list(range(len(ids))) if kind=='G' else list(reversed(range(len(ids))))
            lines=['CHROM\tPOS\tID\tREF\tALT\t'+'\t'.join(ids[i] for i in order)]
            for j,m in enumerate(meta):
                values=[]
                for i in order:
                    value=G[a][j][i] if kind=='G' else H[a][j][i]
                    if kind=='G' and ref and math.isfinite(value):value=H[a][j][i]-value
                    if fractional and a==0 and j==0 and i==2 and kind=='G':value=.25
                    values.append('NA' if math.isnan(value) else format(value,'.17g'))
                lines.append('\t'.join(m)+'\t'+'\t'.join(values))
            file.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        files.append((label,*paths))
    manifest=folder/'manifest.tsv'
    manifest.write_text('ancestry\tgenotype\thapcount\tcounted_allele\n'+''.join(f'{label}\t{g.name}\t{h.name}\t'+('REF' if ref else 'ALT')+'\n'for label,g,h in files),encoding='utf-8')
    return manifest

manifest=write_input(OUT/'source')
ref=write_input(OUT/'ref_source',ref=True)
fractional=write_input(OUT/'fractional_source',fractional=True)
base=['--pack-local','--local-manifest',str(manifest),'--local-block-variants','2']
first=run('default_pack',base)

def inspect(prefix,order,block):
    data=Path(str(prefix)+'.localbed').read_bytes()
    magic,n,m,k,b,footer,direction=struct.unpack_from('<8sQQIIQB',data)
    check(prefix.name+':header',(magic,n,m,k,b,direction)==(b'SQRLOC01',7,3,2,block,0))
    cursor=41;labels=[]
    for a in range(k):
        size=struct.unpack_from('<H',data,cursor)[0];cursor+=2;labels.append(data[cursor:cursor+size].decode());cursor+=size
    check(prefix.name+':labels',labels==['A','B'])
    count=struct.unpack_from('<Q',data,footer)[0]
    offsets=struct.unpack_from('<'+'Q'*count,data,footer+8)
    check(prefix.name+':footer_length',footer+16+8*count==len(data))
    bimcrc,famcrc=struct.unpack_from('<II',data,footer+8+8*count)
    check(prefix.name+':bim_crc',bimcrc==zlib.crc32(Path(str(prefix)+'.localbim').read_bytes()))
    check(prefix.name+':fam_crc',famcrc==zlib.crc32(Path(str(prefix)+'.localfam').read_bytes()))
    for q,offset in enumerate(offsets):
        nr,rawbytes,compressed,crc=struct.unpack_from('<IQQI',data,offset)
        expected=bytearray()
        for j in range(q*block,min(3,(q+1)*block)):
            for a in range(2):
                for kind in ['G','H']:
                    plane=bytearray(2)
                    for target,i in enumerate(order):
                        value=G[a][j][i] if kind=='G' else H[a][j][i]
                        if kind=='G' and H[a][j][i]==0 and math.isnan(value):value=0
                        code=3 if math.isnan(value) else int(value)
                        plane[target//4]|=code<<(2*(target%4))
                    expected.extend(plane)
        check(prefix.name+':raw_dimensions_'+str(q),nr==min(block,3-q*block) and rawbytes==len(expected))
        check(prefix.name+':raw_crc_'+str(q),crc==zlib.crc32(expected))
        nextoff=offsets[q+1] if q+1<count else footer
        check(prefix.name+':block_extent_'+str(q),offset+24+compressed==nextoff)

inspect(first,list(range(7)),2)
# Preserve arbitrary valid six-column sidecar fields and original CRLF bytes.
order=list(reversed(range(7)))
fam=OUT/'original.fam';fam.write_bytes(''.join(f'family{i}  {ids[i]}  0 0 2 {10+i}\r\n' for i in order).encode())
bim=OUT/'original.bim';bim.write_bytes(''.join(f'{m[0]}  {m[2]}  {j+.25} {m[1]} {m[4]} {m[3]}\r\n' for j,m in enumerate(meta)).encode())
reordered=run('sidecar_pack',base+['--local-fam',str(fam),'--local-bim',str(bim)])
check('sidecars:FAM_bytes',Path(str(reordered)+'.localfam').read_bytes()==fam.read_bytes())
check('sidecars:BIM_bytes',Path(str(reordered)+'.localbim').read_bytes()==bim.read_bytes())
inspect(reordered,order,2)
# REF source normalizes once to canonical ALT, preserving the expected planes.
canonical=run('ref_pack',['--pack-local','--local-manifest',str(ref),'--local-block-variants','2'])
inspect(canonical,list(range(7)),2)
run('repack',['--pack-local','--local-prefix',str(first),'--local-block-variants','1','--local-zstd-level','-5'])
inspect(OUT/'repack',list(range(7)),1)
# A known-H point-posterior collection is an eligible lossless hardcall source.
from make_local_bgen_fixtures import write_bgen
point_dir=OUT/'point_bgen';point_dir.mkdir(exist_ok=True)
for a,label in enumerate(['A','B']):
    sample_order=list(range(7)) if a==0 else list(reversed(range(7)))
    records=[]
    for j,m in enumerate(meta):
        probabilities=[]
        for i in sample_order:
            if math.isnan(G[a][j][i]):probabilities.append(None)
            else:
                vector=[0.]*(H[a][j][i]+1);vector[int(G[a][j][i])]=1.;probabilities.append(vector)
        records.append({'metadata':m,'H':[H[a][j][i]for i in sample_order],'P':probabilities})
    write_bgen(point_dir/(label+'.bgen'),[ids[i]for i in sample_order],records,16,'zlib')
point_manifest=point_dir/'manifest.tsv'
point_manifest.write_text('ancestry\tbgen\tcounted_allele\nA\tA.bgen\tALT\nB\tB.bgen\tALT\n',encoding='utf-8')
point=run('point_posterior_pack',['--pack-local','--local-bgen-manifest',str(point_manifest),'--local-block-variants','2'])
inspect(point,list(range(7)),2)

run('fractional_rejected',['--pack-local','--local-manifest',str(fractional)],success=False,fragment='not hard-called')
# An uncertain posterior is rejected even when its expected G may be integer.
run('posterior_rejected',['--pack-local','--local-bgen-manifest','examples/local/bgen_demo/bits16_zlib.localbgen.tsv'],success=False,fragment='not hard-called')
subset=OUT/'subset.fam';subset.write_text('\n'.join(fam.read_text().splitlines()[:-1])+'\n')
run('subset_rejected',base+['--local-fam',str(subset)],success=False,fragment='same complete sample set')
wrong=OUT/'wrong.bim';wrong.write_text('\n'.join(reversed(bim.read_text().splitlines()))+'\n')
run('order_rejected',base+['--local-bim',str(wrong)],success=False,fragment='variant order')
# Existing complete outputs and a lone existing sidecar must survive failures.
original={ext:Path(str(first)+ext).read_bytes()for ext in ['.localbed','.localbim','.localfam']}
cp=subprocess.run([str(BINARY),*base,'--out',str(first)],cwd=CPP,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30);commands+=1
check('existing_complete:exit',cp.returncode!=0)
check('existing_complete:preserved',all(Path(str(first)+ext).read_bytes()==content for ext,content in original.items()))
lone=OUT/'existing_sidecar';Path(str(lone)+'.localbim').write_bytes(b'preserve\x00sentinel')
cp=subprocess.run([str(BINARY),*base,'--out',str(lone)],cwd=CPP,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30);commands+=1
check('existing_sidecar:exit',cp.returncode!=0)
check('existing_sidecar:preserved',Path(str(lone)+'.localbim').read_bytes()==b'preserve\x00sentinel' and not Path(str(lone)+'.localbed').exists()and not Path(str(lone)+'.localfam').exists())
# The same statistical reader sees identical matrices after sample reordering.
ph=OUT/'pheno.tsv';ph.write_text('IID\tE\tR\n'+''.join(f'{id}\t{.2+i*.3}\t{math.sin(i+.7)}\n'for i,id in enumerate(ids)))
common=['--method','SPAGxEmix_local','--pheno',str(ph),'--resid-name','R','--envir-name','E','--mac','0','--maf','0','--geno','1']
text=run('read_text',common+['--local-manifest',str(manifest)])
packed=run('read_packed',common+['--local-prefix',str(reordered)])
for ext in ['.R.SPAGxE_local','.R.SPAGxE_local_CCT']:
    def rows(file):
        with file.open(encoding='utf-8',newline='') as stream:return list(csv.DictReader(stream,delimiter='\t'))
    expected=rows(Path(str(text)+ext));actual=rows(Path(str(packed)+ext))
    key=lambda row:(row['ID'],row.get('ANCESTRY',''))
    left={key(row):row for row in expected};right={key(row):row for row in actual}
    check('roundtrip:'+ext+':keys',left.keys()==right.keys())
    for ident,row in left.items():
        peer=right.get(ident,{})
        for name,value in row.items():
            observed=peer.get(name,'MISSING');ok=value==observed
            if not ok:
                try:
                    a=float(value);b=float(observed)
                    ok=(math.isnan(a)and math.isnan(b))or (math.isfinite(a)and math.isfinite(b)and abs(a-b)<=1e-10*max(1,abs(a),abs(b)))
                except ValueError:pass
            check('roundtrip:'+ext+':'+str(ident)+':'+name,ok,value+' vs '+observed if not ok else '')
report={'binary_sha256':hashlib.sha256(BINARY.read_bytes()).hexdigest(),'commands':commands,'checks':len(checks),'failure_count':sum(not x['passed']for x in checks),'records':checks}
(OUT/'local_pack_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:report[k]for k in ['commands','checks','failure_count']},indent=2))
if report['failure_count']:
    print('\n'.join(x['name']+': '+x['detail']for x in checks if not x['passed']))
sys.exit(bool(report['failure_count']))
