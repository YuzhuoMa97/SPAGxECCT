"""Repeat native four-trait analyses across methods, CGFs and thread counts.
Run from cpp: python tests/stress_native_models.py [binary] [output-dir] [repetitions]
Default: 48 processes (3 methods x 2 CGFs x 2 thread counts x 4 repetitions).
"""
import collections,hashlib,json,os,pathlib,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
BIN=pathlib.Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/"build/bin"/("SPAGxE.exe" if os.name=="nt" else "SPAGxE")
OUT=pathlib.Path(sys.argv[2]).resolve() if len(sys.argv)>2 else ROOT.parent.parent/"validation/native_stability"
REPS=int(sys.argv[3]) if len(sys.argv)>3 else 4
OUT.mkdir(parents=True,exist_ok=True)
rows=[];baseline={};failures=[]
for rep in range(REPS):
 for mode in ("base","mix","plus"):
  for exact in (True,):
   for threads in (1,4):
    tag=f"{mode}_{'exact' if exact else 'fast'}_t{threads}_r{rep}"
    prefix=OUT/tag
    cmd=[str(BIN),"--bfile","examples/synthetic/synthetic","--pheno","examples/synthetic/synthetic.pheno",
         "--pheno-name","Quantitative,Binary,Time:Event,Ordinal","--covar-name","X1,E,PC1,PC2","--envir-name","E",
         "--threads",str(threads),"--out",os.path.relpath(prefix,ROOT)]
    if mode=="mix":cmd+=["--method","SPAGxEmix","--pc-cols","PC1,PC2"]
    if mode=="plus":cmd+=["--grm-file","examples/synthetic/synthetic.grm.tsv"]
    started=time.monotonic()
    try:
     p=subprocess.run(cmd,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,encoding="utf-8",errors="replace",timeout=60)
     code=p.returncode
     if code:
      failures.append(f"{tag}: process exit {code}")
      (OUT/(tag+".console.txt")).write_text(p.stdout,encoding="utf-8")
    except subprocess.TimeoutExpired as err:
     code="TIMEOUT";failures.append(tag+": timeout")
    elapsed=time.monotonic()-started
    hashes={}
    if code==0:
     meth="SPAGxEmix" if mode=="mix" else "SPAGxE"
     for trait in ("Quantitative","Binary","Time_Event","Ordinal"):
      path=pathlib.Path(str(prefix)+"."+trait+"."+meth)
      sha=hashlib.sha256(path.read_bytes()).hexdigest();hashes[trait]=sha
      key=(mode,exact,trait)
      if key in baseline and baseline[key]!=sha:failures.append(f"{tag}: nondeterministic {trait}")
      baseline.setdefault(key,sha)
    rows.append({"run":tag,"exit_code":code,"elapsed_seconds":elapsed,"sha256":hashes})
    summary={"passed":not failures,"runs":len(rows),"failures":failures,"details":rows}
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps({"passed":not failures,"runs":len(rows),"failures":failures},indent=2))
sys.exit(1 if failures else 0)
