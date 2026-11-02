"""Run SPAGxE and compare independent R normal/SPA/Wald/CCT references.
Run from cpp: python tests/verify_oracle.py [binary] [validation-dir].
All paths written by this harness remain in the chosen validation directory.
"""
import csv,json,math,os,pathlib,subprocess,sys
ROOT=pathlib.Path(__file__).resolve().parents[1]
BIN=pathlib.Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/"build"/"bin"/("SPAGxE.exe" if os.name=="nt" else "SPAGxE")
OUT=pathlib.Path(sys.argv[2]).resolve() if len(sys.argv)>2 else ROOT.parent.parent/"validation"
OUT.mkdir(parents=True,exist_ok=True)
def tab(p):
    with open(p,encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
REF=tab(ROOT/"examples/synthetic/oracle_reference.tsv")
TRAITS={"quantitative":"Quantitative","binary":"Binary","survival":"Time_Event","ordinal":"Ordinal"}
summaries=[];diffrows=[];failures=[];null_checks=[]
RESIDUAL = os.environ.get("SPAGXE_ORACLE_RESIDUAL_MODE") == "1"
if RESIDUAL: TRAITS={"quantitative":"R_quant","binary":"R_binary","survival":"R_survival","ordinal":"R_ordinal"}
def num(x):
    try:return float(x)
    except (TypeError,ValueError):return math.nan
for mode in ("base","mix","plus"):
  for exact in (False,):
    tag=("oracle_residual_" if RESIDUAL else "oracle_")+mode+("_exact" if exact else "_guarded_partial")
    prefix=OUT/tag
    args=[str(BIN),"--bfile","examples/synthetic/synthetic","--pheno","examples/synthetic/synthetic.pheno",
          "--pheno-name","Quantitative,Binary,Time:Event,Ordinal","--covar-name","X1,E,PC1,PC2","--envir-name","E",
          "--regression-model","auto","--threads",os.environ.get("SPAGXE_TEST_THREADS","1"),"--save-resid","--out",str(prefix)]
    if RESIDUAL:
        args.remove("--save-resid")
        args[args.index("--pheno-name")]="--resid-name"
        args[args.index("Quantitative,Binary,Time:Event,Ordinal")]="R_quant,R_binary,R_survival,R_ordinal"
    if mode=="mix":args+=["--method","SPAGxEmix","--pc-cols","PC1,PC2"]
    if mode=="plus":args+=["--grm-file","examples/synthetic/synthetic.grm.tsv"]
    # Keep outputs within the explicit validation directory.
    args[args.index("--out")+1]=os.path.relpath(prefix,ROOT)
    proc=subprocess.run(args,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,encoding="utf-8",errors="replace",timeout=60)
    (OUT/(tag+".console.txt")).write_text(proc.stdout,encoding="utf-8")
    if proc.returncode:failures.append(tag+": process exit "+str(proc.returncode));continue
    if not RESIDUAL:
      native_res={r["IID"]:r for r in tab(str(prefix)+".null.resid")}
      source_res={r["IID"]:r for r in tab(ROOT/"examples/synthetic/synthetic.pheno")}
      for native_name,ref_name in [("Quantitative","R_quant"),("Binary","R_binary"),("Time_Event","R_survival"),("Ordinal","R_ordinal")]:
        err=max(abs(num(native_res[i][native_name])-num(source_res[i][ref_name])) for i in source_res)
        ok=math.isfinite(err) and err<1e-5
        null_checks.append({"run":tag,"trait":native_name,"subjects":len(source_res),"max_abs_error":err,"pass":ok})
        if not ok:failures.append(f"{tag}: native residual {native_name} error {err}")
    for tr,cn in TRAITS.items():
      meth="SPAGxEmix" if mode=="mix" else "SPAGxE"
      raw_result=tab(str(prefix)+"."+cn+"."+meth)
      result={r["ID"]:r for r in raw_result}
      if len(raw_result)!=48 or len(result)!=48:
        failures.append(tag+":"+tr+": expected exactly 48 unique input variant rows")
      mono=result.get("rsOracle48",{})
      if not math.isnan(num(mono.get("LOG10P_GxE"))):
        failures.append(tag+":"+tr+": monomorphic variant unexpectedly tested")
      refs=[r for r in REF if r["MODE"]==mode and r["TRAIT"]==tr]
      maxerr={};compared=0;bad=0
      for r in refs:
        rr=result.get(r["ID"])
        if rr is None:failures.append(tag+":"+tr+":"+r["ID"]+": absent");bad+=1;continue
        fields={"LOG10P_G":"LOG10P_G","Z_G":"Z_G","Z_Norm":"Z_Norm_GxE",
                ("LOG10P_EXACT_SPA" if exact else "LOG10P_SPA"):"LOG10P_SPA_GxE",
                ("LOG10P_EXACT_FINAL" if exact else "LOG10P_FINAL"):"LOG10P_GxE",
                "LOG10P_WALD":"LOG10P_Wald_GxE"}
        if RESIDUAL:
          fields.pop("LOG10P_WALD",None)
          fields.pop("LOG10P_EXACT_FINAL" if exact else "LOG10P_FINAL",None)
        for oracle,actual in fields.items():
          a=num(rr[actual]);b=num(r[oracle])
          # Oracle is minor oriented. Executable reports ALT direction.
          if oracle in ("Z_G","Z_Norm") and num(rr["ALT_FREQ"])>.5:b=-b
          if math.isnan(a) and math.isnan(b):continue
          error=abs(a-b);tol=1e-5 if oracle=="LOG10P_WALD" else 2e-5
          if not math.isfinite(error):error=math.inf
          maxerr[oracle]=max(maxerr.get(oracle,0),error)
          ok=error<=tol
          if not ok:bad+=1;failures.append(f"{tag}:{tr}:{r['ID']}:{oracle}: {a} vs {b}, err={error}")
          diffrows.append({"run":tag,"trait":tr,"ID":r["ID"],"field":oracle,"expected":b,"actual":a,"abs_error":error,"pass":ok})
          compared+=1
        br="B" if num(rr["BRANCH_GxE"])==2 else "A"
        if br!=r["BRANCH"]:bad+=1;failures.append(tag+":"+tr+":"+r["ID"]+": branch")
      summaries.append({"run":tag,"trait":tr,"variants":len(refs),"comparisons":compared,"failed":bad,"max_abs_error":maxerr})
with open(OUT/("oracle_residual_cpp_comparison.tsv" if RESIDUAL else "oracle_cpp_comparison.tsv"),"w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["run","trait","ID","field","expected","actual","abs_error","pass"],delimiter="\t");w.writeheader();w.writerows(diffrows)
summary={"cgf_policy":"guarded_partial_normal","passed":not failures,"reference_rows":len(REF),"native_residual_checks":null_checks,"summaries":summaries,"failure_count":len(failures),"failures":failures}
(OUT/("oracle_residual_cpp_summary.json" if RESIDUAL else "oracle_cpp_summary.json")).write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps({"passed":not failures,"runs":len(summaries),"comparisons":len(diffrows),"failure_count":len(failures),"first_failures":failures[:10]},indent=2))
sys.exit(0 if not failures else 1)
