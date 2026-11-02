"""Compare current outputs with preserved R using strict invariants and SPA deltas.
Run from cpp after reference_original_R.R:
python tests/verify_original_R.py BINARY REFERENCE_DIR OUTPUT_DIR
R is a validation dependency only. All tested analyses run in C++.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import os
import subprocess
import sys

CPP=Path(__file__).resolve().parents[1]
BIN=Path(sys.argv[1]).resolve()
REFERENCE=Path(sys.argv[2]).resolve()
OUT=Path(sys.argv[3]).resolve()
OUT.mkdir(parents=True,exist_ok=True)

def tab(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as stream:return list(csv.DictReader(stream,delimiter="\t"))

def number(value):
    try:return float(value)
    except (ValueError,TypeError):return math.nan

def relative(path):return os.path.relpath(path,CPP)

def cct(ls):
    ps=[.999 if x==0 else 10**(-x) for x in ls]
    t=sum(1/(math.pi*p) if p<1e-16 else 1/math.tan(math.pi*p) for p in ps)/len(ps)
    return math.log10(math.pi)+math.log10(t) if t>1e15 else -math.log10(math.atan(1/t)/math.pi if t>0 else .5-math.atan(t)/math.pi)

R=tab(REFERENCE/"original_R_components.tsv")
F=tab(REFERENCE/"formula_profiles_components.tsv")
key=lambda r:tuple(r[x] for x in ("TRAIT","ENV","MODE","ID"))
source={(key(r),r["ROOT"]):r for r in R if r["SCENARIO"]=="native_residual"}
formula={(key(r),r["PROFILE"]):r for r in F if r["SCENARIO"]=="native_residual"}
trait_cols={"quantitative":("Quantitative","R_quant"),"binary":("Binary","R_binary"),"survival":("Time_Event","R_survival"),"ordinal":("Ordinal","R_ordinal")}
checks=[];failures=[];runs=[];legacy_default_deltas=[]

def check(run,detail,wanted,actual,tolerance=2e-5,required=True):
    if math.isnan(wanted) and math.isnan(actual):error=0.0;passed=True
    elif wanted==actual:error=0.0;passed=True
    else:error=abs(actual-wanted);passed=math.isfinite(error) and error<=tolerance
    checks.append({"run":run,"detail":detail,"expected":wanted,"actual":actual,"abs_error":error,"pass":passed,"required":required})
    if required and not passed:failures.append(f"{run}:{detail}: {actual} vs {wanted}, error {error}")

fixture=CPP/"examples/synthetic"
for input_mode in ("native","residual"):
  for en in ("E","Eshift"):
    for mode in ("base","plus","mix"):
      for cgftype in ("guarded_partial_normal",):
        tag="_".join((input_mode,en,mode,cgftype));prefix=OUT/tag
        args=[str(BIN),"--bfile",relative(fixture/"synthetic"),"--pheno",relative(fixture/"synthetic.pheno"),"--envir-name",en,"--threads","1","--mac","0","--out",relative(prefix)]
        if input_mode=="native":args += ["--pheno-name","Quantitative,Binary,Time:Event,Ordinal","--covar-name","X1,"+en+",PC1,PC2","--regression-model","auto"]
        else:args += ["--resid-name","R_quant,R_binary,R_survival,R_ordinal"]
        if mode=="plus":args += ["--grm-file",relative(fixture/"synthetic.grm.tsv")]
        if mode=="mix":args += ["--method","SPAGxEmix","--pc-cols","PC1,PC2"]
        proc=subprocess.run(args,cwd=CPP,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding="utf-8",errors="replace",timeout=90)
        (OUT/(tag+".console.txt")).write_text(proc.stdout,encoding="utf-8")
        if proc.returncode:failures.append(f"{tag}: process exited {proc.returncode}: {proc.stdout[-1000:]}");continue
        start=len(checks)
        for trait,(native_col,resid_col) in trait_cols.items():
            column=native_col if input_mode=="native" else resid_col
            result={r["ID"]:r for r in tab(str(prefix)+"."+column+(".SPAGxEmix" if mode=="mix" else ".SPAGxE"))}
            if len(result)!=48:failures.append(f"{tag}:{trait}: expected 48 variants, got {len(result)}")
            for ident,row in result.items():
                if ident=="rsOracle48":continue
                kr=(trait,en,mode,ident);rr=source[kr,"tight"];fr=formula[kr,"paper_full"]
                wanted_spa=number(rr["LOG10P_SPA"])
                check(tag,f"{trait}:{ident}:SPA_full_reference_delta",wanted_spa,number(row["LOG10P_SPA_Gx"+en]),required=False)
                check(tag,f"{trait}:{ident}:marginal",number(rr["LOG10P_G"]),number(row["LOG10P_G"]))
                check(tag,f"{trait}:{ident}:normal_Z_magnitude",abs(number(fr["Z_NORM"])),abs(number(row["Z_Norm_Gx"+en])))
                check(tag,f"{trait}:{ident}:branch",2 if rr["BRANCH"]=="B" else 1,number(row["BRANCH_Gx"+en]),0)
                if input_mode=="native" and mode!="plus" and rr["BRANCH"]=="B":
                    wanted_wald=number(rr["LOG10P_WALD"])
                    check(tag,f"{trait}:{ident}:Wald",wanted_wald,number(row["LOG10P_Wald_Gx"+en]))
                    wanted_final=number(rr["LOG10P_FINAL"])
                else:wanted_final=wanted_spa
                check(tag,f"{trait}:{ident}:final_full_reference_delta",wanted_final,number(row["LOG10P_Gx"+en]),required=False)
                if cgftype=="guarded_partial_normal":
                    old=source[kr,"default"]
                    legacy_default_deltas.append({"run":tag,"trait":trait,"ID":ident,"original_default_log10p_SPA":number(old["LOG10P_SPA"]),"native_log10p_SPA":number(row["LOG10P_SPA_Gx"+en]),"difference":abs(number(old["LOG10P_SPA"])-number(row["LOG10P_SPA_Gx"+en]))})
        runs.append({"run":tag,"comparison_count":len(checks)-start})

for name,rows in (("paper_cli_vs_original_R.tsv",checks),("native_vs_original_default_solver.tsv",legacy_default_deltas)):
    if rows:
        with (OUT/name).open("w",encoding="utf-8",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter="\t");writer.writeheader();writer.writerows(rows)
summary={"comparison_role":"strict model invariants; original-R SPA/final deltas are descriptive","cgf_policy":"guarded_partial_normal","passed":not failures,"binary":str(BIN),"binary_sha256":hashlib.sha256(BIN.read_bytes()).hexdigest(),"run_count":len(runs),"comparison_count":len(checks),"strict_invariant_checks":sum(x["required"] for x in checks),"cgf_delta_comparisons":sum(not x["required"] for x in checks),"max_abs_error":max((x["abs_error"] for x in checks if x["required"] and math.isfinite(x["abs_error"])),default=0),"max_original_default_solver_difference":max((x["difference"] for x in legacy_default_deltas if math.isfinite(x["difference"])),default=0),"runs":runs,"failure_count":len(failures),"failures":failures}
(OUT/"paper_cli_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps({k:summary[k] for k in ("passed","run_count","comparison_count","max_abs_error","max_original_default_solver_difference","failure_count")},indent=2))
if failures:print("\n".join(failures[:25]))
sys.exit(0 if not failures else 1)
