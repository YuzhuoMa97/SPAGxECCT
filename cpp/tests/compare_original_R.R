# Cross-check preserved original R routines against an independent oracle.
# Run from cpp: Rscript tests/compare_original_R.R
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ordinal))
suppressPackageStartupMessages(library(dplyr))
source("../R/SPAGxECCT.R"); source("../R/SPAGxEPlus_functions_MYZ.R")
# Optional tighter root tolerance diagnoses legacy Newton stopping error.
# Preserved source files are never edited.
refined <- identical(Sys.getenv("SPAGXE_REFINE_LEGACY_ROOT"),"1")
if(refined){
 legacy_root <- fastgetroot_H1_new
 fastgetroot_H1_new <- function(init.t,R,s,MAFVec,tol=1e-12,maxiter=100) legacy_root(init.t,R,s,MAFVec,tol=1e-12,maxiter=maxiter)
}
ph<-read.delim("examples/synthetic/synthetic.pheno",check.names=FALSE)
G<-as.matrix(read.delim("examples/synthetic/synthetic.genotypes.tsv",row.names=1))
ref<-read.delim("examples/synthetic/oracle_reference.tsv",check.names=FALSE)
gr<-read.delim("examples/synthetic/synthetic.grm.tsv")
n<-nrow(ph);phi<-diag(n);phi[cbind(seq(1,n,2),seq(2,n,2))]<-.125;phi[cbind(seq(2,n,2),seq(1,n,2))]<-.125
Cova.mtx<-as.matrix(ph[,c("X1","PC1","PC2")]);E<-ph$E
rows<-list();k<-0L
for(trait in c("quantitative","binary","survival","ordinal")){
 R<-switch(trait,quantitative=ph$R_quant,binary=ph$R_binary,survival=ph$R_survival,ordinal=ph$R_ordinal)
 Phen.mtx<-data.frame(IID=ph$IID,Y=if(trait=="binary")ph$Binary else if(trait=="ordinal")ph$Ordinal else ph$Quantitative,surv.time=ph$Time,event=ph$Event)
 la<-as.numeric(crossprod(R,phi%*%(E*R))/crossprod(R,phi%*%R));rn<-(E-la)*R
 obj<-list(ResidMat=data.frame(SubjID=ph$IID,Resid=R,RE=E*R,R.new=rn),R=R,R_GRM_R=as.numeric(crossprod(R,phi%*%R)),R_GRM_RE=as.numeric(crossprod(R,phi%*%(E*R))),R.new=rn,R.new_GRM_R.new=as.numeric(crossprod(rn,phi%*%rn)),sparseGRM=gr)
 for(mode in c("base","mix","plus"))for(j in seq_len(47)){
  g<-G[,j]
  ans<-tryCatch({
   if(mode=="base")SPAGxE_CCT_one_SNP(traits=if(trait=="ordinal")"categorical" else trait,g=g,R=R,E=E,Phen.mtx=Phen.mtx,Cova.mtx=Cova.mtx)
   else if(mode=="mix")SPAGxEmix_CCT_one_SNP(traits=if(trait=="ordinal")"categorical" else trait,g=g,R=R,E=E,Phen.mtx=Phen.mtx,Cova.mtx=Cova.mtx,topPCs=as.matrix(ph[,c("PC1","PC2")]))
   else SPAGxE_Plus_one_SNP(g=g,E=E,Phen.mtx=Phen.mtx,Cova.mtx=Cova.mtx,obj.SPAGxE_Plus_Nullmodel=obj)
  },error=function(e)rep(NA_real_,10))
  rr<-ref[ref$ID==sprintf("rsOracle%02d",j)&ref$MODE==mode&ref$TRAIT==trait,]
  observed<-if(mode=="plus")-log10(ans[3]) else -log10(ans[5])
  expected<-if(mode=="plus")rr$LOG10P_EXACT_SPA else rr$LOG10P_EXACT_FINAL
  k<-k+1L;rows[[k]]<-data.frame(ID=sprintf("rsOracle%02d",j),MODE=mode,TRAIT=trait,original_log10p=observed,independent_log10p=expected,abs_error=abs(observed-expected),pass=is.finite(observed)&&abs(observed-expected)<if(refined)1e-6 else 5e-4)
 }
}
out<-do.call(rbind,rows)
write.table(out,if(refined)"../../validation/original_R_refined_comparison.tsv" else "../../validation/original_R_comparison.tsv",sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c(paste("Rows:",nrow(out)),paste("Absolute log10p tolerance:",if(refined)1e-6 else 5e-4),paste("Passed:",sum(out$pass)),paste("Failed:",sum(!out$pass)),paste("Maximum absolute -log10p error:",max(out$abs_error))),if(refined)"../../validation/original_R_refined_summary.txt" else "../../validation/original_R_summary.txt")
cat(paste("Original R comparisons:",nrow(out),"passed",sum(out$pass),"max error",max(out$abs_error),"\n"))
if(any(!out$pass))quit(status=1)
