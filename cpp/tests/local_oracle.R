# Independent Bin(H,q) local-ancestry oracle, scalar R + uniroot.
# Run from cpp: Rscript tests/local_oracle.R [fixture-root]
# Native raw-trait model comparisons use R lm/glm/coxph/clm, never C++ fits.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))args[1] else "examples/local"
options(digits=17)
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ordinal))
ph<-read.delim(file.path(root,"samples.pheno"),check.names=FALSE)
native<-list()
native$quantitative<-resid(lm(Quantitative~X+E,ph))
native$binary<-resid(glm(Binary~X+E,ph,family=binomial()),type="response")
native$survival<-resid(coxph(Surv(Time,Event)~X+E,ph,ties="breslow"),type="martingale")
ord<-clm(as.factor(Ordinal)~X+E,data=ph)
eta<-as.numeric(as.matrix(ph[,c("X","E")])%*%ord$beta)
nu<-cbind(0,plogis(outer(eta,ord$alpha,function(a,b)b-a)),1)
yy<-as.integer(factor(ph$Ordinal));native$ordinal<-nu[cbind(seq_len(nrow(ph)),yy)]+nu[cbind(seq_len(nrow(ph)),yy+1)]-1
residualNames<-c(quantitative="R_quant",binary="R_binary",survival="R_survival",ordinal="R_ordinal",masked="R_masked")
logadd<-function(a,b){z<-pmax(a,b);out<-z+log1p(exp(-abs(a-b)));out[is.infinite(z)&z<0]<--Inf;out}
logsum<-function(x){if(all(x==-Inf))return(-Inf);z<-max(x);z+log(sum(exp(x-z)))}
source("tests/partial_reference.R")
cctlog<-function(lp){
 lp<-lp[!is.na(lp)];if(!length(lp))return(NA_real_);if(any(lp==Inf))return(Inf)
 # Match the documented legacy P=1 boundary convention, using a stable
 # cotangent expansion for tails too small for ordinary tan subtraction.
 lp[lp==0]<--log10(.999);p<-10^(-lp)
 tt<-mean(ifelse(p<1e-16,1/(pi*p),tan(pi*(.5-p))))
 out<-if(tt>1e15)log10(pi)+log10(tt) else -log10(pcauchy(tt,lower.tail=FALSE))
 max(0,out)
}
localscore<-function(g,h,w,q,hybrid=FALSE,cutoff=2){
 mu<-q*sum(h*w);variance<-q*(1-q)*sum(h*w*w);dev<-sum(g*w)-mu;z<-dev/sqrt(variance)
 if(!is.finite(z)||variance<=0)return(c(LOG10P=NA_real_,Z=NA_real_,VAR=variance))
 normal<-pnorm(-abs(z),log.p=TRUE)+log(2)
 if(abs(z)<=cutoff)return(c(LOG10P=-normal/log(10),Z=z,VAR=variance))
 value<-reference_guarded_spa(w,h,q,mu,abs(dev),variance,partial=hybrid)
 c(LOG10P=value,Z=z,VAR=variance)
}
wald<-function(dat,g,E,allH,trait){
 design<-cbind(1,dat$X,dat$E);rank<-qr(design,tol=1e-10)$rank
 for(a in seq_len(ncol(allH))){candidate<-cbind(design,allH[,a]);rk<-qr(candidate,tol=1e-10)$rank;if(rk>rank){design<-candidate;rank<-rk}}
 d<-as.data.frame(design[,-1,drop=FALSE]);names(d)<-paste0("C",seq_len(ncol(d)));d$g<-g;d$ge<-g*E
 rhs<-paste(c(names(d)[grepl("^C",names(d))],"g","ge"),collapse="+")
 ans<-tryCatch({
  if(trait=="quantitative"){d$y<-dat$Quantitative;p<-coef(summary(lm(as.formula(paste("y~",rhs)),d)))["ge","Pr(>|t|)"]}
  else if(trait=="binary"){d$y<-dat$Binary;p<-coef(summary(glm(as.formula(paste("y~",rhs)),d,family=binomial())))["ge","Pr(>|z|)"]}
  else if(trait=="survival"){d$time<-dat$Time;d$event<-dat$Event;p<-coef(summary(coxph(as.formula(paste("Surv(time,event)~",rhs)),d,ties="breslow")))["ge","Pr(>|z|)"]}
  else {d$y<-factor(dat$Ordinal);p<-coef(summary(clm(as.formula(paste("y~",rhs)),data=d)))["ge","Pr(>|z|)"]}
  -log10(p)
 },error=function(e)NA_real_)
 as.numeric(ans)
}
readLocal<-function(directory){
 manifest<-read.delim(file.path(directory,"cohort.local.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
 gl<-hl<-list();ids<-NULL;meta<-NULL
 for(a in seq_len(nrow(manifest))){
  gg<-read.delim(file.path(directory,manifest$genotype[a]),check.names=FALSE);hh<-read.delim(file.path(directory,manifest$hapcount[a]),check.names=FALSE)
  if(is.null(ids)){ids<-names(gg)[-(1:5)];meta<-gg[,1:5,drop=FALSE]}
  gl[[a]]<-as.matrix(gg[,-(1:5)]);hl[[a]]<-as.matrix(hh[,-(1:5)])
 }
 list(ancestry=manifest$ancestry,G=gl,H=hl,IID=ids,meta=meta)
}
rows<-list();ki<-0L;combined<-list();kc<-0L;moments<-list();km<-0L
for(dataset in c("two","three")){
 source<-readLocal(if(dataset=="two")root else file.path(root,"three_ancestry"))
 ph0<-ph[match(source$IID,ph$IID),,drop=FALSE]
 for(mode in c("native","residual"))for(trait in if(mode=="native")names(native)else names(residualNames)){
  residual<-if(mode=="native")native[[trait]][match(source$IID,ph$IID)]else ph0[[residualNames[trait]]]
  baseKeep<-which(is.finite(residual)&is.finite(ph0$E)&is.finite(ph0$Eshift));dat0<-ph0[baseKeep,,drop=FALSE];r0<-residual[baseKeep]
  for(j in seq_len(nrow(source$meta))){
   allH0<-sapply(source$H,function(x)x[j,baseKeep]);keep<-which(complete.cases(allH0));allH<-allH0[keep,,drop=FALSE];R<-r0[keep];dat<-dat0[keep,,drop=FALSE]
   for(a in seq_along(source$ancestry)){
    h<-allH[,a];g<-source$G[[a]][j,baseKeep][keep];g[is.na(g)&h==0]<-0
    observed<-is.finite(g);af<-if(sum(h[observed])>0)sum(g[observed])/sum(h[observed])else NA_real_;q<-min(af,1-af)
    missing<-mean(!observed);g[!observed]<-h[!observed]*af;flipped<-is.finite(af)&&af>.5
    if(flipped)g<-h-g
    valid<-is.finite(q)&&q>0&&q<1&&sum(h*R^2)>0
    zg<-if(valid)sum((g-q*h)*R)/sqrt(q*(1-q)*sum(h*R^2))else NA_real_
    marginal<-if(valid)-(pnorm(-abs(zg),log.p=TRUE)+log(2))/log(10)else NA_real_
    branch<-if(valid&&2*pnorm(-abs(zg))<=.001)"B"else "A"
    for(en in if(mode=="native")"E"else c("E","Eshift")){
     E<-dat[[en]];lambda<-if(valid)sum(h*E*R^2)/sum(h*R^2)else NA_real_
     w<-if(valid){if(branch=="A")(E-lambda)*R else E*resid(lm(R~g))}else rep(0,length(R))
     pw<-if(valid&&branch=="B"&&mode=="native")wald(dat,g,E,allH,trait)else NA_real_
     for(cgf in c("exact","guarded_partial_normal")){
      test<-if(valid)localscore(g,h,w,q,cgf=="guarded_partial_normal")else c(LOG10P=NA_real_,Z=NA_real_,VAR=NA_real_)
      final<-if(is.finite(pw))cctlog(c(test["LOG10P"],pw))else test["LOG10P"]
      ki<-ki+1L;rows[[ki]]<-data.frame(DATASET=dataset,MODE=mode,TRAIT=trait,ENV=en,CGF=cgf,ID=source$meta$ID[j],ANCESTRY=source$ancestry[a],N=length(R),N_H_MISSING=length(baseKeep)-length(R),AF=af,MAC=sum(h)*q,MISS_RATE=missing,FLIPPED=as.integer(flipped),BRANCH=branch,LOG10P_G=marginal,Z_G=zg,Z_Norm=test["Z"],VAR=test["VAR"],LOG10P_SPA=test["LOG10P"],LOG10P_WALD=pw,LOG10P_FINAL=final,check.names=FALSE)
     }
     km<-km+1L;moments[[km]]<-data.frame(DATASET=dataset,MODE=mode,TRAIT=trait,ENV=en,ID=source$meta$ID[j],ANCESTRY=source$ancestry[a],LAMBDA=lambda,MEAN=q*sum(h*w),SCORE=sum(g*w),VAR=q*(1-q)*sum(h*w*w))
    }
   }
  }
 }
}
ref<-do.call(rbind,rows);mom<-do.call(rbind,moments)
groups<-unique(ref[,c("DATASET","MODE","TRAIT","ENV","CGF","ID")])
for(i in seq_len(nrow(groups))){key<-groups[i,];take<-rep(TRUE,nrow(ref));for(nm in names(key))take<-take&ref[[nm]]==key[[nm]];values<-ref$LOG10P_FINAL[take];kc<-kc+1L;combined[[kc]]<-cbind(key,LOG10P_ANCESTRY_CCT=cctlog(values),N_ANCESTRIES=sum(!is.na(values)))}
wr<-function(x,file)write.table(x,file.path(root,file),sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")
write_partial_reference_audit(file.path(root,"partial_reference_audit.tsv"))
wr(ref,"local_reference.tsv");wr(mom,"local_moments.tsv");wr(do.call(rbind,combined),"local_combined_reference.tsv")
wr(data.frame(IID=ph$IID,R_quant_native=native$quantitative,R_binary_native=native$binary,R_survival_native=native$survival,R_ordinal_native=native$ordinal),"native_R_residuals.tsv")
# Compare the preserved original function on its valid input domain. Missing G/H
# and high-AF variants are excluded because the source has known sum(NA), 2*q
# imputation and 2-G instead of H-G defects, which the C++ port corrects.
leg<-new.env(parent=globalenv());sys.source("../R/SPAGxECCT.R",envir=leg)
legacyRoot<-leg$fastgetroot_H1_new_local
leg$fastgetroot_H1_new_local<-function(init.t,R,s,MAFVec,haplo_numVec,tol=1e-12,maxiter=100)legacyRoot(init.t,R,s,MAFVec,haplo_numVec,tol=1e-12,maxiter=maxiter)
# Set the original Cox Wald's tie convention explicitly in memory. On this
# untied fixture the estimate is unchanged, and original source is preserved.
leg$coxph<-function(...)survival::coxph(...,ties="breslow")
legacy<-list();kl<-0L
for(dataset in c("two","three")){
 src<-readLocal(if(dataset=="two")root else file.path(root,"three_ancestry"))
 for(trait in names(native))for(j in seq_len(nrow(src$meta)))for(a in seq_along(src$ancestry)){
  g<-src$G[[a]][j,];h<-src$H[[a]][j,];allH<-sapply(src$H,function(x)x[j,]);af<-sum(g)/sum(h)
  if(anyNA(g)||anyNA(allH)||!is.finite(af)||af<=0||af>.5)next
  expected<-ref[ref$DATASET==dataset&ref$MODE=="native"&ref$TRAIT==trait&ref$ENV=="E"&ref$CGF=="exact"&ref$ID==src$meta$ID[j]&ref$ANCESTRY==src$ancestry[a],]
  phen<-data.frame(Y=if(trait=="binary")ph$Binary else if(trait=="ordinal")ph$Ordinal else ph$Quantitative,surv.time=ph$Time,event=ph$Event)
  ans<-tryCatch(suppressWarnings(leg$SPAGxEmixCCT_localance_one_SNP(traits=if(trait=="ordinal")"categorical"else trait,g=g,haplo_numVec=h,R=native[[trait]],E=ph$E,Phen.mtx=phen,Cova.mtx=as.matrix(ph["X"]),Cova.haplo.mtx=allH[,seq_len(ncol(allH)-1),drop=FALSE],min.maf=1e-5)),error=function(e)e)
  for(field in c("SPA","FINAL",if(expected$BRANCH=="B")"WALD")){
   observed<-if(inherits(ans,"error"))NA_real_ else -log10(ans[switch(field,SPA=3,FINAL=5,WALD=4)])
   oracle<-expected[[paste0("LOG10P_",field)]];error<-abs(observed-oracle);kl<-kl+1L
   legacy[[kl]]<-data.frame(DATASET=dataset,TRAIT=trait,ID=src$meta$ID[j],ANCESTRY=src$ancestry[a],FIELD=field,ORIGINAL=observed,INDEPENDENT=oracle,ERROR=error,PASS=is.finite(error)&&error<1e-6,NOTE=if(inherits(ans,"error"))conditionMessage(ans)else "Original source unchanged; root tolerance 1e-12 in memory; valid low-AF complete-input domain")
  }
 }
}
legacyComparison<-do.call(rbind,legacy);wr(legacyComparison,"original_local_comparison.tsv")
writeLines(c(paste("Independent local reference rows:",nrow(ref)),paste("Dense moment rows:",nrow(mom)),paste("Ancestry CCT rows:",length(combined)),paste("Original R valid-domain comparisons:",nrow(legacyComparison)),paste("Original R passed:",sum(legacyComparison$PASS)),paste("Original R max log10p error:",max(legacyComparison$ERROR,na.rm=TRUE)),"R lm/glm/coxph Breslow/clm native residuals and full-model Wald; H rank-basis ancestry adjustment; scalar full and guarded partial-normal Bin(H,q) CGFs and independent bracketed root solving.","Missing H: common complete-H mask across ancestries; missing G: H*q; H=0 structural zero; high-AF minor dosage: H-G.","All data are deterministic examples; R is needed only for independent validation, not C++ analysis."),file.path(root,"local_oracle_summary.txt"))
cat(paste("Wrote",nrow(ref),"local reference rows and",length(combined),"combined rows.\n"))
