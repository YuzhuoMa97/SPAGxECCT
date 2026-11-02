# Validate local-BGEN using the independent guarded SPA policy and preserved-R deltas.
# Run from cpp: Rscript tests/verify_local_bgen_R.R [binary] [output-directory]
# Requires R packages survival and ordinal for the four-phenotype checks only.
# The paired text G/H are computed from stored quantized probabilities by the
# fixture generator, independently of the native BGEN reader.
options(digits=17)
source("tests/partial_reference.R")
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ordinal))
args<-commandArgs(trailingOnly=TRUE)
binary<-if(length(args))args[1]else if(.Platform$OS.type=="windows")"build/bin/SPAGxE.exe"else "build/bin/SPAGxE"
out<-if(length(args)>1)args[2]else "results/local_bgen_R"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
base<-"examples/local/bgen_demo"
stored<-read.delim(file.path(base,"samples.pheno"),check.names=FALSE)
ph<-read.delim("examples/local/samples.pheno",check.names=FALSE)
stopifnot(identical(stored$IID,ph$IID))
residuals<-list(quantitative=resid(lm(Quantitative~X+E,ph)),binary=resid(glm(Binary~X+E,ph,family=binomial()),type="response"),survival=resid(coxph(Surv(Time,Event)~X+E,ph,ties="breslow"),type="martingale"))
ordinalFit<-clm(as.factor(Ordinal)~X+E,data=ph)
eta<-as.numeric(as.matrix(ph[,c("X","E")])%*%ordinalFit$beta)
nu<-cbind(0,plogis(outer(eta,ordinalFit$alpha,function(a,b)b-a)),1)
yy<-as.integer(factor(ph$Ordinal))
residuals$ordinal<-nu[cbind(seq_len(nrow(ph)),yy)]+nu[cbind(seq_len(nrow(ph)),yy+1)]-1
loadOriginal<-function(tight){
 e<-new.env(parent=globalenv());sys.source("../R/SPAGxECCT.R",envir=e)
 if(tight){r<-e$fastgetroot_H1_new_local;e$fastgetroot_H1_new_local<-function(init.t,R,s,MAFVec,haplo_numVec,tol=1e-12,maxiter=100)r(init.t,R,s,MAFVec,haplo_numVec,tol=1e-12,maxiter=maxiter)}
 e
}
original<-list(default=loadOriginal(FALSE),tight=loadOriginal(TRUE))
readPaired<-function(stem){
 mf<-read.delim(file.path(base,paste0(stem,".local.tsv")),check.names=FALSE)
 G<-H<-list()
 for(a in seq_len(nrow(mf))){
  gg<-read.delim(file.path(base,mf$genotype[a]),check.names=FALSE)
  hh<-read.delim(file.path(base,mf$hapcount[a]),check.names=FALSE)
  stopifnot(identical(names(gg)[-(1:5)],stored$IID),identical(gg[,1:5],hh[,1:5]))
  G[[a]]<-as.matrix(gg[,-(1:5)]);H[[a]]<-as.matrix(hh[,-(1:5)])
 }
 list(ancestry=mf$ancestry,meta=gg[,1:5],G=G,H=H)
}
rows<-excluded<-checks<-list();kr<-ke<-kc<-0L
waldDiagnostics<-list();kw<-0L
canonicalWald<-function(g,allH,trait,E){
 design<-cbind(1,ph$X,E);rk<-qr(design,tol=1e-10)$rank
 for(a in seq_len(ncol(allH))){candidate<-cbind(design,allH[,a]);newRank<-qr(candidate,tol=1e-10)$rank;if(newRank>rk){design<-candidate;rk<-newRank}}
 design<-cbind(design,g,g*E);rank<-qr(design,tol=1e-10)$rank
 nuisance<-design[,-ncol(design),drop=FALSE];nuisanceRank<-qr(nuisance,tol=1e-10)$rank
 base<-list(RANK=rank,COLUMNS=ncol(design),CONDITION=kappa(design,exact=TRUE),NUISANCE_RANK=nuisanceRank,INTERACTION_IDENTIFIABLE=rank==nuisanceRank+1,REDUCED_COLUMNS=NA_integer_,CONVERGED=FALSE,TIGHT_LOG10P=NA_real_,DEFAULT_LOG10P=NA_real_,DEFAULT_REEVAL_LOG10P=NA_real_,DEFAULT_SCORE_NORM=NA_real_,TIGHT_SCORE_NORM=NA_real_)
 if(!base$INTERACTION_IDENTIFIABLE)return(base)
 # Keep a full-rank nuisance basis while always reserving the final column
 # for the interaction. A redundant main G does not invalidate estimable G:E.
 basis<-matrix(numeric(0),nrow(design),0);basisRank<-0L
 for(column in seq_len(ncol(nuisance))){candidate<-cbind(basis,nuisance[,column]);newRank<-qr(candidate,tol=1e-10)$rank;if(newRank>basisRank){basis<-candidate;basisRank<-newRank}}
 design<-cbind(basis,design[,ncol(design)]);base$REDUCED_COLUMNS<-ncol(design)
 d<-as.data.frame(design[,-1,drop=FALSE]);names(d)<-paste0("C",seq_len(ncol(d)))
 col<-names(d)[ncol(d)]
 tryCatch({
  if(trait=="quantitative"){
   d$y<-ph$Quantitative;fit<-lm(y~.,d);p<-coef(summary(fit))[col,"Pr(>|t|)"];pd<-p;base$CONVERGED<-TRUE
  }else if(trait=="binary"){
   d$y<-ph$Binary;default<-glm(y~.,d,family=binomial());pd<-coef(summary(default))[col,"Pr(>|z|)"]
   defaultMu<-fitted(default);defaultInfo<-crossprod(design,design*as.numeric(defaultMu*(1-defaultMu)));defaultSE<-sqrt(solve(defaultInfo)[ncol(design),ncol(design)])
   base$DEFAULT_REEVAL_LOG10P<--log10(2*pnorm(-abs(tail(coef(default),1)/defaultSE)))
   base$DEFAULT_SCORE_NORM<-max(abs(crossprod(design,ph$Binary-defaultMu)))
   fit<-glm(y~.,d,family=binomial(),control=glm.control(epsilon=1e-12,maxit=200));mu<-fitted(fit);info<-crossprod(design,design*as.numeric(mu*(1-mu)))
   se<-sqrt(solve(info)[ncol(design),ncol(design)]);z<-tail(coef(fit),1)/se;p<-2*pnorm(-abs(z));base$CONVERGED<-fit$converged
   base$TIGHT_SCORE_NORM<-max(abs(crossprod(design,ph$Binary-mu)))
  }else if(trait=="survival"){
   d$time<-ph$Time;d$event<-ph$Event;rhs<-paste(names(d)[grepl("^C",names(d))],collapse="+");form<-as.formula(paste("Surv(time,event)~",rhs))
   fit<-coxph(form,d,ties="breslow",control=coxph.control(eps=1e-12,iter.max=100));p<-coef(summary(fit))[col,"Pr(>|z|)"];pd<-p;base$CONVERGED<-TRUE
  }else{
   d$y<-factor(ph$Ordinal);fit<-clm(y~.,data=d,control=clm.control(gradTol=1e-10,maxIter=1000));p<-coef(summary(fit))[col,"Pr(>|z|)"];pd<-p;base$CONVERGED<-all(is.finite(coef(fit)))&&max(abs(fit$gradient))<1e-6
   base$TIGHT_SCORE_NORM<-max(abs(fit$gradient))
  }
  base$TIGHT_LOG10P<--log10(p);base$DEFAULT_LOG10P<--log10(pd);base
 },error=function(e)base)
}
combine<-function(logs){p<-10^-logs;p[p==1]<-.999;t<-mean(ifelse(p<1e-16,1/(pi*p),tan(pi*(.5-p))));-log10(if(t>1e15)1/(pi*t)else pcauchy(t,lower.tail=FALSE))}
bracketLocalSPA<-function(g,h,R,E,q,cutoff){
 marginal<-(sum(g*R)-q*sum(h*R))/sqrt(q*(1-q)*sum(h*R^2))
 branch<-if(2*pnorm(-abs(marginal))<=cutoff)2 else 1
 w<-if(branch==2)E*resid(lm(R~g))else (E-sum(h*E*R^2)/sum(h*R^2))*R
 mean<-q*sum(h*w);variance<-q*(1-q)*sum(h*w*w);dev<-abs(sum(g*w)-mean);z<-dev/sqrt(variance)
 if(z<=2)return(-(pnorm(-z,log.p=TRUE)+log(2))/log(10))
 reference_guarded_spa(w,h,q,mean,dev,variance,partial=TRUE)
}
check<-function(tag,field,expected,actual,tol=2e-5){
 required<-!field%in%c("SPA","FINAL_A","FINAL_IDENTIFIABLE_ORIGINAL_SPA_LEG")
 same<-is.na(expected)&&is.na(actual)||identical(as.numeric(expected),as.numeric(actual))
 error<-if(same)0 else abs(actual-expected)
 passed<-same||is.finite(error)&&error<tol
 kc<<-kc+1L;checks[[kc]]<<-data.frame(RUN=tag,FIELD=field,EXPECTED=expected,NATIVE=actual,ABS_DELTA=error,PASS=passed,REQUIRED=required)
}
run<-function(stem,input,cutoff){
 src<-readPaired(stem);tag<-paste(stem,input,if(cutoff==1)"force_B"else "default_branch",sep="_")
 prefix<-file.path(out,tag)
 cmd<-c("--method","SPAGxEmix_local","--local-bgen-manifest",file.path(base,paste0(stem,".localbgen.tsv")),"--envir-name","E","--maf","0","--mac","0","--geno","1","--spagxe-marginal-cutoff",format(cutoff,scientific=FALSE),"--out",prefix)
 isResidual<-input!="native"
 supplied<-stored
 if(input=="tail_residual"){
  target<-src$G[[1]][match("point_calls",src$meta$ID),]
  supplied$R<-resid(lm((15*target*stored$E+sin(seq_along(target)*.37))~stored$E+ph$X))
  suppliedPath<-paste0(prefix,".input.pheno");write.table(supplied,suppliedPath,sep="\t",quote=FALSE,row.names=FALSE)
 }else suppliedPath<-file.path(base,"samples.pheno")
 if(isResidual)cmd<-c(cmd,"--pheno",suppliedPath,"--resid-name","R") else cmd<-c(cmd,"--pheno","examples/local/samples.pheno","--pheno-name","Quantitative,Binary,Ordinal,Time:Event","--covar-name","X,E","--regression-model","auto")
 log<-paste0(prefix,".console.txt");status<-system2(binary,shQuote(cmd),stdout=log,stderr=log)
 if(status!=0)stop(paste("Native local-BGEN process failed:",tag,"; see",log))
 traits<-if(isResidual)"residual"else c("quantitative","binary","ordinal","survival")
 nativeName<-c(residual="R",quantitative="Quantitative",binary="Binary",ordinal="Ordinal",survival="Time_Event")
 for(trait in traits){
  native<-read.delim(paste0(prefix,".",nativeName[trait],".SPAGxE_local"),check.names=FALSE)
  R<-if(trait=="residual")supplied$R else residuals[[trait]]
  E<-if(trait=="residual")stored$E else ph$E
  phen<-data.frame(Y=switch(trait,residual=supplied$R,quantitative=ph$Quantitative,binary=ph$Binary,ordinal=ph$Ordinal,survival=ph$Quantitative),surv.time=ph$Time,event=ph$Event)
  for(j in seq_len(nrow(src$meta)))for(a in seq_along(src$ancestry)){
   g<-src$G[[a]][j,];h<-src$H[[a]][j,];allH<-sapply(src$H,function(x)x[j,]);q<-sum(g)/sum(h)
   reason<-if(anyNA(allH))"unknown_H"else if(anyNA(g))"missing_G"else if(!is.finite(q)||sum(h)==0)"zero_trials"else if(q>.5)"high_AF"else if(q<=0||q>=1)"boundary_AF"else if(qr(cbind(1,g))$rank<2)"constant_G"else "safe_complete_minor"
   if(reason!="safe_complete_minor"){
    ke<<-ke+1L;excluded[[ke]]<<-data.frame(RUN=tag,TRAIT=trait,ID=src$meta$ID[j],ANCESTRY=src$ancestry[a],REASON=reason)
    next
   }
   nn<-native[native$ID==src$meta$ID[j]&native$ANCESTRY==src$ancestry[a],,drop=FALSE]
   if(nrow(nn)!=1)stop("Expected exactly one native ancestry result")
   canonical<-if(!isResidual)canonicalWald(g,allH,trait,E)else NULL
   independentSPA<-bracketLocalSPA(g,h,R,E,q,cutoff)
   for(rootMode in names(original)){
    ans<-tryCatch(suppressWarnings(original[[rootMode]]$SPAGxEmixCCT_localance_one_SNP(traits=if(trait%in%c("quantitative","residual"))"quantitative"else if(trait=="ordinal")"categorical"else trait,g=g,haplo_numVec=h,R=R,E=E,Phen.mtx=phen,Cova.mtx=as.matrix(ph["X"]),Cova.haplo.mtx=allH[,seq_len(ncol(allH)-1),drop=FALSE],epsilon=cutoff,min.maf=1e-8)),error=function(e)e)
    val<-if(inherits(ans,"error"))rep(NA_real_,10)else ans
    branch<-if(is.finite(val[7])&&val[7]<=cutoff)2 else 1
    kr<<-kr+1L;rows[[kr]]<<-data.frame(FORMAT=stem,INPUT=input,CUTOFF=cutoff,TRAIT=trait,ID=src$meta$ID[j],ANCESTRY=src$ancestry[a],ROOT=rootMode,ALT_AF=q,ORIGINAL_STATUS=if(inherits(ans,"error"))conditionMessage(ans)else "OK",ORIGINAL_BRANCH=branch,NATIVE_BRANCH=nn$BRANCH_GxE,ORIGINAL_LOG10P_SPA=-log10(val[3]),INDEPENDENT_LOG10P_SPA=independentSPA,NATIVE_LOG10P_SPA=nn$LOG10P_SPA_GxE,NATIVE_Z_NORM=nn$Z_Norm_GxE,NATIVE_SPA_STATUS=nn$SPA_STATUS_GxE,ABS_SPA_DELTA=abs(-log10(val[3])-nn$LOG10P_SPA_GxE),ORIGINAL_LOG10P_WALD=-log10(val[4]),NATIVE_LOG10P_WALD=nn$LOG10P_Wald_GxE,ORIGINAL_LOG10P_FINAL=-log10(val[5]),NATIVE_LOG10P_FINAL=nn$LOG10P_GxE,check.names=FALSE)
    if(rootMode=="tight"){
     key<-paste(tag,trait,src$meta$ID[j],src$ancestry[a],sep=":")
     check(key,"SPA",-log10(val[3]),nn$LOG10P_SPA_GxE)
     check(key,"SPA_INDEPENDENT_BRACKET",independentSPA,nn$LOG10P_SPA_GxE)
     check(key,"Z_G",val[10],nn$Z_G,2e-7)
     check(key,"AF",val[1],nn$ALT_FREQ,1e-12)
     check(key,"BRANCH",branch,nn$BRANCH_GxE,1e-12)
     if(!isResidual){
      if(branch==1)check(key,"FINAL_A",-log10(val[5]),nn$LOG10P_GxE)
      if(branch==2){
       wd<-abs(-log10(val[4])-nn$LOG10P_Wald_GxE);fd<-abs(-log10(val[5])-nn$LOG10P_GxE)
       tightFinalOriginal<-if(is.finite(canonical$TIGHT_LOG10P))combine(c(-log10(val[3]),canonical$TIGHT_LOG10P))else NA_real_
       tightFinalNative<-if(is.finite(canonical$TIGHT_LOG10P))combine(c(nn$LOG10P_SPA_GxE,canonical$TIGHT_LOG10P))else NA_real_
       boundaryDifference<-(val[3]==1)!=(nn$LOG10P_SPA_GxE==0)
       kw<<-kw+1L;waldDiagnostics[[kw]]<<-data.frame(RUN=key,ORIGINAL_WALD=-log10(val[4]),NATIVE_WALD=nn$LOG10P_Wald_GxE,ABS_WALD_DELTA=wd,ORIGINAL_FINAL=-log10(val[5]),NATIVE_FINAL=nn$LOG10P_GxE,ABS_FINAL_DELTA=fd,ORIGINAL_WITHIN_2e_5=is.finite(wd)&&wd<2e-5&&is.finite(fd)&&fd<2e-5,EXACT_P1_ROUNDING_DIFFERENCE=boundaryDifference,TIGHT_FINAL_WITH_ORIGINAL_SPA=tightFinalOriginal,TIGHT_FINAL_WITH_NATIVE_SPA=tightFinalNative,ABS_TIGHT_WALD_DELTA=abs(canonical$TIGHT_LOG10P-nn$LOG10P_Wald_GxE),ABS_SAME_LEG_CCT_DELTA=abs(tightFinalNative-nn$LOG10P_GxE),as.data.frame(canonical),check.names=FALSE)
       if(canonical$INTERACTION_IDENTIFIABLE&&canonical$CONVERGED){
        check(key,"WALD_IDENTIFIABLE_TIGHT_FIT",canonical$TIGHT_LOG10P,nn$LOG10P_Wald_GxE)
        check(key,"FINAL_IDENTIFIABLE_SAME_SPA_LEG",tightFinalNative,nn$LOG10P_GxE)
        if(!boundaryDifference)check(key,"FINAL_IDENTIFIABLE_ORIGINAL_SPA_LEG",tightFinalOriginal,nn$LOG10P_GxE)
       }
      }
     }
    }
   }
  }
 }
}
formats<-sub(".localbgen.tsv$","",list.files(base,pattern="^bits[0-9]+_(raw|zlib|zstd)\\.localbgen\\.tsv$"))
stopifnot(length(formats)==15)
for(stem in formats)for(cutoff in c(.001,1))run(stem,"residual",cutoff)
for(cutoff in c(.001,1))run("bits16_zlib","tail_residual",cutoff)
for(cutoff in c(.001,1))run("bits16_zlib","native",cutoff)
reference<-do.call(rbind,rows);checkTable<-do.call(rbind,checks)
write_partial_reference_audit(file.path(out,"partial_reference_audit.tsv"))
write.table(reference,file.path(out,"localbgen_source_reference.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write.table(do.call(rbind,excluded),file.path(out,"original_R_excluded_domains.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write.table(checkTable,file.path(out,"localbgen_R_checks.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
wd<-do.call(rbind,waldDiagnostics)
write.table(wd,file.path(out,"localbgen_Wald_diagnostics.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write.table(wd[wd$EXACT_P1_ROUNDING_DIFFERENCE,],file.path(out,"localbgen_exact_P1_rounding.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
waldSummary<-data.frame(BRANCH_B_GROUPS=nrow(wd),IDENTIFIABLE_INTERACTIONS=sum(wd$INTERACTION_IDENTIFIABLE),REDUNDANT_NUISANCE_GROUPS=sum(wd$RANK<wd$COLUMNS),CONVERGED_CANONICAL_FITS=sum(wd$CONVERGED),MAX_CANONICAL_WALD_LOG10P_DELTA=max(wd$ABS_TIGHT_WALD_DELTA,na.rm=TRUE),MAX_SAME_SPA_LEG_CCT_LOG10P_DELTA=max(wd$ABS_SAME_LEG_CCT_DELTA,na.rm=TRUE),ORIGINAL_WALD_DIFFERENCES_OVER_2e_5=sum(!is.finite(wd$ABS_WALD_DELTA)|wd$ABS_WALD_DELTA>2e-5),ORIGINAL_FINAL_DIFFERENCES_OVER_2e_5=sum(!is.finite(wd$ABS_FINAL_DELTA)|wd$ABS_FINAL_DELTA>2e-5),ORIGINAL_EXACT_P1_ROUNDING_GROUPS=sum(wd$EXACT_P1_ROUNDING_DIFFERENCE),MAX_ORIGINAL_WALD_LOG10P_DELTA=max(wd$ABS_WALD_DELTA,na.rm=TRUE),MAX_ORIGINAL_FINAL_LOG10P_DELTA=max(wd$ABS_FINAL_DELTA,na.rm=TRUE))
write.table(waldSummary,file.path(out,"localbgen_Wald_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
summaries<-list();kk<-0L
for(input in unique(reference$INPUT))for(rootMode in names(original)){
 r<-reference[reference$INPUT==input&reference$ROOT==rootMode,]
 kk<-kk+1L;summaries[[kk]]<-data.frame(INPUT=input,ROOT=rootMode,SOURCE_CALLS=nrow(r),SOURCE_ERRORS=sum(r$ORIGINAL_STATUS!="OK"),N_BEYOND_NORMAL_CUTOFF=sum(abs(r$NATIVE_Z_NORM)>2),MAX_ABS_SPA_LOG10P_DELTA=max(r$ABS_SPA_DELTA,na.rm=TRUE))
}
write.table(do.call(rbind,summaries),file.path(out,"localbgen_R_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("Actual preserved original-R local function; source files unchanged.","Original-R SPA/final differences are descriptive; the independent guarded partial policy is the current SPA acceptance reference.","default: original Newton tolerance; tight: only root tolerance overridden in memory to 1e-12.","Paired text G/H independently computed from stored quantized P; native input is BGEN.","15 bit-depth/codec formats, residual SPA A/B, strong-tail residuals, and representative native four-phenotype fits.","Default cutoff .001 and forced Branch B cutoff 1 are separate test conditions.","High AF, missing G/H, zero-trial and degenerate-G domains are explicitly excluded from original-R equivalence.","All identifiable native Wald targets are independently checked after reducing nuisance rank; all their CCT results are checked using the same native SPA leg.","Original default-Wald and exact-P1/near-P1 final-CCT differences remain in separate diagnostics; passing does not claim universal original-final-P equivalence.",paste("Actual original-source calls:",nrow(reference)),paste("Tight native/reference checks:",nrow(checkTable)),paste("Failed checks:",sum(checkTable$REQUIRED&!checkTable$PASS)),paste("Original final differences over 2e-5:",waldSummary$ORIGINAL_FINAL_DIFFERENCES_OVER_2e_5),paste("Original exact-P1/near-P1 rounding groups:",waldSummary$ORIGINAL_EXACT_P1_ROUNDING_GROUPS)),file.path(out,"README.txt"))
cat(nrow(reference)," actual original-R calls; ",nrow(checkTable)," checks; failures: ",sum(checkTable$REQUIRED&!checkTable$PASS),"\n",sep="")
if(any(checkTable$REQUIRED&!checkTable$PASS))stop("Local-BGEN/original-R comparison failed; inspect localbgen_R_checks.tsv")
