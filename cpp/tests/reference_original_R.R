# Historical full-CGF and unguarded IQR formula profiles are reference/delta
# models, not runtime selectors or the current guarded_partial_normal policy.
# Execute the preserved original R functions and compare full/accelerated CGFs.
# Run from cpp: Rscript tests/reference_original_R.R [output-directory]
# R is a validation dependency only. Preserved source is never overwritten.
options(digits=17)
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ordinal))
suppressPackageStartupMessages(library(dplyr))
args<-commandArgs(trailingOnly=TRUE)
base<-"examples/synthetic";out<-if(length(args))args[1]else "results/original_R_reference"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ph<-read.delim(file.path(base,"synthetic.pheno"),check.names=FALSE)
G<-as.matrix(read.delim(file.path(base,"synthetic.genotypes.tsv"),row.names=1))
colnames(G)<-read.table(file.path(base,"synthetic.bim"),stringsAsFactors=FALSE)[,2]
gr<-read.delim(file.path(base,"synthetic.grm.tsv"));n<-nrow(ph);phi<-diag(n)
ii<-match(gr$ID1,ph$IID);jj<-match(gr$ID2,ph$IID);phi[cbind(ii,jj)]<-gr$Value;phi[cbind(jj,ii)]<-gr$Value
wr<-function(x,f)write.table(x,file.path(out,f),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
wr(data.frame(IID=ph$IID,G,check.names=FALSE),"genotypes_named.tsv")
newLegacy<-function(refined=FALSE){
 e<-new.env(parent=globalenv());sys.source("../R/SPAGxECCT.R",envir=e);sys.source("../R/SPAGxEPlus_functions_MYZ.R",envir=e)
 if(refined){r<-e$fastgetroot_H1_new;e$fastgetroot_H1_new<-function(init.t,R,s,MAFVec,tol=1e-12,maxiter=100)r(init.t,R,s,MAFVec,tol=1e-12,maxiter=maxiter)}
 e
}
e0<-newLegacy(FALSE);e1<-newLegacy(TRUE)
residNames<-c(quantitative="R_quant",binary="R_binary",survival="R_survival",ordinal="R_ordinal")
traitPhen<-function(trait) data.frame(IID=ph$IID,Y=switch(trait,quantitative=ph$Quantitative,binary=ph$Binary,survival=ph$Quantitative,ordinal=ph$Ordinal),surv.time=ph$Time,event=ph$Event)
objFromR<-function(R,E){la<-as.numeric(crossprod(R,phi%*%(E*R))/crossprod(R,phi%*%R));rn<-(E-la)*R
 list(ResidMat=data.frame(SubjID=ph$IID,Resid=R,RE=E*R,R.new=rn),R=R,R_GRM_R=as.numeric(crossprod(R,phi%*%R)),R_GRM_RE=as.numeric(crossprod(R,phi%*%(E*R))),R.new=rn,R.new_GRM_R.new=as.numeric(crossprod(rn,phi%*%rn)),sparseGRM=gr)}
# Explicitly exercise the real original Plus_Nullmodel interface. Earlier
# Existing single-SNP checks construct this object from dense formulas.
interface<-list();ki<-0L
record<-function(name,pass,detail){ki<<-ki+1L;interface[[ki]]<<-data.frame(CHECK=name,PASS=pass,DETAIL=detail)}
callNull<-function(env,trait,data=ph,E=ph$E,gIDs=ph$IID){
 form<-switch(trait,quantitative=Quantitative~X1+E+PC1+PC2,binary=Binary~X1+E+PC1+PC2,survival=Surv(Time,Event)~X1+E+PC1+PC2)
 args<-list(traits=trait,formula=form,data=data,pIDs=data$IID,gIDs=gIDs,sparseGRM=gr,E=E)
 if(trait=="binary")args$family<-binomial();if(trait=="survival")args$ties<-"breslow"
 do.call(env$SPAGxE_Plus_Nullmodel,args)
}
clean<-newLegacy(TRUE)
err<-tryCatch({callNull(clean,"quantitative");"NO_ERROR"},error=function(e)conditionMessage(e))
record("original_nullmodel_without_global_Phen.mtx",grepl("Phen.mtx",err),err)
originalObjects<-list()
for(trait in c("quantitative","binary","survival")){
 env<-newLegacy(TRUE);env$Phen.mtx<-traitPhen(trait);o<-callNull(env,trait);originalObjects[[trait]]<-o
 dense<-objFromR(ph[[residNames[trait]]],ph$E)
 for(nm in c("R","R_GRM_R","R_GRM_RE","R.new","R.new_GRM_R.new")){
  error<-max(abs(o[[nm]]-dense[[nm]]));record(paste("original_real_nullmodel",trait,nm,sep=":"),is.finite(error)&&error<1e-6,paste("max_abs_error",error))
 }
 # Original check_input returns p2g but the original null function never uses
 # it to reorder R into the requested genotype order.
 perm<-rev(seq_len(n));env$Phen.mtx<-traitPhen(trait)[perm,,drop=FALSE]
 oo<-callNull(env,trait,ph[perm,,drop=FALSE],ph$E[perm],ph$IID)
 reorderError<-max(abs(oo$R-o$R));idReorderedError<-max(abs(oo$R[match(ph$IID,oo$ResidMat$SubjID)]-o$R))
 record(paste("original_p2g_alignment",trait,sep=":"),reorderError>1e-6&&idReorderedError<1e-6,paste("raw_vector_error",reorderError,"after_IID_reorder_error",idReorderedError))
}
wr(do.call(rbind,interface),"original_plus_nullmodel_interface.tsv")

# Original R p values, including default vs tight Newton tolerance, split into
# their actual SPA, Wald, normal, marginal and final-CCT fields.
rows<-list();k<-0L
for(rootMode in c("default","tight")){
 env<-if(rootMode=="default")e0 else e1
 for(scenario in c("native_residual","external_offset"))for(trait in if(scenario=="external_offset")"quantitative"else names(residNames))for(en in c("E","Eshift")){
  R<-if(scenario=="external_offset")ph$R_offset else ph[[residNames[trait]]];E<-ph[[en]];phen<-traitPhen(trait);obj<-objFromR(R,E)
  for(mode in c("base","plus","mix"))for(j in seq_len(ncol(G)-1)){
   ans<-tryCatch(suppressWarnings({
    if(mode=="base")env$SPAGxE_CCT_one_SNP(traits=if(trait=="ordinal")"categorical"else trait,g=G[,j],R=R,E=E,Phen.mtx=phen,Cova.mtx=as.matrix(ph[,c("X1","PC1","PC2")]))
    else if(mode=="mix")env$SPAGxEmix_CCT_one_SNP(traits=if(trait=="ordinal")"categorical"else trait,g=G[,j],R=R,E=E,Phen.mtx=phen,Cova.mtx=as.matrix(ph[,c("X1","PC1","PC2")]),topPCs=as.matrix(ph[,c("PC1","PC2")]))
    else env$SPAGxE_Plus_one_SNP(g=G[,j],E=E,Phen.mtx=phen,Cova.mtx=as.matrix(ph[,c("X1","PC1","PC2")]),obj.SPAGxE_Plus_Nullmodel=obj)
   }),error=function(e)e)
   values<-if(inherits(ans,"error"))rep(NA_real_,12)else ans
   if(mode=="plus"){spa<-values[3];pw<-NA_real_;fin<-spa;pn<-values[4];pg<-values[5];zg<-values[8]}else{spa<-values[3];pw<-values[4];fin<-values[5];pn<-values[6];pg<-values[7];zg<-values[10]}
   k<-k+1L;rows[[k]]<-data.frame(SCENARIO=scenario,ROOT=rootMode,TRAIT=trait,ENV=en,MODE=mode,ID=colnames(G)[j],AF_INPUT=mean(G[,j],na.rm=TRUE)/2,BRANCH=if(is.finite(pg)&&pg<=.001)"B"else "A",LOG10P_G=-log10(pg),Z_G=zg,LOG10P_SPA=-log10(spa),LOG10P_WALD=if(mode=="plus")NA_real_ else -log10(pw),LOG10P_NORMAL=-log10(pn),LOG10P_FINAL=-log10(fin),STATUS=if(inherits(ans,"error"))conditionMessage(ans)else "OK")
  }
 }
}
original<-do.call(rbind,rows);wr(original,"original_R_components.tsv")

logadd<-function(a,b){z<-pmax(a,b);x<-z+log1p(exp(-abs(a-b)));x[is.infinite(z)&z<0]<--Inf;x}
logsum<-function(x){if(all(x==-Inf))return(-Inf);z<-max(x);z+log(sum(exp(x-z)))}
partition<-function(w){qu<-quantile(w,c(.25,.75),names=FALSE,type=7);iq<-diff(qu);r<-1.5;pick<-function(z)which(w<qu[1]-z*iq|w>qu[2]+z*iq);sel<-pick(r);it<-0;if(iq>0)while(!length(sel)&&it<200){r<-r*.8;sel<-pick(r);it<-it+1};sel}
fitAF<-function(g,pcs){
 q<-mean(g)/2;nn<-length(g);if(q*2*nn<=20)return(rep(q,nn))
 fit<-lm(g~pcs);qi<-fitted(fit)/2
 bad<-mean(qi<0)>.1
 if(!bad)return(pmax(0,pmin(1,qi)))
 selected<-which(coef(summary(fit))[-1,4]<.05);if(!length(selected))return(rep(q,nn))
 carrier<-as.numeric(round(g)!=0)
 if(sum(round(g))==0||sum(2-round(g))==0)return(rep(q,nn))
 fit<-suppressWarnings(glm(carrier~pcs[,selected,drop=FALSE],family=binomial(),control=glm.control(maxit=100,epsilon=1e-12)))
 1-sqrt(1-fitted(fit))
}
spa<-function(g,w,qi,phi=NULL,hybrid=FALSE){
 mu<-sum(2*qi*w);v0<-sum(2*qi*(1-qi)*w*w)
 variance<-if(is.null(phi))v0 else as.numeric(2*qi[1]*(1-qi[1])*crossprod(w,phi%*%w))
 score<-sum(g*w);z<-(score-mu)/sqrt(variance)
 if(!is.finite(z)||variance<=0)return(c(LOG10P=NA_real_,Z=z,VAR=variance,MEAN=mu,SCORE=score))
 normal<-pnorm(-abs(z),log.p=TRUE)+log(2)
 if(abs(z)<=2)return(c(LOG10P=-normal/log(10),Z=z,VAR=variance,MEAN=mu,SCORE=score))
 io<-if(hybrid)partition(w)else seq_along(w);ic<-setdiff(seq_along(w),io)
 m0<-sum(2*qi[ic]*w[ic]);vN<-sum(2*qi[ic]*(1-qi[ic])*w[ic]^2);wi<-w[io];q<-qi[io]
 cgf<-function(t){den<-logadd(log1p(-q),log(q)+t*wi);pp<-plogis(qlogis(q)+t*wi);c(K=sum(2*den)+t*m0+.5*t*t*vN,K1=sum(2*wi*pp)+m0+t*vN,K2=sum(2*wi*wi*pp*(1-pp))+vN)}
 d<-abs(score-mu)*sqrt(v0/variance)
 one<-function(target,upper){f<-function(t)cgf(t)["K1"]-target;b<-1/sqrt(v0);it<-0
  if(upper){lo<-0;hi<-b;while(f(hi)<0&&it<100){hi<-hi*2;it<-it+1}}else{hi<-0;lo<--b;while(f(lo)>0&&it<100){lo<-lo*2;it<-it+1}}
  if(f(lo)*f(hi)>0)return(NA_real_)
  rt<-uniroot(f,c(lo,hi),tol=1e-12,maxiter=300)$root;cv<-cgf(rt);zz<-sign(rt)*sqrt(max(0,2*(rt*target-cv["K"])));vv<-rt*sqrt(cv["K2"]);a<-zz+log(vv/zz)/zz
  pnorm(a,lower.tail=!upper,log.p=TRUE)
 }
 pp<-c(one(mu+d,TRUE),one(mu-d,FALSE));lp<-if(anyNA(pp))normal else logsum(pp)
 c(LOG10P=-min(lp,0)/log(10),Z=z,VAR=variance,MEAN=mu,SCORE=score)
}
theory<-list();kt<-0L
for(scenario in c("native_residual","external_offset"))for(trait in if(scenario=="external_offset")"quantitative"else names(residNames))for(en in c("E","Eshift"))for(mode in c("base","plus","mix"))for(j in seq_len(ncol(G)-1)){
 R<-if(scenario=="external_offset")ph$R_offset else ph[[residNames[trait]]];E<-ph[[en]]
 rawg<-G[,j];altq<-mean(rawg,na.rm=TRUE)/2;rawg[is.na(rawg)]<-2*altq
 for(profile in c("paper_full","paper_fast")){
  g<-if(altq>.5)2-rawg else rawg
  q<-mean(g)/2;qi<-if(mode=="mix")fitAF(g,as.matrix(ph[,c("PC1","PC2")]))else rep(q,n)
  varianceG<-if(mode=="plus")as.numeric(2*q*(1-q)*crossprod(R,phi%*%R))else sum(2*qi*(1-qi)*R*R)
  meanG<-sum(2*qi*R)
  zg<-(sum(g*R)-meanG)/sqrt(varianceG);pg<-2*pnorm(-abs(zg));branch<-if(is.finite(pg)&&pg<=.001)"B"else "A"
  if(branch=="A"){
   lambda<-if(mode=="plus")as.numeric(crossprod(R,phi%*%(E*R))/crossprod(R,phi%*%R))else sum(2*qi*(1-qi)*E*R^2)/sum(2*qi*(1-qi)*R^2)
   w<-(E-lambda)*R;extra<-0
  }else{
   rt<-resid(lm(R~g));w<-E*rt;extra<-0
  }
  a<-spa(g,w,qi,if(mode=="plus")phi else NULL,profile!="paper_full")
  kt<-kt+1L;theory[[kt]]<-data.frame(SCENARIO=scenario,TRAIT=trait,ENV=en,MODE=mode,ID=colnames(G)[j],PROFILE=profile,AF_INPUT=altq,BRANCH=branch,LOG10P_G=-(pnorm(-abs(zg),log.p=TRUE)+log(2))/log(10),Z_G=zg,EXTRA_BRANCH_B_LAMBDA=extra,LOG10P_SPA=a["LOG10P"],Z_NORM=a["Z"],VAR=a["VAR"],MEAN=a["MEAN"],SCORE=a["SCORE"],check.names=FALSE)
 }
}
predicted<-do.call(rbind,theory);wr(predicted,"formula_profiles_components.tsv")

# Counterexample inventory: numerical differences are evidence of
# non-equivalence, not evidence by themselves that one source is correct.
keys<-c("SCENARIO","TRAIT","ENV","MODE","ID")
wide<-merge(predicted[predicted$PROFILE=="paper_full",],predicted[predicted$PROFILE=="paper_fast",],by=keys,suffixes=c("_full","_fast"))
wide$CGF_ABS_DELTA<-abs(wide$LOG10P_SPA_full-wide$LOG10P_SPA_fast)
wr(wide[order(-wide$CGF_ABS_DELTA),],"full_vs_fast_counterexamples.tsv")
tol<-merge(original[original$ROOT=="default",],original[original$ROOT=="tight",],by=keys,suffixes=c("_default","_tight"))
tol$TOL_ABS_DELTA<-abs(tol$LOG10P_SPA_default-tol$LOG10P_SPA_tight)
wr(tol[order(-tol$TOL_ABS_DELTA),],"legacy_root_counterexamples.tsv")
writeLines(c(paste("Original R component rows:",nrow(original)),paste("Formula profile rows:",nrow(predicted)),paste("Real Plus_Nullmodel/interface checks:",length(interface)),paste("Largest full vs fast log10p difference:",max(wide$CGF_ABS_DELTA,na.rm=TRUE)),"Original source files were not modified. Root=tight changes only Newton tolerance in memory.","The default-root and tight-root results are separate; failure-domain corrections must not be described as original-R equivalence."),file.path(out,"summary.txt"))
cat(paste("Wrote",nrow(original),"original-R component rows and",nrow(predicted),"formula-profile rows.\n"))
