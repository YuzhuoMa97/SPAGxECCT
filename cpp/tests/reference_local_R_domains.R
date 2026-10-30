# Execute actual preserved local functions, including their failure domains.
# Run from cpp: Rscript tests/reference_local_R_domains.R [output-directory]
# Only this package's synthetic two/three-ancestry input is read.
options(digits=17)
suppressPackageStartupMessages(library(survival));suppressPackageStartupMessages(library(ordinal))
args<-commandArgs(trailingOnly=TRUE)
root<-"examples/local";out<-if(length(args))args[1]else "results/local_original_R_domains"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ph<-read.delim(file.path(root,"samples.pheno"),check.names=FALSE)
rs<-list(quantitative=resid(lm(Quantitative~X+E,ph)),binary=resid(glm(Binary~X+E,ph,family=binomial()),type="response"),survival=resid(coxph(Surv(Time,Event)~X+E,ph,ties="breslow"),type="martingale"))
ord<-clm(as.factor(Ordinal)~X+E,data=ph);eta<-as.numeric(as.matrix(ph[,c("X","E")])%*%ord$beta);nu<-cbind(0,plogis(outer(eta,ord$alpha,function(a,b)b-a)),1);yy<-as.integer(factor(ph$Ordinal));rs$ordinal<-nu[cbind(seq_len(nrow(ph)),yy)]+nu[cbind(seq_len(nrow(ph)),yy+1)]-1
loadLegacy<-function(tight){e<-new.env(parent=globalenv());sys.source("../R/SPAGxECCT.R",envir=e)
 if(tight){r<-e$fastgetroot_H1_new_local;e$fastgetroot_H1_new_local<-function(init.t,R,s,MAFVec,haplo_numVec,tol=1e-12,maxiter=100)r(init.t,R,s,MAFVec,haplo_numVec,tol=1e-12,maxiter=maxiter)}
 e}
legacy<-list(default=loadLegacy(FALSE),tight=loadLegacy(TRUE))
logsum<-function(x){if(all(x==-Inf))return(-Inf);z<-max(x);z+log(sum(exp(x-z)))}
logadd<-function(a,b){z<-pmax(a,b);x<-z+log1p(exp(-abs(a-b)));x[is.infinite(z)&z<0]<--Inf;x}
independentSPA<-function(g,h,w,q){mu<-q*sum(h*w);v<-q*(1-q)*sum(h*w*w);s<-sum(g*w);z<-(s-mu)/sqrt(v)
 if(!is.finite(z)||v<=0)return(c(LOG10P=NA_real_,Z=NA_real_,VAR=v))
 if(abs(z)<=2)return(c(LOG10P=-(pnorm(-abs(z),log.p=TRUE)+log(2))/log(10),Z=z,VAR=v))
 use<-h>0;hi<-h[use];wi<-w[use]
 cgf<-function(t){den<-logadd(rep(log1p(-q),length(wi)),log(q)+t*wi);pii<-plogis(qlogis(q)+t*wi);c(K=sum(hi*den),K1=sum(hi*wi*pii),K2=sum(hi*wi*wi*pii*(1-pii)))}
 one<-function(target,upper){f<-function(t)cgf(t)["K1"]-target;b<-1/sqrt(v);it<-0
  if(upper){lo<-0;up<-b;while(f(up)<0&&it<100){up<-up*2;it<-it+1}}else{up<-0;lo<--b;while(f(lo)>0&&it<100){lo<-lo*2;it<-it+1}}
  if(f(lo)*f(up)>0)return(NA_real_);rt<-uniroot(f,c(lo,up),tol=1e-12,maxiter=300)$root;cv<-cgf(rt);r<-sign(rt)*sqrt(max(0,2*(rt*target-cv["K"])));u<-rt*sqrt(cv["K2"]);a<-r+log(u/r)/r;pnorm(a,lower.tail=!upper,log.p=TRUE)}
 lp<-c(one(mu+abs(s-mu),TRUE),one(mu-abs(s-mu),FALSE));if(anyNA(lp))lp<-pnorm(-abs(z),log.p=TRUE)+log(2)else lp<-logsum(lp)
 c(LOG10P=-min(lp,0)/log(10),Z=z,VAR=v)}
independentWald<-function(dat,g,allH,trait){design<-cbind(1,dat$X,dat$E);rk<-qr(design,tol=1e-10)$rank
 for(a in seq_len(ncol(allH))){candidate<-cbind(design,allH[,a]);kk<-qr(candidate,tol=1e-10)$rank;if(kk>rk){design<-candidate;rk<-kk}}
 d<-as.data.frame(design[,-1,drop=FALSE]);names(d)<-paste0("C",seq_len(ncol(d)));d$g<-g;d$ge<-g*dat$E;rhs<-paste(c(names(d)[grepl("^C",names(d))],"g","ge"),collapse="+")
 tryCatch({if(trait=="quantitative"){d$y<-dat$Quantitative;p<-coef(summary(lm(as.formula(paste("y~",rhs)),d)))["ge","Pr(>|t|)"]}else if(trait=="binary"){d$y<-dat$Binary;p<-coef(summary(glm(as.formula(paste("y~",rhs)),data=d,family=binomial())))["ge","Pr(>|z|)"]}else if(trait=="survival"){d$time<-dat$Time;d$event<-dat$Event;p<-coef(summary(coxph(as.formula(paste("Surv(time,event)~",rhs)),data=d,ties="breslow")))["ge","Pr(>|z|)"]}else{d$y<-factor(dat$Ordinal);p<-coef(summary(clm(as.formula(paste("y~",rhs)),data=d)))["ge","Pr(>|z|)"]};-log10(p)},error=function(e)NA_real_)}
cct<-function(lp){lp<-lp[!is.na(lp)];if(!length(lp))return(NA_real_);if(any(lp==Inf))return(Inf);lp[lp==0]<--log10(.999);p<-10^-lp;t<-mean(ifelse(p<1e-16,1/(pi*p),tan(pi*(.5-p))));if(t>1e15)log10(pi)+log10(t)else -log10(pcauchy(t,lower.tail=FALSE))}
readSynthetic<-function(dir){mf<-read.delim(file.path(dir,"cohort.local.tsv"),check.names=FALSE);G<-H<-list();for(a in seq_len(nrow(mf))){g<-read.delim(file.path(dir,mf$genotype[a]),check.names=FALSE);h<-read.delim(file.path(dir,mf$hapcount[a]),check.names=FALSE);G[[a]]<-as.matrix(g[,-(1:5)]);H[[a]]<-as.matrix(h[,-(1:5)])};list(labels=mf$ancestry,meta=g[,1:5],G=G,H=H)}
rows<-list();k<-0L
for(dataset in c("two","three")){
 src<-readSynthetic(if(dataset=="two")root else file.path(root,"three_ancestry"))
 for(trait in names(rs))for(j in seq_len(nrow(src$meta)))for(a in seq_along(src$labels)){
  rawG<-src$G[[a]][j,];rawH<-src$H[[a]][j,];allH<-sapply(src$H,function(x)x[j,]);af<-sum(rawG)/sum(rawH)
  domain<-if(anyNA(allH))"unknown_H"else if(anyNA(rawG))if(all(is.na(rawG)<= (rawH==0)))"structural_zero_missing_G"else "missing_G"else if(!is.finite(af)||sum(rawH)==0)"zero_trial"else if(af>.5)"high_AF_variable_H"else "safe_complete_minor"
  keep<-which(complete.cases(allH));hh<-allH[keep,,drop=FALSE];h<-rawH[keep];g<-rawG[keep];R<-rs[[trait]][keep];dat<-ph[keep,,drop=FALSE]
  g[is.na(g)&h==0]<-0;obs<-is.finite(g);qAlt<-if(sum(h[obs])>0)sum(g[obs])/sum(h[obs])else NA_real_;q<-min(qAlt,1-qAlt);g[!obs]<-h[!obs]*qAlt;if(is.finite(qAlt)&&qAlt>.5)g<-h-g
  valid<-is.finite(q)&&q>0&&q<1&&sum(h*R^2)>0;zg<-if(valid)sum((g-q*h)*R)/sqrt(q*(1-q)*sum(h*R^2))else NA_real_;branch<-if(valid&&2*pnorm(-abs(zg))<=.001)"B"else "A"
  w<-if(valid){if(branch=="B")dat$E*resid(lm(R~g))else(dat$E-sum(h*dat$E*R^2)/sum(h*R^2))*R}else rep(0,length(R))
  calc<-if(valid)independentSPA(g,h,w,q)else c(LOG10P=NA_real_,Z=NA_real_,VAR=NA_real_);pw<-if(valid&&branch=="B")independentWald(dat,g,hh,trait)else NA_real_;final<-if(is.finite(pw))cct(c(calc["LOG10P"],pw))else calc["LOG10P"]
  for(rootMode in names(legacy)){
   phen<-data.frame(Y=if(trait=="binary")ph$Binary else if(trait=="ordinal")ph$Ordinal else ph$Quantitative,surv.time=ph$Time,event=ph$Event)
   ans<-tryCatch(suppressWarnings(legacy[[rootMode]]$SPAGxEmixCCT_localance_one_SNP(traits=if(trait=="ordinal")"categorical"else trait,g=rawG,haplo_numVec=rawH,R=rs[[trait]],E=ph$E,Phen.mtx=phen,Cova.mtx=as.matrix(ph["X"]),Cova.haplo.mtx=allH[,seq_len(ncol(allH)-1),drop=FALSE],min.maf=1e-5)),error=function(e)e)
   ov<-if(inherits(ans,"error"))rep(NA_real_,10)else ans
   k<-k+1L;rows[[k]]<-data.frame(DATASET=dataset,TRAIT=trait,ID=src$meta$ID[j],ANCESTRY=src$labels[a],DOMAIN=domain,ROOT=rootMode,ORIGINAL_STATUS=if(inherits(ans,"error"))conditionMessage(ans)else if(is.na(ov[3]))"RETURNED_NA"else "OK",ORIGINAL_AF=ov[1],CORRECTED_ALT_AF=qAlt,N_ORIGINAL=length(rawG),N_CORRECTED=length(R),N_H_REMOVED=length(rawG)-length(R),ORIGINAL_LOG10P_SPA=-log10(ov[3]),CORRECTED_LOG10P_SPA=calc["LOG10P"],ORIGINAL_LOG10P_FINAL=-log10(ov[5]),CORRECTED_LOG10P_FINAL=final,ORIGINAL_LOG10P_WALD=-log10(ov[4]),CORRECTED_LOG10P_WALD=pw,ABS_SPA_DELTA=abs(-log10(ov[3])-calc["LOG10P"]),check.names=FALSE)
  }
 }
}
results<-do.call(rbind,rows)
write.table(results,file.path(out,"original_local_domains.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
safe<-results[results$DOMAIN=="safe_complete_minor",];summary<-list();kk<-0L
for(rootMode in names(legacy))for(domain in unique(results$DOMAIN)){
 r<-results[results$ROOT==rootMode&results$DOMAIN==domain,];finite<-is.finite(r$ABS_SPA_DELTA);fd<-abs(r$ORIGINAL_LOG10P_FINAL-r$CORRECTED_LOG10P_FINAL);wd<-abs(r$ORIGINAL_LOG10P_WALD-r$CORRECTED_LOG10P_WALD);kk<-kk+1L;summary[[kk]]<-data.frame(ROOT=rootMode,DOMAIN=domain,ROWS=nrow(r),ORIGINAL_ERRORS=sum(!r$ORIGINAL_STATUS%in%c("OK","RETURNED_NA")),ORIGINAL_NA=sum(r$ORIGINAL_STATUS=="RETURNED_NA"),FINITE_COMPARISONS=sum(finite),MAX_ABS_SPA_LOG10P_DELTA=if(any(finite))max(r$ABS_SPA_DELTA[finite])else NA_real_,MAX_ABS_FINAL_LOG10P_DELTA=if(any(is.finite(fd)))max(fd[is.finite(fd)])else NA_real_,BRANCH_B_WALD_COMPARISONS=sum(is.finite(wd)),MAX_ABS_WALD_LOG10P_DELTA=if(any(is.finite(wd)))max(wd[is.finite(wd)])else NA_real_,CORRECTED_VALID=sum(is.finite(r$CORRECTED_LOG10P_SPA)))}
write.table(do.call(rbind,summary),file.path(out,"domain_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
cat(paste("Actual unmodified original-local rows:",nrow(results),"; safe complete rows:",nrow(safe),".\n"))
