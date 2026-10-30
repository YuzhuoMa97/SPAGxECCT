# Original-R and native AF-cascade counterexamples using valid hardcalls.
# Run from cpp: Rscript tests/verify_af_cascade.R [binary] [output-directory]
# The executable has one published-method implementation; historical GRAB
# AF/CGF calculations below are independent comparisons, not product modes.
# No original source file is modified. R is a validation dependency only.
options(digits=17)
args<-commandArgs(trailingOnly=TRUE)
binary<-if(length(args))args[1]else if(.Platform$OS.type=="windows")"build/bin/SPAGxE.exe"else "build/bin/SPAGxE"
out<-if(length(args)>1)args[2]else "results/af_cascade"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# Load only scalar-CGF helper definitions, never the other script's workflow.
for(expr in parse("tests/reference_original_R.R"))if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&is.symbol(expr[[2]])&&as.character(expr[[2]])%in%c("logadd","logsum","partition","spa"))eval(expr)
legacy<-new.env(parent=globalenv());sys.source("../R/SPAGxECCT.R",envir=legacy)
originalRoot<-legacy$fastgetroot_H1_new
legacy$fastgetroot_H1_new<-function(init.t,R,s,MAFVec,tol=1e-12,maxiter=100)originalRoot(init.t,R,s,MAFVec,tol=1e-12,maxiter=maxiter)
cases<-list(
 upper_overshoot=list(x=c(rep(3,60),rep(1,40),rep(0,300)),g=c(rep(2,58),rep(0,2),rep(2,38),rep(0,2),rep(1,10),rep(0,290))),
 negative_boundary=list(x=c(rep(-1,40),rep(0,320),rep(1,40)),g=c(rep(1,3),rep(0,37),rep(1,10),rep(0,310),rep(2,35),rep(0,5))),
 low_carrier_count=list(x=seq(-1,1,length.out=400),g=rep(0,400)))
cases$low_carrier_count$g[c(1,seq(373,399,2))]<-2
results<-checks<-list();k<-kc<-0L
for(name in names(cases)){
 g<-cases[[name]]$g;PC<-as.numeric(scale(cases[[name]]$x));nn<-length(g);E<-sin(seq_len(nn)*.31)+.2*cos(seq_len(nn)*.11)
 R<-resid(lm((2*g*E+sin(seq_len(nn)*.7))~PC+E));q<-mean(g)/2
 fit<-lm(g~PC);ols<-fitted(fit)/2;negative<-mean(ols<0);outside<-mean(ols<0|ols>1)
 significant<-coef(summary(fit))["PC","Pr(>|t|)"]<.05;carrier<-as.numeric(round(g)!=0)
 logistic<-glm(carrier~PC,family=binomial(),control=glm.control(epsilon=1e-12,maxit=100));lq<-1-sqrt(1-fitted(logistic))
 stopifnot(q<=.5,q>0,all(g%in%0:2),logistic$converged,significant)
 qr<-if(negative>.1)lq else pmax(0,pmin(1,ols))
 qp<-if(outside>=.1)lq else pmax(0,pmin(1,ols))
 qg<-if(outside>=.1&&sum(carrier)<=20)rep(q,nn)else qp
 references<-list(original_R=qr,paper=qp,grab2=qg)
 record<-function(label,value,qi=NULL){k<<-k+1L;results[[k]]<<-data.frame(CASE=name,AF_STRATEGY=label,N=nn,POOLED_AF=q,MAC=sum(g),CARRIERS=sum(carrier),NEGATIVE_OLS_RATIO=negative,UPPER_OLS_RATIO=mean(ols>1),OUTSIDE_OLS_RATIO=outside,OLS_MIN=min(ols),OLS_MAX=max(ols),LOGISTIC_CONVERGED=logistic$converged,AF_MIN=if(is.null(qi))NA_real_ else min(qi),AF_MAX=if(is.null(qi))NA_real_ else max(qi),LOG10P_SPA=value,check.names=FALSE)}
 expected<-list()
 for(flavour in names(references)){
  qi<-references[[flavour]];v<-2*qi*(1-qi);lambda<-sum(v*E*R^2)/sum(v*R^2)
  full<-unname(spa(g,(E-lambda)*R,qi,NULL,FALSE)["LOG10P"])
  record(paste0("independent_",flavour,"_full"),full,qi)
  expected[[flavour]]<-if(flavour=="grab2")unname(spa(g,(E-lambda)*R,qi,NULL,TRUE)["LOG10P"])else full
  if(flavour=="grab2")record("independent_grab2_default_hybrid",expected[[flavour]],qi)
 }
 ans<-legacy$SPAGxEmix_CCT_one_SNP(traits="quantitative",g=g,R=R,E=E,Phen.mtx=data.frame(Y=R),Cova.mtx=as.matrix(PC),topPCs=as.matrix(PC),epsilon=0)
 record("actual_original_R_tight",-log10(ans[3]))
 id<-sprintf("af%03d",seq_len(nn));prefix<-file.path(out,name)
 write.table(data.frame(IID=id,R=R,E=E,PC=PC,G=g),paste0(prefix,".pheno"),sep="\t",quote=FALSE,row.names=FALSE)
 write.table(data.frame(id,id,0,0,0,-9),paste0(prefix,".fam"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
 write.table(data.frame(1,name,0,1,"A","G"),paste0(prefix,".bim"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
 bed<-file(paste0(prefix,".bed"),"wb");writeBin(as.raw(c(0x6c,0x1b,1)),bed);co<-c(3L,2L,0L)[g+1L];writeBin(as.raw(colSums(matrix(co,4)*c(1,4,16,64))),bed);close(bed)
 for(profile in "native"){
  output<-paste0(prefix,".",profile)
  command<-c("--method","SPAGxEmix","--bfile",prefix,"--pheno",paste0(prefix,".pheno"),"--resid-name","R","--envir-name","E","--pc-cols","PC","--spagxe-marginal-cutoff","0","--maf","0","--mac","0","--out",output)
  log<-paste0(output,".console.txt");status<-system2(binary,shQuote(command),stdout=log,stderr=log)
  native<-if(status==0)read.delim(paste0(output,".R.SPAGxEmix"),check.names=FALSE)$LOG10P_SPA_GxE else NA_real_
  record("actual_SPAGxE",native)
  kc<-kc+1L;checks[[kc]]<-data.frame(CASE=name,TARGET="SPAGxE",EXPECTED_LOG10P=expected$paper,NATIVE_LOG10P=native,ABS_DELTA=abs(native-expected$paper),PASS=status==0&&is.finite(native)&&abs(native-expected$paper)<2e-5)
 }
 kc<-kc+1L;checks[[kc]]<-data.frame(CASE=name,TARGET="original_R_tight",EXPECTED_LOG10P=expected$original_R,NATIVE_LOG10P=-log10(ans[3]),ABS_DELTA=abs(-log10(ans[3])-expected$original_R),PASS=is.finite(ans[3])&&abs(-log10(ans[3])-expected$original_R)<2e-5)
}
write.table(do.call(rbind,results),file.path(out,"af_counterexamples.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
checks<-do.call(rbind,checks);write.table(checks,file.path(out,"af_regression_checks.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
cat(nrow(checks)," AF checks; failures: ",sum(!checks$PASS),"; maximum -log10(P) delta: ",max(checks$ABS_DELTA),"\n",sep="")
if(any(!checks$PASS))stop("AF-cascade reference comparison failed; see output tables")
