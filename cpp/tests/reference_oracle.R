# Independent SPAGxE numerical oracle. Uses only R model fitting and distributions.
args<-commandArgs(trailingOnly=TRUE); od<-if(length(args))args[1] else "examples/synthetic"
dir.create(od,recursive=TRUE,showWarnings=FALSE); set.seed(20261024)
n<-800L;m<-48L;IID<-sprintf("S%04d",seq_len(n))
PC1<-as.numeric(scale(rnorm(n)));PC2<-as.numeric(scale(rnorm(n)));X1<-rnorm(n);E<-rnorm(n)+.4*PC1
G<-sapply(seq_len(m),function(j)rbinom(n,2,plogis(qlogis(.02+.32*((j-1)%%12)/11)+if(j%%3==0)1.1*PC1 else 0)))
G[,2]<-rbinom(n,2,.3);G[,47]<-2-G[,3];G[,48]<-0L
Quantitative<-.3*X1+.65*E+2.4*G[,2]+.7*G[,4]*E+rt(n,5)
Binary<-rbinom(n,1,plogis(-2+.25*X1+.3*E+2*G[,2]+.4*G[,4]*E))
et<-rexp(n,exp(.2*X1+.3*E+.8*G[,2]+.25*G[,4]*E));ct<-rexp(n,.65)
Time<-pmin(et,ct);Event<-as.integer(et<=ct)
Ordinal<-as.integer(cut(.3*X1+.5*E+1.1*G[,2]+.5*G[,4]*E+rlogis(n),c(-Inf,-1,0,1.2,Inf),labels=FALSE))
G[seq(3,n,47),11]<-NA_integer_
R<-resid(lm(Quantitative~X1+E+PC1+PC2))
Rb<-resid(glm(Binary~X1+E+PC1+PC2,family=binomial()),type="response")
Rs<-resid(survival::coxph(survival::Surv(Time,Event)~X1+E+PC1+PC2,ties="breslow"),type="martingale")
ofit<-ordinal::clm(as.factor(Ordinal)~X1+E+PC1+PC2)
oeta<-as.numeric(cbind(X1,E,PC1,PC2)%*%ofit$beta)
nu<-cbind(0,plogis(outer(oeta,ofit$alpha,function(a,b)b-a)),1)
Ro<-nu[cbind(seq_len(n),Ordinal)]+nu[cbind(seq_len(n),Ordinal+1L)]-1
ph<-data.frame(IID,X1,E,Eshift=E+3,PC1,PC2,Quantitative,Binary,Time,Event,Ordinal,R_quant=R,R_binary=Rb,R_survival=Rs,R_offset=R+.35,R_ordinal=Ro)
wr<-function(x,f,header=TRUE)write.table(x,file.path(od,f),sep="\t",quote=FALSE,row.names=FALSE,col.names=header)
wr(ph,"synthetic.pheno");wr(data.frame(IID,IID,0,0,0,-9),"synthetic.fam",FALSE)
snps<-sprintf("rsOracle%02d",seq_len(m));wr(data.frame(1,snps,0,seq_len(m)*100,"A","G"),"synthetic.bim",FALSE)
# BED SNP-major A1 dosage: 00=2, 10=1, 11=0, 01=missing.
bed<-file(file.path(od,"synthetic.bed"),"wb");writeBin(as.raw(c(0x6c,0x1b,0x01)),bed)
for(j in seq_len(m)){co<-ifelse(is.na(G[,j]),1L,c(3L,2L,0L)[G[,j]+1L]);writeBin(as.raw(colSums(matrix(c(co,rep(0L,(-n)%%4L)),4)*c(1,4,16,64))),bed)}
close(bed);wr(data.frame(IID,G),"synthetic.genotypes.tsv")
gr<-rbind(data.frame(ID1=IID,ID2=IID,Value=1),data.frame(ID1=IID[seq(1,n,2)],ID2=IID[seq(2,n,2)],Value=.125))
wr(gr,"synthetic.grm.tsv");phi<-diag(n);phi[cbind(seq(1,n,2),seq(2,n,2))]<-.125;phi[cbind(seq(2,n,2),seq(1,n,2))]<-.125
logadd<-function(a,b){z<-pmax(a,b);z+log1p(exp(-abs(a-b)))}
logsum<-function(x){if(all(x==-Inf))return(-Inf);z<-max(x);z+log(sum(exp(x-z)))}
source("tests/partial_reference.R")
# Full reference and guarded partial policy use independently solved scalar R CGFs.
spa<-function(g,w,q,grm=NULL,cutoff=2,hybrid=FALSE){
 q<-rep(q,length.out=length(w));mu<-sum(2*q*w);v0<-sum(2*q*(1-q)*w*w)
 variance<-if(is.null(grm))v0 else as.numeric(2*q[1]*(1-q[1])*crossprod(w,grm%*%w))
 dev<-sum(g*w)-mu;z<-dev/sqrt(variance);ln<-pnorm(-abs(z),log.p=TRUE)+log(2)
 if(!is.finite(z)||variance<=0)return(c(log10p=NA,z=z,variance=variance))
 if(abs(z)<=cutoff)return(c(log10p=-ln/log(10),z=z,variance=variance))
 d<-abs(dev)*sqrt(v0/variance)
 value<-reference_guarded_spa(w,2,q,mu,d,v0,partial=hybrid)
 c(log10p=value,z=z,variance=variance)
}
cct<-function(p){if(any(is.na(p)))return(NA_real_);if(any(p==0))return(0);p[p==1]<-.999;tt<-mean(ifelse(p<1e-16,1/(pi*p),tan(pi*(.5-p))));if(tt>1e15)1/(pi*tt) else pcauchy(tt,lower.tail=FALSE)}
individualAF<-function(g,pcs){
 q<-mean(g)/2;nn<-length(g);if(2*nn*q<=20)return(rep(q,nn));fit<-lm(g~pcs);af<-fitted(fit)/2
 if(mean(af<0|af>1)<.1)return(pmax(0,pmin(1,af)))
 ps<-coef(summary(fit))[-1,4];selected<-which(ps<.05);if(!length(selected))return(rep(q,nn));carrier<-as.numeric(g>.5)
 if(sum(carrier)==0||sum(carrier)==nn)return(rep(q,nn));lf<-suppressWarnings(glm(carrier~pcs[,selected,drop=FALSE],family=binomial()));eta<-predict(lf,type="link");sp<-pmax(eta,0)+log1p(exp(-abs(eta)));af<--expm1(-sp/2)
 if(!isTRUE(lf$converged)||any(!is.finite(af))||all(af==0)||all(af==1))return(rep(q,nn));af
}
wald<-function(g,trait){
 dat<-ph;dat$g<-g
 if(trait=="quantitative")return(coef(summary(lm(Quantitative~X1+E+PC1+PC2+g+g:E,dat)))["E:g","Pr(>|t|)"])
 if(trait=="binary")return(coef(summary(glm(Binary~X1+E+PC1+PC2+g+g:E,dat,family=binomial())))["E:g","Pr(>|z|)"])
 if(trait=="survival")return(coef(summary(survival::coxph(survival::Surv(Time,Event)~X1+E+PC1+PC2+g+g:E,dat,ties="breslow")))["E:g","Pr(>|z|)"])
 if(trait=="ordinal")return(coef(summary(ordinal::clm(as.factor(Ordinal)~X1+E+PC1+PC2+g+g:E,data=dat)))["E:g","Pr(>|z|)"])
}
testone<-function(g,R,E,mode="base",trait="quantitative",cutoff=2){
 af<-mean(g,na.rm=TRUE)/2;g[is.na(g)]<-2*af;if(af>.5){g<-2-g;af<-1-af};if(af==0)return(NULL)
 q<-if(mode=="mix")individualAF(g,cbind(PC1,PC2)) else rep(af,n);vr<-2*q*(1-q);sg<-sum(g*R);mug<-sum(2*q*R)
 vg<-if(mode=="plus")vr[1]*as.numeric(crossprod(R,phi%*%R)) else sum(vr*R^2);zg<-(sg-mug)/sqrt(vg);pg<-2*pnorm(-abs(zg));branch<-if(pg<=.001)"B" else "A"
 if(branch=="A"){la<-if(mode=="plus")as.numeric(crossprod(R,phi%*%(E*R))/crossprod(R,phi%*%R)) else sum(vr*E*R^2)/sum(vr*R^2);w<-(E-la)*R}else{rt<-resid(lm(R~g));w<-E*rt}
 a<-spa(g,w,q,if(mode=="plus")phi else NULL,cutoff,FALSE);h<-spa(g,w,q,if(mode=="plus")phi else NULL,cutoff,TRUE)
 pw<-if(branch=="B"&&mode!="plus")wald(g,trait) else NA_real_
 final<-if(is.finite(pw))-log10(cct(c(10^(-h[1]),pw))) else h[1];efinal<-if(is.finite(pw))-log10(cct(c(10^(-a[1]),pw))) else a[1]
 c(branch=branch,af=af,zG=zg,log10pG=-log10(pg),zNorm=h[2],variance=h[3],log10pSPA=h[1],log10pExactSPA=a[1],log10pWald=-log10(pw),log10pFinal=final,log10pExactFinal=efinal)
}
rows<-list();ri<-0L
for(mode in c("base","mix","plus"))for(trait in c("quantitative","binary","survival","ordinal")){
 rv<-switch(trait,quantitative=R,binary=Rb,survival=Rs,ordinal=Ro)
 for(j in seq_len(m-1L)){r<-testone(G[,j],rv,E,mode,trait);if(!is.null(r)){ri<-ri+1L;rows[[ri]]<-data.frame(ID=snps[j],MODE=mode,TRAIT=trait,t(unname(r)),check.names=FALSE)}}
}
reference<-do.call(rbind,rows);names(reference)<-c("ID","MODE","TRAIT","BRANCH","AF","Z_G","LOG10P_G","Z_Norm","VAR","LOG10P_SPA","LOG10P_EXACT_SPA","LOG10P_WALD","LOG10P_FINAL","LOG10P_EXACT_FINAL")
wr(reference,"oracle_reference.tsv")
write_partial_reference_audit(file.path(od,"partial_reference_audit.tsv"))
writeLines(c("Deterministic fixture seed 20261024.","800 subjects and 48 variants.","Marker 2: Branch B; marker 11: missing; marker 47: allele flip of marker 3; marker 48: monomorphic.","Independent reference includes base, mix, GRM, normal/SPA/CCT, full and guarded partial-normal CGF.","Cox uses Breslow ties. R is required for validation only, never for SPAGxE runtime."),file.path(od,"README.txt"))
cat(sprintf("Wrote %d independent reference rows to %s\n",nrow(reference),od))

