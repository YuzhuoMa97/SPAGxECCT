# Independent R tan()/pcauchy() oracle for native weighted global/local CCT.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
x <- read.delim(args[1],check.names=FALSE)
one <- function(a,b,wg,wl) {
  if(is.na(a)||is.na(b)) return(NA_real_)
  p <- 10^(-c(a,b));p[p==1] <- .999
  w <- c(wg,wl);w <- w/sum(w)
  small <- p<1e-16
  terms <- numeric(2)
  terms[small] <- w[small]/p[small]/pi
  # cot(pi*p) equals tan(pi*(.5-p)) but avoids cancellation around pi/2
  # when p is small. Near one use 1-p from expm1 rather than subtraction.
  positive <- !small & p<=.5
  negative <- !small & p>.5
  terms[positive] <- w[positive]/tan(pi*p[positive])
  u <- -expm1(-c(a,b)[negative]*log(10));u[c(a,b)[negative]==0] <- .001
  terms[negative] <- -w[negative]/tan(pi*u)
  T <- sum(terms)
  # Stable upper-tail form; avoid 1-pcauchy cancellation for small P.
  logp <- pcauchy(T,lower.tail=FALSE,log.p=TRUE)
  -logp/log(10)
}
ans <- mapply(one,x$GLOBAL,x$LOCAL,x$WG,x$WL)
write.table(data.frame(EXPECTED=ans),args[2],quote=FALSE,sep='\t',row.names=FALSE)
