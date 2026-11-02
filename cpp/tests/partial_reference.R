# Independent scalar R oracle for the declared guarded partial-normal CGF.
# Uses R quantile, uniroot and pnorm; no production C++ implementation is called.
partial_reference_audit <- list()

reference_cgf <- function(w, h, q, central = integer()) {
    exact <- setdiff(seq_along(w), central)
    deterministic <- exact[q[exact] == 1]
    shift <- sum(h[deterministic] * w[deterministic])
    random <- exact[h[exact] > 0 & w[exact] != 0 & q[exact] > 0 & q[exact] < 1]
    a <- w[random]; b <- h[random]; p <- q[random]
    mu <- shift + sum(h[central] * q[central] * w[central])
    variance <- sum(h[central] * q[central] * (1-q[central]) * w[central]^2)
    function(t) {
        u <- t*a
        # Log-space Bernoulli moment, with positive/negative tilts separated.
        logmoment <- numeric(length(u))
        positive <- u >= 0
        logmoment[positive] <- u[positive] + log(p[positive] + (1-p[positive])*exp(-u[positive]))
        logmoment[!positive] <- log1p(p[!positive]*expm1(u[!positive]))
        tilted <- plogis(qlogis(p)+u)
        c(K=sum(b*logmoment)+t*mu+.5*t*t*variance,
          K1=sum(b*a*tilted)+mu+t*variance,
          K2=sum(b*a*a*tilted*(1-tilted))+variance)
    }
}

reference_support <- function(w, h, q) {
    active <- h > 0 & w != 0
    random <- active & q > 0 & q < 1
    shift <- sum(h[active & q == 1]*w[active & q == 1])
    w0 <- w[random]; h0 <- h[random]; q0 <- q[random]
    low <- shift + sum(h0*pmin(w0,0))
    high <- shift + sum(h0*pmax(w0,0))
    lower_mass <- sum(h0*ifelse(w0>0,log1p(-q0),log(q0)))
    upper_mass <- sum(h0*ifelse(w0>0,log(q0),log1p(-q0)))
    scale <- sum(abs(h[active]*w[active]))
    function(target, upper) {
        tol <- 8*.Machine$double.eps*max(scale,abs(target))
        if(target < low-tol) return(if(upper) 0 else -Inf)
        if(target > high+tol) return(if(upper) -Inf else 0)
        if(abs(target-low) <= tol) return(if(upper) 0 else lower_mass)
        if(abs(target-high) <= tol) return(if(upper) upper_mass else 0)
        NA_real_
    }
}

reference_tail <- function(cgf, target, upper, variance, support=NULL) {
    if(!is.null(support)) {
        bound <- support(target, upper)
        if(!is.na(bound)) return(list(logp=bound, root=NA_real_, value=NULL))
    }
    f <- function(t) unname(cgf(t)["K1"]-target)
    width <- 1/sqrt(variance)
    if(upper) {
        lo <- 0; hi <- width
        for(i in seq_len(150)) {if(f(hi)>=0) break; hi <- hi*2}
    } else {
        lo <- -width; hi <- 0
        for(i in seq_len(150)) {if(f(lo)<=0) break; lo <- lo*2}
    }
    root <- uniroot(f,c(lo,hi),tol=1e-12,maxiter=400)$root
    value <- cgf(root)
    r <- sign(root)*sqrt(2*(root*target-value["K"]))
    v <- root*sqrt(value["K2"])
    corrected <- r+log(v/r)/r
    list(logp=unname(pnorm(corrected,lower.tail=!upper,log.p=TRUE)),root=root,value=value)
}

reference_guarded_spa <- function(w, h, q, mu, deviation, variance, partial=TRUE) {
    h <- rep(h,length.out=length(w)); q <- rep(q,length.out=length(w))
    full <- reference_cgf(w,h,q)
    support <- reference_support(w,h,q)
    targets <- c(mu+abs(deviation),mu-abs(deviation))
    eligible <- which(h>0 & w!=0 & q>0 & q<1)
    outliers <- integer()
    central <- seq_along(w)
    if(length(eligible)) {
        quartile <- quantile(w[eligible],c(.25,.75),names=FALSE,type=7)
        fence <- 1.5*diff(quartile)
        outliers <- eligible[w[eligible]<quartile[1]-fence | w[eligible]>quartile[2]+fence]
        central <- setdiff(seq_along(w),outliers)
    }
    central_variance <- sum(h[central]*q[central]*(1-q[central])*w[central]^2)
    use_partial <- partial && length(outliers)>0 && central_variance>0
    policy <- if(partial) "full_partition" else "full_reference"
    errors <- c(K0=NA_real_,ROOT=NA_real_,K2=NA_real_)
    if(use_partial && ( !is.na(support(targets[1],TRUE)) || !is.na(support(targets[2],FALSE)))) {
        use_partial <- FALSE; policy <- "full_support"
    }
    tails <- NULL
    if(use_partial) {
        approximation <- reference_cgf(w,h,q,central)
        candidate <- tryCatch(lapply(seq_along(targets),function(j)
            reference_tail(approximation,targets[j],j==1,variance)),error=function(e)NULL)
        accepted <- !is.null(candidate)
        if(accepted) {
            diagnostics <- lapply(seq_along(targets),function(j) {
                result <- candidate[[j]]; actual <- full(result$root); approx <- result$value
                k0 <- abs(actual["K"]-approx["K"])
                root_error <- abs(actual["K1"]-targets[j])/sqrt(actual["K2"])
                curvature <- abs(actual["K2"]-approx["K2"])/actual["K2"]
                eps <- 64*.Machine$double.eps
                floor0 <- eps*(1+abs(actual["K"])+abs(approx["K"])+abs(result$root*targets[j]))
                floor1 <- eps*(1+abs(targets[j])+abs(actual["K1"]))/sqrt(actual["K2"])
                okay <- all(is.finite(c(k0,root_error,curvature,result$logp))) && actual["K2"]>0 &&
                        k0<=1e-3+floor0 && root_error<=1e-3+floor1 && curvature<=1e-3+eps
                list(errors=unname(c(k0,root_error,curvature)),accepted=isTRUE(okay))
            })
            errors[] <- apply(do.call(rbind,lapply(diagnostics,`[[`,"errors")),2,max)
            accepted <- all(vapply(diagnostics,`[[`,logical(1),"accepted"))
        }
        if(accepted) {tails <- candidate; policy <- "partial"} else policy <- "full_guard"
    }
    if(is.null(tails)) tails <- lapply(seq_along(targets),function(j)
        reference_tail(full,targets[j],j==1,variance,support))
    lp <- vapply(tails,`[[`,numeric(1),"logp")
    largest <- max(lp)
    combined <- if(largest == -Inf) -Inf else largest+log(sum(exp(lp-largest)))
    if(partial) {
        partial_reference_audit[[length(partial_reference_audit)+1L]] <<- data.frame(
            POLICY=policy,N=length(w),N_ELIGIBLE=length(eligible),N_EXACT=length(outliers),
            N_CENTRAL=length(central),CENTRAL_VARIANCE=central_variance,
            K0_ERROR=errors[1],ROOT_ERROR=errors[2],K2_ERROR=errors[3],LOG10P=-min(combined,0)/log(10))
    }
    -min(combined,0)/log(10)
}

write_partial_reference_audit <- function(path) {
    if(length(partial_reference_audit)) write.table(do.call(rbind,partial_reference_audit),
        path,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
