##### SPAmix functions with input of local ancestry information

# Gets NULL model residuals for SPA-G
SPA_G_Get_Resid = function(traits="survival/binary",
                           formula=NULL,
                           data=NULL,
                           pIDs=NULL,
                           gIDs=NULL,
                           range=c(-100,100),
                           length.out = 10000,
                           ...)
  
{
  if(traits=="survival"){
    Call = match.call()
    ### Fit a Cox model
    obj.coxph = coxph(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.coxph, range)
    
    ### Get the covariate matrix to adjust for genotype
    resid = obj.coxph$residuals
    
    re = resid
  }
  else if(traits=="binary"){
    Call = match.call()
    ### Fit a logistic model
    obj.logistic = glm(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.logistic, range)
    
    ### Get the covariate matrix to adjust for genotype
    mu = obj.logistic$fitted.values
    resid = obj.logistic$y - mu
    re = resid
  }
  return(re)
}

# MAF is the minor allele frequency of the tested marker in the index ancestry
# MAFVec is the vector (length = sample size) of minor allele frequency of the tested marker in the index ancestry
# haplo_numVec is the vector (length = sample size) of the number haplotype of the tested marker in the index ancestry

# The MGF of G (genotype from the index ancestry)
M_G0 = function(t, MAF, haplo_num){
  re = (1 - MAF + MAF * exp(t))^haplo_num
  return(re)
}

# The first derivative of the MGF of G (genotype from the index ancestry)
M_G1 = function(t, MAF, haplo_num){
  re = haplo_num*(MAF * exp(t))*(1 - MAF + MAF * exp(t))^(haplo_num - 1)
  return(re)                           
}

# The second derivative of the MGF of G (genotype from the index ancestry)
M_G2 = function(t, MAF, haplo_num){
  re = haplo_num * (haplo_num-1) * (MAF * exp(t))^2 * (1 - MAF + MAF * exp(t))^(haplo_num - 2) +
    haplo_num * (MAF * exp(t)) * (1 - MAF + MAF * exp(t)) ^ (haplo_num - 1)
  return(re)
}


# The CGF of G (genotype from the index ancestry)
K_G0 = function(t, MAF, haplo_num){
  re = log(M_G0(t, MAF, haplo_num))
  return(re)
}

# The first derivative of the CGF of G (genotype from the index ancestry)
K_G1 = function(t, MAF, haplo_num){
  re = M_G1(t, MAF, haplo_num)/M_G0(t, MAF, haplo_num)
  return(re)
}

# The second derivative of the CGF of G (genotype from the index ancestry)
K_G2 = function(t, MAF, haplo_num){
  re = (M_G0(t, MAF, haplo_num)*M_G2(t, MAF, haplo_num)-M_G1(t, MAF, haplo_num)^2)/M_G0(t, MAF, haplo_num)^2
  return(re)
}

# The CGF of score test statistic for the index ancestry
H_org = function(t, R, MAFVec, haplo_numVec){
  
  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(K_G0(t1*R, MAFVec, haplo_numVec))
  }
  return(out)
}

# The first derivative of the CGF of score test statistic for the index ancestry
H1_adj = function(t, R, s, MAFVec, haplo_numVec)
{
  n.t = length(t)
  out = rep(0,n.t)
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(R*K_G1(t1*R, MAFVec, haplo_numVec)) - s
  }
  return(out)
}


# The second derivative of the CGF of score test statistic for the index ancestry
H2 = function(t, R, MAFVec, haplo_numVec)
{
  n.t = length(t)
  out = rep(0,n.t)
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(R^2*K_G2(t1*R, MAFVec, haplo_numVec))
  }
  return(out)
}


GetProb_SPA_G = function(MAFVec, R, haplo_numVec, s, lower.tail){
  
  out = uniroot(H1_adj, c(-20,20), extendInt = "upX",
                R=R, s=s, MAFVec=MAFVec, haplo_numVec=haplo_numVec)
  zeta = out$root
  
  k1 = H_org(t=zeta, R=R, MAFVec=MAFVec, haplo_numVec=haplo_numVec)
  k2 = H2(t=zeta, R=R, MAFVec=MAFVec, haplo_numVec=haplo_numVec)
  
  temp1 = zeta * s - k1
  
  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  
  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  re = pval
  return(re)
}

# updataed on 2023-10-24
fastgetroot_H1_new = function(init.t,
                              R,          # residuals
                              s,          # observed test statistic
                              MAFVec,
                              haplo_numVec,
                              tol = .Machine$double.eps^0.25,
                              maxiter = 100)
{
  t = init.t;
  H1_adj_eval = 0
  diff.t = Inf
  converge = T
  for(iter in 1:maxiter){
    old.t = t
    old.diff.t = diff.t
    old.H1 = H1_adj_eval
    H1_adj_eval = H1_adj(t, R, s, MAFVec, haplo_numVec)
    H2_eval = H2(t, R, MAFVec, haplo_numVec)
    
    diff.t = -1 * H1_adj_eval / H2_eval
    
    # updated on 2023-10-24
    if(is.na(diff.t)){
      diff.t =  5
    }
    
    if(is.na(H1_adj_eval) | is.infinite(H1_adj_eval)){
      # checked it on 07/05:
      # if the solution 't' tends to infinity, 'K2_eval' tends to 0, and 'K1_eval' tends to 0 very slowly.
      # then we can set the one sided p value as 0 (instead of setting converge = F)
      t = sign(s)*Inf
      H2_eval = 0;
      break;
    }
    if(sign(H1_adj_eval) != sign(old.H1)){
      while(abs(diff.t) > abs(old.diff.t) - tol){
        diff.t = diff.t/2
      }
    }
    if(is.na(diff.t)){
      diff.t =  5
    }
    if(abs(diff.t) < tol) break;
    t = old.t + diff.t
  }
  
  # updated on 2023-10-24
  
  if(iter == maxiter) {
    converge = F
    t = NA}
  
  return(list(root = t,
              iter = iter,
              converge = converge,
              H2_eval = H2_eval))
}


# updataed on 2023-10-24

GetProb_SPA_G_new = function(init.t = 0, MAFVec, R, haplo_numVec, s, lower.tail){
  
  # out = uniroot(H1_adj, c(-10,10), extendInt = "upX",
  #               R=R, s=s, MAFVec=MAFVec)
  out = fastgetroot_H1_new(init.t, R = R, s = s, haplo_numVec = haplo_numVec, MAFVec = MAFVec)
  
  zeta = out$root
  
  k1 = H_org(zeta, R=R, MAFVec=MAFVec, haplo_numVec = haplo_numVec)
  k2 = H2(zeta, R=R, MAFVec=MAFVec, haplo_numVec = haplo_numVec)
  
  temp1 = zeta * s - k1
  
  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  
  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  
  if(is.na(pval)){
    pval =  0
  }
  
  re = pval
  return(re)
}





#### SPAmix #########################################################################

SPAmix_localance_one_SNP = function(g,
                                    R,
                                    haplo_numVec,
                                    # obj.null,
                                    Cutoff = 2,
                                    impute.method = "fixed",
                                    missing.cutoff = 0.15,
                                    min.maf = 0.00001,          # update on 2022-08-16 : replace 0.0001 by 0.000001
                                    G.model = "Add")
{
  ## calculate MAF and update genotype vector
  # MAF = mean(g, na.rm=T)/2     # MAF in the index ancestry
  MAF = sum(g)/sum(haplo_numVec)   # MAF in the index ancestry
  
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N
  
  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }
  ####################################################### AF, not MAF  
  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }
  
  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)
  
  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA, NA, NA))
  
  ## Score statistic
  S = sum(g * R)
  
  ## estimated variance without adjusting for covariates
  # N1set = 1:N
  # N0 = 0
  
  N1set = which(g!=0)  # position of non-zero genotypes
  N0 = N-length(N1set)
  
  MAC = sum(g)    
  MAF.est.Vec = c(rep(MAF, N))
  
  S.mean = sum(MAF.est.Vec * haplo_numVec * R)
  
  g.var.est.Vec = haplo_numVec * MAF.est.Vec * (1 - MAF.est.Vec)
  S.var = sum(R^2 * g.var.est.Vec)
  
  z = (S - S.mean)/sqrt(S.var)
  
  if(abs(z) < Cutoff){
    pval.norm = pnorm(abs(z), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.mean, S.var, z, MAC))  
  }
  
  pval1 = GetProb_SPA_G_new(MAFVec=MAF.est.Vec, R=R, haplo_numVec=haplo_numVec, s=max(S, (2*S.mean-S)), lower.tail = FALSE) # SPA-G p value 
  pval2 = GetProb_SPA_G_new(MAFVec=MAF.est.Vec, R=R, haplo_numVec=haplo_numVec, s=min(S, (2*S.mean-S)), lower.tail = TRUE) # SPA-G p value 
  
  pval3 = pnorm(abs(z), lower.tail = FALSE) # Normal
  pval4 = pnorm(-abs(z), lower.tail = TRUE) # Normal
  
  pval.spa.G = pval1 + pval2
  pval.norm = pval3 + pval4
  
  pval = c(pval.spa.G, pval.norm) # 2 elements: element 1 is from SPA, element 2 is from Normal
  
  return(c(MAF, missing.rate, pval, S, S.mean, S.var, z, MAC))  
}


SPAmix_localance = function(Geno.mtx,
                            R,
                            haplo.mtx,
                            # obj.null,
                            Cutoff = 2,
                            impute.method = "fixed",
                            missing.cutoff = 0.15,
                            min.maf = 0.00001,          
                            G.model = "Add")
  
  
{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  G.model=G.model)
  
  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))
  
  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 9) 
  colnames(output) = c("MAF","missing.rate","Pvalue.SPAmix.index.ance","Pvalue.norm.index.ance",
                       "Stat","Mean","Var","z", "MAC") # update on 2022-10-05 : MAF.est.negative.num 
  rownames(output) = colnames(Geno.mtx)
  
  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())
  
  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    # print(i)
    g = Geno.mtx[,i]
    haplo_numVec = haplo.mtx[,i]
    output.one.SNP = SPAmix_localance_one_SNP(g,
                                              R,
                                              haplo_numVec,
                                              Cutoff,
                                              impute.method,
                                              missing.cutoff,
                                              min.maf,         
                                              G.model)
    
    output[i,] = output.one.SNP
  }
  
  output = as.data.frame(cbind(rsID = colnames(Geno.mtx), output))
  
  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}




check_input = function(pIDs, gIDs, obj, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = obj$residuals
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")
  
  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)
  
  return(p2g)
}

check_input_Resid = function(pIDs, gIDs, R, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = R
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")
  
  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)
  
  return(p2g)
}



################## 2024-04-25 GxE --------------------------------------------------------------------------------------------------------------------






SPAGxE_CCT_mix_localance_one_SNP = function(traits="survival/binary/quantitative/categorical",
                                            g,                     # ancestry-specific genotype vector
                                            haplo_numVec,
                                            R,                     # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                                            E,                     # environmental factor
                                            Phen.mtx,              # include surv.time, event
                                            Cova.mtx,              #  other covariates (such as age, gender, and top PCs) excluding E
                                            epsilon = 0.001,       # a fixed value
                                            Cutoff = 2,
                                            impute.method = "fixed",
                                            missing.cutoff = 0.15,
                                            min.maf = 0.00001,
                                            G.model = "Add")

{
  MAF = sum(g)/sum(haplo_numVec)   # MAF in the index ancestry
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N
  
  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }
  
  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }
  
  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)
  
  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA, NA, NA, NA))
  
  S1 = sum(g*R)                    # test statistic for marginal genetic effect
  
  MAC = sum(g)    
  MAF.est.Vec = c(rep(MAF, N))
  
  S1.mean = sum(MAF.est.Vec * haplo_numVec * R)
  
  g.var.est.Vec = haplo_numVec * MAF.est.Vec * (1 - MAF.est.Vec)
  
  VarS1 = sum(R^2 * g.var.est.Vec)         # estimated variance of S1
  
  Z1 = (S1 - S1.mean) / sqrt(VarS1)            # standardize S1
  pval.norm1 = pnorm(-abs(Z1))*2   # p value for marginal genetic effect from normal approximation
  
  print(pval.norm1)
  
  if(pval.norm1 > epsilon){
    
    S2 = sum(g*E*R)                # test statistic for marginal GxE effect
    lambda = sum(haplo_numVec * E * R^2) / sum(haplo_numVec * R^2)   # lambda = Cov/VarS1
    S_GxE = S2 - lambda * S1       # new test statistic for marginal GxE effect
    R.new = (E - lambda) * R       # new residuals
    
    ################### SPA
    
    res = SPAmix_localance_one_SNP(g=g, R=R.new, haplo_numVec=haplo_numVec, min.maf=min.maf)        # the third element is p-value from SPA-G, the fourth element is p-value from normal approximation
    
    
    pval.spaG = res[3]                                # p value from SPA-G
    pval.norm = res[4]                                # p-value from normal approximation
    
    # updated on 2023-12-27
    pval.output = c(pval.spaG, pval.spaG, pval.spaG, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
  }else{
    
    W = cbind(1, g)
    
    # updated on 2023-12-27
    R0 = R - (W)%*%(solve(t(W)%*%W))%*%(t(W)%*%R) # null model residuals adjusting for G
    
    S_GxE0 = sum(g*E*R0)  # test statistic for marginal GxE effect
    R.new0 = E*R0
    
    # updated on 2023-10-24
    
    res = SPAmix_localance_one_SNP(g=g, R=R.new0, haplo_numVec=haplo_numVec, min.maf=min.maf)              # the third element is p-value from SPA-G, the fourth element is p-value from normal approximation
    
    pval.spaGxE = res[3]                                # p value from SPA-G
    pval.norm = res[4]                                # p-value from normal approximation
    
    # Wald test
    
    if(traits=="survival"){
      
      # data0 = coxph(formula = Surv(Phen.mtx$surv.time, Phen.mtx$event) ~ as.matrix(Cova.mtx)+E+g+g*E, iter.max = 1000)
      data0 = coxph(formula = Surv(Phen.mtx$surv.time, Phen.mtx$event) ~ as.matrix(Cova.mtx)+haplo_numVec+E+g+g*E, iter.max = 1000) # updated  
      
      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect
      
      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))
      
      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }
    
    if(traits=="binary"){
      
      # data0 = glm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E, family = "binomial") # updated on 2023-12-27
      data0 = glm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+haplo_numVec+E+g+g*E, family = "binomial") # updated 
      
      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect
      
      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))
      
      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }
    
    if(traits=="quantitative"){
      
      # data0 = lm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E) # updated on 2023-12-27
      data0 = lm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+haplo_numVec+E+g+g*E) # updated  
      
      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|t|)"]     # Wald test p value for marginal GxE effect
      
      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))
      
      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }
    
    if(traits=="categorical"){
      
      # data0 = clm(formula = as.factor(Phen.mtx$Y) ~ as.matrix(Cova.mtx)+E+g+g*E)
      data0 = clm(formula = as.factor(Phen.mtx$Y) ~ as.matrix(Cova.mtx)+haplo_numVec+E+g+g*E) # updated  
      
      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect
      
      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))
      
      
      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }
  }
  
  output.one.snp = c(MAF, missing.rate, pval.output, S1, VarS1, Z1)
  
  return(output.one.snp)
}







SPAGxE_CCT_mix_localance = function(traits="survival/binary/quantitative/categorical",
                                    Geno.mtx,
                                    haplo.mtx,
                                    R,                     # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                                    Cutoff = 2,
                                    E,                     # environmental factor
                                    Phen.mtx,              # include surv.time, event
                                    Cova.mtx,              #  other covariates (such as age, gender, and top PCs) excluding E
                                    epsilon = 0.001,       # a fixed value
                                    impute.method = "fixed",
                                    missing.cutoff = 0.15,
                                    min.maf = 0.00001,          
                                    G.model = "Add")


{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  G.model=G.model)
  
  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))
  
  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 10) 
  
  colnames(output) = c("MAF","missing.rate",
                       "p.value.spaGxE.index.ance","p.value.spaGxE.Wald.index.ance","p.value.spaGxE.CCT.Wald.index.ance","p.value.normGxE.index.ance",
                       "p.value.betaG.index.ance", "Stat.betaG","Var.betaG","z.betaG") 
  
  rownames(output) = colnames(Geno.mtx)
  
  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())
  
  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    
    g = Geno.mtx[,i]
    haplo_numVec = haplo.mtx[,i]

    output.one.SNP = SPAGxE_CCT_mix_localance_one_SNP(traits=traits,
                                                      g,                     # ancestry-specific genotype vector
                                                      haplo_numVec,
                                                      R,                     # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                                                      E,                     # environmental factor
                                                      Phen.mtx,              # include surv.time, event
                                                      Cova.mtx,              #  other covariates (such as age, gender, and top PCs) excluding E
                                                      epsilon,               # a fixed value
                                                      Cutoff,
                                                      impute.method,
                                                      missing.cutoff,
                                                      min.maf,
                                                      G.model)
    
    
    output[i,] = output.one.SNP
  }
  
  output = as.data.frame(cbind(rsID = colnames(Geno.mtx), output))
  
  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}



































