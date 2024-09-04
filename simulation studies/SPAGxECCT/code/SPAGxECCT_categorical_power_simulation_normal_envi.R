# cd /gdata01/user/yuzhuoma/GxE/code_ocean/
# sbatch -J categorical_typeIerror --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript SPAGxECCT_categorical_power_simulation_normal_envi.R $SLURM_ARRAY_TASK_ID'

# args=commandArgs(TRUE)
# print(args)
# print(sessionInfo())
# n.cpu=as.numeric(args)

N = 50000  # sample size
nSNP = 10  # number of SNPs with marginal GxE effects

idx.ratiosVec = c(1, 2, 3)
idx.ratios = idx.ratiosVec[3] # phenotypic distribution (ratio)
MAFVec = c(0.3, 0.05, 0.005)
GammaVec = seq(0, 2, 0.2)

if(idx.ratios == 1){
  ratios = c(1,1,1,1)
}
if(idx.ratios == 2){
  ratios = c(10,1,1,1)
}
if(idx.ratios == 3){
  ratios = c(100,1,1,1)
}

MAF = MAFVec[2]      # minor allele frequencies of SNPs with marginal GxE effects
Gamma = GammaVec[3]  # marginal GxE effect size


### sub-functions to simulate ordinal categorical phenotype
#### lower function to estimate epsilons given a ratio (phenotypic distribution). Will be used in data.simu.categorical().
getProb = function(eps,
                   eta.true,
                   prob,
                   seed)
{
  set.seed(seed)
  n = length(eta.true)
  mu = exp(eps-eta.true) / (1+exp(eps-eta.true))
  Y.latent = runif(n)
  diffprob = mean(Y.latent < mu) - prob
  return(diffprob)
}

getEps = function(ratios,
                  eta.true,
                  seed)
{
  sumR = sum(ratios)
  cumR = 0
  J = length(ratios)
  Eps = c()
  for(i in 1:(J-1)){
    cumR = cumR + ratios[i]
    eps = uniroot(getProb, c(-100,100), eta.true = eta.true, prob = cumR/sumR, seed = seed)
    Eps = c(Eps, eps$root)
  }
  return(Eps)
}

#### Main function to simulate categorical phenotype
data.simu.categorical = function(N,                            # Sample size
                                 ratios,                       # Phenotypic distribution
                                 gamma,
                                 g,                            # Genotype
                                 beta.true = c(0.5, 0.5, 0.5), # Fixed non-genetic coefficients
                                 bVec = 0)                     # Random effect
{
  seed = sample(1e9, 1);
  cat("random number seed is", seed)
  set.seed(seed)
  
  # covariate matrix
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  X = cbind(X1, X2, E)
  
  # ordinal categorical outcome
  eta.true = X %*% beta.true + g*E*gamma + bVec
  Eps = getEps(ratios, eta.true, seed)
  Y.latent = runif(N)
  Y = rep(0, N)
  for(g in 1:length(Eps)){
    mu = exp(Eps[g]-eta.true)/(1+exp(Eps[g]-eta.true))
    Y[Y.latent > mu] = g
  }
  Y = Y + 1 # added because levels are 1, ..., J
  out = data.frame(Y, X1, X2, E, bVec) 
  return(out)
}

#### Function to get residuals from genotype-independent (covariates-only) model for categorical phenotype analyses
library(ordinal)

get_cate_resid = function(formula,
                          data)     # without intercept
{
  obj.polmm = clm(formula, data = data)
  coeff = obj.polmm$coefficients
  eps = obj.polmm$alpha
  beta = obj.polmm$beta
  Cov.mtx = as.matrix(data[,c("X1","X2","E")])
  eta.est = Cov.mtx %*% beta
  J = length(eps) + 1
  
  nu.est = c()
  for (j in 1:length(eps)) {
    l = eps[j] - eta.est
    v = exp(l)/(1 + exp(l))
    nu.est = cbind(nu.est, v)
  }
  nu.est = cbind(0, nu.est, 1)
  
  mu.est = c()
  for (j in 1:J) {
    diff = nu.est[,j+1] - nu.est[,j]
    mu.est = cbind(mu.est, diff)
  }
  
  R = c()
  for (j in 1:J) {
    Rvec = ((nu.est[,j]*(1-nu.est[,j]) - nu.est[,j+1]*(1-nu.est[,j+1]))/mu.est[,j]) - 
      ((nu.est[,J]*(1-nu.est[,J]) - nu.est[,J+1]*(1-nu.est[,J+1]))/mu.est[,J])
    
    R = cbind(R, Rvec)
  }
  
  Y_tilde = c()
  for (i in 1:length(Phen.mtx$Y)) {
    y_tilde = rep(0,J)
    y_tilde[Phen.mtx$Y[i]] = 1
    Y_tilde = rbind(Y_tilde, y_tilde)
  }
  
  resid = rowSums(R * (Y_tilde - mu.est))
  # re = resid
  return(resid)
}


########## generate genotype matrix with marginal GxE effect beta_GxE!=0

MAFVec = runif(nSNP, MAF, MAF)
GenoMat = t(matrix(rbinom(N*nSNP, 2, MAFVec), nSNP, N))
rownames(GenoMat) = paste0("IID-",1:N)
colnames(GenoMat) = paste0("SNP-",1:nSNP)


#### calculate p values to evaluate powers of SPAGxECCT

#### install SPAGxECCT package
# library(devtools)  # author version: 2.4.5
# install_github("YuzhuoMa97/SPAGxECCT")
library(SPAGxECCT)
# ?SPAGxECCT  # manual of SPAGxECCT package

categorical_res = c()
for (i in 1:nSNP) {
  print(i)
  G = GenoMat[,i]  # genotype to test
  
  data = data.simu.categorical(N = N, ratios = ratios, gamma = Gamma, g = G) # simulate categorical phenotype
  
  ### phenotype matrix
  Phen.mtx = data.frame(ID = paste0("IID-",1:N),
                        Y = data$Y,
                        X1 = data$X1,
                        X2 = data$X2,
                        E = data$E)
  
  ### fit genotype-independent (covariate-only) model
  ### null model residuals
  R = get_cate_resid(as.factor(Y) ~ X1+X2+E, data = Phen.mtx)
  
  #### calculate p values for SNPs with marginal GxE effect 
  categorical_res_oneSNP = SPAGxE_CCT(traits = "categorical",             # ordinal categorical trait analysis
                                      Geno.mtx = as.matrix(G),            # genotype vector
                                      R = R,                              # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                                      E = Phen.mtx$E,                     # environmental factor
                                      Phen.mtx = Phen.mtx,                # include surv.time, event
                                      Cova.mtx = Phen.mtx[,c("X1","X2")]) # covariate matrix excluding the environmental factor E
  
  categorical_res = rbind(categorical_res, categorical_res_oneSNP)
}

# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with ordinal categorical phenotypes
head(categorical_res)




