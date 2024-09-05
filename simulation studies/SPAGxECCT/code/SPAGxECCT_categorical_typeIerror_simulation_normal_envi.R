# cd /gdata01/user/yuzhuoma/GxE/code_ocean/
# sbatch -J categorical_typeIerror --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript SPAGxECCT_categorical_typeIerror_simulation_normal_envi.R $SLURM_ARRAY_TASK_ID'

# args=commandArgs(TRUE)
# print(args)
# print(sessionInfo())
# n.cpu=as.numeric(args)

N = 10000  # sample size
#### simulate genotypes of 1000 SNPs with genetic main effects
nSNP = 1000                                                    # number of SNPs with genetic main effects
MAF = runif(nSNP, 0.05, 0.5)                                   # minor allele frequencies of SNPs with genetic main effects
GenoMat_1 = t(matrix(rbinom(N*nSNP, 2, MAF), nSNP, N))         # Genotype matrix of SNPs with genetic main effects
rownames(GenoMat_1) = paste0("IID-",1:N)
colnames(GenoMat_1) = paste0("SNP-",1:nSNP)

idx.ratiosVec = c(1, 2, 3)
idx.ratios = idx.ratiosVec[1] # phenotypic distribution (ratio)

if(idx.ratios == 1){
  ratios = c(1,1,1,1)
}
if(idx.ratios == 2){
  ratios = c(10,1,1,1)
}
if(idx.ratios == 3){
  ratios = c(100,1,1,1)
}

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
                                 nSNP,                         # number of SNPs
                                 ratios,                       # Phenotypic distribution
                                 GMat,                         # Genotype matrix
                                 beta.true = c(0.5, 0.5, 0.5), # Fixed non-genetic coefficients
                                 bVec = 0)                     # Random effect
{
  seed = sample(1e9, 1);
  cat("random number seed is", seed)
  set.seed(seed)
  gamma = runif(nSNP, -0.4, 0.4) # Genetic effect vector
  
  # covariate matrix
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  X = cbind(X1, X2, E)
  
  # ordinal categorical outcome
  eta.true = X %*% beta.true + GMat%*%gamma + bVec
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



#### simulate ordinal categorical phenotype with certain prevalence
data = data.simu.categorical(N = N, nSNP = nSNP, ratios = ratios, GMat = GenoMat_1)
Phen.mtx = data.frame(ID = paste0("IID-",1:N),
                      Y = data$Y,
                      X1 = data$X1,
                      X2 = data$X2,
                      E = data$E)



#### simulate genotypes of 10000 SNPs without genetic main effects
nSNP_0 = 10000                                                         # number of SNPs without genetic main effects
MAF_0 = 0.01                                                           # minor allele frequencies of SNPs without genetic main effects
MAFVec_0 = runif(nSNP_0, MAF_0, MAF_0)                                   
GenoMat_0 = t(matrix(rbinom(N*nSNP_0, 2, MAFVec_0), nSNP_0, N))        # Genotype matrix of SNPs without genetic main effects
rownames(GenoMat_0) = paste0("IID-",1:N)
colnames(GenoMat_0) = paste0("SNP-",1:nSNP_0)


#### evaluate type one error rates using SPAGxECCT

#### install and library SPAGxECCT package
library(remotes)                            # remotes library requires less dependency packages than devtools
install_github("YuzhuoMa97/SPAGxECCT")      # The INSTALL_opts is required in Windows OS.
library(SPAGxECCT)
# ?SPAGxECCT                                # manual of SPAGxECCT package
# or source("~/capsule/code/R/SPAGxECCT.R")


X1 = Phen.mtx$X1                   # Covariate 1
X2 = Phen.mtx$X2                   # Covariate 2
E = Phen.mtx$E                     # environmental factor (normal distributed)
Cova.mtx = Phen.mtx[,c("X1","X2")] # covariate matrix excluding the environmental factor E

### fit genotype-independent (covariate-only) model
### null model residuals

R = get_cate_resid(as.factor(Y) ~ X1+X2+E, data = Phen.mtx)

#### calculate p values for SNPs without marginal genetic effect 
categorical_res_0 = SPAGxE_CCT(traits = "categorical",           # categorical trait analysis
                               Geno.mtx = GenoMat_0,             # genotype vector
                               R = R,                            # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                               E = E,                            # environmental factor
                               Phen.mtx = Phen.mtx,              # include surv.time, event
                               Cova.mtx = Cova.mtx)              # covariate matrix excluding the environmental factor E

# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with categorical phenotypes
head(categorical_res_0)

#### calculate p values for SNPs with marginal genetic effect 
categorical_res_1 = SPAGxE_CCT(traits = "categorical",           # categorical trait analysis
                               Geno.mtx = GenoMat_1,             # genotype vector
                               R = R,                            # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                               E = E,                            # environmental factor
                               Phen.mtx = Phen.mtx,              # include surv.time, event
                               Cova.mtx = Cova.mtx)              # covariate matrix excluding the environmental factor E

# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with categorical phenotypes
head(categorical_res_1)

#### save results
write.csv(categorical_res_0, file = "~/capsule/results/SPAGxECCT_categorical_typeIerror_simulation_normal_envi_example_results_0.csv")
write.csv(categorical_res_1, file = "~/capsule/results/SPAGxECCT_categorical_typeIerror_simulation_normal_envi_example_results_1.csv")
