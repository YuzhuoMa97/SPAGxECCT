# cd /gdata01/user/yuzhuoma/GxE/code_ocean/
# sbatch -J binary_typeIerror --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript SPAGxECCT_binary_typeIerror_simulation_normal_envi.R $SLURM_ARRAY_TASK_ID'

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

prevalenceVec = c(0.05, 0.2)
prevalence = prevalenceVec[1] # prevalence

#### lower function to estimate beta0 given an prevalence. Will be used in data.simu.surv().
f.binary.new = function(N,                  # Sample size
                        nSNP,               # SNPs with genetic main effects
                        prevalence,         # Prevalence
                        beta0,              # Intercept
                        # gamma,            # Genetic effect
                        GMat,               # Genotype matrix
                        seed,               # random seed
                        bVec = 0)           # Additional effect, could be random effect. If bVec = 0 (default), then no additional effect is included.
{
  set.seed(seed)
  gamma = runif(nSNP, -0.4, 0.4) # Genetic effect vector
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  betas = c(0.5, 0.5, 0.5)
  
  eta = beta0 + betas[1]*X1 + betas[2]*X2 + betas[3]*E + GMat%*%gamma + bVec
  mu = exp(eta) / (1 + exp(eta))   # The probability being a case given the covariates, genotypes, and addition effect
  Y = rbinom(N, 1, mu)    # Case-control status
  re = sum(Y) - N * prevalence
  return(re)
}

#### Main function to simulate binary phenotype
data.simu.binary.new = function(N,              # Sample size
                                nSNP,           # SNPs with genetic main effects
                                prevalence,     # Prevalence
                                # gamma,        # Genetic effect
                                GMat,           # Genotype matrix
                                bVec = 0)       # Additional effect, could be random effect. If bVec = 0 (default), then no additional effect is included.
{ 
  beta0 = 0
  
  while (beta0 == 0) {
    seed = sample(1e9, 1);
    cat("random number seed is", seed)
    beta0 = uniroot(f.binary.new, c(-30,30), N = N, nSNP = nSNP,
                    prevalence = prevalence, GMat = GMat, seed = seed, bVec = bVec)
    beta0 = beta0$root
  }
  
  set.seed(seed)
  gamma = runif(nSNP, -0.4, 0.4) # Genetic effect vector
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  betas = c(0.5, 0.5, 0.5)
  
  eta = beta0 + betas[1]*X1 + betas[2]*X2 + betas[3]*E + GMat%*%gamma + bVec
  mu = exp(eta) / (1 + exp(eta))   # The probability being a case given the covariates, genotypes, and addition effect
  Y = rbinom(N, 1, mu)    # Case-control status
  out = data.frame(Y, X1, X2, E, bVec)
  return(out)
}



#### simulate binary phenotype with certain prevalence
data = data.simu.binary.new(N = N, nSNP = nSNP, prevalence = prevalence, GMat = GenoMat_1)
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

R = SPA_G_Get_Resid("binary",
                    Y~X1+X2+E,
                    data=Phen.mtx,
                    pIDs=Phen.mtx$ID,
                    gIDs=Phen.mtx$ID)

#### calculate p values for SNPs without marginal genetic effect 
binary_res_0 = SPAGxE_CCT(traits = "binary",                # binary trait analysis
                          Geno.mtx = GenoMat_0,             # genotype vector
                          R = R,                            # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                          E = E,                            # environmental factor
                          Phen.mtx = Phen.mtx,              # include surv.time, event
                          Cova.mtx = Cova.mtx)              # covariate matrix excluding the environmental factor E

# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with binary phenotypes
head(binary_res_0)

#### calculate p values for SNPs with marginal genetic effect 
binary_res_1 = SPAGxE_CCT(traits = "binary",                # binary trait analysis
                          Geno.mtx = GenoMat_1,             # genotype vector
                          R = R,                            # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                          E = E,                            # environmental factor
                          Phen.mtx = Phen.mtx,              # include surv.time, event
                          Cova.mtx = Cova.mtx)              # covariate matrix excluding the environmental factor E

# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with binary phenotypes
head(binary_res_1)


#### save results
write.csv(binary_res_0, file = "~/capsule/results/SPAGxECCT_binary_typeIerror_simulation_normal_envi_example_results_0.csv")
write.csv(binary_res_1, file = "~/capsule/results/SPAGxECCT_binary_typeIerror_simulation_normal_envi_example_results_1.csv")

