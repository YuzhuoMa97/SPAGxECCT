# cd /gdata01/user/yuzhuoma/GxE/code_ocean/
# sbatch -J binary_admixed --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript SPAGxEmixCCT_AdmixedPopulationAnalysis_binary_simulation_normal_envi.R $SLURM_ARRAY_TASK_ID'

# args=commandArgs(TRUE)
# print(args)
# print(sessionInfo())
# n.cpu=as.numeric(args)

library(dplyr)
library(tidyr)
library(data.table)

#### install and library SPAGxECCT package
library(remotes)                            # remotes library requires less dependency packages than devtools
install_github("YuzhuoMa97/SPAGxECCT")
library(SPAGxECCT)
# ?SPAGxECCT                                  # manual of SPAGxECCT package
# or source("~/capsule/code/R/SPAGxECCT.R")

N = 10000     # sample size of an admixed population from European (EUR) population and East Asia (EAS) population
N1 = N/2
N2 = N/2
nSNP = 1000   # number of SNPs to test

load("~/capsule/code/simulation studies/SPAGxEmixCCT/data/a.population.structure.6.RData") # true ancestry proportion vector
load("~/capsule/code/simulation studies/SPAGxEmixCCT/data/top10PCs_MAFsumover0.1.RData")   # SNP-derived top 10 PCs corresponding to the above ancestry proportion vector  

ancestry_Vec = a6 # true ancestry proportion vector
colnames(ancestry_Vec) = c("AnceProportion_EUR", "AnceProportion_EAS") # AnceProportion_EUR + AnceProportion_EAS = 1
rownames(ancestry_Vec) = paste0("IID-",1:N)

prevalence_EUR = 0.1  # prevalence in EUR
prevalence_EAS = 0.01 # prevalence in EAS

MAF_EUR = 0.3  # minor allele frequency in EUR
MAF_EAS = 0.01 # minor allele frequency in EAS

MAF_admixed = ancestry_Vec %*% c(prevalence_EUR, prevalence_EAS)  # minor allele frequency in the admixed population
colnames(MAF_admixed) = "MAF"

GenoMat_admixed = matrix(rbinom(N*nSNP, 2, MAF_admixed), N, nSNP) # genotype matrix of the admixed population
colnames(GenoMat_admixed) = paste0("SNP-", 1:nSNP)
rownames(GenoMat_admixed) = paste0("IID-", 1:N)

top4PCs = top10PCs[,1:4] # top 4 PCs (used in SPAGxEmixCCT)


#### sub-functions to simulate binary phenotype
#### simulate binary phenotype for admixed population analysis (e.g. EUR and EAS)

#### lower function to estimate beta0 given a prevalence. Will be used in data.simu.binary().
f1.binary = function(N1,                # Sample size
                     prevalence1,       # Prevalence 
                     beta0,             # Intercept
                     gamma,             # Genetic effect
                     g1,                # Genotype vector
                     gamma_GxE,         # Marginal GxE effect
                     bVec = 0)          # Additional effect, could be random effect. If bVec = 0 (default), then no additional effect is included.
{
  #set.seed(54769)
  X21 = rnorm(N)           # Covariate vector 1
  X31 = rbinom(N, 1, 0.5)  # Covariate vector 2
  E1 = rnorm(N1)
  eta = beta0 + 0.5 * X21 + 0.5 * X31 + 0.5 * E1 + gamma * g1 + gamma_GxE * g1 * E1 + bVec
  mu = exp(eta) / (1 + exp(eta))   # The probability being a case given the covariates, genotypes, and addition effect
  Y = rbinom(N, 1, mu)    # Case-control status
  re = sum(Y) - N * prevalence1
  return(re)
}

f2.binary = function(N1,
                     N2,
                     prevalence1,
                     prevalence2,
                     beta1,
                     gamma,
                     g1,
                     g2,
                     gamma_GxE,
                     bVec = 0)
  # seed1,
  # seed2)
{
  beta0 = uniroot(f1.binary, c(-100,100), N1 = N1, prevalence1 = prevalence1, gamma = gamma, g1 = g1, gamma_GxE = gamma_GxE, bVec = bVec)
  beta0 = beta0$root
  # set.seed(seed2)
  X1 = c(rep(0, N1), rep(1, N2))
  X11 = X1[c(1 : N1)]
  X12 = X1[c((N1 + 1) : N)]
  X22 = rnorm(N2)
  X32 = rbinom(N2, 1, 0.5)
  E2 = rnorm(N2)
  
  eta = beta0 + beta1 * X12 + 0.5 * X22 + 0.5 * X32 + 0.5 * E2 + gamma * g2 + gamma_GxE * g2 * E2 + bVec
  mu = exp(eta) / (1 + exp(eta))   # The probability being a case given the covariates, genotypes, and addition effect
  Y = rbinom(N, 1, mu)    # Case-control status
  re = sum(Y) - N * prevalence2
  return(re)
}

#### Main function to simulate binary phenotype (with individual-specific prevalence)
data.simu.binary.admix = function(N,
                                  N1,
                                  N2,
                                  prevalence1,
                                  prevalence2,
                                  gamma,
                                  gamma_GxE,
                                  g,               # genotype vector for N samples from an admixed population
                                  g1,
                                  g2,
                                  # seed1,
                                  # seed2,
                                  bVec = 0,                  
                                  ance.data)       # ancestry proportion vector 
{
  beta0 = uniroot(f1.binary, c(-100,100), N1 = N1, prevalence1 = prevalence1, gamma = gamma, g1 = g1, gamma_GxE = gamma_GxE, bVec = bVec)
  beta0 = beta0$root
  
  beta1 =  uniroot(f2.binary, c(-100,100), N1 = N1, N2 = N2, prevalence1 = prevalence1, 
                   prevalence2 = prevalence2, gamma = gamma, g1 = g1, g2 = g2, gamma_GxE = gamma_GxE, bVec = bVec)
  beta1 = beta1$root
  
  X1 = c(rep(0, N1), rep(1, N2))
  X2 = rnorm(N)
  X3 = rbinom(N, 1, 0.5)
  E = rnorm(N)              # environmental factor
  
  # betas = c(0.5, 0.5)     # Coefficient vector of fixed effects
  eta = beta0 + ance.data[,2]*beta1 + X2*0.5 + X3*0.5 + E*0.5 + gamma * g + gamma_GxE * g * E + bVec
  mu = exp(eta) / (1 + exp(eta))   # The probability being a case given the covariates, genotypes, and addition effect
  Y = rbinom(N, 1, mu)    # Case-control status
  out = data.frame(Y, ance.data[,2], X2, X3, E)
  colnames(out)= c("Y", "X1.new", "X2", "X3", "E")
  return(out)
}

#### simulate binary phenotype in the admixed population

data = data.simu.binary.admix(N, N1, N2, prevalence1 = prevalence_EUR, prevalence2 = prevalence_EAS, 
                              gamma = 0, gamma_GxE = 0, g = 0, g1 = 0, g2 = 0, ance.data = ancestry_Vec)

Phen.mat = data.frame(Prevalence = paste0("Prevalence in EUR = ", prevalence_EUR, sep=", ", "Prevalence in EAS = ", prevalence_EAS),
                      ID = paste0("IID-",1:N),
                      Y = data$Y,
                      X1.new = data$X1.new,
                      X2 = data$X2,
                      X3 = data$X3,
                      E = data$E,
                      PC1 = top10PCs[,1],
                      PC2 = top10PCs[,2],
                      PC3 = top10PCs[,3],
                      PC4 = top10PCs[,4])

Cova.mtx = Phen.mat %>% select(PC1, PC2, PC3, PC4, X2, X3) # a covariate matrix excluding the environmental factor E

#### fit null model
#### null model residuals

R = SPA_G_Get_Resid("binary",
                    Y~PC1+PC2+PC3+PC4+X2+X3+E,
                    data=Phen.mat,
                    pIDs=Phen.mat$ID,
                    gIDs=Phen.mat$ID)


#### calculate p-values using SPAGxEmixCCT
#### output p-values of SPAGxEmixCCT, SPAGxEmix, SPAGxEmix-Wald, and NormGxEmix
binary_admixed_res = SPAGxEmix_CCT(traits = "binary",
                                   Geno.mtx = GenoMat_admixed,
                                   R = R,
                                   E = Phen.mat$E,           # environmental factor
                                   Phen.mtx =  Phen.mat,     # include Y 
                                   Cova.mtx = Cova.mtx,      # other covariates (such as age, gender, and top PCs) excluding E
                                   topPCs = top4PCs)  


# we recommand using column of 'p.value.spaGxE.CCT.Wald.index.ance' to associate genotype with phenotypes
head(binary_admixed_res)

#### save results
write.csv(binary_admixed_res, file = "~/capsule/results/SPAGxEmixCCT_AdmixedPopulationAnalysis_binary_simulation_normal_envi_example_results.csv")







