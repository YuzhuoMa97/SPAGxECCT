# cd /gdata01/user/yuzhuoma/GxE/code_ocean/
# sbatch -J surv_typeIerror --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript SPAGxECCT_survival_typeIerror_simulation_normal_envi.R $SLURM_ARRAY_TASK_ID'

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

ERVec = c(0.01, 0.1, 0.5)
ER = ERVec[1] # event rate

############################################################

library(survival)
#### lower function to estimate beta0 given an event rate. Will be used in data.simu.surv().
f.surv.new = function(N,                  # Sample size
                      nSNP,               # number of SNPs with genetic main effects
                      event.rate,         # Event rate
                      lambda,             # Scale parameter
                      # gamma,            # Genetic effect
                      GMat,               # Genotype matrix
                      seed,               # random seed
                      bVec = 0)           # Additional effect, could be random effect. If bVec = 0 (default), then no additional effect is included.
{
  scale0 = lambda
  shape0 = 2
  set.seed(seed)
  gamma = runif(nSNP, -0.4, 0.4) # Genetic effect vector
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  betas = c(0.5, 0.5, 0.5)
  cens = rweibull(N, shape = 1, scale = 0.15)
  eps <- runif(N, 0, 1)
  time = (-log(eps) * exp(- betas[1]*X1 - betas[2]*X2 - betas[3]*E - GMat%*%gamma - bVec)) ^ (1/shape0) * scale0
  surv.time = pmin(time, cens)
  event = ifelse(time < cens, 1, 0)
  re = mean(event) - event.rate
  return(re)
}

#### Main function to simulate time-to-event phenotype
data.simu.surv.new = function(N,              # Sample size
                              nSNP,           # SNPs with genetic main effects
                              event.rate,     # Event rate
                              # gamma,        # Genetic effect
                              GMat,           # Genotype matrix
                              bVec = 0)       # Additional effect, could be random effect. If bVec = 0 (default), then no additional effect is included.
{ 
  scale0 = 0
  
  while (scale0 == 0) {
    seed = sample(1e9, 1);
    cat("random number seed is", seed)
    scale0 = uniroot(f.surv.new, c(-100000,100000), N = N, nSNP = nSNP,
                     event.rate = event.rate, GMat = GMat, seed = seed, bVec = bVec)
    scale0 = scale0$root
  }
  
  shape0 = 2
  set.seed(seed)
  gamma = runif(nSNP, -0.4, 0.4) # Genetic effect vector
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  betas = c(0.5, 0.5, 0.5)
  cens = rweibull(N, shape = 1, scale = 0.15) 
  eps <- runif(N, 0, 1)
  time = (-log(eps) * exp(- betas[1]*X1 - betas[2]*X2 - betas[3]*E - GMat%*%gamma - bVec)) ^ (1/shape0) * scale0
  surv.time = pmin(time, cens)
  event = ifelse(time < cens, 1, 0)
  out = data.frame(surv.time, event, X1, X2, E, bVec)
  return(out)
}



#### simulate time-to-event phenotype with certain event rate
data = data.simu.surv.new(N = N, nSNP = nSNP, event.rate = ER, GMat = GenoMat_1)
Phen.mtx = data.frame(ID = paste0("IID-",1:N),
                      event = data$event,
                      surv.time = data$surv.time,
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

R = SPA_G_Get_Resid("survival",
                    Surv(surv.time,event)~X1+X2+E,
                    data=Phen.mtx,
                    pIDs=Phen.mtx$ID,
                    gIDs=Phen.mtx$ID)

#### calculate p values for SNPs without marginal genetic effect 
survival_res_0 = SPAGxE_CCT(traits = "survival",              # time-to-event trait analysis
                            Geno.mtx = GenoMat_0,             # genotype vector
                            R = R,                            # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                            E = E,                            # environmental factor
                            Phen.mtx = Phen.mtx,              # include surv.time, event
                            Cova.mtx = Cova.mtx)              # covariate matrix excluding the environmental factor E

# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
head(survival_res_0)

#### calculate p values for SNPs with marginal genetic effect 
survival_res_1 = SPAGxE_CCT(traits = "survival",              # time-to-event trait analysis
                            Geno.mtx = GenoMat_1,             # genotype vector
                            R = R,                            # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                            E = E,                            # environmental factor
                            Phen.mtx = Phen.mtx,              # include surv.time, event
                            Cova.mtx = Cova.mtx)              # covariate matrix excluding the environmental factor E

# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
head(survival_res_1)

#### save results
write.csv(survival_res_0, file = "~/capsule/results/SPAGxECCT_survival_typeIerror_simulation_normal_envi_example_results_0.csv")
write.csv(survival_res_1, file = "~/capsule/results/SPAGxECCT_survival_typeIerror_simulation_normal_envi_example_results_1.csv")
