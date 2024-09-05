# cd /gdata01/user/yuzhuoma/GxE/code_ocean/
# sbatch -J surv_typeIerror --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript SPAGxECCT_survival_power_simulation_normal_envi.R $SLURM_ARRAY_TASK_ID'

# args=commandArgs(TRUE)
# print(args)
# print(sessionInfo())
# n.cpu=as.numeric(args)

N = 50000  # sample size
nSNP = 10  # number of SNPs with marginal GxE effects

ERVec = c(0.01, 0.1, 0.5)
MAFVec = c(0.3, 0.05, 0.005)
GammaVec = seq(0, 2, 0.2)

ER = ERVec[2]        # event rate
MAF = MAFVec[2]      # minor allele frequencies of SNPs with marginal GxE effects
Gamma = GammaVec[3]  # marginal GxE effect size


### sub-functions to simulate time-to-event phenotype
library(survival)
#### lower function to estimate beta0 given an event rate. Will be used in data.simu.surv().
f = function(N,
             event.rate,
             lambda,
             gamma,
             g,
             seed,
             bVec = 0)
{
  scale0 = lambda
  shape0 = 2
  set.seed(seed)
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  cens = rweibull(N, shape=1, scale=0.15)
  eps <- runif(N, 0, 1)
  time = (-log(eps)*exp(- X1*0.5 - X2*0.5 - E*0.5 - g*E*gamma - bVec))^(1/shape0)*scale0
  surv.time = pmin(time, cens)
  event = ifelse(time<cens, 1, 0)
  re = mean(event) - event.rate
  return(re)
}

data.simu.surv = function(N,
                          event.rate,
                          gamma,
                          g,
                          bVec = 0)
{
  seed = sample(1e9, 1);
  cat("random number seed is", seed)
  scale0 = uniroot(f, c(-100000,100000), N = N, event.rate = event.rate,
                   gamma = gamma, g = g, seed = seed, bVec = bVec)
  scale0 = scale0$root
  shape0 = 2
  set.seed(seed)
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  cens = rweibull(N, shape=1, scale=0.15)
  eps <- runif(N, 0, 1)
  time = (-log(eps)*exp(- X1*0.5 - X2*0.5 - E*0.5 - g*E*gamma - bVec))^(1/shape0)*scale0
  surv.time = pmin(time, cens)
  event = ifelse(time<cens, 1, 0)
  out = data.frame(surv.time, event, X1, X2, E, bVec)
  return(out)
}



########## generate genotype matrix with marginal GxE effect beta_GxE!=0

MAFVec = runif(nSNP, MAF, MAF)
GenoMat = t(matrix(rbinom(N*nSNP, 2, MAFVec), nSNP, N))
rownames(GenoMat) = paste0("IID-",1:N)
colnames(GenoMat) = paste0("SNP-",1:nSNP)


#### calculate p values to evaluate powers of SPAGxECCT

#### install and library SPAGxECCT package
library(remotes)                            # remotes library requires less dependency packages than devtools
install_github("YuzhuoMa97/SPAGxECCT")      # The INSTALL_opts is required in Windows OS.
library(SPAGxECCT)
# ?SPAGxECCT                                # manual of SPAGxECCT package
# or source("~/capsule/code/R/SPAGxECCT.R")

survival_res = c()
for (i in 1:nSNP) {
  print(i)
  G = GenoMat[,i]  # genotype to test
  
  data = data.simu.surv(N = N, event.rate = ER, gamma = Gamma, g = G) # simulate time-to-event phenotype
  
  ### phenotype matrix
  Phen.mtx = data.frame(ID = paste0("IID-",1:N),
                        event = data$event,
                        surv.time = data$surv.time,
                        X1 = data$X1,
                        X2 = data$X2,
                        E = data$E)
  
  ### fit genotype-independent (covariate-only) model
  ### null model residuals
  R = SPA_G_Get_Resid("survival",
                      Surv(surv.time,event)~X1+X2+E,
                      data=Phen.mtx,
                      pIDs=Phen.mtx$ID,
                      gIDs=Phen.mtx$ID)
  
  #### calculate p values for SNPs with marginal GxE effect 
  survival_res_oneSNP = SPAGxE_CCT(traits = "survival",                # time-to-event trait analysis
                                   Geno.mtx = as.matrix(G),            # genotype vector
                                   R = R,                              # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                                   E = Phen.mtx$E,                     # environmental factor
                                   Phen.mtx = Phen.mtx,                # include surv.time, event
                                   Cova.mtx = Phen.mtx[,c("X1","X2")]) # covariate matrix excluding the environmental factor E
  
  survival_res = rbind(survival_res, survival_res_oneSNP)
}

# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
head(survival_res)

#### save results
write.csv(survival_res, file = "~/capsule/results/SPAGxECCT_survival_power_simulation_normal_envi_example_results.csv")


