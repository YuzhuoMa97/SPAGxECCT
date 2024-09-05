# cd /gdata01/user/yuzhuoma/GxE/code_ocean/
# sbatch -J surv_admixed --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript SPAGxEmixCCT_AdmixedPopulationAnalysis_survival_simulation_normal_envi.R $SLURM_ARRAY_TASK_ID'

# args=commandArgs(TRUE)
# print(args)
# print(sessionInfo())
# n.cpu=as.numeric(args)

library(dplyr)
library(tidyr)
library(data.table)
library(survival)

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
colnames(ancestry_Vec) = c("AnceProportion_EUR", "AnceProportion_EAS")
rownames(ancestry_Vec) = paste0("IID-",1:N)

ER_EUR = 0.1  # event rate in EUR
ER_EAS = 0.01 # event rate in EAS

MAF_EUR = 0.3  # minor allele frequency in EUR
MAF_EAS = 0.01 # minor allele frequency in EAS

MAF_admixed = ancestry_Vec %*% c(ER_EUR, ER_EAS) # minor allele frequency in the admixed population
colnames(MAF_admixed) = "MAF"

GenoMat_admixed = matrix(rbinom(N*nSNP, 2, MAF_admixed), N, nSNP) # genotype matrix of the admixed population
colnames(GenoMat_admixed) = paste0("SNP-", 1:nSNP)
rownames(GenoMat_admixed) = paste0("IID-", 1:N)

top4PCs = top10PCs[,1:4] # top 4 PCs (used in SPAGxEmixCCT)


### sub-functions to simulate time-to-event phenotype
library(survival)
#### lower function to estimate beta0 given an event rate. Will be used in data.simu.surv().
#### simulate time-to-event phenotype for admixed population (e.g. EUR and EAS)

f1 = function(N1,
              event.rate1,
              lamda,
              gamma,
              g1,
              gamma_GxE,         # Marginal GxE effect
              bVec = 0)               
{
  scale0 = lamda
  shape0 = 2
  # set.seed(seed1)
  X21 = rnorm(N1)
  X31 = rbinom(N1, 1, 0.5)
  E1 = rnorm(N1)
  
  cens = rweibull(N1, shape=1, scale=0.15)
  eps <- runif(N1, 0, 1)
  time = (-log(eps)*exp(-X21*0.5-X31*0.5-E1*0.5-g1*gamma-g1*E1*gamma_GxE-bVec))^(1/shape0)*scale0
  surv.time = pmin(time, cens)
  event = ifelse(time<cens, 1, 0)
  re = mean(event) - event.rate1
  return(re)
}

f2 = function(N1,
              N2,
              event.rate1,
              event.rate2,
              b,
              gamma,
              g1,
              g2,
              gamma_GxE,
              bVec = 0)
  # seed1,
  # seed2)
{
  scale0 = uniroot(f1, c(-100,100), N1 = N1, event.rate1 = event.rate1, gamma = gamma, g1 = g1, gamma_GxE = gamma_GxE, bVec = bVec)
  scale0 = scale0$root
  shape0 = 2
  # set.seed(seed2)
  X1 = c(rep(0, N1), rep(1, N2))
  X11 = X1[c(1 : N1)]
  X12 = X1[c((N1 + 1) : N)]
  X22 = rnorm(N2)
  X32 = rbinom(N2, 1, 0.5)
  E2 = rnorm(N2)
  
  cens = rweibull(N2, shape=1, scale=0.15)
  eps <- runif(N2, 0, 1)
  time = (-log(eps)*exp(-X12*b-X22*0.5-X32*0.5-E2*0.5-g2*gamma-g2*E2*gamma_GxE-bVec))^(1/shape0)*scale0
  surv.time = pmin(time, cens)
  event = ifelse(time<cens, 1, 0)
  re = mean(event) - event.rate2
  return(re)
}

####  simulate time-to-event phenotype for population admixed with two subpopulations
####  event rates are different for different subjects (individual-specific event rates)
data.simu.surv.admix = function(N,
                                N1,
                                N2,
                                event.rate1,
                                event.rate2,
                                gamma,
                                gamma_GxE,
                                g,               
                                g1,
                                g2,
                                # seed1,
                                # seed2,
                                bVec = 0,                  
                                a)           # ancestry proportion vector
{
  scale0 = uniroot(f1, c(-100,100), N1 = N1, event.rate1 = event.rate1, gamma = gamma, g1 = g1, gamma_GxE = gamma_GxE, bVec = bVec)
  scale0 = scale0$root
  shape0 = 2
  b = uniroot(f2, c(-100,100), N1 = N1, N2 = N2, event.rate1 = event.rate1, 
              event.rate2 = event.rate2, gamma = gamma, g1 = g1, g2 = g2, gamma_GxE = gamma_GxE, bVec = bVec)
  b = b$root
  
  X1 = c(rep(0, N1), rep(1, N2))
  X2 = rnorm(N)
  X3 = rbinom(N, 1, 0.5)
  E = rnorm(N)                    # environmental factor
  
  # set.seed(1)
  cens = rweibull(N, shape=1, scale=0.15) 
  eps <- runif(N, 0, 1)
  time = (-log(eps)*exp(-a[,2]*b-X2*0.5-X3*0.5-E*0.5-g*gamma-g*E*gamma_GxE))^(1/shape0)*scale0
  surv.time = pmin(time, cens)
  event = ifelse(time<cens, 1, 0)
  out = data.frame(surv.time, event, cens, a[,2], X2, X3, E)
  colnames(out)= c("surv.time","event", "cens", "X1.new", "X2", "X3", "E")
  return(out)
}

#### simulate time-to-event phenotype in the admixed population

data = data.simu.surv.admix(N, N1, N2, event.rate1 = ER_EUR, event.rate2 = ER_EAS, 
                            gamma = 0, gamma_GxE = 0, g = 0, g1 = 0, g2 = 0, a = ancestry_Vec)

Phen.mat = data.frame(ER = paste0("ER in EUR = ", ER_EUR, sep=", ", "ER in EAS = ", ER_EAS),
                      ID = paste0("IID-",1:N),
                      event = data$event,
                      surv.time = data$surv.time,
                      X1.new = data$X1.new,
                      X2 = data$X2,
                      X3 = data$X3,
                      E = data$E,
                      PC1 = top10PCs[,1],
                      PC2 = top10PCs[,2],
                      PC3 = top10PCs[,3],
                      PC4 = top10PCs[,4])

Cova.mtx = Phen.mat %>% select(PC1, PC2, PC3, PC4, X2, X3) # a covariate matrix excluding the environmental factor E

### fit null model
### null model residuals

R = SPA_G_Get_Resid("survival",
                    Surv(surv.time,event)~PC1+PC2+PC3+PC4+X2+X3+E,
                    data=Phen.mat,
                    pIDs=Phen.mat$ID,
                    gIDs=Phen.mat$ID)


#### calculate p-values using SPAGxEmixCCT
#### output p-values of SPAGxEmixCCT, SPAGxEmix, SPAGxEmix-Wald, and NormGxEmix
survival_admixed_res = SPAGxEmix_CCT(traits = "survival",
                                     Geno.mtx = GenoMat_admixed,
                                     R = R,
                                     E = Phen.mat$E,                    # environmental factor
                                     Phen.mtx =  Phen.mat,     # include surv.time, event 
                                     Cova.mtx = Cova.mtx,      # other covariates (such as age, gender, and top PCs) excluding E
                                     topPCs = top4PCs)  


# we recommand using column of 'p.value.spaGxE.CCT.Wald.index.ance' to associate genotype with phenotypes
head(survival_admixed_res)

#### save results
write.csv(survival_admixed_res, file = "~/capsule/results/SPAGxEmixCCT_AdmixedPopulationAnalysis_survival_simulation_normal_envi_example_results.csv")







