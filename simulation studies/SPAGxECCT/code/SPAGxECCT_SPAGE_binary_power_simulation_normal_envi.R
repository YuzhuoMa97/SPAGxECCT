# cd /gdata01/user/yuzhuoma/GxE/code_ocean/
# sbatch -J binary_typeIerror --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript SPAGxECCT_SPAGE_binary_power_simulation_normal_envi.R $SLURM_ARRAY_TASK_ID'

# args=commandArgs(TRUE)
# print(args)
# print(sessionInfo())
# n.cpu=as.numeric(args)

N = 50000  # sample size
nSNP = 10  # number of SNPs with marginal GxE effects

prevalenceVec = c(0.01, 0.05, 0.2) 
MAFVec = c(0.3, 0.05, 0.005)
GammaVec = seq(0, 2, 0.2)

prevalence = prevalenceVec[1]        # prevalence
MAF = MAFVec[3]                      # minor allele frequencies of SNPs with marginal GxE effects
Gamma = GammaVec[8]                  # marginal GxE effect size


### sub-functions to simulate binary phenotype
#### lower function to estimate beta0 given an prevalence. Will be used in data.simu.binary().
f = function(N,
             prevalence,
             beta0,             # Intercept
             gamma,
             g,
             seed,
             bVec = 0)
{
  
  set.seed(seed)
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  
  eta = beta0 + X1*0.5 + X2*0.5 + E*0.5 + g*E*gamma + bVec
  mu = exp(eta) / (1 + exp(eta))   # The probability being a case given the covariates, genotypes, and addition effect
  Y = rbinom(N, 1, mu)    # Case-control status
  re = sum(Y) - N * prevalence
  
  return(re)
}

data.simu.binary = function(N,
                            prevalence,
                            gamma,
                            g,
                            bVec = 0)
{
  seed = sample(1e9, 1);
  cat("random number seed is", seed)
  beta0 = uniroot(f, c(-30,30), N = N, prevalence = prevalence,
                  gamma = gamma, g = g, seed = seed, bVec = bVec)
  beta0 = beta0$root
  set.seed(seed)
  X1 = rnorm(N)
  X2 = rbinom(N, 1, 0.5)
  E = rnorm(N)
  eta = beta0 + X1*0.5 + X2*0.5 + E*0.5 + g*E*gamma + bVec
  mu = exp(eta) / (1 + exp(eta))   # The probability being a case given the covariates, genotypes, and addition effect
  Y = rbinom(N, 1, mu)             # Case-control status
  out = data.frame(Y, X1, X2, E, bVec)
  return(out)
}



########## generate genotype matrix with marginal GxE effect beta_GxE!=0

MAFVec = runif(nSNP, MAF, MAF)
GenoMat = t(matrix(rbinom(N*nSNP, 2, MAFVec), nSNP, N))
rownames(GenoMat) = paste0("IID-",1:N)
colnames(GenoMat) = paste0("SNP-",1:nSNP)


#### calculate p values to evaluate powers of SPAGxECCT and SPAGE

#### install SPAGxECCT package
# library(devtools)  # author version: 2.4.5
# install_github("YuzhuoMa97/SPAGxECCT")
library(SPAGxECCT)
# ?SPAGxECCT  # manual of SPAGxECCT package

#### SPAGE functions
source("/gdata01/user/yuzhuoma/GxE/code/SPAGE/main.R")                            # SPAGE
source("/gdata01/user/yuzhuoma/GxE/code/SPAGE/fit-null.R")
source("/gdata01/user/yuzhuoma/GxE/code/SPAGE/data-simu.R")


binary_res_SPAGxECCT = c()
binary_res_SPAGE = c()

for (i in 1:nSNP) {
  print(i)
  G = GenoMat[,i]  # genotype to test
  
  data = data.simu.binary(N = N, prevalence = prevalence, gamma = Gamma, g = G) # simulate binary phenotype
  
  ### phenotype matrix
  Phen.mtx = data.frame(ID = paste0("IID-",1:N),
                        Y = data$Y,
                        X1 = data$X1,
                        X2 = data$X2,
                        E = data$E)
  
  ### fit genotype-independent (covariate-only) model
  ### null model residuals
  R = SPA_G_Get_Resid("binary",
                      Y~X1+X2+E,
                      data=Phen.mtx,
                      pIDs=Phen.mtx$ID,
                      gIDs=Phen.mtx$ID)
  
  #### calculate p values for SNPs with marginal GxE effect using SPAGxECCT
  binary_res_SPAGxECCT_oneSNP = SPAGxE_CCT(traits = "binary",                  # binary trait analysis
                                           Geno.mtx = as.matrix(G),            # genotype vector
                                           R = R,                              # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                                           E = Phen.mtx$E,                     # environmental factor
                                           Phen.mtx = Phen.mtx,                # include surv.time, event
                                           Cova.mtx = Phen.mtx[,c("X1","X2")]) # covariate matrix excluding the environmental factor E
  
  binary_res_SPAGxECCT = rbind(binary_res_SPAGxECCT, binary_res_SPAGxECCT_oneSNP)
  
  
  
  #### calculate p values for SNPs with marginal GxE effect using SPAGE
  
  obj.null.SPAGE = SPAGE_Null_Model(Y~X1+X2+E,
                                    subjectID = paste0("IID-",1:dim(Phen.mtx)[1]), 
                                    data = Phen.mtx, 
                                    out_type = "D")
  
  subjectID = paste0("IID-",1:dim(Phen.mtx)[1])
  Envn.mtx = as.matrix(Phen.mtx)[,"E",drop=FALSE]
  Envn.mtx = as.numeric(Envn.mtx)
  Envn.mtx = as.matrix(Envn.mtx)
  colnames(Envn.mtx) = "E"
  
  binary_res_SPAGE_oneSNP = SPAGE.one.SNP(g = G, obj.null = obj.null.SPAGE, Envn.mtx = Envn.mtx)
  binary_res_SPAGE = rbind(binary_res_SPAGE, binary_res_SPAGE_oneSNP)
}


rownames(binary_res_SPAGE) = c(1 : length(binary_res_SPAGE[,1]))
colnames(binary_res_SPAGE) = c("MAF","missing.rate","p.value.BetaG","p.value.spa-E", 
                               "p.value.norm-E", "p.value.Firth-E","Stat-E","Var-E", "z-E")

head(binary_res_SPAGE)

# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with binary phenotypes
head(binary_res_SPAGxECCT)




