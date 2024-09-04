# cd /gdata01/user/yuzhuoma/GxE/code_ocean/
# sbatch -J quantitative_localance --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript SPAGxEmixCCT_local_ancestry_quantitative_power_simulation_normal_envi.R $SLURM_ARRAY_TASK_ID'

# args=commandArgs(TRUE)
# print(args)
# print(sessionInfo())
# n.cpu=as.numeric(args)

library(dplyr)
library(data.table)
library(SPAGxECCT)

N = 10000 # sample size
nSNP = 10 # number of SNPs with marginal GxE effects

prevalence = 0.2 # prevalence in the admixed population
MAF_ance1 = 0.1  # minor allele frequencies of SNPs with marginal GxE effects in ancestry 1
MAF_ance2 = 0.05 # minor allele frequencies of SNPs with marginal GxE effects in ancestry 2

Gamma_ance1 = 0.2 # marginal GxE effect size of ancestry 1
Gamma_ance2 = 0.5   # marginal GxE effect size of ancestry 2

#### Main function to simulate quantitative phenotype
data.simu.quantitative= function(N,                # Sample size
                                 beta0,             # intercept term
                                 gamma1,            # Genetic effect
                                 gamma2,
                                 g1,
                                 g2,                # Genotype vector
                                 bVec = 0,
                                 seed)         # Additional effect, could be random effect. If bVec = 0 (default), then no additional effect is included.
{
  set.seed(seed)
  Cov1 = rnorm(N)           # Covariate 1
  Cov2 = rbinom(N, 1, 0.5)  # Covariate 2
  
  E = rnorm(N)            # Environmental factor
  
  betas = c(0.5, 0.5, 0.5)     # Coefficient vector of fixed effects
  eta = beta0 + betas[1] * Cov1 + betas[2] * Cov2 + betas[3] * E + gamma1 * g1 * E + gamma2 * g2 * E + bVec
  epsilon = rnorm(N)
  y = eta + epsilon
  out = data.frame(y, Cov1, Cov2, E)
  return(out)
}



#### generate global ancestry

theta = 0.5                                 # expected proportion of ancestry 
sigma = 0.125                               # sd 
MAFVec1 = runif(nSNP, MAF_ance1, MAF_ance1) # MAF vector of ancestry 1
MAFVec2 = runif(nSNP, MAF_ance2, MAF_ance2) # MAF vector of ancestry 2

set.seed(528)

alpha = rnorm(N, theta, sigma)  # individual global ancestry
alpha = pmin(1, pmax(alpha, 0)) # forced to be in [0,1]

# For each individual, draw a local ancestry count of ancestry 2
l = lapply(1:N, function(i){
  rbinom(nSNP, 2, alpha[i]) 
}) %>% do.call("rbind",.) %>% as.matrix()

G1 = lapply(1:N, function(i){
  g = rbinom(nSNP, 2-l[i,] , MAFVec1)
}) %>% do.call("rbind",.) %>% as.matrix() # risk allele from ancestry 1

G2 = lapply(1:N, function(i){
  g = rbinom(nSNP, l[i,] , MAFVec2)
})%>% do.call("rbind",.) %>% as.matrix() # risk allele from ancestry 2

rownames(G1) =rownames(G2) = rownames(l) = paste0("IID-",1:N)
colnames(G1) =colnames(G2) = colnames(l) =  paste0("rs",1:nSNP)
G = G1 + G2

# combine genotypes
combG1 = apply(t(G1)-2 * MAFVec1, 2, sum) 
combG2 = apply(t(G2)-2 * MAFVec2, 2, sum)


Geno.mtx = G
Geno.mtx.ance1 = G1
Geno.mtx.ance2 = G2
Comb.Geno.mtx.ance1 = combG1
Comb.Geno.mtx.ance2 = combG2

# save as hapcount and dosage files separately
anc1.hapcount =  data.frame(CHROM = 1, POS = 1:nSNP, ID = colnames(l), REF = 1:nSNP, ALT = 1:nSNP+1) %>% cbind(2-t(l))
anc2.hapcount =  data.frame(CHROM = 1, POS = 1:nSNP, ID = colnames(l), REF = 1:nSNP, ALT = 1:nSNP+1) %>% cbind(t(l))


haplo.mtx.ance1 = t(anc1.hapcount[,-1:-5]) # matrix of number of haplotypes (local ancestry counts) at the locus in ancestry 1
haplo.mtx.ance2 = t(anc2.hapcount[,-1:-5]) # matrix of number of haplotypes (local ancestry counts) at the locus in ancestry 2

colnames(haplo.mtx.ance1) = colnames(Geno.mtx.ance1)
colnames(haplo.mtx.ance2) = colnames(Geno.mtx.ance2)


#### simulate quantitative phenotypes with heterogeneous GxE effects across ancestries

seed = 10000 # set seed

pheno = data.simu.quantitative(N=N, gamma1 = Gamma_ance1, gamma2 = Gamma_ance2
                               , g1 = Comb.Geno.mtx.ance1, g2 = Comb.Geno.mtx.ance2
                               , beta0 = 0, seed = seed) %>% mutate(ID = paste0("IID-",1:N)) %>% select(ID, y, Cov1, Cov2, E) %>% rename(IID=ID)


load("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/topPCs/top10PCs.RData") # load in top 10 PCs corresponding to the global ancestry alpha

Pheno.mtx = cbind(pheno, top10PCs[,1:4]) %>% rename(PC1="Comp.1", PC2="Comp.2", PC3="Comp.3", PC4="Comp.4") %>% mutate(Y = y)
Cova.mtx = Pheno.mtx %>% select(PC1, PC2, PC3, PC4, Cov1, Cov2) # a covariate matrix excluding the environmental factor E and local ancestry counts (the number of haplotypes).
E = Pheno.mtx$E # a numeric environmental factor with each element as an environmental factor value of an individual.

#### Cova.haplo.mtx.list
#### a list with each element as a matrix of local ancestry counts (the number of haplotypes) of an ancestry (i.e., haplo.mtx of an ancestry).
Cova.haplo.mtx.list = list(haplo.mtx.ance1 = haplo.mtx.ance1,
                           haplo.mtx.ance2 = haplo.mtx.ance2)

### fit null model
obj.linear = lm(y ~ Cov1 + Cov2  + E + PC1 + PC2 + PC3 + PC4, data = Pheno.mtx)
resid = obj.linear$residuals

### calculate p values for ancestry 1 using SPAGxEmixCCT_local
quantitative_res_ance1 = SPAGxEmixCCT_localance(traits = "quantitative",
                                                Geno.mtx = Geno.mtx.ance1,
                                                R = resid,
                                                haplo.mtx = haplo.mtx.ance1,
                                                E = E,
                                                Phen.mtx = Pheno.mtx,
                                                Cova.mtx = Cova.mtx,
                                                Cova.haplo.mtx.list = Cova.haplo.mtx.list)


colnames(quantitative_res_ance1) = c("Marker", "MAF.ance1","missing.rate.ance1",
                                     "Pvalue.spaGxE.ance1","Pvalue.spaGxE.Wald.ance1", "Pvalue.spaGxE.CCT.Wald.ance1",
                                     "Pvalue.normGxE.ance1", "Pvalue.betaG.ance1",
                                     "Stat.betaG.ance1","Var.betaG.ance1","z.betaG.ance1")

### calculate p values for ancestry 2 using SPAGxEmixCCT_local
quantitative_res_ance2 = SPAGxEmixCCT_localance(traits = "quantitative",
                                                Geno.mtx = Geno.mtx.ance2,
                                                R = resid,
                                                haplo.mtx = haplo.mtx.ance2,
                                                E = E,
                                                Phen.mtx = Pheno.mtx,
                                                Cova.mtx = Cova.mtx,
                                                Cova.haplo.mtx.list = Cova.haplo.mtx.list)


colnames(quantitative_res_ance2) = c("Marker", "MAF.ance2","missing.rate.ance2",
                                     "Pvalue.spaGxE.ance2","Pvalue.spaGxE.Wald.ance2", "Pvalue.spaGxE.CCT.Wald.ance2",
                                     "Pvalue.normGxE.ance2", "Pvalue.betaG.ance2",
                                     "Stat.betaG.ance2","Var.betaG.ance2","z.betaG.ance2")

#### merge data frame
quantitative.res = merge(quantitative_res_ance1, quantitative_res_ance2)

# we recommand using column of 'p.value.spaGxE.CCT.Wald.index.ance' to associate genotype with phenotypes
head(quantitative.res)



