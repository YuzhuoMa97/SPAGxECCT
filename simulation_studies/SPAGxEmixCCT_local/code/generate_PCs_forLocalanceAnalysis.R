# cd /gdata01/user/yuzhuoma/SPA-G/tractor/code/typeIerror/
# sbatch -J PCA --mem=4000M -t 1-0:0 --array=1-1 -o log/%A_%a.log --wrap='Rscript generate_PCs_forLocalanceAnalysis.R $SLURM_ARRAY_TASK_ID'

args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
##########################################
library(dplyr)
library(data.table)



theta = 0.5 # expected proportion of ancestry 
sigma = 0.125 #sd 
nSNP = 100000 # num of SNP
N = 10000 # sample size
MAFVec1 = runif(nSNP, 0.2, 0.5) # MAF vector of ancestry 1
MAFVec2 = runif(nSNP, 0.2, 0.5) # MAF vector of ancestry 2


# load("/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/global_ancestry/global_ancestry.RData")

theta = 0.5 #expected proportion of ancestry 
sigma = 0.125 #sd 
nSNP = 1000 # num of SNP
N = 10000 # sample size

set.seed(528)

alpha = rnorm(N, theta, sigma)  # individual global ancestry
alpha = pmin(1, pmax(alpha, 0)) # forced to be in [0,1]

# For each individual, draw a local ancestry count of ancestry 2
l_PC = lapply(1:N, function(i){
  rbinom(nSNP, 2, alpha[i]) 
}) %>% do.call("rbind",.) %>% as.matrix()

G1_PC = lapply(1:N, function(i){
  g = rbinom(nSNP, 2-l_PC[i,] , MAFVec1)
}) %>% do.call("rbind",.) %>% as.matrix() # risk allele from ancestry 1

G2_PC = lapply(1:N, function(i){
  g = rbinom(nSNP, l_PC[i,] , MAFVec2)
})%>% do.call("rbind",.) %>% as.matrix() # risk allele from ancestry 2

rownames(G1_PC) =rownames(G2_PC) = rownames(l_PC) = paste0("IID-",1:N)
colnames(G1_PC) =colnames(G2_PC) = colnames(l_PC) =  paste0("rs",1:nSNP)
G_PC = G1_PC + G2_PC


##########################################################
GRM_PCA = function(Geno.mtx, N, nSNP)
{
  p.est0 = (1/2) * colMeans(Geno.mtx)
  p.est0.mtx = matrix(0, N, nSNP)
  for (i in 1:N) {
    p.est0.mtx[i,] = p.est0
  }
  Z = (Geno.mtx - 2 * p.est0.mtx)/sqrt(2 * p.est0.mtx * (1 - p.est0.mtx))
  GRM = Z %*% t(Z) / nSNP
  PC = princomp(covmat = GRM)
  top10_PC = PC$loadings[,1:10]
  return(top10_PC)
}



top10PCs = GRM_PCA(Geno.mtx = G_PC, N = N, nSNP = nSNP)
save(top10PCs, file = "/gdata01/user/yuzhuoma/SPA-G/tractor/data/typeIerror_v1/topPCs/top10PCs.RData")
