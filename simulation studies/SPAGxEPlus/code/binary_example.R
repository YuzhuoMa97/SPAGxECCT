# cd /gdata01/user/yuzhuoma/GxE/F/GxE/code/binary/typeIerror/4mem/

setwd("/gdata01/user/yuzhuoma/GxE/F/GxE/SPAGxEPlus_test_data/")

library(data.table)
library(lme4)
source("/gdata01/user/yuzhuoma/SPAGxEmixPlus/simulation_admix_relateness_test/code/typeIerror/binary/pheno_hetero/4mem/20240916_SPAGxECCT_pcakage.R")


prevalence = 0.01
True.MAF = 0.01
Pheno.ID = 1

#### load genotype data (10000 independent SNPs)
load(paste0("/gdata01/user/yuzhuoma/GxE/F/GxE/data/binary/typeIerror/4mem/genotype/MAF-", True.MAF, "/GenoMat.Rdata"))

#### load phenotype data
load(paste0("/gdata01/user/yuzhuoma/GxE/F/GxE/data/binary/typeIerror/4mem/phenotype/tau1/prevalence-", prevalence, "/Phenotype-", Pheno.ID,".Rdata"))

load("/gdata01/user/yuzhuoma/GxE/F/GxE/data/binary/typeIerror/4mem/GRM/GRM_nSNP100000.RData")
sparseGRM = data.table::fread("/gdata01/user/yuzhuoma/GxE/F/GxE/data/binary/typeIerror/4mem/GRM/SparseGRM.txt")

library(Matrix)
GRM[GRM < 0.05] <- 0



####################

N = dim(Phen.mtx)[1]    # sample size
# nSNP = 100            # number of SNPs
Cov1 = Phen.mtx$Cov1        # Covariate 1
Cov2 = Phen.mtx$Cov2        # Covariate 2
E = Phen.mtx$E        # E
Cova.mtx = Phen.mtx[,c("Cov1","Cov2")]


# Gets NULL model residuals for SPA-G
SPA_G_Get_Resid = function(traits="survival/binary",
                           formula=NULL,
                           data=NULL,
                           pIDs=NULL,
                           gIDs=NULL,
                           range=c(-100,100),
                           length.out = 10000,
                           ...)

{
  if(traits=="survival"){
    Call = match.call()
    ### Fit a Cox model
    obj.coxph = coxph(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.coxph, range)

    ### Get the covariate matrix to adjust for genotype
    resid = obj.coxph$residuals

    re = resid
  }
  else if(traits=="binary"){
    Call = match.call()
    ### Fit a logistic model
    obj.logistic = glm(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.logistic, range)

    ### Get the covariate matrix to adjust for genotype
    mu = obj.logistic$fitted.values
    resid = obj.logistic$y - mu
    re = resid
  }
  return(re)
}

check_input = function(pIDs, gIDs, obj, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = obj$residuals
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")

  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)

  return(p2g)
}

check_input_Resid = function(pIDs, gIDs, R, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = R
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")

  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)

  return(p2g)
}


###########################################################################

R = SPA_G_Get_Resid(traits = "binary",
                    Y~Cov1+Cov2+E,family=binomial(link="logit"),
                    data=Phen.mtx,
                    pIDs=Phen.mtx$IID,
                    gIDs=rownames(GenoMat))

ResidMat = data.frame(SubjID = Phen.mtx$IID, Resid = R)


# S1 = sum(g*R)                    # test statistic for marginal genetic effect
# S2 = sum(g*E*R)

# RE = R*E
#
# R_GRM_R = as.numeric(t(R) %*% GRM %*% R)
# R_GRM_RE = as.numeric(t(R) %*% GRM %*% RE)
# lambda = R_GRM_RE/R_GRM_R   # lambda = Cov/VarS1
# R.new = (E - lambda) * R       # new residuals
# R.new_GRM_R.new = as.numeric(t(R.new) %*% GRM %*% R.new)

library(dplyr)


obj.SPAGxE_Plus_Nullmodel = SPAGxE_Plus_Nullmodel(traits = "binary",
                                                  Y~Cov1+Cov2+E,family=binomial(link="logit"),
                                                  data=Phen.mtx,
                                                  pIDs=Phen.mtx$IID,
                                                  gIDs=rownames(GenoMat),
                                                  sparseGRM = sparseGRM,
                                                  E = E)



# g = GenoMat[,1]

binary.res = SPAGxE_Plus(Geno.mtx = GenoMat,
                         # ResidMat = ResidMat,
                         E = E,
                         Phen.mtx = Phen.mtx,
                         Cova.mtx = Cova.mtx,
                         sparseGRM = sparseGRM,
                         obj.SPAGxE_Plus_Nullmodel = obj.SPAGxE_Plus_Nullmodel)






summary(binary.res)

binary.res = as.data.frame(binary.res)



data_fig = binary.res %>% mutate(Pvalue = p.value.spaGxE)
# data_fig = binary.res %>% mutate(Pvalue = p.value.normGxE)

data_fig$est.Pvalue = rank(data_fig$Pvalue)/(nrow(data_fig)+1)



library(ggplot2)


p0 = ggplot(data_fig, aes(-log10(est.Pvalue), -log10(Pvalue))) +
  geom_point(alpha = 0.7) + geom_abline(slope = 1, intercept = 0)  # + facet_grid(Prev~MAF, scales="fixed") +
  labs(x="Expected(-log10P)",y="Observed(-log10P)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "bottom") +

  theme(
    plot.title = element_text(hjust = 0.5),
    # legend.title = element_blank(),
    # legend.position="none",
    # panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    legend.text.align = 0  # Left-align the legend text
    # panel.grid.minor = element_blank()
    # legend.title=element_blank()
  )


ggsave(p0, filename = paste0("test2-prevalence-",prevalence,"-MAF-", True.MAF,"2024-09-21.png"), width = 6, height = 6,
       path = "/gdata01/user/yuzhuoma/SPAGxEmixPlus/")
















