#' Fits a genotype-independent model
#'
#' Fits a null linear regression model for quantitative traits, a null logistic regression model for binary traits ,or a null Cox proportional hazards model for time-to-event traits and then calculates residuals under a genotype-independent model.
#' @param traits a character value corresponding to phenotype. It should be "binary" for binary phenotype, "survival" for time-to-event phenotype, and "others" for other types of phenotype.
#' @param formula a formula to be passed to function lm(), glm(), or coxph(). For more details, please refer to package survival.
#' @param data a data.frame in which to interpret the variables named in the formula
#' @param pIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order in the formula.
#' @param gIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order of the Geno.mtx (i.e. the input of the function SPAGxE_CCT()).
#' @param range a two-element numeric vector (default: c(-100,100)) to specify the domain of the empirical CGF.
#' @param length.out a positive integer (default: 9999) for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.
#' @param ... Other arguments passed to function glm() or coxph(). For more details, please refer to package survival.
#' @return Residuals after fitting a genotype-independent (covariate-only) model.
#' @examples
#' # example 1  time-to-event phenotype
#' # Simulation phenotype and genotype
#' N = 10000
#' nSNP = 100
#' MAF = 0.1
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       event=rbinom(N,1,0.5),
#'                       surv.time=runif(N),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5),
#'                       E = rnorm(N))
#'
#' Cova.mtx = Phen.mtx[,c("Cov1","Cov2")]
#'
#' E = Phen.mtx$E
#'
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#'
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#'
#' R = SPA_G_Get_Resid("survival",
#'                     Surv(surv.time,event)~Cov1+Cov2+E,
#'                     data=Phen.mtx,
#'                     pIDs=Phen.mtx$ID,
#'                     gIDs=paste0("IID-",1:N))
#' @export
#' @import survival

SPA_G_Get_Resid = function(traits="survival/binary/quantitative",
                           formula=NULL,
                           data=NULL,
                           pIDs=NULL,
                           gIDs=NULL,
                           range=c(-100,100),
                           length.out = 10000,
                           ...)

{
  if(traits=="quantitative"){
    Call = match.call()
    ### Fit a linear model
    obj.linear = lm(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.linear, range)

    ### Get the covariate matrix to adjust for genotype
    resid = obj.linear$residuals

    re = resid
  }

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




#' A scalable and accurate framework for large-scale genome-wide gene-environment interaction (GxE) analysis.
#'
#' A scalable and accurate analysis framework for a large-scale genome-wide SPAGxECCT implementation of gene-environmental interaction (GxE) analyses of quantitative traits, binary traits, time-to-event traits, and ordinal categorical traits.
#' @param GenoFile a character of genotype file. Currently, we support two formats of genotype input including PLINK and BGEN. Other formats such as VCF will be added later. Please refer to \code{?GRAB.ReadGeno} for more details.
#' @param GenoFileIndex additional index file(s) corresponding to GenoFile. Please refer to \code{?GRAB.ReadGeno} for more details.
#' @param control a list of parameters to decide which markers to extract. Please refer to \code{?GRAB.ReadGeno} for more details.
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missing genotypes should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param R model residuals after fitting a genotype-independent model (i.e., a covariate-only model in which marginal genetic effect and GxE effect are 0)
#' @param E a numeric environmental factor with each element as an environmental factor value of an individual.
#' @param Phen.mtx phenotype dataframe at least including three columns of ID, surv.time and event for time-to-event trait analysis, two columns of ID and linear phenotype Y for linear trait analysis, two columns of ID and binary phenotype Y for binary trait analysis, or two columns of ID and ordinal categorical phenotype Y for ordinal categorical trait analysis.
#' @param Cova.mtx a covariate matrix excluding the environmental factor E.
#' @param epsilon a numeric value (default: 0.001) to specify the p-value cutoff for betaG estimation. Please see details for more information.
#' @param Cutoff a numeric value (default: 2) to specify the standard deviation cutoff to be used.
#'               If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation,
#'               otherwise, its p value is calculated based on a saddlepoint approximation.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e., 2MAF).
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#'                       Any variant with missing rate higher than this cutoff will be excluded from the analysis.
#' @param min.maf a numeric value (default: 0.001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param G.model model type
#' @details To run SPAGxECCT, the following two steps are required:
#' \itemize{
#'   \item Step 1: Use function SPA_G_Null_Model() or other functions to fit a genotype-independent (covariate-only) model to get residuals under a genotype-independent (covariate-only) model.
#'   \item Step 2: Use function SPAGxE_CCT() to calculate p value for each genetic variant to conduct a GxE analysis.
#' }
#'
#'
#'
#' SPAGxE_CCT() uses a hybrid strategy with both saddlepoint approximation and normal distribution approximation.
#' Generally speaking, saddlepoint approximation is more accurate than, but a little slower than, the traditional normal distribution approximation.
#' Hence, when the score statistic is close to 0 (i.e., p-values are not small), we use the normal distribution approximation.
#' And when the score statistic is far away from 0 (i.e., p-values are small), we use the saddlepoint approximation.
#' Argument 'Cutoff' is to specify the standard deviation cutoff.
#'
#' To calibrate the score statistics, SPAGxE_CCT() uses martingale residuals which are calculated via R package survival for time-to-event trait analysis, residuals from glm() for binary trait analysis, resuduals from lm() for quantitative trait analysis, and residuals from clm() for ordinal categorical trait analysis via R package ordinal.
#' All extentions (such as strata, ties, left-censoring) supported by package survival could also be used in SPAGxE_CCT().
#' Time-varying covariates are also supported by splitting each subject into several observations.
#'
#' Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
#' To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx.
#'
#' @return an R matrix with the following columns
#' \item{MAF}{Minor allele frequencies}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.spaGxE}{p value from SPAGxE method}
#' \item{p.value.spaGxE.Wald}{p value from SPAGxE_Wald method}
#' \item{p.value.spaGxE.CCT.Wald}{p value (recommanded) from SPAGxE_CCT method}
#' \item{p.value.normGxE}{p value from NormGxE method (based on a normal distribution approximation)}
#' \item{p.value.betaG}{p value of the marginal genetic effect based on a normal distribution approximation}
#' \item{Stat.betaG}{score statistics testing for marginal genetic effect}
#' \item{Var.betaG}{estimated variances of the score statistics testing for marginal genetic effect}
#' \item{z.betaG}{z values (using Var1) corresponding to the score statistics testing for marginal genetic effect}
#' @examples
#'
#' # example 1  time-to-event phenotype (genotype data input provided in the R matrix format)
#' library(SPAGxECCT)
#' # Simulate phenotype and genotype
#' N = 10000   # sample size
#' nSNP = 100  # number of SNPs
#' MAF = 0.1   # minor allele frequency
#'
#' # phenotype data
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       event=rbinom(N,1,0.5),
#'                       surv.time=runif(N),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5),
#'                       E = rnorm(N))
#'
#' Cova.mtx = Phen.mtx[,c("Cov1","Cov2")]         # covariates dataframe excluding environmental factor E
#' E = Phen.mtx$E                                 # environmental factor E
#'
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP) # genotype matrix
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#'
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#'
#' # fit a genotype-independent model
#'
#' R = SPA_G_Get_Resid("survival",
#'                     Surv(surv.time,event)~Cov1+Cov2+E,
#'                     data=Phen.mtx,
#'                     pIDs=Phen.mtx$ID,
#'                     gIDs=paste0("IID-",1:N))
#'
#' # calculate p values
#'
#' survival.res = SPAGxE_CCT(traits = "survival",                     # trait type
#'                           Geno.mtx = Geno.mtx,                     # genotype matrix
#'                           R = R,                                   # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
#'                           E = E,                                   # environmental factor
#'                           Phen.mtx = Phen.mtx,                     # phenotype dataframe including surv.time
#'                           Cova.mtx = Cova.mtx)                     # a covariate matrix excluding the environmental factor E
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
#' head(survival.res)
#'
#'
#' # example 2  binary phenotype (genotype data input provided in the R matrix format)
#' library(SPAGxECCT)
#' # Simulate phenotype and genotype
#' N = 10000   # sample size
#' nSNP = 100  # number of SNPs
#' MAF = 0.1   # minor allele frequency
#'
#' # phenotype data
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       Y = rbinom(N,1,0.5),
#'                       Cov1 = rnorm(N),
#'                       Cov2 = rbinom(N,1,0.5),
#'                       E = rnorm(N))
#'
#' Cova.mtx = Phen.mtx[,c("Cov1","Cov2")]           # covariates dataframe excluding environmental factor E
#' E = Phen.mtx$E                                   # environmental factor E
#'
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)   # genotype matrix
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#'
#' # fit a genotype-independent model
#'
#' R = SPA_G_Get_Resid(traits = "binary",
#'                     glm(formula = Y ~ Cov1+Cov2+E, data = Phen.mtx, family = "binomial"),
#'                     data = Phen.mtx,
#'                     pIDs = Phen.mtx$ID,
#'                     gIDs = paste0("IID-",1:N))
#'
#' # calculate p values
#'
#' binary.res = SPAGxE_CCT(traits = "binary",                       # trait type
#'                         Geno.mtx = Geno.mtx,                     # genotype matrix
#'                         R = R,                                   # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
#'                         E = E,                                   # environmental factor
#'                         Phen.mtx = Phen.mtx,                     # phenotype dataframe including Y
#'                         Cova.mtx = Cova.mtx)                     # a covariate matrix excluding the environmental factor E
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with binary phenotypes
#' head(binary.res)
#'
#'
#' # example 3  quantitative phenotype (genotype data input provided in the R matrix format)
#' library(SPAGxECCT)
#' # Simulate phenotype and genotype
#' N = 10000   # sample size
#' nSNP = 100  # number of SNPs
#' MAF = 0.1   # minor allele frequency
#'
#' # phenotype data
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       Y=rnorm(N,1,0.5),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5),
#'                       E = rnorm(N))
#'
#' Cova.mtx = Phen.mtx[,c("Cov1","Cov2")]            # covariates dataframe excluding environmental factor E
#' E = Phen.mtx$E                                    # environmental factor E
#'
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)    # genotype matrix
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#'
#' # fit a genotype-independent model
#'
#' R = SPA_G_Get_Resid("quantitative",
#'                     lm(formula = Y ~ Cov1+Cov2+E, data = Phen.mtx),
#'                     data=Phen.mtx,
#'                     pIDs=Phen.mtx$ID,
#'                     gIDs=paste0("IID-",1:N))
#'
#' # calculate p values
#'
#' quantitative.res = SPAGxE_CCT(traits = "quantitative",          # trait type
#'                               Geno.mtx = Geno.mtx,              # genotype matrix
#'                               R = R,                            # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
#'                               E = E,                            # environmental factor
#'                               Phen.mtx = Phen.mtx,              # phenotype dataframe including Y
#'                               Cova.mtx = Cova.mtx)              # a covariate matrix excluding the environmental factor E
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with binary phenotypes
#' head(quantitative.res)
#'
#' # example 4  analysis of time-to-event phenotype (genotype data input provided in PLINK file format)
#' library(SPAGxECCT)
#' # Simulate phenotype and genotype
#' N = 10000  # sample size
#'
#' # phenotype data
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       event=rbinom(N,1,0.5),
#'                       surv.time=runif(N),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5),
#'                       E = rnorm(N))
#'
#' Cova.mtx = Phen.mtx[,c("Cov1","Cov2")]                                     # covariates dataframe excluding environmental factor E
#' E = Phen.mtx$E                                                             # environmental factor E
#' GenoFile = system.file("", "GenoMat_SPAGxE.bed", package = "SPAGxECCT")    # PLINK format for genotype data
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#'
#' # fit a genotype-independent model
#'
#' R = SPA_G_Get_Resid("survival",
#'                     Surv(surv.time,event)~Cov1+Cov2+E,
#'                     data=Phen.mtx,
#'                     pIDs=Phen.mtx$ID,
#'                     gIDs=paste0("IID-",1:N))
#'
#' # calculate p values
#'
#' survival.res = SPAGxE_CCT(traits = "survival",                     # trait type
#'                           GenoFile = GenoFile,                     # a character of genotype file
#'                           R = R,                                   # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
#'                           E = E,                                   # environmental factor
#'                           Phen.mtx = Phen.mtx,                     # phenotype dataframe
#'                           Cova.mtx = Cova.mtx)                     # a covariate matrix excluding the environmental factor E
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
#' head(survival.res)
#'
#' # example 5  analysis of time-to-event phenotype (genotype input using BGEN file format)
#' library(SPAGxECCT)
#' # Simulate phenotype and genotype
#' N = 10000  # sample size
#'
#' # phenotype data
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       event=rbinom(N,1,0.5),
#'                       surv.time=runif(N),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5),
#'                       E = rnorm(N))
#'
#' Cova.mtx = Phen.mtx[,c("Cov1","Cov2")]                                     # covariates dataframe excluding environmental factor E
#' E = Phen.mtx$E                                                             # environmental factor E
#'
#' # BGEN format for genotype data
#' GenoFile = system.file("", "GenoMat_SPAGxE.bgen", package = "SPAGxECCT")
#' GenoFileIndex = c(system.file("", "GenoMat_SPAGxE.bgen.bgi", package = "SPAGxECCT"),
#'                   system.file("", "GenoMat_SPAGxE.sample", package = "SPAGxECCT"))
#'
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#'
#' # fit a genotype-independent model
#'
#' R = SPA_G_Get_Resid("survival",
#'                     Surv(surv.time,event)~Cov1+Cov2+E,
#'                     data=Phen.mtx,
#'                     pIDs=Phen.mtx$ID,
#'                     gIDs=paste0("IID-",1:N))
#'
#' # calculate p values
#'
#' survival.res = SPAGxE_CCT(traits = "survival",                     # trait type
#'                           GenoFile = GenoFile,                     # a character of genotype file
#'                           GenoFileIndex = GenoFileIndex,           # additional index file(s) corresponding to GenoFile.
#'                           R = R,                                   # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
#'                           E = E,                                   # environmental factor
#'                           Phen.mtx = Phen.mtx,                     # phenotype dataframe
#'                           Cova.mtx = Cova.mtx)                     # a covariate matrix excluding the environmental factor E
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
#' head(survival.res)
#'
#' @export
#' @import survival
#' @import GRAB

SPAGxE_CCT = function(traits="survival/binary/quantitative/categorical",
                      GenoFile = NULL,
                      GenoFileIndex = NULL,
                      control = list(AllMarkers = TRUE),
                      Geno.mtx = NULL,
                      R,
                      E,
                      Phen.mtx,
                      Cova.mtx,
                      epsilon = 0.001,
                      Cutoff = 2,
                      impute.method = "fixed",
                      missing.cutoff = 0.15,
                      min.maf = 0.001,
                      G.model = "Add")
{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  G.model=G.model)


  # update on 2024-12
  # suppressPackageStartupMessages(library("GRAB",quietly = T))

  if(is.null(Geno.mtx)){
    Geno.mtx = GRAB::GRAB.ReadGeno(GenoFile = GenoFile,
                                   GenoFileIndex = GenoFileIndex,
                                   SampleIDs = Phen.mtx$ID,
                                   control = control)$GenoMat
  }



  check_input_Resid(pIDs = Phen.mtx$ID, gIDs = rownames(Geno.mtx), R = R)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 10) # update on 2023-12-27
  colnames(output) = c("MAF","missing.rate",
                       "p.value.spaGxE","p.value.spaGxE.Wald","p.value.spaGxE.CCT.Wald","p.value.normGxE",
                       "p.value.betaG", "Stat.betaG","Var.betaG","z.betaG") # update on 2023-12-27
  rownames(output) = colnames(Geno.mtx)

  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    g = Geno.mtx[,i]

    output.one.SNP = SPAGxE_CCT_one_SNP(traits=traits,
                                        g,                     # genotype vector
                                        R,                     # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                                        E,                     # environmental factor
                                        Phen.mtx,              # phenotype dataframe
                                        Cova.mtx,              # other covariates (such as age, sex, and top PCs) excluding E
                                        epsilon,               # a fixed value
                                        Cutoff,
                                        impute.method,
                                        missing.cutoff,
                                        min.maf,
                                        G.model)
    output[i,] = output.one.SNP
  }

  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}




#' A scalable and accurate framework for large-scale genome-wide gene-environment interaction (GxE) analysis (One-SNP-version).
#'
#' One-SNP-version SPAGxE_CCT function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function SPAGxE_CCT(). NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function SPAGxE_CCT().
#' @export

SPAGxE_CCT_one_SNP = function(traits="survival/binary/quantitative/categorical",
                              g,                     # genotype vector
                              R,                     # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                              E,                     # environmental factor
                              Phen.mtx,              # phenotype dataframe
                              Cova.mtx,              # other covariates (such as age, sex, and top PCs) excluding E
                              epsilon = 0.001,       # a fixed value
                              Cutoff = 2,
                              impute.method = "fixed",
                              missing.cutoff = 0.15,
                              min.maf = 0.001,
                              G.model = "Add")

{
  MAF = mean(g, na.rm=T)/2
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA, NA, NA, NA))

  VarG = 2 * MAF * (1 - MAF)


  S1 = sum(g*R)                    # test statistic for marginal genetic effect
  VarS1 = sum(R^2) * VarG          # estimated variance of S1
  Z1 = S1 / sqrt(VarS1)            # standardize S1
  pval.norm1 = pnorm(-abs(Z1))*2   # p value for marginal genetic effect from normal approximation

  # print(pval.norm1)

  if(pval.norm1 > epsilon){

    S2 = sum(g*E*R)                # test statistic for marginal GxE effect
    lambda = sum(E*R^2)/sum(R^2)   # lambda = Cov/VarS1
    S_GxE = S2 - lambda * S1       # new test statistic for marginal GxE effect
    R.new = (E - lambda) * R       # new residuals

    ################### SPA
    ### cutoff = 2
    res = SPA_G_one_SNP_homo_new(g=g, R=R.new)        # the third element is p-value from SPA-G, the fourth element is p-value from normal approximation
    pval.spaG = res[3]                                # p value from SPA
    pval.norm = res[4]                                # p-value from normal approximation

    # update on 2023-12-27
    pval.output = c(pval.spaG, pval.spaG, pval.spaG, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
  }else{

    W = cbind(1, g)

    # update on 2023-12-27
    R0 = R - (W)%*%(solve(t(W)%*%W))%*%(t(W)%*%R) # null model residuals adjusting for G

    S_GxE0 = sum(g*E*R0)  # test statistic for marginal GxE effect
    R.new0 = E*R0

    # update on 2023-10-24

    res = SPA_G_one_SNP_homo_new(g=g, R=R.new0, min.maf=min.maf)              # the third element is p-value from SPA-G, the fourth element is p-value from normal approximation
    pval.spaGxE = res[3]                                                      # p value from SPA
    pval.norm = res[4]                                                        # p-value from normal approximation

    # Wald test

    if(traits=="survival"){

      data0 = coxph(formula = Surv(Phen.mtx$surv.time, Phen.mtx$event) ~ as.matrix(Cova.mtx)+E+g+g*E, iter.max = 1000)

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }

    if(traits=="binary"){

      data0 = glm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E, family = "binomial") # update on 2023-12-27

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }

    if(traits=="quantitative"){

      data0 = lm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E) # update on 2023-12-27

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|t|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }

    if(traits=="categorical"){

      data0 = clm(formula = as.factor(Phen.mtx$Y) ~ as.matrix(Cova.mtx)+E+g+g*E)

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))


      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }
  }

  output.one.snp = c(MAF, missing.rate, pval.output, S1, VarS1, Z1)

  return(output.one.snp)
}



#' A scalable and accurate framework for large-scale genome-wide gene-environment interaction (GxE) analysis in admixed populations.
#'
#' A scalable and accurate analysis framework for a large-scale genome-wide SPAGxECCT implementation of gene-environmental interaction (GxE) analyses of quantitative traits, binary traits, time-to-event traits, and ordinal categorical traits.
#' @param GenoFile a character of genotype file. Currently, we support two formats of genotype input including PLINK and BGEN. Other formats such as VCF will be added later. Please refer to \code{?GRAB.ReadGeno} for more details.
#' @param GenoFileIndex additional index file(s) corresponding to GenoFile. Please refer to \code{?GRAB.ReadGeno} for more details.
#' @param control a list of parameters to decide which markers to extract. Please refer to \code{?GRAB.ReadGeno} for more details.
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missing genotypes should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param R model residuals after fitting a genotype-independent model (i.e., a covariate-only model in which marginal genetic effect and GxE effect are 0)
#' @param E a numeric environmental factor with each element as an environmental factor value of an individual.
#' @param Phen.mtx phenotype dataframe at least including three columns of ID, surv.time and event for time-to-event trait analysis, two columns of ID and linear phenotype Y for linear trait analysis, two columns of ID and binary phenotype Y for binary trait analysis, or two columns of ID and ordinal categorical phenotype Y for ordinal categorical trait analysis.
#' @param Cova.mtx a covariate matrix excluding the environmental factor E.
#' @param topPCs a covariate matrix including the SNP-derived principle components (PCs) containing all information of population structure.
#' @param epsilon a numeric value (default: 0.001) to specify the p-value cutoff for betaG estimation. Please see details for more information.
#' @param Cutoff a numeric value (default: 2) to specify the standard deviation cutoff to be used.
#'               If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation,
#'               otherwise, its p value is calculated based on a saddlepoint approximation.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e., 2MAF).
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#'                       Any variant with missing rate higher than this cutoff will be excluded from the analysis.
#' @param min.maf a numeric value (default: 0.0001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param G.model model type
#' @details To run SPAGxEmixCCT, the following two steps are required:
#' \itemize{
#'   \item Step 1: Use function SPA_G_Null_Model() or other functions to fit a genotype-independent (covariate-only) model to get residuals under a genotype-independent (covariate-only) model.
#'   \item Step 2: Use function SPAGxEmix_CCT() to calculate p value for each genetic variant to conduct a GxE analysis.
#' }
#'
#' SPAGxEmix_CCT() is an extension of SPAGxE_CCT() which uses individual-level allele frequencies to account for population stratification in an admixed popuklation.
#' SPAGxEmix_CCT() uses a hybrid strategy with both saddlepoint approximation and normal distribution approximation.
#' Generally speaking, saddlepoint approximation is more accurate than, but a little slower than, the traditional normal distribution approximation.
#' Hence, when the score statistic is close to 0 (i.e. p-values are not small), we use the normal distribution approximation.
#' And when the score statistic is far away from 0 (i.e. p-values are small), we use the saddlepoint approximation.
#' Argument 'Cutoff' is to specify the standard deviation cutoff.
#'
#' To calibrate the score statistics, SPAGxEmix_CCT() uses martingale residuals which are calculated via R package survival for time-to-event trait analysis, residuals from glm() for binary trait analysis, resuduals from lm() for quantitative trait analysis, and residuals from clm() for ordinal categorical trait analysis via R package ordinal.
#' All extentions (such as strata, ties, left-censoring) supported by package survival could also be used in SPAGxEmix_CCT().
#' Time-varying covariates are also supported by splitting each subject into several observations.
#'
#' Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
#' To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx.
#'
#' @return an R matrix with the following columns
#' \item{MAF}{Minor allele frequencies calculated with half of mean value of genotypes}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.spaGxE}{p value from SPAGxE method}
#' \item{p.value.spaGxE.Wald}{p value from SPAGxE_Wald method}
#' \item{p.value.spaGxE.CCT.Wald}{p value (recommanded) from SPAGxE_CCT method}
#' \item{p.value.normGxE}{p value from NormGxE method (based on a normal distribution approximation)}
#' \item{p.value.betaG}{p value of the marginal genetic effect based on a normal distribution approximation}
#' \item{Stat.betaG}{score statistics testing for marginal genetic effect}
#' \item{Var.betaG}{estimated variances of the score statistics testing for marginal genetic effect}
#' \item{z.betaG}{z values (using Var1) corresponding to the score statistics testing for marginal genetic effect}
#' \item{MAF.est.negative.num}{numbers of negative individual-level MAF estimated using linear regression model}
#' \item{MAC}{minor allele counts}
#' @examples
#'
#'# example 1  time-to-event phenotype (genotype input using R matrix format)
#'library(SPAGxECCT)
#'# Simulate phenotype and genotype
#'N = 10000
#'N.population1 = N/2
#'N.population2 = N/2
#'
#'nSNP = 100
#'MAF.population1 = 0.1
#'MAF.population2 = 0.3
#'
#'Phen.mtx.population1 = data.frame(ID = paste0("IID-",1:N.population1),
#'                                  event=rbinom(N.population1,1,0.5),
#'                                  surv.time=runif(N.population1),
#'                                  Cov1=rnorm(N.population1),
#'                                  Cov2=rbinom(N.population1,1,0.5),
#'                                  E = rnorm(N.population1),
#'                                  PC1 = 1)
#'
#'Phen.mtx.population2 = data.frame(ID = paste0("IID-",(N.population1+1):N),
#'                                  event=rbinom(N.population2,1,0.5),
#'                                  surv.time=runif(N.population2),
#'                                  Cov1=rnorm(N.population2),
#'                                  Cov2=rbinom(N.population2,1,0.5),
#'                                  E = rnorm(N.population2),
#'                                  PC1 = 0)
#'
#'Phen.mtx = rbind.data.frame(Phen.mtx.population1,
#'                            Phen.mtx.population2)     # phenotype dataframe
#'Cova.mtx = Phen.mtx[,c("Cov1","Cov2", "PC1")]         # covariate matrix excluding environmental factor
#'E = Phen.mtx$E                                        # environmental factor
#'
#'Geno.mtx.population1 = matrix(rbinom(N.population1*nSNP,2,MAF.population1),N.population1,nSNP)
#'Geno.mtx.population2 = matrix(rbinom(N.population2*nSNP,2,MAF.population2),N.population2,nSNP)
#'Geno.mtx = rbind(Geno.mtx.population1, Geno.mtx.population2)
#'
#'# NOTE: The row and column names of genotype matrix are required.
#'rownames(Geno.mtx) = paste0("IID-",1:N)
#'colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#'
#'# Attach the survival package so that we can use its function Surv()
#'library(survival)
#'
#'### fit a genotype-independent model for the SPAGxEmix_CCT analysis
#'
#'R = SPA_G_Get_Resid("survival",
#'                    Surv(surv.time,event)~Cov1+Cov2+PC1+E,
#'                    data=Phen.mtx,
#'                    pIDs=Phen.mtx$ID,
#'                    gIDs=paste0("IID-",1:N))
#'
#'### calculate p values
#'
#'survival.res = SPAGxEmix_CCT(traits = "survival",                     # trait type
#'                             Geno.mtx = Geno.mtx,                     # a character of genotype file
#'                             R = R,                                   # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
#'                             E = E,                                   # environmental factor
#'                             Phen.mtx = Phen.mtx,                     # phenotype dataframe
#'                             Cova.mtx = Cova.mtx,                     # a covariate matrix excluding the environmental factor E
#'                             topPCs = Cova.mtx[,"PC1"])               # PCs
#'
#'# we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
#'head(survival.res)
#'
#' # example 2  analysis of time-to-event phenotype (genotype input using PLINK file format)
#' # Simulate phenotype
#' N = 10000
#' N.population1 = N/2
#' N.population2 = N/2
#'
#' Phen.mtx.population1 = data.frame(ID = paste0("IID-",1:N.population1),
#'                                   event=rbinom(N.population1,1,0.5),
#'                                   surv.time=runif(N.population1),
#'                                   Cov1=rnorm(N.population1),
#'                                   Cov2=rbinom(N.population1,1,0.5),
#'                                   E = rnorm(N.population1),
#'                                   PC1 = 1)
#'
#' Phen.mtx.population2 = data.frame(ID = paste0("IID-",(N.population1+1):N),
#'                                   event=rbinom(N.population2,1,0.5),
#'                                   surv.time=runif(N.population2),
#'                                   Cov1=rnorm(N.population2),
#'                                   Cov2=rbinom(N.population2,1,0.5),
#'                                   E = rnorm(N.population2),
#'                                   PC1 = 0)
#'
#' Phen.mtx = rbind.data.frame(Phen.mtx.population1,
#'                             Phen.mtx.population2)  # phenotype dataframe
#' Cova.mtx = Phen.mtx[,c("Cov1","Cov2", "PC1")]      # covariate matrix excluding environmental factor
#' E = Phen.mtx$E                                     # environmental factor
#'
#' # PLINK format
#' GenoFile = system.file("", "GenoMat_SPAGxEmix.bed", package = "SPAGxECCT")
#'
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#'
#' R = SPA_G_Get_Resid("survival",
#'                     Surv(surv.time,event)~Cov1+Cov2+PC1+E,
#'                     data=Phen.mtx,
#'                     pIDs=Phen.mtx$ID,
#'                     gIDs=paste0("IID-",1:N))
#'
#' survival.res = SPAGxEmix_CCT(traits = "survival",                     # trait type
#'                              GenoFile = GenoFile,                     # a character of genotype file
#'                              R = R,                                   # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
#'                              E = E,                                   # environmental factor
#'                              Phen.mtx = Phen.mtx,                     # phenotype dataframe
#'                              Cova.mtx = Cova.mtx,                     # a covariate matrix excluding the environmental factor E
#'                              topPCs = Cova.mtx[,"PC1"])               # PCs
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
#' head(survival.res)
#'
#' # example 3  analysis of time-to-event phenotype (genotype input using BGEN file format)
#' library(SPAGxECCT)
#' # Simulate phenotype
#' N = 10000
#' N.population1 = N/2
#' N.population2 = N/2
#'
#' Phen.mtx.population1 = data.frame(ID = paste0("IID-",1:N.population1),
#'                                   event = rbinom(N.population1,1,0.5),
#'                                   surv.time = runif(N.population1),
#'                                   Cov1 = rnorm(N.population1),
#'                                   Cov2 = rbinom(N.population1,1,0.5),
#'                                   E = rnorm(N.population1),
#'                                   PC1 = 1)
#'
#' Phen.mtx.population2 = data.frame(ID = paste0("IID-",(N.population1+1):N),
#'                                   event = rbinom(N.population2,1,0.5),
#'                                   surv.time = runif(N.population2),
#'                                   Cov1 = rnorm(N.population2),
#'                                   Cov2 = rbinom(N.population2,1,0.5),
#'                                   E = rnorm(N.population2),
#'                                   PC1 = 0)
#'
#' Phen.mtx = rbind.data.frame(Phen.mtx.population1,
#'                             Phen.mtx.population2)  # phenotype dataframe
#' Cova.mtx = Phen.mtx[,c("Cov1","Cov2", "PC1")]      # covariate matrix excluding environmental factor
#' E = Phen.mtx$E                                     # environmental factor
#'
#' # BGEN format for genotype data
#' GenoFile = system.file("", "GenoMat_SPAGxEmix.bgen", package = "SPAGxECCT")
#' GenoFileIndex = c(system.file("", "GenoMat_SPAGxEmix.bgen.bgi", package = "SPAGxECCT"),
#'                   system.file("", "GenoMat_SPAGxEmix.sample", package = "SPAGxECCT"))
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#'
#' # fit a genotype-independent model for the SPAGxEmix_CCT analysis
#'
#' R = SPA_G_Get_Resid("survival",
#'                     Surv(surv.time,event)~Cov1+Cov2+PC1+E,
#'                     data=Phen.mtx,
#'                     pIDs=Phen.mtx$ID,
#'                     gIDs=paste0("IID-",1:N))
#'
#' # calculate p values
#'
#' survival.res = SPAGxEmix_CCT(traits = "survival",                     # trait type
#'                              GenoFile = GenoFile,                     # genotype file
#'                              GenoFileIndex = GenoFileIndex,           # additional index file(s) corresponding to GenoFile.
#'                              R = R,                                   # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
#'                              E = E,                                   # environmental factor
#'                              Phen.mtx = Phen.mtx,                     # phenotype dataframe
#'                              Cova.mtx = Cova.mtx,                     # a covariate matrix excluding the environmental factor E
#'                              topPCs = Cova.mtx[,"PC1"])               # PCs
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
#' head(survival.res)
#'
#' @export
#' @import survival

SPAGxEmix_CCT = function(traits="survival/binary/quantitative/categorical",
                         GenoFile = NULL,
                         GenoFileIndex = NULL,
                         control = list(AllMarkers = TRUE),
                         Geno.mtx = NULL,
                         R,                                                 # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                         E,                                                 # environmental factor
                         Phen.mtx,                                          # phenotype dataframe
                         Cova.mtx,                                          # other covariates (such as age, sex, and top PCs) excluding E
                         topPCs,                                            # a covariate matrix including the SNP-derived principle components (PCs) containing all information of population structure
                         epsilon = 0.001,                                   # a fixed value
                         Cutoff = 2,
                         impute.method = "fixed",
                         missing.cutoff = 0.15,
                         min.maf = 0.001,
                         G.model = "Add",
                         topPCs.pvalue.cutoff = 0.05,
                         MAF.est.negative.ratio.cutoff = 0.1)


{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  G.model=G.model)

  # update on 2024-12
  # suppressPackageStartupMessages(library("GRAB",quietly = T))

  if(is.null(Geno.mtx)){
    Geno.mtx = GRAB::GRAB.ReadGeno(GenoFile = GenoFile,
                                   GenoFileIndex = GenoFileIndex,
                                   SampleIDs = Phen.mtx$ID,
                                   control = control)$GenoMat
  }

  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 12)

  colnames(output) = c("MAF","missing.rate",
                       "p.value.spaGxE","p.value.spaGxE.Wald","p.value.spaGxE.CCT.Wald","p.value.normGxE",
                       "p.value.betaG", "Stat.betaG","Var.betaG","z.betaG",
                       "MAF.est.negative.num" ,"MAC") # update on 2024-04-16

  rownames(output) = colnames(Geno.mtx)

  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    g = Geno.mtx[,i]

    # print(i)

    output.one.SNP = SPAGxEmix_CCT_one_SNP(traits=traits,
                                           g,                     # genotype vector
                                           R,                     # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                                           E,                     # environmental factor
                                           Phen.mtx,              # phenotype dataframe
                                           Cova.mtx,              # other covariates (such as age, sex, and top PCs) excluding E
                                           topPCs,                # a covariate matrix including the SNP-derived principle components (PCs) containing all information of population structure
                                           epsilon,               # a fixed value
                                           Cutoff,
                                           impute.method,
                                           missing.cutoff,
                                           min.maf,
                                           G.model,
                                           topPCs.pvalue.cutoff,
                                           MAF.est.negative.ratio.cutoff)

    output[i,] = output.one.SNP
  }

  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}


#' A scalable and accurate framework for large-scale genome-wide gene-environment interaction (GxE) analysis in admixed populations (One-SNP-version).
#'
#' One-SNP-version SPAGxEmix_CCT() function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function SPAGxEmix_CCT(). NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function SPAGxEmix_CCT().
#' @export

SPAGxEmix_CCT_one_SNP = function(traits="survival/binary/quantitative/categorical",
                                 g,                     # genotype vector
                                 R,                     # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                                 E,                     # environmental factor
                                 Phen.mtx,              # phenotype dataframe
                                 Cova.mtx,              # other covariates (such as age, sex, and top PCs) excluding E
                                 topPCs,
                                 epsilon = 0.001,       # a fixed value
                                 Cutoff = 2,
                                 impute.method = "fixed",
                                 missing.cutoff = 0.15,
                                 min.maf = 0.001,
                                 G.model = "Add",
                                 topPCs.pvalue.cutoff = 0.05,
                                 MAF.est.negative.ratio.cutoff = 0.1)

{
  MAF = mean(g, na.rm=T)/2
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA, NA, NA, NA,
             NA, NA))



  N1set = which(g!=0)  # position of non-zero genotypes
  N0 = N-length(N1set)

  MAC = MAF*2*N      #  update on 2022-10-05


  if(MAC <= 20){
    MAF.est = c(rep(MAF, N))
    MAF.est.negative.num = 0
    topPCs.pvalueVec = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  }else{
    Fit.lm = lm(g~topPCs)
    MAF.est = Fit.lm$fitted.values/2

    topPCs.pvalueVec = summary(Fit.lm)$coefficients[,4][-1]

    MAF.est.negative.num = length(MAF.est[which(MAF.est < 0)])
    MAF.est.negative.ratio = MAF.est.negative.num/N

    if(MAF.est.negative.ratio > MAF.est.negative.ratio.cutoff){
      if(length(which(topPCs.pvalueVec < topPCs.pvalue.cutoff))==0){
        MAF.est = c(rep(MAF, N)) # MAF.all
      }else{
        selected.PCs = topPCs[,which(topPCs.pvalueVec < topPCs.pvalue.cutoff)]

        # update on 2023-02-13: MAF.est1
        g.tilde = round(g)

        # update on 2023-03-23: some variants sum(g.tilde)==0, thus MAF.est=0, and result in inflated type I error rates
        if(sum(g.tilde)==0 | sum(2-g.tilde)==0){
          MAF.est = c(rep(MAF, N)) # MAF.all
        }else{
          g.tilde[which(g.tilde != 0)] = 1
          Fit = glm(g.tilde~selected.PCs, family = binomial(link="logit"))
          MAF.est = 1- sqrt(1 - Fit$fitted.values)

        }
      }
    }else{
      MAF0 = 0          #  update on 2022-10-28
      MAF.est = ifelse(MAF.est <= 0, MAF0, MAF.est)
      MAF.est = ifelse(MAF.est >= 1, 1 - MAF0, MAF.est)
    }
  }



  g.var.est = 2 * MAF.est * (1 - MAF.est)

  S1 = sum(g * R)                             # test statistic for marginal genetic effect
  meanS1 = 2 * sum(R * MAF.est)
  VarS1 = sum(R^2 * g.var.est)                # estimated variance of S1
  Z1 = (S1 - meanS1) / sqrt(VarS1)            # standardize S1
  pval.norm1 = pnorm(-abs(Z1))*2              # p value for marginal genetic effect from normal approximation

  # print(pval.norm1)

  if(pval.norm1 > epsilon){

    S2 = sum(g*E*R)                                         # test statistic for marginal GxE effect
    lambda = sum(g.var.est * E* R^2)/sum(g.var.est * R^2)   # lambda = Cov/VarS1
    S_GxE = S2 - lambda * S1                                # new test statistic for marginal GxE effect
    R.new = (E - lambda) * R                                # new residuals

    S_GxE.mean = 2 * sum(R.new * MAF.est)
    S_GxE.var = sum(R.new^2 * g.var.est)

    z_GxE = (S_GxE - S_GxE.mean)/sqrt(S_GxE.var)

    if(abs(z_GxE) < Cutoff){
      pval.norm = pnorm(abs(z_GxE), lower.tail = FALSE)*2
      pval.output = c(pval.norm, pval.norm, pval.norm, pval.norm, pval.norm1) # 5 elements: SPAGxEmix, SPAGxEmixCCT_Wald, SPAGxEmixCCT_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }else{

      # update on 2023-10-24
      pval1 = GetProb_SPA_G_new(init.t = 0, MAF.est, R.new, max(S_GxE, (2*S_GxE.mean-S_GxE)), lower.tail = FALSE) # SPA-G p value
      pval2 = GetProb_SPA_G_new(init.t = 0, MAF.est, R.new, min(S_GxE, (2*S_GxE.mean-S_GxE)), lower.tail = TRUE)  # SPA-G p value

      # update on 2023-10-24
      if(is.na(pval2)){
        pval2 = pval1
      }

      pval3 = pnorm(abs(z_GxE), lower.tail = FALSE) # Normal
      pval4 = pnorm(-abs(z_GxE), lower.tail = TRUE) # Normal

      pval.spaG = pval1 + pval2
      pval.norm = pval3 + pval4


      # update on 2023-12-27
      pval.output = c(pval.spaG, pval.spaG, pval.spaG, pval.norm, pval.norm1)
    } # 5 elements: SPAGxEmix, SPAGxEmixCCT_Wald, SPAGxEmixCCT_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
  }else{

    W = cbind(1, g)

    # update on 2023-12-27
    R0 = R - (W)%*%(solve(t(W)%*%W))%*%(t(W)%*%R) # null model residuals adjusting for G

    S_GxE0 = sum(g*E*R0)  # test statistic for marginal GxE effect
    R.new0 = E*R0


    S_GxE0.mean = 2 * sum(R.new0 * MAF.est)
    S_GxE0.var = sum(R.new0^2 * g.var.est)

    z_GxE0 = (S_GxE0 - S_GxE0.mean)/sqrt(S_GxE0.var)

    if(abs(z_GxE0) < Cutoff){
      pval.norm = pnorm(abs(z_GxE0), lower.tail = FALSE)*2
      pval.spaGxE = pval.norm
    }else{

      # update on 2023-10-24
      pval1 = GetProb_SPA_G_new(init.t = 0, MAF.est, R.new0, max(S_GxE0, (2*S_GxE0.mean-S_GxE0)), lower.tail = FALSE) # SPA-G p value
      pval2 = GetProb_SPA_G_new(init.t = 0, MAF.est, R.new0, min(S_GxE0, (2*S_GxE0.mean-S_GxE0)), lower.tail = TRUE)  # SPA-G p value

      # update on 2023-10-24
      if(is.na(pval2)){
        pval2 = pval1
      }

      pval3 = pnorm(abs(z_GxE0), lower.tail = FALSE) # Normal
      pval4 = pnorm(-abs(z_GxE0), lower.tail = TRUE) # Normal

      pval.spaGxE = pval1 + pval2
      pval.norm = pval3 + pval4
    }

    # Wald test

    if(traits=="survival"){

      data0 = coxph(formula = Surv(Phen.mtx$surv.time, Phen.mtx$event) ~ as.matrix(Cova.mtx)+E+g+g*E, iter.max = 1000)

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # update on 2024-04-21
      if(is.na(pval.wald)){
        pval.wald = pval.spaGxE
      }

      # use CCT to conbine p values from SPAGxEmix and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix, SPAGxEmixCCT_Wald, SPAGxEmixCCT_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }

    if(traits=="binary"){

      data0 = glm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E, family = "binomial") # update on 2023-12-27

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # update on 2024-04-21
      if(is.na(pval.wald)){
        pval.wald = pval.spaGxE
      }

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix, SPAGxEmixCCT_Wald, SPAGxEmixCCT_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }

    if(traits=="quantitative"){

      data0 = lm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E) # update on 2023-12-27

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|t|)"]     # Wald test p value for marginal GxE effect

      # update on 2024-04-21
      if(is.na(pval.wald)){
        pval.wald = pval.spaGxE
      }

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix, SPAGxEmixCCT_Wald, SPAGxEmixCCT_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }

    if(traits=="categorical"){

      data0 = clm(formula = as.factor(Phen.mtx$Y) ~ as.matrix(Cova.mtx)+E+g+g*E)

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # update on 2024-04-21
      if(is.na(pval.wald)){
        pval.wald = pval.spaGxE
      }

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))


      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix, SPAGxEmixCCT_Wald, SPAGxEmixCCT_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }
  }

  output.one.snp = c(MAF, missing.rate, pval.output, S1, VarS1, Z1,
                     MAF.est.negative.num, MAC)

  return(output.one.snp)
}







#' A scalable and accurate framework to identify ancestry-specific GxE effects by incorporating local ancestry for large-scale genome-wide gene-environment interaction (GxE) analyses in admixed populations.
#'
#' A scalable and accurate analysis framework to efficiently identify ancestry-specific GxE effects by incorporating local ancestry for a large-scale genome-wide gene-environmental interaction (GxE) analyses of quantitative traits, binary traits, time-to-event traits, and ordinal categorical traits in admixed populations.
#' @param Geno.mtx a numeric ancestry-specific genotype matrix with each row as an individual and each column as a genetic variant from an ancestry.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missing genotypes should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param haplo.mtx matrix of local ancestry counts (the number of haplotypes) of an ancestry to analyze.
#'                  Each row represents an individual and each column represents local ancestry counts (the number of haplotypes) of the ancestry to analyze at a genetic variant.
#' @param R model residuals after fitting a genotype-independent model (i.e., a covariate-only model in which marginal genetic effect and GxE effect are 0)
#' @param E a numeric environmental factor with each element as an environmental factor value of an individual.
#' @param Phen.mtx phenotype dataframe at least including three columns of ID, surv.time and event for time-to-event trait analysis, two columns of ID and linear phenotype Y for linear trait analysis, two columns of ID and binary phenotype Y for binary trait analysis, or two columns of ID and ordinal categorical phenotype Y for ordinal categorical trait analysis.
#' @param Cova.mtx a covariate matrix excluding the environmental factor E and local ancestry counts (the number of haplotypes).
#' @param Cova.haplo.mtx.list a list with each element as a matrix of local ancestry counts (the number of haplotypes) of an ancestry (i.e., haplo.mtx of an ancestry).
#'                            If all samples are from K ancestries, then Cova.haplo.mtx.list with (K-1) elements from (K-1) ancestries is enough.
#'                            Names of all elements of Cova.haplo.mtx.list are needed.
#' @param epsilon a numeric value (default: 0.001) to specify the p-value cutoff for betaG estimation. Please see details for more information.
#' @param Cutoff a numeric value (Default: 2) to specify the standard deviation cutoff to be used.
#'               If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation,
#'               otherwise, its p value is calculated based on a saddlepoint approximation.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. haplo.num * MAF).
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#'                       Any variant with missing rate higher than this cutoff will be excluded from the analysis.
#' @param min.maf a numeric value (default: 0.0001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param G.model model type
#' @details To run SPAGxEmixCCT_local, the following two steps are required:
#' \itemize{
#'   \item Step 1: Use function SPA_G_Null_Model() or other functions to fit a genotype-independent (covariate-only) model to get residuals under a genotype-independent model.
#'   \item Step 2: Use function SPAGxEmixCCT_localance() to calculate p value for each genetic variant to conduct a ancestry-specific GxE analysis by incorporating local ancestry.
#' }
#'
#' SPAGxEmixCCT_localance() is an extension of SPAGxEmix_CCT() which identifies ancestry-specific GxE effects by incorporating local ancestry in an admixed population analysis.
#' SPAGxEmixCCT_localance() uses a hybrid strategy with both saddlepoint approximation and normal distribution approximation.
#' Generally speaking, saddlepoint approximation is more accurate than, but a little slower than, the traditional normal distribution approximation.
#' Hence, when the score statistic is close to 0 (i.e. p-values are not small), we use the normal distribution approximation.
#' And when the score statistic is far away from 0 (i.e. p-values are small), we use the saddlepoint approximation.
#' Argument 'Cutoff' is to specify the standard deviation cutoff.
#'
#' To calibrate the score statistics, SPAGxEmixCCT_localance() uses martingale residuals which are calculated via R package survival for time-to-event trait analysis, residuals from glm() for binary trait analysis, resuduals from lm() for quantitative trait analysis, and residuals from clm() for ordinal categorical trait analysis via R package ordinal.
#' All extentions (such as strata, ties, left-censoring) supported by package survival could also be used in SPAGxEmixCCT_localance().
#' Time-varying covariates are also supported by splitting each subject into several observations.
#'
#' Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
#' To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx.
#'
#' @return an R matrix with the following columns
#' \item{MAF}{ancestry-specific minor allele frequencies calculated with ancestry-specific genotypes}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.spaGxE.index.ance}{p value from SPAGxEmix_local method}
#' \item{p.value.spaGxE.Wald.index.ance}{p value from SPAGxEmix_Wald_local method}
#' \item{p.value.spaGxE.CCT.Wald.index.ance}{p value (recommanded) from SPAGxEmixCCT_local method}
#' \item{p.value.normGxE.index.ance}{p value from NormmixGxE_local method (based on a normal distribution approximation)}
#' \item{p.value.betaG.index.ance}{p value of the ancestry-specific marginal genetic effect based on a normal distribution approximation}
#' \item{Stat.betaG.index.ance}{ancestry-specific score statistics testing for ancestry-specific marginal genetic effect}
#' \item{Var.betaG.index.ance}{estimated variances of the ancestry-specific score statistics testing for ancestry-specific marginal genetic effect}
#' \item{z.betaG.index.ance}{z values (using Var1) corresponding to the ancestry-specific score statistics testing for ancestry-specific marginal genetic effect}
#' @examples
#'
#'# example 1  binary phenotype
#'library(SPAGxECCT)
#'# load in phenotype and genotype
#'
#'data("Pheno.mtx")           # phenotype data
#'data("Geno.mtx")            # genotype data
#'data("Geno.mtx.ance1")      # ancestry-specific genotype data for ancestry 1
#'data("Geno.mtx.ance2")      # ancestry-specific genotype data for ancestry 2
#'data("haplo.mtx.ance1")     # local ancestry counts for ancestry 1
#'data("haplo.mtx.ance2")     # local ancestry counts for ancestry 2
#'
#'Cova.mtx = Pheno.mtx[,c("PC1", "PC2", "PC3", "PC4", "Cov1", "Cov2")]    # Covariate matrix excluding environmental factor
#'E = Pheno.mtx$E                                                         # environmental factor
#'
#'### fit a genotype-independent model for the SPAGxEmix_CCT_local analysis
#'
#'resid  = SPA_G_Get_Resid(traits = "binary",
#'                         y ~ Cov1 + Cov2  + E + PC1 + PC2 + PC3 + PC4,family=binomial(link="logit"),
#'                         data=Pheno.mtx,
#'                         pIDs=Pheno.mtx$IID,
#'                         gIDs=rownames(Geno.mtx))
#'
#'### Cova.haplo.mtx.list
#'
#'Cova.haplo.mtx.list = list(haplo.mtx.ance1 = haplo.mtx.ance1,
#'                           haplo.mtx.ance2 = haplo.mtx.ance2)
#'
#'### calculate p values for ancestry 1
#'
#'binary_res_ance1 = SPAGxEmixCCT_localance(traits = "binary",
#'                                          Geno.mtx = Geno.mtx.ance1,
#'                                          R = resid,
#'                                          haplo.mtx = haplo.mtx.ance1,
#'                                          E = E,
#'                                          Phen.mtx = Pheno.mtx,
#'                                          Cova.mtx = Cova.mtx,
#'                                          Cova.haplo.mtx.list = Cova.haplo.mtx.list)
#'
#'
#'colnames(binary_res_ance1) = c("Marker", "MAF.ance1","missing.rate.ance1",
#'                               "Pvalue.spaGxE.ance1","Pvalue.spaGxE.Wald.ance1", "Pvalue.spaGxE.CCT.Wald.ance1",
#'                               "Pvalue.normGxE.ance1", "Pvalue.betaG.ance1",
#'                               "Stat.betaG.ance1","Var.betaG.ance1","z.betaG.ance1")
#'
#'### calculate p values for ancestry 2
#'
#'binary_res_ance2 = SPAGxEmixCCT_localance(traits = "binary",
#'                                          Geno.mtx = Geno.mtx.ance2,
#'                                          R = resid,
#'                                          haplo.mtx = haplo.mtx.ance2,
#'                                          E = E,
#'                                          Phen.mtx = Pheno.mtx,
#'                                          Cova.mtx = Cova.mtx,
#'                                          Cova.haplo.mtx.list = Cova.haplo.mtx.list)
#'
#'
#'colnames(binary_res_ance2) = c("Marker", "MAF.ance2","missing.rate.ance2",
#'                               "Pvalue.spaGxE.ance2","Pvalue.spaGxE.Wald.ance2", "Pvalue.spaGxE.CCT.Wald.ance2",
#'                               "Pvalue.normGxE.ance2", "Pvalue.betaG.ance2",
#'                               "Stat.betaG.ance2","Var.betaG.ance2","z.betaG.ance2")
#'
#'### merge data frame
#'binary.res = merge(binary_res_ance1, binary_res_ance2)
#'
#'# we recommand using column of 'p.value.spaGxE.CCT.Wald.index.ance' to associate genotype with phenotypes
#'head(binary.res)
#'
#' @export
#' @import survival


SPAGxEmixCCT_localance = function(traits="survival/binary/quantitative/categorical",
                                  Geno.mtx,
                                  haplo.mtx,
                                  R,                     # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                                  Cutoff = 2,
                                  E,                     # environmental factor
                                  Phen.mtx,              # phenotype dataframe
                                  Cova.mtx,              # other covariates (such as age, sex, and top PCs) excluding E
                                  Cova.haplo.mtx.list,
                                  epsilon = 0.001,       # a fixed value
                                  impute.method = "fixed",
                                  missing.cutoff = 0.15,
                                  min.maf = 0.0001,
                                  G.model = "Add")


{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  G.model=G.model)

  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 10)

  colnames(output) = c("MAF","missing.rate",
                       "p.value.spaGxE.index.ance","p.value.spaGxE.Wald.index.ance","p.value.spaGxE.CCT.Wald.index.ance","p.value.normGxE.index.ance",
                       "p.value.betaG.index.ance", "Stat.betaG.index.ance","Var.betaG.index.ance","z.betaG.index.ance")

  rownames(output) = colnames(Geno.mtx)

  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){

    g = Geno.mtx[,i]
    haplo_numVec = haplo.mtx[,i]

    Cova.haplo.mtx = c()

    for(name in names(Cova.haplo.mtx.list)) {

      Cova.haplo.mtx = cbind(Cova.haplo.mtx,
                             Cova.haplo.mtx.list[[name]][,i])
    }

    output.one.SNP = SPAGxEmixCCT_localance_one_SNP(traits=traits,
                                                    g,                     # ancestry-specific genotype vector
                                                    haplo_numVec,          # a vector of local ancestry counts (the number of haplotypes) of the ancestry to analyze at a genetic variant
                                                    R,                     # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                                                    E,                     # environmental factor
                                                    Phen.mtx,              # phenotype dataframe
                                                    Cova.mtx,              # other covariates (such as age, sex, and top PCs) excluding E
                                                    Cova.haplo.mtx,        # matrix of local ancestry haplotypes used in Wald test (for two-way admixed populations, only index haplo_numVec is needed instead of haplo.mtx)
                                                    epsilon,               # a fixed value
                                                    Cutoff,
                                                    impute.method,
                                                    missing.cutoff,
                                                    min.maf,
                                                    G.model)


    output[i,] = output.one.SNP
  }

  output = as.data.frame(cbind(rsID = colnames(Geno.mtx), output))

  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}




#' A scalable and accurate analysis framework for a large-scale genome-wide gene-environmental interaction (GxE) analyses of quantitative traits, binary traits, time-to-event traits, and ordinal categorical traits to identify ancestry-specific GxE effects by incorporating local ancestry (One-SNP-version).
#'
#' One-SNP-version SPAGxEmixCCT_localance() function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric ancestry-specific genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param haplo_numVec a vector of local ancestry counts (the number of haplotypes) of the ancestry to analyze at a genetic variant.
#' @param others the same as function SPAGxEmixCCT_localance(). NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function SPAGxEmixCCT_localance().
#' @export

SPAGxEmixCCT_localance_one_SNP = function(traits="survival/binary/quantitative/categorical",
                                          g,                     # ancestry-specific genotype vector
                                          haplo_numVec,          # a vector of local ancestry counts
                                          R,                     # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                                          E,                     # environmental factor
                                          Phen.mtx,              # phenotype dataframe
                                          Cova.mtx,              # other covariates (such as age, sex, and top PCs) excluding E
                                          Cova.haplo.mtx,        # matrix of local ancestry haplotypes used in Wald test (for two-way admixed populations, only index haplo_numVec is needed instead of haplo.mtx)
                                          epsilon = 0.001,       # a fixed value
                                          Cutoff = 2,
                                          impute.method = "fixed",
                                          missing.cutoff = 0.15,
                                          min.maf = 0.0001,
                                          G.model = "Add")

{
  MAF = sum(g)/sum(haplo_numVec)   # MAF in the ancestry to analyze
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA, NA, NA, NA))

  S1 = sum(g*R)                    # test statistic for ancestry-specific marginal genetic effect

  MAC = sum(g)
  MAF.est.Vec = c(rep(MAF, N))

  S1.mean = sum(MAF.est.Vec * haplo_numVec * R)

  g.var.est.Vec = haplo_numVec * MAF.est.Vec * (1 - MAF.est.Vec)

  VarS1 = sum(R^2 * g.var.est.Vec)         # estimated variance of S1

  Z1 = (S1 - S1.mean) / sqrt(VarS1)        # standardize S1
  pval.norm1 = pnorm(-abs(Z1))*2           # p value for ancestry-specific marginal genetic effect from normal approximation

  # print(pval.norm1)

  if(pval.norm1 > epsilon){

    S2 = sum(g*E*R)                                                  # ancestry-specific test statistic for ancestry-specific marginal GxE effect
    lambda = sum(haplo_numVec * E * R^2) / sum(haplo_numVec * R^2)   # lambda = Cov/VarS1
    S_GxE = S2 - lambda * S1                                         # new test statistic for ancestry-specific marginal GxE effect
    R.new = (E - lambda) * R                                         # new residuals

    ################### SPA

    res = SPAmix_localance_one_SNP(g=g, R=R.new, haplo_numVec=haplo_numVec, min.maf=min.maf)  # the third element is p-value from SPA-G, the fourth element is p-value from normal approximation


    pval.spaG = res[3]                                                                        # p value from SPA
    pval.norm = res[4]                                                                        # p value from normal approximation

    # update on 2023-12-27
    pval.output = c(pval.spaG, pval.spaG, pval.spaG, pval.norm, pval.norm1) # 5 elements: SPAGxEmix_local, SPAGxEmixCCT_local_Wald, SPAGxEmixCCT_local_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
  }else{

    W = cbind(1, g)

    # update on 2023-12-27
    R0 = R - (W)%*%(solve(t(W)%*%W))%*%(t(W)%*%R) # null model residuals adjusting for G

    S_GxE0 = sum(g*E*R0)  # test statistic for marginal GxE effect
    R.new0 = E*R0

    # update on 2023-10-24

    res = SPAmix_localance_one_SNP(g=g, R=R.new0, haplo_numVec=haplo_numVec, min.maf=min.maf)       # the third element is p-value from SPA-G, the fourth element is p-value from normal approximation

    pval.spaGxE = res[3]                                                                            # p value from SPA
    pval.norm = res[4]                                                                              # p-value from normal approximation

    # Wald test

    if(traits=="survival"){

      data0 = coxph(formula = Surv(Phen.mtx$surv.time, Phen.mtx$event) ~ as.matrix(Cova.haplo.mtx)+as.matrix(Cova.mtx)+E+g+g*E, iter.max = 1000) # update

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # use CCT to combine p values from SPAGxEmixCCT_local and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix_local, SPAGxEmixCCT_local_Wald, SPAGxEmixCCT_local_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }

    if(traits=="binary"){

      data0 = glm(formula = Phen.mtx$Y ~ as.matrix(Cova.haplo.mtx)+as.matrix(Cova.mtx)+E+g+g*E, family = "binomial") # update

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxEmixCCT_local and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix_local, SPAGxEmixCCT_local_Wald, SPAGxEmixCCT_local_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }

    if(traits=="quantitative"){

      data0 = lm(formula = Phen.mtx$Y ~ as.matrix(Cova.haplo.mtx)+as.matrix(Cova.mtx)+haplo_numVec+E+g+g*E) # update

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|t|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxEmixCCT_local and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix_local, SPAGxEmixCCT_local_Wald, SPAGxEmixCCT_local_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }

    if(traits=="categorical"){

      data0 = clm(formula = as.factor(Phen.mtx$Y) ~ as.matrix(Cova.haplo.mtx)+as.matrix(Cova.mtx)+haplo_numVec+E+g+g*E) # update

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxEmixCCT_local and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))


      # update on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix_local, SPAGxEmixCCT_local_Wald, SPAGxEmixCCT_local_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }
  }

  output.one.snp = c(MAF, missing.rate, pval.output, S1, VarS1, Z1)

  return(output.one.snp)
}


#### other functions

M_G0 = function(t, MAF){
  re = (1 - MAF + MAF * exp(t))^2
  return(re)
}

M_G1 = function(t, MAF){
  re = 2*(MAF * exp(t))*(1 - MAF + MAF * exp(t))
  return(re)
}

M_G2 = function(t, MAF){
  re = 2*(MAF * exp(t))^2 + 2*(MAF * exp(t))*(1 - MAF + MAF * exp(t))
  return(re)
}

K_G0 = function(t_R,    # t x R
                MAF){
  re = log(M_G0(t_R, MAF))
  return(re)
}

K_G1 = function(t, MAF){
  re = M_G1(t, MAF)/M_G0(t, MAF)
  return(re)
}

K_G2 = function(t, MAF){
  re = (M_G0(t, MAF)*M_G2(t, MAF)-M_G1(t, MAF)^2)/M_G0(t, MAF)^2
  return(re)
}

H_org = function(t, R, MAFVec){

  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(K_G0(t1*R, MAFVec))
  }
  return(out)
}

H1_adj = function(t, R, s, MAFVec)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(R*K_G1(t1*R, MAFVec)) - s
  }
  return(out)
}

H2 = function(t, R, MAFVec)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(R^2*K_G2(t1*R, MAFVec))
  }
  return(out)
}



fastgetroot_H1 = function(init.t,
                          R,          # residuals
                          s,          # observed test statistic
                          MAFVec,
                          tol = .Machine$double.eps^0.25,
                          maxiter = 100)
{
  t = init.t;
  H1_adj_eval = 0
  diff.t = Inf
  converge = T
  for(iter in 1:maxiter){
    old.t = t
    old.diff.t = diff.t
    old.H1 = H1_adj_eval
    H1_adj_eval = H1_adj(t, R, s, MAFVec)
    H2_eval = H2(t, R, MAFVec)

    diff.t = -1 * H1_adj_eval / H2_eval
    if(is.na(H1_adj_eval) | is.infinite(H1_adj_eval)){
      # checked it on 07/05:
      # if the solution 't' tends to infinity, 'K2_eval' tends to 0, and 'K1_eval' tends to 0 very slowly.
      # then we can set the one sided p value as 0 (instead of setting converge = F)
      t = sign(s)*Inf
      H2_eval = 0;
      break;
    }
    if(sign(H1_adj_eval) != sign(old.H1)){
      while(abs(diff.t) > abs(old.diff.t) - tol){
        diff.t = diff.t/2
      }
    }
    if(is.na(diff.t)){
      diff.t =  5
    }
    if(abs(diff.t) < tol) break;
    t = old.t + diff.t
  }
  if(iter == maxiter) converge = F
  return(list(root = t,
              iter = iter,
              converge = converge,
              H2_eval = H2_eval))
}



GetProb_SPA_G = function(init.t = 0, MAFVec, R, s, lower.tail){

  # out = uniroot(H1_adj, c(-10,10), extendInt = "upX",
  #               R=R, s=s, MAFVec=MAFVec)
  out = fastgetroot_H1(init.t, R = R, s = s, MAFVec = MAFVec)

  zeta = out$root

  k1 = H_org(zeta, R=R, MAFVec=MAFVec)
  k2 = H2(zeta, R=R, MAFVec=MAFVec)

  temp1 = zeta * s - k1

  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}

  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)

  if(is.na(pval)){
    pval =  0
  }

  re = pval
  return(re)
}


# updataed on 2023-10-24
fastgetroot_H1_new = function(init.t,
                              R,          # residuals
                              s,          # observed test statistic
                              MAFVec,
                              tol = .Machine$double.eps^0.25,
                              maxiter = 100)
{
  t = init.t;
  H1_adj_eval = 0
  diff.t = Inf
  converge = T
  for(iter in 1:maxiter){
    old.t = t
    old.diff.t = diff.t
    old.H1 = H1_adj_eval
    H1_adj_eval = H1_adj(t, R, s, MAFVec)
    H2_eval = H2(t, R, MAFVec)

    diff.t = -1 * H1_adj_eval / H2_eval

    # update on 2023-10-24
    if(is.na(diff.t)){
      diff.t =  5
    }

    if(is.na(H1_adj_eval) | is.infinite(H1_adj_eval)){
      # checked it on 07/05:
      # if the solution 't' tends to infinity, 'K2_eval' tends to 0, and 'K1_eval' tends to 0 very slowly.
      # then we can set the one sided p value as 0 (instead of setting converge = F)
      t = sign(s)*Inf
      H2_eval = 0;
      break;
    }
    if(sign(H1_adj_eval) != sign(old.H1)){
      while(abs(diff.t) > abs(old.diff.t) - tol){
        diff.t = diff.t/2
      }
    }
    if(is.na(diff.t)){
      diff.t =  5
    }
    if(abs(diff.t) < tol) break;
    t = old.t + diff.t
  }

  # update on 2023-10-24

  if(iter == maxiter) {
    converge = F
    t = NA}

  return(list(root = t,
              iter = iter,
              converge = converge,
              H2_eval = H2_eval))
}


# updataed on 2023-10-24

GetProb_SPA_G_new = function(init.t = 0, MAFVec, R, s, lower.tail){

  # out = uniroot(H1_adj, c(-10,10), extendInt = "upX",
  #               R=R, s=s, MAFVec=MAFVec)
  out = fastgetroot_H1_new(init.t, R = R, s = s, MAFVec = MAFVec)

  zeta = out$root

  k1 = H_org(zeta, R=R, MAFVec=MAFVec)
  k2 = H2(zeta, R=R, MAFVec=MAFVec)

  temp1 = zeta * s - k1

  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}

  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)

  if(is.na(pval)){
    pval =  0
  }

  re = pval
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

check_input_Resid = function(pIDs, gIDs, R, range=c(-100,100))
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


############# SPA-G function for homogeneous population
SPA_G_one_SNP_homo = function(g,
                              R,
                              Cutoff = 2,
                              impute.method = "fixed",
                              missing.cutoff = 0.15,
                              min.maf = 0.00001,       # update on 2023-12-27
                              G.model = "Add")
{
  ## calculate MAF and update genotype vector
  MAF = mean(g, na.rm=T)/2
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA))

  #if(!is.null(obj.null$p2g))
  # g = g[obj.null$p2g]

  ## Score statistic
  S = sum(g * R)

  ## estimated variance without adjusting for covariates
  N1set = 1:N
  N0 = 0
  MAF.est = MAF

  g.var.est = 2 * MAF.est * (1 - MAF.est)
  S.var = sum(R^2 * g.var.est)
  z = S/sqrt(S.var)

  if(abs(z) < Cutoff){
    pval.norm = pnorm(abs(z), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.var, z))
  }

  pval1 = GetProb_SPA_G(init.t = 0, MAF.est, R, abs(S), lower.tail = FALSE)
  pval2 = GetProb_SPA_G(init.t = 0, MAF.est, R, -abs(S), lower.tail = TRUE)

  pval3 = pnorm(abs(z), lower.tail = FALSE) # Normal
  pval4 = pnorm(-abs(z), lower.tail = TRUE) # Normal

  pval.spa.G = pval1 + pval2
  pval.norm = pval3 + pval4

  # if(abs(z) < Cutoff){
  #   pval.spa.G = pval.norm
  # }

  pval = c(pval.spa.G, pval.norm) # 2 elements: element 1 is from EmpSPA-G, element 2 is from Normal

  return(c(MAF, missing.rate, pval, S, S.var, z))
}





############# SPA-G function for homogeneous population -- 2023-10-24
SPA_G_one_SNP_homo_new = function(g,
                                  R,
                                  Cutoff = 2,
                                  impute.method = "fixed",
                                  missing.cutoff = 0.15,
                                  min.maf = 0.00001,       # update on 2023-12-27
                                  G.model = "Add")
{
  ## calculate MAF and update genotype vector
  MAF = mean(g, na.rm=T)/2
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA))

  #if(!is.null(obj.null$p2g))
  # g = g[obj.null$p2g]

  ## Score statistic
  S = sum(g * R)

  ## estimated variance without adjusting for covariates
  N1set = 1:N
  N0 = 0
  MAF.est = MAF

  g.var.est = 2 * MAF.est * (1 - MAF.est)
  S.var = sum(R^2 * g.var.est)

  # update on 2023-10-24
  S.mean = 2 * sum(R * MAF.est)

  z = (S - S.mean)/sqrt(S.var)


  if(abs(z) < Cutoff){
    pval.norm = pnorm(abs(z), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.var, z))
  }

  # update on 2023-10-24
  pval1 = GetProb_SPA_G_new(init.t = 0, MAF.est, R, max(S, (2*S.mean-S)), lower.tail = FALSE) # SPA-G p value
  pval2 = GetProb_SPA_G_new(init.t = 0, MAF.est, R, min(S, (2*S.mean-S)), lower.tail = TRUE)  # SPA-G p value

  # update on 2023-10-24
  if(is.na(pval2)){
    pval2 = pval1
  }

  pval3 = pnorm(abs(z), lower.tail = FALSE) # Normal
  pval4 = pnorm(-abs(z), lower.tail = TRUE) # Normal

  pval.spa.G = pval1 + pval2
  pval.norm = pval3 + pval4

  pval = c(pval.spa.G, pval.norm) # 2 elements: element 1 is from SPA-G, element 2 is from Normal

  return(c(MAF, missing.rate, pval, S, S.var, z))
}




CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  which.one = which(pvals==1)
  is.one <- (length(which.one)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    pvals[which.one] = 0.999
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}



#### other functions with input of local ancestry information

# MAF is the minor allele frequency of the tested marker in the ancestry to analyze
# MAFVec is the vector (length = sample size) of minor allele frequency of the tested marker in the ancestry to analyze
# haplo_numVec is the vector (length = sample size) of the number haplotype of the tested marker in the ancestry to analyze

# The MGF of G (genotype from the ancestry to analyze)
M_G0_local = function(t, MAF, haplo_num){
  re = (1 - MAF + MAF * exp(t))^haplo_num
  return(re)
}

# The first derivative of the MGF of G (genotype from the ancestry to analyze)
M_G1_local = function(t, MAF, haplo_num){
  re = haplo_num*(MAF * exp(t))*(1 - MAF + MAF * exp(t))^(haplo_num - 1)
  return(re)
}

# The second derivative of the MGF of G (genotype from the ancestry to analyze)
M_G2_local = function(t, MAF, haplo_num){
  re = haplo_num * (haplo_num-1) * (MAF * exp(t))^2 * (1 - MAF + MAF * exp(t))^(haplo_num - 2) +
    haplo_num * (MAF * exp(t)) * (1 - MAF + MAF * exp(t)) ^ (haplo_num - 1)
  return(re)
}


# The CGF of G (genotype from the ancestry to analyze)
K_G0_local = function(t, MAF, haplo_num){
  re = log(M_G0_local(t, MAF, haplo_num))
  return(re)
}

# The first derivative of the CGF of G (genotype from the ancestry to analyze)
K_G1_local = function(t, MAF, haplo_num){
  re = M_G1_local(t, MAF, haplo_num)/M_G0_local(t, MAF, haplo_num)
  return(re)
}

# The second derivative of the CGF of G (genotype from the ancestry to analyze)
K_G2_local = function(t, MAF, haplo_num){
  re = (M_G0_local(t, MAF, haplo_num)*M_G2_local(t, MAF, haplo_num)-M_G1_local(t, MAF, haplo_num)^2)/M_G0_local(t, MAF, haplo_num)^2
  return(re)
}

# The CGF of score test statistic for the ancestry to analyze
H_org_local = function(t, R, MAFVec, haplo_numVec){

  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(K_G0_local(t1*R, MAFVec, haplo_numVec))
  }
  return(out)
}

# The first derivative of the CGF of score test statistic for the ancestry to analyze
H1_adj_local = function(t, R, s, MAFVec, haplo_numVec)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(R*K_G1_local(t1*R, MAFVec, haplo_numVec)) - s
  }
  return(out)
}


# The second derivative of the CGF of score test statistic for the ancestry to analyze
H2_local = function(t, R, MAFVec, haplo_numVec)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(R^2*K_G2_local(t1*R, MAFVec, haplo_numVec))
  }
  return(out)
}


GetProb_SPA_G_local = function(MAFVec, R, haplo_numVec, s, lower.tail){

  out = uniroot(H1_adj_local, c(-20,20), extendInt = "upX",
                R=R, s=s, MAFVec=MAFVec, haplo_numVec=haplo_numVec)
  zeta = out$root

  k1 = H_org_local(t=zeta, R=R, MAFVec=MAFVec, haplo_numVec=haplo_numVec)
  k2 = H2_local(t=zeta, R=R, MAFVec=MAFVec, haplo_numVec=haplo_numVec)

  temp1 = zeta * s - k1

  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}

  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  re = pval
  return(re)
}

# update on 2023-10-24
fastgetroot_H1_new_local = function(init.t,
                                    R,          # residuals
                                    s,          # observed test statistic
                                    MAFVec,
                                    haplo_numVec,
                                    tol = .Machine$double.eps^0.25,
                                    maxiter = 100)
{
  t = init.t;
  H1_adj_eval = 0
  diff.t = Inf
  converge = T
  for(iter in 1:maxiter){
    old.t = t
    old.diff.t = diff.t
    old.H1 = H1_adj_eval
    H1_adj_eval = H1_adj_local(t, R, s, MAFVec, haplo_numVec)
    H2_eval = H2_local(t, R, MAFVec, haplo_numVec)

    diff.t = -1 * H1_adj_eval / H2_eval

    # update on 2023-10-24
    if(is.na(diff.t)){
      diff.t =  5
    }

    if(is.na(H1_adj_eval) | is.infinite(H1_adj_eval)){
      # checked it on 07/05:
      # if the solution 't' tends to infinity, 'K2_eval' tends to 0, and 'K1_eval' tends to 0 very slowly.
      # then we can set the one sided p value as 0 (instead of setting converge = F)
      t = sign(s)*Inf
      H2_eval = 0;
      break;
    }
    if(sign(H1_adj_eval) != sign(old.H1)){
      while(abs(diff.t) > abs(old.diff.t) - tol){
        diff.t = diff.t/2
      }
    }
    if(is.na(diff.t)){
      diff.t =  5
    }
    if(abs(diff.t) < tol) break;
    t = old.t + diff.t
  }

  # update on 2023-10-24

  if(iter == maxiter) {
    converge = F
    t = NA}

  return(list(root = t,
              iter = iter,
              converge = converge,
              H2_eval = H2_eval))
}


# update on 2023-10-24

GetProb_SPA_G_new_local = function(init.t = 0, MAFVec, R, haplo_numVec, s, lower.tail){

  # out = uniroot(H1_adj, c(-10,10), extendInt = "upX",
  #               R=R, s=s, MAFVec=MAFVec)
  out = fastgetroot_H1_new_local(init.t, R = R, s = s, haplo_numVec = haplo_numVec, MAFVec = MAFVec)

  zeta = out$root

  k1 = H_org_local(zeta, R=R, MAFVec=MAFVec, haplo_numVec = haplo_numVec)
  k2 = H2_local(zeta, R=R, MAFVec=MAFVec, haplo_numVec = haplo_numVec)

  temp1 = zeta * s - k1

  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}

  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)

  if(is.na(pval)){
    pval =  0
  }

  re = pval
  return(re)
}





#### SPAmix #########################################################################

SPAmix_localance_one_SNP = function(g,
                                    R,
                                    haplo_numVec,
                                    Cutoff = 2,
                                    impute.method = "fixed",
                                    missing.cutoff = 0.15,
                                    min.maf = 0.00001,          # update on 2022-08-16 : replace 0.0001 by 0.000001
                                    G.model = "Add")
{
  ## calculate MAF and update genotype vector
  # MAF = mean(g, na.rm=T)/2       # MAF
  MAF = sum(g)/sum(haplo_numVec)   # MAF in the ancestry to analyze

  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }
  ####################################################### AF, not MAF
  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA, NA, NA))

  ## Score statistic
  S = sum(g * R)

  ## estimated variance without adjusting for covariates
  # N1set = 1:N
  # N0 = 0

  N1set = which(g!=0)  # position of non-zero genotypes
  N0 = N-length(N1set)

  MAC = sum(g)
  MAF.est.Vec = c(rep(MAF, N))

  S.mean = sum(MAF.est.Vec * haplo_numVec * R)

  g.var.est.Vec = haplo_numVec * MAF.est.Vec * (1 - MAF.est.Vec)
  S.var = sum(R^2 * g.var.est.Vec)

  z = (S - S.mean)/sqrt(S.var)

  if(abs(z) < Cutoff){
    pval.norm = pnorm(abs(z), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.mean, S.var, z, MAC))
  }

  pval1 = GetProb_SPA_G_new_local(MAFVec=MAF.est.Vec, R=R, haplo_numVec=haplo_numVec, s=max(S, (2*S.mean-S)), lower.tail = FALSE) # SPA
  pval2 = GetProb_SPA_G_new_local(MAFVec=MAF.est.Vec, R=R, haplo_numVec=haplo_numVec, s=min(S, (2*S.mean-S)), lower.tail = TRUE)  # SPA

  pval3 = pnorm(abs(z), lower.tail = FALSE) # Normal distribution approximation
  pval4 = pnorm(-abs(z), lower.tail = TRUE) # Normal distribution approximation

  pval.spa.G = pval1 + pval2
  pval.norm = pval3 + pval4

  pval = c(pval.spa.G, pval.norm) # 2 elements: element 1 is from SPA, element 2 is from Normal distribution approximation

  return(c(MAF, missing.rate, pval, S, S.mean, S.var, z, MAC))
}


SPAmix_localance = function(Geno.mtx,
                            R,
                            haplo.mtx,
                            Cutoff = 2,
                            impute.method = "fixed",
                            missing.cutoff = 0.15,
                            min.maf = 0.00001,
                            G.model = "Add")
{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  G.model=G.model)

  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 9)
  colnames(output) = c("MAF","missing.rate","Pvalue.SPAmix.index.ance","Pvalue.norm.index.ance",
                       "Stat","Mean","Var","z", "MAC") # update on 2022-10-05 : MAF.est.negative.num
  rownames(output) = colnames(Geno.mtx)

  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    g = Geno.mtx[,i]
    haplo_numVec = haplo.mtx[,i]
    output.one.SNP = SPAmix_localance_one_SNP(g,
                                              R,
                                              haplo_numVec,
                                              Cutoff,
                                              impute.method,
                                              missing.cutoff,
                                              min.maf,
                                              G.model)

    output[i,] = output.one.SNP
  }

  output = as.data.frame(cbind(rsID = colnames(Geno.mtx), output))

  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}













