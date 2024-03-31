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
#' @import ordinal

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




#' SaddlePoint Approximation implementation of gene-environmental interaction (GxE) analysis of complex traits.
#'
#' A scalable and accurate analysis framework for a large-scale genome-wide SPAGxECCT implementation of gene-environmental interaction (GxE) analyses of quantitative traits, binary traits, time-to-event traits, and ordinal categorical traits.
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missing genotypes should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param R model residuals after fitting a genotype-independent model (i.e., a covariate-only model in which marginal genetic effect and GxE effect are 0)
#' @param E a numeric environmental factor with each element as a evnironmental factor value of an individual.
#' @param Phen.mtx phenotype dataframe at least including three columns of ID, surv.time and event for time-to-event trait analysis, two columns of ID and linear phenotype Y for linear trait analysis, two columns of ID and binary phenotype Y for binary trait analysis, or two columns of ID and ordinal categorical phenotype Y for ordinal categorical trait analysis.
#' @param Cova.mtx a covariate matrix excluding the environmental factor E.
#' @param epsilon a numeric value (default: 0.001) to specify the p-value cutoff for betaG estimation. Please see details for more information.
#' @param Cutoff a numeric value (Default: 2) to specify the standard deviation cutoff to be used.
#'               If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation,
#'               otherwise, its p value is calculated based on a saddlepoint approximation.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#'                       Any variant with missing rate higher than this cutoff will be excluded from the analysis.
#' @param min.maf a numeric value (default: 0.00001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param G.model model type
#' @details To run SPAGxECCT, the following two steps are required:
#' \itemize{
#'   \item Step 1: Use function SPA_G_Null_Model() to fit a genotype-independent (covariate-only) model.
#'   \item Step 2: Use function SPAGxE_CCT() to calculate p value for each genetic variant to conduct a GxE analysis.
#' }
#'
#' SPAGxE_CCT() uses a hybrid strategy with both saddlepoint approximation and normal distribution approximation.
#' Generally speaking, saddlepoint approximation is more accurate than, but a little slower than, the traditional normal distribution approximation.
#' Hence, when the score statistic is close to 0 (i.e. p-values are not small), we use the normal distribution approximation.
#' And when the score statistic is far away from 0 (i.e. p-values are small), we use the saddlepoint approximation.
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
#'
#' survival.res = SPAGxE_CCT("survival",
#'                           Geno.mtx,                     # genotype vector
#'                           R,                     # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
#'                           E,                     # environmental factor
#'                           Phen.mtx,              # include surv.time, event
#'                           Cova.mtx)
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with time-to-event phenotypes
#' head(survival.res)
#'
#'
#' # example 2  binary phenotype
#' # Simulation phenotype and genotype
#' N = 10000
#' nSNP = 100
#' MAF = 0.1
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       Y=rbinom(N,1,0.5),
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
#' R = SPA_G_Get_Resid("binary",
#'                     glm(formula = Y ~ Cov1+Cov2+E, data = Phen.mtx, family = "binomial"),
#'                     data=Phen.mtx,
#'                     pIDs=Phen.mtx$ID,
#'                     gIDs=paste0("IID-",1:N))
#'
#' binary.res = SPAGxE_CCT("binary",
#'                         Geno.mtx,                     # genotype vector
#'                         R,                     # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
#'                         E,                     # environmental factor
#'                         Phen.mtx,              # include surv.time, event
#'                         Cova.mtx)
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with binary phenotypes
#' head(binary.res)
#'
#'
#' # example 3  quantitative phenotype
#' # Simulation phenotype and genotype
#' N = 10000
#' nSNP = 100
#' MAF = 0.1
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       Y=rnorm(N,1,0.5),
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
#' R = SPA_G_Get_Resid("quantitative",
#'                     lm(formula = Y ~ Cov1+Cov2+E, data = Phen.mtx),
#'                     data=Phen.mtx,
#'                     pIDs=Phen.mtx$ID,
#'                     gIDs=paste0("IID-",1:N))
#'
#' quantitative.res = SPAGxE_CCT("quantitative",
#'                                Geno.mtx,                     # genotype vector
#'                                R,                     # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
#'                                E,                     # environmental factor
#'                                Phen.mtx,              # include surv.time, event
#'                                Cova.mtx)
#'
#' # we recommand using column of 'p.value.spaGxE.CCT.Wald' to associate genotype with binary phenotypes
#' head(quantitative.res)
#'
#' @export
#' @import survival
#' @import ordinal

SPAGxE_CCT = function(traits="survival/binary/quantitative/categorical",
                      Geno.mtx,              # genotype vector
                      R,                     # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                      E,                     # environmental factor
                      Phen.mtx,              # include surv.time, event
                      Cova.mtx,              #  other covariates (such as age, gender, and top PCs) excluding E
                      epsilon = 0.001,       # a fixed value
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
  check_input_Resid(pIDs = Phen.mtx$ID, gIDs = rownames(Geno.mtx), R = R)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 10) # updated on 2023-12-27
  colnames(output) = c("MAF","missing.rate",
                       "p.value.spaGxE","p.value.spaGxE.Wald","p.value.spaGxE.CCT.Wald","p.value.normGxE",
                       "p.value.betaG", "Stat.betaG","Var.betaG","z.betaG") # updated on 2023-12-27
  rownames(output) = colnames(Geno.mtx)

  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    g = Geno.mtx[,i]
    print(i)
    output.one.SNP = SPAGxE_CCT_one_SNP(traits=traits,
                                        g,                     # genotype vector
                                        R,                     # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                                        E,                     # environmental factor
                                        Phen.mtx,              # include surv.time, event, X1, X2, and E
                                        Cova.mtx,              #  other covariates (such as age, gender, and top PCs) excluding E
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




#' SaddlePoint Approximation implementation of gene-environmental interaction (GxE) analysis of complex traits (One-SNP-version).
#'
#' One-SNP-version SPAGxE_CCT function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function SPAGxE_CCT_surv. NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function SPAGxE_CCT.
#' @export

SPAGxE_CCT_one_SNP = function(traits="survival/binary/quantitative/categorical",
                              g,                     # genotype vector
                              R,                     # null model residuals (null model in which marginal genetic effect and GxE effect are 0)
                              E,                     # environmental factor
                              Phen.mtx,              # include surv.time, event
                              Cova.mtx,              #  other covariates (such as age, gender, and top PCs) excluding E
                              epsilon = 0.001,       # a fixed value
                              Cutoff = 2,
                              impute.method = "fixed",
                              missing.cutoff = 0.15,
                              min.maf = 0.00001,
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

  print(pval.norm1)

  if(pval.norm1 > epsilon){

    S2 = sum(g*E*R)                # test statistic for marginal GxE effect
    lambda = sum(E*R^2)/sum(R^2)   # lambda = Cov/VarS1
    S_GxE = S2 - lambda * S1       # new test statistic for marginal GxE effect
    R.new = (E - lambda) * R       # new residuals

    ################### SPA
    ### cutoff = 2
    res = SPA_G_one_SNP_homo_new(g=g, R=R.new)        # the third element is p-value from SPA-G, the fourth element is p-value from normal approximation
    pval.spaG = res[3]                                # p value from SPA-G
    pval.norm = res[4]                                # p-value from normal approximation

    # updated on 2023-12-27
    pval.output = c(pval.spaG, pval.spaG, pval.spaG, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
  }else{

    W = cbind(1, g)

    # updated on 2023-12-27
    R0 = R - (W)%*%(solve(t(W)%*%W))%*%(t(W)%*%R) # null model residuals adjusting for G

    S_GxE0 = sum(g*E*R0)  # test statistic for marginal GxE effect
    R.new0 = E*R0

    # updated on 2023-10-24

    res = SPA_G_one_SNP_homo_new(g=g, R=R.new0, min.maf=min.maf)              # the third element is p-value from SPA-G, the fourth element is p-value from normal approximation
    pval.spaGxE = res[3]                                # p value from SPA-G
    pval.norm = res[4]                                # p-value from normal approximation

    # Wald test

    if(traits=="survival"){

      data0 = coxph(formula = Surv(Phen.mtx$surv.time, Phen.mtx$event) ~ as.matrix(Cova.mtx)+E+g+g*E, iter.max = 1000)

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }

    if(traits=="binary"){

      data0 = glm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E, family = "binomial") # updated on 2023-12-27

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }

    if(traits=="quantitative"){

      data0 = lm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E) # updated on 2023-12-27

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|t|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))

      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
    }

    if(traits=="categorical"){

      data0 = clm(formula = as.factor(Phen.mtx$Y) ~ as.matrix(Cova.mtx)+E+g+g*E)

      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect

      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))


      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxE, SPAGxE_Wald, SPAGxE_CCT, SPAGxE_normal, and normal p-value for marginal GxE effect, and normal p value for marginal genetic effect
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

    # updated on 2023-10-24
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

  # updated on 2023-10-24

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
                              min.maf = 0.00001,       # updated on 2023-12-27
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
                                  min.maf = 0.00001,       # updated on 2023-12-27
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

  # updated on 2023-10-24
  S.mean = 2 * sum(R * MAF.est)

  z = (S - S.mean)/sqrt(S.var)


  if(abs(z) < Cutoff){
    pval.norm = pnorm(abs(z), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.var, z))
  }

  # updated on 2023-10-24
  pval1 = GetProb_SPA_G_new(init.t = 0, MAF.est, R, max(S, (2*S.mean-S)), lower.tail = FALSE) # SPA-G p value
  pval2 = GetProb_SPA_G_new(init.t = 0, MAF.est, R, min(S, (2*S.mean-S)), lower.tail = TRUE) # SPA-G p value

  # updated on 2023-10-24
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


















