#' A scalable and accurate analytical framework to account for familial relatedness in large-scale genome-wide gene-environment interaction (GxE) analysis.
#'
#' A scalable and accurate analytical framework to account for familial relatedness in large-scale genome-wide gene-environmental interaction (GxE) analyses of quantitative, binary, time-to-event, ordinal categorical, and longitudinal traits.
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missing genotypes should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param E a numeric environmental factor with each element as an environmental factor value of an individual.
#' @param Phen.mtx phenotype dataframe at least including three columns of ID, surv.time and event for time-to-event trait analysis, two columns of ID and linear phenotype Y for linear trait analysis, two columns of ID and binary phenotype Y for binary trait analysis, or two columns of ID and ordinal categorical phenotype Y for ordinal categorical trait analysis.
#' @param Cova.mtx a covariate matrix excluding the environmental factor E.
#' @param sparseGRM a three-column sparse GRM file with the first column as "ID1",the second column as "ID2", and the last column as "Value" (i.e., two times of kinship coefficient) without information of distant genetic relatedness (such as population structure).
#' @param obj.SPAGxE_Plus_Nullmodel a null model object from SPAGxE_Plus_Nullmodel()
#' @param epsilon a numeric value (default: 0.001) to specify the p-value cutoff for betaG estimation. Please see details for more information.
#' @param Cutoff a numeric value (Default: 2) to specify the standard deviation cutoff to be used.
#'               If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation,
#'               otherwise, its p value is calculated based on a saddlepoint approximation.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value.
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#'                       Any variant with missing rate higher than this cutoff will be excluded from the analysis.
#' @param min.maf a numeric value (default: 0.001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param G.model model type
#' @details To run SPAGxE_Plus, the following two steps are required:
#' \itemize{
#'   \item Step 1: Use function SPAGxE_Plus_Nullmodel() to fit a genotype-independent (covariate-only) model to get residuals under a genotype-independent (covariate-only) model.
#'   \item Step 2: Use function SPAGxE_Plus() to calculate p value for each genetic variant to conduct a GxE analysis.
#' }
#'
#' SPAGxE_Plus() is an extension of SPAGxE_CCT() which additionally uses sparse GRM to account for familial relatedness.
#' SPAGxE_Plus() uses a hybrid strategy with both saddlepoint approximation and normal distribution approximation.
#' Generally speaking, saddlepoint approximation is more accurate than, but a little slower than, the traditional normal distribution approximation.
#' Hence, when the score statistic is close to 0 (i.e., p-values are not small), we use the normal distribution approximation.
#' And when the score statistic is far away from 0 (i.e., p-values are small), we use the saddlepoint approximation.
#' Argument 'Cutoff' is to specify the standard deviation cutoff.
#'
#' Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
#' To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx.
#'
#' @return an R matrix with the following columns
#' \item{MAF}{Minor allele frequencies calculated with half of mean value of genotypes}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.spaGxE.plus}{p value from SPAGxE method}
#' \item{p.value.spaGxE.plus.Wald}{p value from SPAGxE_Wald method}
#' \item{p.value.spaGxE.plus.CCT.Wald}{p value (recommanded) from SPAGxE_CCT method}
#' \item{p.value.normGxE.plus}{p value from NormGxE method (based on a normal distribution approximation)}
#' \item{p.value.betaG}{p value of the marginal genetic effect based on a normal distribution approximation}
#' \item{Stat.betaG}{score statistics testing for marginal genetic effect}
#' \item{Var.betaG}{estimated variances of the score statistics testing for marginal genetic effect}
#' \item{z.betaG}{z values (using Var1) corresponding to the score statistics testing for marginal genetic effect}
#'
#' @export
#' @import survival
#' @import lme4
#' @import dplyr

SPAGxE_Plus = function(Geno.mtx,                                          # genotype vector
                       E,                                                 # environmental factor
                       Phen.mtx,                                          # phenotype dataframe
                       Cova.mtx,                                          # other covariates (such as age, gender, and top PCs) excluding E
                       sparseGRM,
                       obj.SPAGxE_Plus_Nullmodel,
                       epsilon = 0.001,                                   # a fixed value
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

  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 8) # updated on 2023-12-27
  colnames(output) = c("MAF","missing.rate",
                       "p.value.spaGxE.plus", "p.value.normGxE.plus",
                       "p.value.betaG", "Stat.betaG","Var.betaG","z.betaG")
  rownames(output) = colnames(Geno.mtx)


  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    g = Geno.mtx[,i]

    print(i)

    output.one.SNP = SPAGxE_Plus_one_SNP(g,                     # genotype vector
                                         E,                     # environmental factor
                                         Phen.mtx,              # phenotype dataframe
                                         Cova.mtx,              # other covariates (such as age, gender, and top PCs) excluding E
                                         sparseGRM,
                                         obj.SPAGxE_Plus_Nullmodel,
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

#' A scalable and accurate framework to account for familial relatedness for large-scale genome-wide gene-environment interaction (GxE) analysis (One-SNP-version).
#'
#' One-SNP-version SPAGxE_Plus() function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function SPAGxE_Plus(). NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function SPAGxE_Plus().
#' @export

SPAGxE_Plus_one_SNP = function(g,                     # genotype vector
                               E,                     # environmental factor
                               Phen.mtx,              # phenotype dataframe
                               Cova.mtx,              # other covariates (such as age, gender, and top PCs) excluding E
                               sparseGRM,
                               obj.SPAGxE_Plus_Nullmodel,
                               epsilon = 0.001,       # a fixed value
                               Cutoff = 2,
                               impute.method = "fixed",
                               missing.cutoff = 0.15,
                               min.maf = 0.001,
                               G.model = "Add")

{
  ResidMat = obj.SPAGxE_Plus_Nullmodel$ResidMat
  R = obj.SPAGxE_Plus_Nullmodel$R
  R_GRM_R = obj.SPAGxE_Plus_Nullmodel$R_GRM_R
  R_GRM_RE = obj.SPAGxE_Plus_Nullmodel$R_GRM_RE
  R.new = obj.SPAGxE_Plus_Nullmodel$R.new
  R.new_GRM_R.new = obj.SPAGxE_Plus_Nullmodel$R.new_GRM_R.new

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
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA, NA))

  N1set = which(g!=0)  # position of non-zero genotypes
  N0 = N-length(N1set)
  g.var.est = 2 * MAF * (1 - MAF)

  ############################################################################################
  # R = ResidMat$Resid
  S1 = sum(g * R)
  meanS1 = 0
  VarS1 = g.var.est * R_GRM_R
  Z1 = (S1 - meanS1) / sqrt(VarS1)            # standardize S1
  pval.norm1 = pnorm(-abs(Z1))*2              # p value for marginal genetic effect from normal approximation

  if(pval.norm1 > epsilon){

    S2 = sum(g*E*R)                                         # test statistic for marginal GxE effect

    ####

    lambda = R_GRM_RE/R_GRM_R   # lambda = Cov/VarS1

    S_GxE = S2 - lambda * S1                                # new test statistic for marginal GxE effect
    S_GxE.var = g.var.est * R.new_GRM_R.new

    z_GxE = S_GxE/sqrt(S_GxE.var)

    if(abs(z_GxE) < Cutoff){
      pval.norm = pnorm(abs(z_GxE), lower.tail = FALSE)*2
      pval.output = c(pval.norm, pval.norm, pval.norm1) # 3 elements: SPAGxE+, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }else{

      S_GxE.var.SPAGxE = sum(R.new^2 * g.var.est)
      Var.ratio.S_GxE = S_GxE.var.SPAGxE / S_GxE.var

      S_GxE.new = S_GxE * sqrt(Var.ratio.S_GxE)

      pval1 = GetProb_SPA_G_new(init.t = 0, MAF, R.new, abs(S_GxE.new), lower.tail = FALSE)  # SPA-G p value
      pval2 = GetProb_SPA_G_new(init.t = 0, MAF, R.new, -abs(S_GxE.new), lower.tail = TRUE)  # SPA-G p value

      if(is.na(pval2)){
        pval2 = pval1
      }

      pval3 = pnorm(abs(z_GxE), lower.tail = FALSE) # Normal
      pval4 = pnorm(-abs(z_GxE), lower.tail = TRUE) # Normal

      pval.spaG = pval1 + pval2
      pval.norm = pval3 + pval4

      pval.output = c(pval.spaG, pval.norm, pval.norm1)
    } # 3 elements: SPAGxE+, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
  }else{

    W = cbind(1, g)

    R0 = R - (W)%*%(solve(t(W)%*%W))%*%(t(W)%*%R) # null model residuals adjusting for G

    S_GxE0 = sum(g*E*R0)  # test statistic for marginal GxE effect
    R.new0 = E*R0

    ####
    ResidMat = ResidMat %>% mutate(R.new0 = R.new0)

    GRM1 = sparseGRM

    GRM1$R.new0_value1 = ResidMat$R.new0[match(sparseGRM$ID1, ResidMat$SubjID)]
    GRM1$R.new0_value2 = ResidMat$R.new0[match(sparseGRM$ID2, ResidMat$SubjID)]
    GRM1 = GRM1 %>% mutate(Cov_R.new0 = Value * R.new0_value1 * R.new0_value2)

    S_GxE0.mean = 2 * sum(R.new0) * MAF
    S_GxE0.var = 2 * sum(GRM1$Cov_R.new0) * g.var.est - sum(R.new0^2) * g.var.est

    # S_GxE0.var = (GRM1 %>% select(Cov_R.new0_sqrtMAF) %>% sum) * 2 - sum(R.new0^2 * g.var.est)
    # t(R.new0) %*% GRM %*% R.new0 * g.var.est

    z_GxE0 = (S_GxE0 - S_GxE0.mean)/sqrt(S_GxE0.var)

    if(abs(z_GxE0) < Cutoff){
      pval.norm = pnorm(abs(z_GxE0), lower.tail = FALSE)*2
      pval.spaGxE = pval.norm
    }else{

      S_GxE0.var.SPAGxE = sum(R.new0^2) * g.var.est
      Var.ratio.S_GxE0 = S_GxE0.var.SPAGxE / S_GxE0.var

      S_GxE0.new = S_GxE0 * sqrt(Var.ratio.S_GxE0)


      # pval1 = GetProb_SPA_G_new(init.t = 0, MAF.est, R.new0, max(S_GxE0, (2*S_GxE0.mean-S_GxE0)), Var.ratio, lower.tail = FALSE) # SPA-G p value
      # pval2 = GetProb_SPA_G_new_adjust(init.t = 0, MAF.est, R.new0, min(S_GxE0, (2*S_GxE0.mean-S_GxE0)), Var.ratio, lower.tail = TRUE)  # SPA-G p value

      pval1 = GetProb_SPA_G_new(init.t = 0, MAF, R.new0, abs(S_GxE0.new), lower.tail = FALSE)  # SPA-G p value
      pval2 = GetProb_SPA_G_new(init.t = 0, MAF, R.new0, -abs(S_GxE0.new), lower.tail = TRUE)  # SPA-G p value

      # updated on 2023-10-24
      if(is.na(pval2)){
        pval2 = pval1
      }

      pval3 = pnorm(abs(z_GxE0), lower.tail = FALSE) # Normal
      pval4 = pnorm(-abs(z_GxE0), lower.tail = TRUE) # Normal

      pval.spaGxE = pval1 + pval2
      pval.norm = pval3 + pval4
    }

    # updated on 2023-12-27
    pval.output = c(pval.spaGxE, pval.norm, pval.norm1) # SPAGxE+

  }

  output.one.snp = c(MAF, missing.rate, pval.output, S1, VarS1, Z1)

  return(output.one.snp)
}





#' Fits a genotype-independent model of SPAGxE+
#'
#' Fits a null linear regression model for quantitative traits, a null logistic regression model for binary traits ,or a null Cox proportional hazards model for time-to-event traits and then calculates residuals under a genotype-independent model.
#' @param traits a character value corresponding to phenotype. It should be "binary" for binary phenotype, "survival" for time-to-event phenotype, and "others" for other types of phenotype.
#' @param formula a formula to be passed to function lm(), glm(), or coxph(). For more details, please refer to package survival.
#' @param data a data.frame in which to interpret the variables named in the formula
#' @param pIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order in the formula.
#' @param gIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order of the Geno.mtx (i.e. the input of the function SPAGxE_CCT()).
#' @param range a two-element numeric vector (default: c(-100,100)) to specify the domain of the empirical CGF.
#' @param length.out a positive integer (default: 9999) for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.
#' @param sparseGRM a three-column sparse GRM file with the first column as "ID1",the second column as "ID2", and the last column as "Value" (i.e., two times of kinship coefficient) without information of distant genetic relatedness (such as population structure).
#' @param E a numeric environmental factor with each element as an environmental factor value of an individual.
#' @param ... Other arguments passed to function glm() or coxph(). For more details, please refer to package survival.
#' @return Residuals after fitting a genotype-independent (covariate-only) model.

SPAGxE_Plus_Nullmodel = function(traits="survival/binary/quantitative",
                                 formula=NULL,
                                 data=NULL,
                                 pIDs=NULL,
                                 gIDs=NULL,
                                 range=c(-100,100),
                                 length.out = 10000,
                                 sparseGRM,
                                 E,
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


  R = resid
  RE = R*E

  ResidMat = data.frame(SubjID = Phen.mtx$IID, Resid = R, RE = RE)

  #### calculate

  # R_GRM_R,
  # R_GRM_RE,
  # R.new,
  # R.new_GRM_R.new

  #### # update on 2024-09-11 ###############################################################

  sparseGRM$ID1 = as.character(sparseGRM$ID1); sparseGRM$ID2 = as.character(sparseGRM$ID2)

  SubjID.In.Resid = ResidMat$SubjID
  SubjID.In.GRM = unique(c(sparseGRM$ID1, sparseGRM$ID2))

  if(any(!SubjID.In.Resid %in% SubjID.In.GRM))
    stop("At least one subject in residual matrix does not have GRM information.")

  SubjID = SubjID.In.Resid
  sparseGRM = sparseGRM %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)

  GRM1 = sparseGRM
  GRM1$R_value1 = ResidMat$Resid[match(sparseGRM$ID1, ResidMat$SubjID)]
  GRM1$R_value2 = ResidMat$Resid[match(sparseGRM$ID2, ResidMat$SubjID)]
  GRM1 = GRM1 %>% mutate(Cov_R = Value * R_value1 * R_value2)             # update on 2024-09-12

  ####
  R_GRM_R = (GRM1 %>% select(Cov_R) %>% sum) * 2 - sum(R^2)

  ####
  GRM1$RE_value1 = ResidMat$RE[match(sparseGRM$ID1, ResidMat$SubjID)]
  GRM1$RE_value2 = ResidMat$RE[match(sparseGRM$ID2, ResidMat$SubjID)]
  GRM1 = GRM1 %>% mutate(Cov_R_RE_part1 = Value * R_value1 * RE_value2)         # update on 2024-09-12
  GRM1 = GRM1 %>% mutate(Cov_R_RE_part2 = Value * R_value2 * RE_value1)         # update on 2024-09-12
  GRM1 = GRM1 %>% mutate(Cov_RE = Value * RE_value1 * RE_value2)        # update on 2024-09-12

  ####
  RE_GRM_RE = (GRM1 %>% select(Cov_RE) %>% sum) * 2 - sum((RE)^2)

  ####
  R_GRM_RE = (GRM1 %>% select(Cov_R_RE_part1) %>% sum) +
    (GRM1 %>% select(Cov_R_RE_part2) %>% sum) -
    sum(RE * R)
  # Cov_S1_S2 = t(R_sqrtMAF) %*% GRM %*% RE_sqrtMAF

  lambda =  R_GRM_RE/R_GRM_R   # lambda = Cov/VarS1

  ####
  R.new = (E - lambda) * R
  ResidMat = ResidMat %>% mutate(R.new = R.new)

  GRM1$R.new_value1 = ResidMat$R.new[match(sparseGRM$ID1, ResidMat$SubjID)]
  GRM1$R.new_value2 = ResidMat$R.new[match(sparseGRM$ID2, ResidMat$SubjID)]
  GRM1 = GRM1 %>% mutate(Cov_R.new = Value * R.new_value1 * R.new_value2)             # update on 2024-09-12

  ####
  R.new_GRM_R.new = (GRM1 %>% select(Cov_R.new) %>% sum) * 2 - sum(R.new^2)

  obj.SPAGxE_Plus_Nullmodel = list(ResidMat = ResidMat,
                                   R = R,
                                   R_GRM_R = R_GRM_R,
                                   R_GRM_RE = R_GRM_RE,
                                   R.new = R.new,
                                   R.new_GRM_R.new = R.new_GRM_R.new)

  return(obj.SPAGxE_Plus_Nullmodel)
}


