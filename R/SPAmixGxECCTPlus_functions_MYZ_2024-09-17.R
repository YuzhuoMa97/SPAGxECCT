#' A scalable and accurate analytical framework to account for both population structure and family relatedness in large-scale genome-wide gene-environment interaction (GxE) analysis in admixed populations.
#'
#' A scalable and accurate analytical framework to account for both population structure and family relatedness in large-scale genome-wide gene-environmental interaction (GxE) analyses of quantitative traits and binary traits.
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missing genotypes should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param ResidMat a two-column data frame with the first column as "SubjID" and the second column as model residuals after fitting a genotype-independent model (i.e., a covariate-only model in which marginal genetic effect and GxE effect are 0)
#' @param E a numeric environmental factor with each element as an environmental factor value of an individual.
#' @param Phen.mtx phenotype dataframe at least including three columns of ID, surv.time and event for time-to-event trait analysis, two columns of ID and linear phenotype Y for linear trait analysis, two columns of ID and binary phenotype Y for binary trait analysis, or two columns of ID and ordinal categorical phenotype Y for ordinal categorical trait analysis.
#' @param Cova.mtx a covariate matrix excluding the environmental factor E.
#' @param topPCs a covariate matrix including the SNP-derived principle components (PCs) containing all information of population structure (without recent genetic relatedness or kinship information).
#' @param sparseGRM a three-column sparse GRM file with the first column as "ID1",the second column as "ID2", and the last column as "Value" (i.e., two times of kinship coefficient) without information of distant genetic relatedness (such as population structure).
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
#' @details To run SPAGxEmix_CCT_Plus, the following two steps are required:
#' \itemize{
#'   \item Step 1: Use function SPA_G_Null_Model() or other functions to fit a genotype-independent (covariate-only) model to get residuals under a genotype-independent (covariate-only) model.
#'   \item Step 2: Use function SPAGxEmix_CCT_Plus() to calculate p value for each genetic variant to conduct a GxE analysis.
#' }
#'
#' SPAGxEmix_CCT_Plus() is an extension of SPAGxEmix_CCT() which additionally uses sparse GRM to account for family relatedness in an admixed population or multiple populations.
#' SPAGxEmix_CCT_Plus() uses a hybrid strategy with both saddlepoint approximation and normal distribution approximation.
#' Generally speaking, saddlepoint approximation is more accurate than, but a little slower than, the traditional normal distribution approximation.
#' Hence, when the score statistic is close to 0 (i.e. p-values are not small), we use the normal distribution approximation.
#' And when the score statistic is far away from 0 (i.e. p-values are small), we use the saddlepoint approximation.
#' Argument 'Cutoff' is to specify the standard deviation cutoff.
#'
#' Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
#' To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx.
#'
#' @return an R matrix with the following columns
#' \item{MAF}{Minor allele frequencies calculated with half of mean value of genotypes}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.spaGxEmix.plus}{p value from SPAGxE method}
#' \item{p.value.spaGxEmix.plus.Wald}{p value from SPAGxE_Wald method}
#' \item{p.value.spaGxEmix.plus.CCT.Wald}{p value (recommanded) from SPAGxE_CCT method}
#' \item{p.value.normGxEmix.plus}{p value from NormGxE method (based on a normal distribution approximation)}
#' \item{p.value.betaG}{p value of the marginal genetic effect based on a normal distribution approximation}
#' \item{Stat.betaG}{score statistics testing for marginal genetic effect}
#' \item{Var.betaG}{estimated variances of the score statistics testing for marginal genetic effect}
#' \item{z.betaG}{z values (using Var1) corresponding to the score statistics testing for marginal genetic effect}
#' \item{MAF.est.negative.num}{numbers of negative individual-level MAF estimated using linear regression model}
#' \item{MAC}{minor allele counts}
#'
#' @export
#' @import survival
#' @import ordinal
#' @import lme4

SPAGxEmix_CCT_Plus = function(traits="binary/quantitative",
                              Geno.mtx,                                          # genotype vector
                              ResidMat,                                          # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                              E,                                                 # environmental factor
                              Phen.mtx,                                          # phenotype dataframe
                              Cova.mtx,                                          # other covariates (such as age, gender, and top PCs) excluding E
                              topPCs,                                            # a covariate matrix including the SNP-derived principle components (PCs) containing all information of population structure
                              sparseGRM,
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
  
  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))
  
  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 12)
  
  colnames(output) = c("MAF","missing.rate",
                       "p.value.spaGxEmix.plus","p.value.spaGxEmix.plus.Wald","p.value.spaGxEmix.plus.CCT.Wald","p.value.normGxEmix.plus",
                       "p.value.betaG", "Stat.betaG","Var.betaG","z.betaG",
                       "MAF.est.negative.num" ,"MAC") # updated on 2024-04-16
  
  rownames(output) = colnames(Geno.mtx)
  
  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())
  
  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    g = Geno.mtx[,i]
    
    print(i)
    
    output.one.SNP = SPAGxEmix_CCT_Plus_one_SNP(traits=traits,
                                                g,                     # genotype vector
                                                ResidMat,                     # residuals from genotype-independent model (null model in which marginal genetic effect and GxE effect are 0)
                                                E,                     # environmental factor
                                                Phen.mtx,              # phenotype dataframe
                                                Cova.mtx,              # other covariates (such as age, gender, and top PCs) excluding E
                                                topPCs,                # a covariate matrix including the SNP-derived principle components (PCs) containing all information of population structure
                                                sparseGRM,
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


#'A scalable and accurate framework to account for both population structure and family relatedness for large-scale genome-wide gene-environment interaction (GxE) analysis in admixed populations (One-SNP-version).
#'
#' One-SNP-version SPAGxEmix_CCT_Plus() function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function SPAGxEmix_CCT_Plus(). NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function SPAGxEmix_CCT_Plus().
#' @export

SPAGxEmix_CCT_Plus_one_SNP = function(traits="survival/binary/quantitative/categorical",
                                      g,                     # genotype vector
                                      ResidMat,
                                      E,                     # environmental factor
                                      Phen.mtx,              # phenotype dataframe
                                      Cova.mtx,              # other covariates (such as age, gender, and top PCs) excluding E
                                      topPCs,
                                      sparseGRM,
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
  
  MAC = MAF*2*N      #  updated on 2022-10-05
  
  
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
        
        # updated on 2023-02-13: MAF.est1
        g.tilde = round(g)
        
        # updated on 2023-03-23: some variants sum(g.tilde)==0, thus MAF.est=0, and result in inflated type I error rates
        if(sum(g.tilde)==0 | sum(2-g.tilde)==0){
          MAF.est = c(rep(MAF, N)) # MAF.all
        }else{
          g.tilde[which(g.tilde != 0)] = 1
          Fit = glm(g.tilde~selected.PCs, family = binomial(link="logit"))
          MAF.est = 1- sqrt(1 - Fit$fitted.values)
          
        }
      }
    }else{
      MAF0 = 0          #  updated on 2022-10-28
      MAF.est = ifelse(MAF.est <= 0, MAF0, MAF.est)
      MAF.est = ifelse(MAF.est >= 1, 1 - MAF0, MAF.est)
    }
  }
  
  
  
  g.var.est = 2 * MAF.est * (1 - MAF.est)
  
  R_sqrtMAF = R * sqrt(g.var.est)  # update on 2024-09-10
  
  ResidMat_new = ResidMat %>% mutate(Resid_new = R_sqrtMAF) # update on 2024-09-16
  
  #### # update on 2024-09-11 ###############################################################
  
  sparseGRM$ID1 = as.character(sparseGRM$ID1); sparseGRM$ID2 = as.character(sparseGRM$ID2)
  
  SubjID.In.Resid = ResidMat$SubjID
  SubjID.In.GRM = unique(c(sparseGRM$ID1, sparseGRM$ID2))
  
  if(any(!SubjID.In.Resid %in% SubjID.In.GRM))
    stop("At least one subject in residual matrix does not have GRM information.")
  
  SubjID = SubjID.In.Resid
  sparseGRM = sparseGRM %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
  
  GRM1 = sparseGRM
  GRM1$R_sqrtMAF_pos1 = ResidMat_new$Resid_new[match(sparseGRM$ID1, ResidMat_new$SubjID)]
  GRM1$R_sqrtMAF_pos2 = ResidMat_new$Resid_new[match(sparseGRM$ID2, ResidMat_new$SubjID)]
  GRM1 = GRM1 %>% mutate(Cov_R_sqrtMAF = Value * R_sqrtMAF_pos1 * R_sqrtMAF_pos2)         # update on 2024-09-12
  
  ############################################################################################
  S1 = sum(g * R)  
  meanS1 =  2 * sum(R * MAF.est)
  VarS1 = (GRM1 %>% select(Cov_R_sqrtMAF) %>% sum) * 2 - sum(R^2 * g.var.est) # t(R_sqrtMAF) %*% GRM %*% R_sqrtMAF
  Z1 = (S1 - meanS1) / sqrt(VarS1)            # standardize S1
  pval.norm1 = pnorm(-abs(Z1))*2              # p value for marginal genetic effect from normal approximation

  if(pval.norm1 > epsilon){
    
    S2 = sum(g*E*R)                                         # test statistic for marginal GxE effect

    ####
    RE_sqrtMAF =  R * E * sqrt(g.var.est)  # update on 2024-09-16
    
    ResidMat_new = ResidMat_new %>% mutate(RE_sqrtMAF = RE_sqrtMAF) # update on 2024-09-16
    
    GRM1$RE_sqrtMAF_pos1 = ResidMat_new$RE_sqrtMAF[match(sparseGRM$ID1, ResidMat_new$SubjID)]
    GRM1$RE_sqrtMAF_pos2 = ResidMat_new$RE_sqrtMAF[match(sparseGRM$ID2, ResidMat_new$SubjID)]
    GRM1 = GRM1 %>% mutate(Cov_R_sqrtMAF_RE_sqrtMAF_part1 = Value * R_sqrtMAF_pos1 * RE_sqrtMAF_pos2)         # update on 2024-09-12
    GRM1 = GRM1 %>% mutate(Cov_R_sqrtMAF_RE_sqrtMAF_part2 = Value * R_sqrtMAF_pos2 * RE_sqrtMAF_pos1)         # update on 2024-09-12
    GRM1 = GRM1 %>% mutate(Cov_RE_sqrtMAF = Value * RE_sqrtMAF_pos1 * RE_sqrtMAF_pos2)         # update on 2024-09-12
    
    Cov_S1_S2 = (GRM1 %>% select(Cov_R_sqrtMAF_RE_sqrtMAF_part1) %>% sum) +
      (GRM1 %>% select(Cov_R_sqrtMAF_RE_sqrtMAF_part2) %>% sum) - 
      sum(g.var.est * E* R^2)
    # Cov_S1_S2 = t(R_sqrtMAF) %*% GRM %*% RE_sqrtMAF
    
    ####
    lambda = Cov_S1_S2/VarS1   # lambda = Cov/VarS1
    
    S_GxE = S2 - lambda * S1                                # new test statistic for marginal GxE effect
    R.new = (E - lambda) * R                                # new residuals
    R.new_sqrtMAF = R.new * sqrt(g.var.est)  # update on 2024-09-10
    ResidMat_new = ResidMat_new %>% mutate(R.new_sqrtMAF = R.new_sqrtMAF) # update on 2024-09-16
    
    GRM1$R.new_sqrtMAF_pos1 = ResidMat_new$R.new_sqrtMAF[match(sparseGRM$ID1, ResidMat_new$SubjID)]
    GRM1$R.new_sqrtMAF_pos2 = ResidMat_new$R.new_sqrtMAF[match(sparseGRM$ID2, ResidMat_new$SubjID)]
    GRM1 = GRM1 %>% mutate(Cov_R.new_sqrtMAF = Value * R.new_sqrtMAF_pos1 * R.new_sqrtMAF_pos2)         # update on 2024-09-12
    
    S_GxE.mean = 2 * sum(R.new * MAF.est)
    S_GxE.var = (GRM1 %>% select(Cov_R.new_sqrtMAF) %>% sum) * 2 - sum((R.new)^2 * g.var.est)
    
    # S_GxE.var = (GRM1 %>% select(Cov_RE_sqrtMAF) %>% sum) * 2 - sum((R*E)^2 * g.var.est)
    # t(R.new_sqrtMAF) %*% GRM %*% R.new_sqrtMAF
    
    z_GxE = (S_GxE - S_GxE.mean)/sqrt(S_GxE.var)
    
    if(abs(z_GxE) < Cutoff){
      pval.norm = pnorm(abs(z_GxE), lower.tail = FALSE)*2
      pval.output = c(pval.norm, pval.norm, pval.norm, pval.norm, pval.norm1) # 5 elements: SPAGxEmix, SPAGxEmixCCT_Wald, SPAGxEmixCCT_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }else{
      
      # updated on 2024-09-10
      S_GxE.var.SPAGxEmix = sum(R.new^2 * g.var.est)
      Var.ratio = S_GxE.var.SPAGxEmix / S_GxE.var
      
      # updated on 2023-10-24
      pval1 = GetProb_SPA_G_new_adjust(init.t = 0, MAF.est, R.new, max(S_GxE, (2*S_GxE.mean-S_GxE)), Var.ratio, lower.tail = FALSE) # SPA-G p value
      pval2 = GetProb_SPA_G_new_adjust(init.t = 0, MAF.est, R.new, min(S_GxE, (2*S_GxE.mean-S_GxE)), Var.ratio, lower.tail = TRUE)  # SPA-G p value
      
      # updated on 2023-10-24
      if(is.na(pval2)){
        pval2 = pval1
      }
      
      pval3 = pnorm(abs(z_GxE), lower.tail = FALSE) # Normal
      pval4 = pnorm(-abs(z_GxE), lower.tail = TRUE) # Normal
      
      pval.spaG = pval1 + pval2
      pval.norm = pval3 + pval4
      
      pval.output = c(pval.spaG, pval.spaG, pval.spaG, pval.norm, pval.norm1)
    } # 5 elements: SPAGxEmix+, SPAGxEmix_Wald+, SPAGxEmixCCT+, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
  }else{
    
    W = cbind(1, g)
    
    # updated on 2023-12-27
    R0 = R - (W)%*%(solve(t(W)%*%W))%*%(t(W)%*%R) # null model residuals adjusting for G
    
    S_GxE0 = sum(g*E*R0)  # test statistic for marginal GxE effect
    R.new0 = E*R0
    
    ####
    R.new0_sqrtMAF =  R.new0 * sqrt(g.var.est)  # update on 2024-09-16
    ResidMat_new = ResidMat_new %>% mutate(R.new0_sqrtMAF = R.new0_sqrtMAF) # update on 2024-09-16
    
    GRM1$R.new0_sqrtMAF_pos1 = ResidMat_new$R.new0_sqrtMAF[match(sparseGRM$ID1, ResidMat_new$SubjID)]
    GRM1$R.new0_sqrtMAF_pos2 = ResidMat_new$R.new0_sqrtMAF[match(sparseGRM$ID2, ResidMat_new$SubjID)]
    GRM1 = GRM1 %>% mutate(Cov_R.new0_sqrtMAF = Value * R.new0_sqrtMAF_pos1 * R.new0_sqrtMAF_pos2)         # update on 2024-09-12

    S_GxE0.mean = 2 * sum(R.new0 * MAF.est)
    S_GxE0.var = (GRM1 %>% select(Cov_R.new0_sqrtMAF) %>% sum) * 2 - sum(R.new0^2 * g.var.est)
    # t(R.new0_sqrtMAF) %*% GRM %*% R.new0_sqrtMAF
    
    z_GxE0 = (S_GxE0 - S_GxE0.mean)/sqrt(S_GxE0.var)
    
    if(abs(z_GxE0) < Cutoff){
      pval.norm = pnorm(abs(z_GxE0), lower.tail = FALSE)*2
      pval.spaGxE = pval.norm
    }else{
      
      S_GxE0.var.SPAGxEmix = sum(R.new0^2 * g.var.est)
      Var.ratio = S_GxE0.var.SPAGxEmix / S_GxE0.var
      
      # update on 2024-09-16
      pval1 = GetProb_SPA_G_new_adjust(init.t = 0, MAF.est, R.new0, max(S_GxE0, (2*S_GxE0.mean-S_GxE0)), Var.ratio, lower.tail = FALSE) # SPA-G p value
      pval2 = GetProb_SPA_G_new_adjust(init.t = 0, MAF.est, R.new0, min(S_GxE0, (2*S_GxE0.mean-S_GxE0)), Var.ratio, lower.tail = TRUE)  # SPA-G p value
      
      # updated on 2023-10-24
      if(is.na(pval2)){
        pval2 = pval1
      }
      
      pval3 = pnorm(abs(z_GxE0), lower.tail = FALSE) # Normal
      pval4 = pnorm(-abs(z_GxE0), lower.tail = TRUE) # Normal
      
      pval.spaGxE = pval1 + pval2
      pval.norm = pval3 + pval4
    }
    
    # Wald test

    if(traits=="binary"){
      # Phen.mtx.new = Phen.mtx %>% mutate(ID = rownames(Phen.mtx))
      # data00 = glm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E, family = "binomial") # updated on 2023-12-27
      data0 = glmer(Y ~ as.matrix(Cova.mtx)+E+g+g*E + (1 | IID), family = "binomial", data = Phen.mtx) # updated on 2024-09
      
      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|z|)"]     # Wald test p value for marginal GxE effect
      
      # updated on 2024-04-21
      if(is.na(pval.wald)){
        pval.wald = pval.spaGxE
      }
      
      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))
      
      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix, SPAGxEmixCCT_Wald, SPAGxEmixCCT_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }
    
    if(traits=="quantitative"){
      
      # data00 = lm(formula = Phen.mtx$Y ~ as.matrix(Cova.mtx)+E+g+g*E) # updated on 2023-12-27
      data0 = lmer(Y ~ as.matrix(Cova.mtx)+E+g+g*E + (1 | IID), data = Phen.mtx) # updated on 2023-12-27
      
      pval.wald = summary(data0)$coefficients["E:g", "Pr(>|t|)"]     # Wald test p value for marginal GxE effect
      
      # updated on 2024-04-21
      if(is.na(pval.wald)){
        pval.wald = pval.spaGxE
      }
      
      # use CCT to conbine p values from SPAGxE and Wald tests in the case of pval.norm1 < epsilon
      pval_SPAGxE_Wald_CCT = CCT(c(pval.spaGxE,
                                   pval.wald))
      
      # updated on 2023-12-27
      pval.output = c(pval.spaGxE, pval.wald, pval_SPAGxE_Wald_CCT, pval.norm, pval.norm1) # 5 elements: SPAGxEmix, SPAGxEmixCCT_Wald, SPAGxEmixCCT_CCT, normal approximation p-value for marginal GxE effect, and normal approximation p value for marginal genetic effect
    }
    
  }
  
  output.one.snp = c(MAF, missing.rate, pval.output, S1, VarS1, Z1,
                     MAF.est.negative.num, MAC)
  
  return(output.one.snp)
}





#### other functions

GetProb_SPA_G_new_adjust = function(init.t = 0, MAFVec, R, s, Var.ratio, lower.tail){
  
  # out = uniroot(H1_adj, c(-10,10), extendInt = "upX",
  #               R=R, s=s, MAFVec=MAFVec)
  out = fastgetroot_H1_new(init.t, R = R, s = s, MAFVec = MAFVec)
  
  zeta = out$root
  
  k1 = H_org(zeta, R=R, MAFVec=MAFVec)
  k2 = H2(zeta, R=R, MAFVec=MAFVec)
  
  temp1 = zeta * s - k1
  
  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  
  a = w + 1/w * log(v/w)
  
  # pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  pval = pnorm(a * sqrt(Var.ratio), lower.tail = lower.tail)
  
  if(is.na(pval)){
    pval =  0
  }
  
  re = pval
  return(re)
}


# # updataed on 2023-10-24
# fastgetroot_H1_new = function(init.t,
#                               R,          # residuals
#                               s,          # observed test statistic
#                               MAFVec,
#                               tol = .Machine$double.eps^0.25,
#                               maxiter = 100)
# {
#   t = init.t;
#   H1_adj_eval = 0
#   diff.t = Inf
#   converge = T
#   for(iter in 1:maxiter){
#     old.t = t
#     old.diff.t = diff.t
#     old.H1 = H1_adj_eval
#     H1_adj_eval = H1_adj(t, R, s, MAFVec)
#     H2_eval = H2(t, R, MAFVec)
#     
#     diff.t = -1 * H1_adj_eval / H2_eval
#     
#     # updated on 2023-10-24
#     if(is.na(diff.t)){
#       diff.t =  5
#     }
#     
#     if(is.na(H1_adj_eval) | is.infinite(H1_adj_eval)){
#       # checked it on 07/05:
#       # if the solution 't' tends to infinity, 'K2_eval' tends to 0, and 'K1_eval' tends to 0 very slowly.
#       # then we can set the one sided p value as 0 (instead of setting converge = F)
#       t = sign(s)*Inf
#       H2_eval = 0;
#       break;
#     }
#     if(sign(H1_adj_eval) != sign(old.H1)){
#       while(abs(diff.t) > abs(old.diff.t) - tol){
#         diff.t = diff.t/2
#       }
#     }
#     if(is.na(diff.t)){
#       diff.t =  5
#     }
#     if(abs(diff.t) < tol) break;
#     t = old.t + diff.t
#   }
#   
#   # updated on 2023-10-24
#   
#   if(iter == maxiter) {
#     converge = F
#     t = NA}
#   
#   return(list(root = t,
#               iter = iter,
#               converge = converge,
#               H2_eval = H2_eval))
# }

