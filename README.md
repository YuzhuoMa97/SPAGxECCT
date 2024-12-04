# SPAGxE<sub>CCT</sub> 
A scalable and accurate framework for large-scale genome-wide gene-environment interaction (G×E) analysis.
## Software dependencies and operating systems
The package has been tested under linux and windows systems. (The R package will be rewritten with RCPP codes and support genotype data in PLINK format soon!)
## How to install and load this package
```
library(devtools)  # author version: 2.4.5
install_github("YuzhuoMa97/SPAGxECCT")
library(SPAGxECCT)
?SPAGxECCT  # manual of SPAGxECCT package
```
Current version is 1.0.4. For older version and version update information, plesase refer to OldVersions.  

Please do not hesitate to contact me (yuzhuoma@stu.pku.edu.cn) if you meet any problem. Suggestions or comments are also welcome.

## Introduction of SPAGxE<sub>CCT</sub> 

**SPAGxE<sub>CCT</sub> is a G×E analysis framework designed for a wide range of complex traits with intricate structures, including time-to-event, ordinal categorical, binary, quantitative, longitudinal, and other complex traits.** The framework involves two main steps:

- Step 1: SPAGxE<sub>CCT</sub> fits a covariates-only model to calculate model residuals. These covariates include, but are not limited to, confounding factors such as age, sex, SNP-derived principal components (PCs), and environmental factors. The specifics of the model and residuals vary depending on the trait type. Since the covariates-only model is genotype-independent, it only needs to be fitted once across a genome-wide analysis.

- Step 2: SPAGxE<sub>CCT</sub> identifies genetic variants with marginal G×E effects on the trait of interest. First, marginal genetic effects are tested using score statistics. If the marginal genetic effect is not significant, S<sub>G×E</sub> is used as the test statistic to characterize the marginal G×E effect. If significant, S<sub>G×E</sub> is updated to genotype-adjusted test statistics. To balance computational efficiency and accuracy, SPAGxE<sub>CCT</sub> employs a hybrid strategy combining normal distribution approximation and SPA to calculate p-values, as used in previous studies such as [SAIGE](https://saigegit.github.io/SAIGE-doc/) and [SPAGE](https://github.com/WenjianBI/SPAGE). For variants with significant marginal genetic effects, SPAGxE<sub>CCT</sub> additionally calculates p value through Wald test and uses Cauchy combination (CCT) to combine p values from Wald test and the proposed genotype-adjusted test statistics.



![plot](https://github.com/YuzhuoMa97/SPAGxECCT/blob/main/workflow/workflow_SPAGxECCT_MYZ.png)

## Introduction of SPAGxEmix<sub>CCT</sub>

As an extension of SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub> is a G×E analysis framework which is applicable to include individuals from multiple ancestries or multi-way admixed populations. 

SPAGxE<sub>CCT</sub> assumes that genotypes across different individuals follow an identical binomial distribution, a valid assumption in homogeneous populations. However, this assumption may not hold in study cohorts comprising individuals from diverse ancestries. To address this limitation, SPAGxEmix<sub>CCT</sub> modifies the approach to allow for different allele frequencies while maintaining binomial distributions for genotypes. SPAGxEmix<sub>CCT</sub> estimates individual-level allele frequencies using SNP-derived principal components (PCs) and raw genotypes. Similar to SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub> involves two main steps:

- Step 1: SPAGxEmix<sub>CCT</sub> fits a genotype-independent (covariates-only) model and calculates the model residuals. Detailed information is provided in Step 1 of SPAGxE<sub>CCT</sub>.

- Step 2: SPAGxEmix<sub>CCT</sub> identifies genetic variants with marginal G×E effects on the trait of interest. It first estimates the individual-level allele frequencies of the tested variants using SNP-derived PCs and raw genotypes. Next, SPAGxEmix<sub>CCT</sub> evaluates marginal genetic effects using score statistics. If the marginal genetic effect is not significant, S<sub>G×E(mix)</sub> is used as the test statistic to characterize the marginal G×E effect. If significant, S<sub>G×E(mix)</sub> is updated to genotype-adjusted test statistics. The hybrid strategy to balance computational efficiency and accuracy follows the approach used in SPAGxE<sub>CCT</sub>.

Admixed populations are routinely excluded from genomic studies due to concerns over population structure. The admixed population analysis is technically challenging as it requires addressing the complicated patterns of genetic and phenotypic diversities that arise from distinct genetic backgrounds and environmental exposures. As the genetic ancestries are increasingly recognized as continuous rather than discrete, SPAGxEmix<sub>CCT</sub> estimates individual-level allele frequencies to characterize the individual-level genetic ancestries using information from the SNP-derived PCs and raw genotype data in a model-free approach. Thus, it does not necessitate accurate specification of and the availability of appropriate reference population panels for the ancestries contributing to the individual, which might be unknown or not well defined. In addition, SPAGxEmix<sub>CCT</sub> is not sensitive to model misspecification (e.g. missed or biased confounder-trait associations) and trait-based ascertainment.

![plot](https://github.com/YuzhuoMa97/SPAGxECCT/blob/main/workflow/workflow_SPAGxEmixCCT_MYZ.png)

## Introduction of SPAGxEmix<sub>CCT-local</sub>
SPAGxEmix<sub>CCT-local</sub> is a G×E analysis framework designed to efficiently and accurately identify ancestry-specific G×E effects by incorporating local ancestry in multi-way admixed populations.

SPAGxEmix<sub>CCT-local</sub> extends SPAGxEmix<sub>CCT</sub> by integrating local ancestry information to enhance the precision of ancestry-specific G×E effect testing and statistical powers. Furthermore, we introduce SPAGxEmix<sub>CCT-local-global</sub>, which combines p-values from both SPAGxEmix<sub>CCT</sub> and SPAGxEmix<sub>CCT-local</sub>, offering an optimal unified approach for various cross-ancestry genetic architectures. As with SPAGxEmix<sub>CCT</sub>, SPAGxEmix<sub>CCT-local</sub> involves two main steps:

- Step 1: SPAGxEmix<sub>CCT-local</sub> fits a genotype-independent (covariates-only) model and calculates the model residuals. Details are described in Step 1 of SPAGxE<sub>CCT</sub>.
  
- Step 2: SPAGxEmix<sub>CCT-local</sub> identifies genetic variants with marginal G×E effects on the trait of interest. First, SPAGxEmix<sub>CCT</sub> estimates the ancestry-specific allele frequencies for the variants being tested. Next, SPAGxEmix<sub>CCT-local</sub> evaluates ancestry-specific marginal genetic effects using ancestry-specific score statistics. If an ancestry-specific marginal genetic effect is not significant, we use ancestry-specific S<sub>G×E(mix)</sub> as the test statistic to characterize the ancestry-specific marginal G×E effect. If significant, S<sub>G×E(mix)</sub> is updated to ancestry-specific genotype-adjusted test statistics. The hybrid strategy to balance computational efficiency and accuracy follows the approach used in SPAGxE<sub>CCT</sub>.

**Compared to traditional statistical testing methods that account for local ancestry, SPAGxEmix<sub>CCT-local</sub> offers much greater computational efficiency.** SPAGxEmix<sub>CCT-local</sub> leverages local ancestry information to estimate the distribution of ancestry-specific genotypes and the null distribution of test statistics. For most tests (ancestry-specific genetic main effect p-values > 0.001) in a genome-wide G×E analysis, SPAGxEmix<sub>CCT-local</sub> utilizes residuals from a genotype-independent model (fitted only once across the genome-wide analysis) to construct test statistics. In contrast, traditional methods must incorporate local ancestry as covariates and fit a separate model for each variant, which is computationally burdensome for genome-wide analyses.

## Introduction of SPAGxE+

SPAGxE+ is a scalable and accurate G×E analytical framework that controls for sample relatedness and unbalanced phenotypic distribution (e.g., case-control imbalance for binary trait). 
It is applicable to a wide range of complex traits with intricate structures, including binary, quantitative, time-to-event, ordinal categorical, longitudinal, and other complex traits. The framework involves two main steps:

- Step 1: SPAGxE+ fits a covariates-only model to calculate model residuals. These covariates include, but are not limited to, confounding factors such as age, sex, SNP-derived principal components (PCs), and environmental factors. The specifics of the model and residuals vary depending on the trait type. Since the covariates-only model is genotype-independent, it only needs to be fitted once across a genome-wide analysis. Incorporating random effects to account for sample relatedness in null model fitting is optional, rather than required.

- Step 2: SPAGxE+ identifies genetic variants with marginal G×E effects on the trait of interest. First, marginal genetic effects are tested using score statistics. If the marginal genetic effect is not significant, S<sub>G×E</sub> is used as the test statistic to characterize the marginal G×E effect. If significant, S<sub>G×E</sub> is updated to genotype-adjusted test statistics. To balance computational efficiency and accuracy, SPAGxE+ employs a hybrid strategy combining normal distribution approximation and SPA to calculate p-values, as used in previous studies such as [SAIGE](https://saigegit.github.io/SAIGE-doc/) and [SPAGE](https://github.com/WenjianBI/SPAGE). 


Currently, there is still a lack of scalable and accurate gene-environmental interaction analytical frameworks that can control for sample relatedness and be applicable to binary, time-to-event, ordinal categorical, or longitudinal trait analysis. It is usually challenging to detect gene-environmental interaction effects efficiently and accurately while adjusting for sample relateness in a large-scale genome-wide analysis of binary, time-to-event, ordinal categorical, or longitudinal trait. In GWAS for thousands of phenotypes in large biobanks, most binary traits have substantially fewer cases than controls. Existing methods based on logistic mixed model produce inflated type I errors in the presence of unbalanced case-control ratios. 

SPAGxE+ is a scalable and accurate G×E analytical framework that uses saddlepoint approximation to calibrate the null distribution of test statistics while controlling for sample relatedness and case-control imbalance. Compared to SPAGxE, the sparse genetic relationship matrix (GRM) is used for characterizing familial structure in step 2. SPAGxE+ provides accurate p-values even when case-control ratios are extremely unbalanced. Additionally, SPAGxE+ is applicable to other complex traits including time-to-event, ordinal categorical, and longitudinal traits, and it maintains highly accurate even when phenotypic distribution is unbalanced.



## Introduction of SPAGxEmix<sub>CCT</sub>+

As an extension of SPAGxEmix<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub>+ is a scalable and accurate G×E analytical framework. SPAGxEmix<sub>CCT</sub>+ accommodates individuals from multiple ancestries or multi-way admixed populations while controlling for both population structure and familial relatedness. Admixed individuals can be analyzed within a cohort alone or alongside homogeneous groups using this framework.


SPAGxEmix<sub>CCT</sub>+ comprises three main steps:

- Step 0 (a): SPAGxEmix<sub>CCT</sub>+ employs PC-AiR (Conomos et al., 2015, Gen Epi) to compute ancestry-representative principal components (PCs) that capture distant genetic relatedness, such as population structure.
  
- Step 0 (b): SPAGxEmix<sub>CCT</sub>+utilizes PC-Relate (Conomos et al., 2016, AJHG) to estimate an ancestry-adjusted sparse GRM or sparse kinship coefficient matrix, representing recent genetic relatedness.
  
- Step 0 (c): Iterations of Step 0 (a) and Step 0 (b) refine the inference of both population structure (via PC-AiR) and recent genetic relatedness (via PC-Relate).

- Step 1: SPAGxEmix<sub>CCT</sub>+ fits a genotype-independent (covariates-only) model and computes residuals, as detailed in Step 1 of SPAGxE<sub>CCT</sub>. Incorporating random effects to account for sample relatedness during null model fitting is optional.

- Step 2: SPAGxEmix<sub>CCT</sub>+ identifies genetic variants with marginal G×E effects on the trait of interest. It first estimates the individual-level allele frequencies of the tested variants using SNP-derived PCs (from Step 0) and raw genotypes. Next, SPAGxEmix<sub>CCT</sub>+ evaluates marginal genetic effects using score statistics. If the marginal genetic effect is not significant, S<sub>G×E(mix)</sub>+ is used as the test statistic to characterize the marginal G×E effect. If significant, S<sub>G×E(mix)</sub>+ is updated with genotype-adjusted test statistics. The hybrid strategy to balance computational efficiency and accuracy follows the approach used in SPAGxE<sub>CCT</sub>.


The Genetic Relationship Matrix (GRM) contains information on both ancestry and familial relatedness. As introduced in the PC-Relate paper, the empirical GRM has been widely used for inferring population structure (distant genetic relatedness) in samples without close relatives, as well as for estimating recent kinship and heritability of complex traits in single-population studies. However, in multi-ancestry or admixed population studies, a thresholded GRM cannot be sparsified effectively. To address this, we employ PC-AiR and PC-Relate (implemented in the GENESIS R package) to obtain principal components (representing distant genetic relatedness) and sparse kinship coefficients (representing recent genetic relatedness).

PC-AiR is specifically designed for robust inference of population structure in the presence of recent genetic relatedness—whether known or cryptic—among sampled individuals. It performs Principal Component Analysis (PCA) on genome-wide SNP data to detect population structure, while accounting for relatedness to avoid confounding ancestry inference with familial relationships. Unlike standard PCA, PC-AiR adjusts for relatedness, allowing accurate estimation of ancestry that is unaffected by family structure. PC-Relate, on the other hand, uses ancestry-representative principal components to adjust for population structure/ancestry and provide precise estimates of recent genetic relatedness.

To analyze individuals from multi-ancestry or multi-way admixed populations with familial relatedness, SPAGxEmix<sub>CCT</sub>+ first uses PC-AiR to obtain principal components representing population structure. It then applies PC-Relate to regress out the PCs from the genotypes, yielding ancestry-adjusted genotypes and kinship estimates to assess family structure (i.e., recent genetic relatedness). The PC-Relate paper demonstrates that an iterative process alternating between PC-AiR and PC-Relate can enhance inference for both population structure (via PC-AiR) and recent genetic relatedness (via PC-Relate). Typically, two iterations are sufficient to produce accurate ancestry-adjusted principal components and a sparse GRM.

Steps 1 and 2 of the SPAGxEmix<sub>CCT</sub>+ analysis are similar to those in SPAGxEmix<sub>CCT</sub>. In Step 2, the ancestry-adjusted sparse GRM is crucial for characterizing familial structure, and SPA is used to calibrate p-values.

**SPAGxEmix<sub>CCT</sub>+ is a scalable and accurate G×E analytical framework designed to control for both population structure and familial relatedness across diverse population and family structures. It is especially suited for cohorts that include individuals from multi-ancestry or multi-way admixed populations, effectively handling ancestry-specific minor allele frequencies (MAFs) and ancestry-specific case-control ratios (or other ancestry-specific phenotypic distribution such as event rates for time-to-event traits).**

## The computational efficiency and statictical power can be enhanced through incorporating polygenic scores (PGSs) as covariates with fixed effects.

Recent studies have demonstrated that adjusting for polygenic scores (PGSs) can account for polygenic effects and enhance statistical power. In SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub>, and SPAGxEmix<sub>CCT-local</sub>, incorporating PGSs as fixed-effect covariates significantly improves computational efficiency. When accurate PGSs are available, a genotype-independent model can be used for genome-wide analyses, enabling the application of regular score statistics followed by a hybrid test employing normal approximation and SPA. By including PGSs as covariates, we eliminate the need to construct statistics through matrix projection or linear regression with genotype data, further optimizing computational efficiency. Additionally, PGS adjustment boosts statistical power by accounting for polygenic effects, complementing the reduced power typically associated with sparse-GRM methods.





Our approach implements a two-stage strategy, referred to as SPAGxE<sub>CCT</sub>-PGS, SPAGxEmix<sub>CCT</sub>-PGS, and SPAGxEmix<sub>CCT-local</sub>-PGS:

- Stage 1: Perform an initial genome-wide association study (GWAS) using tools like SAIGE, and calculate the Leave-One-Chromosome-Out (LOCO) PGS based on summary statistics. Alternatively, use PGSs from external databases like the PGS Catalog.

- Stage 2: Incorporate the LOCO-PGS as an additional covariate in genome-wide G×E analyses using SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub>, or SPAGxEmix<sub>CCT-local</sub>.

This general strategy has not yet been fully automated into a function. It's important to note that PGS accuracy varies across genetic ancestries. For multi-ancestry or admixed populations, accurate PGSs must be constructed using methods tailored to these groups, such as the PRS<sub>multi</sub> method for multi-ancestry individuals or the GAUDI method for admixed populations. When testing for ancestry-specific G×E effects using SPAGxEmix<sub>CCT-local</sub>, it is essential to use approaches like the pPRS method, which disentangles ancestry mosaics through local ancestry inference before constructing PGSs.

This strategy not only improves computational efficiency but also enhances statistical power. As part of SPAGxE reproducibility efforts, we will demonstrate how to construct LOCO-PGSs using the clumping and thresholding (C+T) method based on summary statistics. We also plan to develop a function to approximate marginal genetic effects, facilitating the construction of PGSs.



## Rank-based inverse normal transformation to enhance statistical power

Model residuals from a fitted null model can be highly unbalanced, especially for longitudinal phenotypes. In such cases, we find that applying a rank-based inverse normal transformation (INT) to the residuals significantly improves empirical power. We propose SPAGxE<sub>CCT(INT)</sub>, SPAGxEmix<sub>CCT(INT)</sub>, and SPAGxEmix<sub>CCT-local(INT)</sub>, where the inverse normal transformed residuals are used as input. The p-values from SPAGxE<sub>CCT(INT)</sub>, SPAGxEmix<sub>CCT(INT)</sub>, or SPAGxEmix<sub>CCT-local(INT)</sub> can then be combined with those from their untransformed counterparts (SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub>, or SPAGxEmix<sub>CCT-local</sub>) using the Cauchy Combination Test (CCT).

# Reproducibility

Scripts to reproduce the experiments performed for the manuscript:

 **A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits (to be updated)**

## UK Biobank data analysis 

In the paper **A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits (to be updated)**, we have applied SPAGxE<sub>CCT</sub> and SPAGxEmix<sub>CCT</sub> to analyze time-to-event traits in UB Biobank. For the SPAGxE<sub>CCT</sub> analyses, 281,299 White British individuals were included. For the SPAGxEmix<sub>CCT</sub> analyses, 338,044 individuals from all ancestries were included. As a universal analysis framework, we also evaluated the performance of SPAGxE<sub>CCT</sub> in ordinal categorical, binary, and quantitative trait analysis.  

**Summary statistics of time-to-event and binary traits in UK Biobank is available [here](https://zenodo.org/records/11571404).**

## Simulation studies

The primary codes used in the SPAGxECCT project's simulation studies are available in the **simulation studies** folder and have also been uploaded to **Code Ocean**.

These codes cover GxE analyses for time-to-event, binary, quantitative, and ordinal categorical traits. With the provided codes, users can reproducibly conduct simulation studies of GxE analyses, including generating genotype data, phenotype data, calculating genetic principal components (PCs), and performing p-value calculations using our proposed methods of SPAGxE<sub>CCT</sub> (for testing GxE effects in homogeneous populations), SPAGxEmix<sub>CCT</sub> (for testing GxE effects in heterogeneous or admixed populations), and SPAGxEmix<sub>CCT-local</sub> (for testing ancestry-specific GxE effects by incorporating local ancestry into analyses of admixed populations).

### 1. For G×E analyses of homogeneous populations, please assess the **simulation studies/SPAGxECCT** folder. 

We simulated a homogeneous population. We simulated binary phenotypes following logistic regression model, quantitative phenotypes following linear regression model, time-to-event phenotypes following Cox PH model, and ordinal categorical phenotypes following proportional odds logistic model.

```
In simulation studies\SPAGxECCT\code\:

SPAGxECCT_binary_typeIerror_simulation_normal_envi.R   # R script to evaluate type I error rates of SPAGxECCT including simulating genotypes, binary phenotypes, and calculating p-values.

SPAGxECCT_survival_typeIerror_simulation_normal_envi.R   # R script to evaluate type I error rates of SPAGxECCT including simulating genotypes, time-to-event phenotypes, and calculate p-values.

SPAGxECCT_categorical_typeIerror_simulation_normal_envi.R   # R script to evaluate type I error rates of SPAGxECCT including simulating genotypes, ordinal categorical phenotypes, and calculate p-values.

SPAGxECCT_binary_power_simulation_normal_envi.R   # R script to evaluate powers of SPAGxECCT including simulating genotypes, binary phenotypes, and calculating p-values.

SPAGxECCT_survival_power_simulation_normal_envi.R   # R script to evaluate powers of SPAGxECCT including simulating genotypes, time-to-event phenotypes, and calculate p-values.

SPAGxECCT_categorical_power_simulation_normal_envi.R   # R script to evaluate powers of SPAGxECCT including simulating genotypes, ordinal categorical phenotypes, and calculate p-values.

```


### 2. For G×E analyses of heterogeneous or admixed populations, please refer to the **simulation studies/SPAGxEmixCCT** folder. 

We simulated a two-way admixed population from EUR ancestry and EAS ancestry with sample size n=10000. We assumed the first 5000 individuals were from a EUR-dominant community and the remaining 5000 individuals were from an EAS-dominant community. In our paper, we used the real MAF values from 1000 Genome Projects to mimic the allele frequency diversity between EUR and EAS. We simulated binary phenotypes with ancestry-specific case-control ratios and time-to-event phenotypes with ancestry-specific event rates.


```
In simulation studies\SPAGxEmixCCT\code\:

SPAGxEmixCCT_AdmixedPopulationAnalysis_binary_simulation_normal_envi.R   # R script to evaluate type I error rates or powers of SPAGxEmixCCT including simulating genotypes, binary phenotypes, and calculating p-values in admixed population analyses.

SPAGxEmixCCT_AdmixedPopulationAnalysis_survival_simulation_normal_envi.R   # R script to evaluate type I error rates or powers of SPAGxEmixCCT including simulating genotypes, time-to-event phenotypes, and calculating p-values in admixed population analyses.

```


### 3. For G×E analyses of using SPAGxEmix<sub>CCT-local</sub> to identify ancestry-specific G×E effects by incorporating local ancestry in analyses of admixed populations, please refer to the **simulation studies/SPAGxEmixCCT-local** folder. 

We simulated a two-way admixed population with sample size n=10000, including ancestry-specific genotypes, local ancestry counts, SNP-derived PCs, and phenotypes.

```
In simulation studies\SPAGxEmixCCT-local\code\:

SPAGxEmixCCT_local_ancestry_binary_power_simulation_normal_envi.R   # R script to evaluate type I error rates or powers of SPAGxEmixCCT-local including simulating local ancestry counts, ancestry-specific genotypes, binary phenotypes, and calculating ancestry-specific p-values corresponding to ancestry-specific marginal GxE effects in admixed population analyses.

SPAGxEmixCCT_local_ancestry_quantitative_power_simulation_normal_envi.R   # R script to evaluate type I error rates or powers of SPAGxEmixCCT-local including simulating local ancestry counts, ancestry-specific genotypes, quantitative phenotypes, and calculating ancestry-specific p-values corresponding to ancestry-specific marginal GxE effects in admixed population analyses.

generate_PCs_forLocalanceAnalysis.R    # R script to generate SNP-derived PCs

```




# 


# Reference
See **A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits** (to be updated) for more details about SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub>, and SPAGxEmix<sub>CCT-local</sub>.







