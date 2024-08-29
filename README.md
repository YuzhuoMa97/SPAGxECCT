# SPAGxE<sub>CCT</sub> 
A scalable and accurate framework for large-scale genome-wide gene-environment interaction (G×E) analysis.
## Software dependencies and operating systems
The package has been tested under linux and windows systems.
## How to install and load this package
```
library(devtools)  # author version: 2.4.5
install_github("YuzhuoMa97/SPAGxECCT")
library(SPAGxECCT)
?SPAGxECCT  # manual of SPAGxECCT package
```
Current version is 0.1.0. For older version and version update information, plesase refer to OldVersions.  

Please do not hesitate to contact me (yuzhuoma@stu.pku.edu.cn) if you meet any problem. Suggestions or comments are also welcome.

## Introduction of SPAGxE<sub>CCT</sub> 
SPAGxE<sub>CCT</sub> is a G×E analysis framework which is applicable to a wide variety of complex traits with intricate structures (e.g. time-to-event, ordinal categorical, binary, quantitative, longitudinal traits, and other complex traits). SPAGxE<sub>CCT</sub> contains two main steps. 

- In step 1, SPAGxE<sub>CCT</sub> fits a covariates-only model and then calculates model residuals. The covariates include, but are not limited to, confounding factors such as age, gender, SNP-derived principal components (PCs), and environmental factors. The model specification and the corresponding model residuals vary depending on the type of trait. As the covariates-only model is genotype-independent, the model fitting and residuals calculation are only required once across a genome-wide analysis.
  
- In step 2, SPAGxE<sub>CCT</sub> identifies genetic variants with marginal G×E effect on the trait of interest. First, SPAGxE<sub>CCT</sub> tests for marginal genetic effect via score statistic. If the marginal genetic effect is not significant, we use S<sub>G×E</sub> as the test statistics to characterize marginal G×E effect. Otherwise, statistics S<sub>G×E</sub> is updated to genotype-adjusted test statistics. To balance the computational efficiency and accuracy, SPAGxE<sub>CCT</sub> employs a hybrid strategy to combine normal distribution approximation and SPA to calculate p values, as in previous studies such as [SAIGE](https://saigegit.github.io/SAIGE-doc/) and [SPAGE](https://github.com/WenjianBI/SPAGE). For variants with significant marginal genetic effect, SPAGxE<sub>CCT</sub> additionally calculates p value through Wald test and then uses Cauchy combination (CCT) to combine p values from Wald test and the proposed genotype-adjusted test statistics.

![plot](https://github.com/YuzhuoMa97/SPAGxECCT/blob/main/workflow/workflow_SPAGxECCT_MYZ.png)

## Introduction of SPAGxEmix<sub>CCT</sub>
As an extension of SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub> is a G×E analysis framework which is applicable to include individuals from multiple ancestries or multi-way admixed populations. 

SPAGxE<sub>CCT</sub> relies on an assumption that the genotypes for different individuals follow an identical binomial distribution. The assumption is usually valid in a homogeneous population. However, if the study cohort consists of individuals from multiple ancestries, this assumption could be violated. To address this issue, we propose SPAGxEmix<sub>CCT</sub> in which genotypes for different individuals still follow binomial distributions but the corresponding allele frequencies could be different. SPAGxEmix<sub>CCT</sub> estimate individual-level allele frequencies using SNP-derived PCs and raw genotypes. Similar as SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub> contains two main steps. 

- In step 1, SPAGxEmix<sub>CCT</sub> fits a covariates-only model and then calculates model residuals. Details can be found in step 1 of SPAGxE<sub>CCT</sub>.
  
- In step 2, SPAGxEmix<sub>CCT</sub> identifies genetic variants with marginal G×E effect on the trait of interest. First, SPAGxEmix<sub>CCT</sub> uses the SNP-derived PCs and raw genotypes to estimate the individual-level allele frequencies of variants being test. Then, SPAGxEmix<sub>CCT</sub> tests for marginal genetic effect via score statistic. If the marginal genetic effect is not significant, we use S<sub>G×E(mix)</sub> as the test statistics to characterize marginal G×E effect. Otherwise, statistics S<sub>G×E(mix)</sub> is updated to genotype-adjusted test statistics. The hybrid strategy to balance the computational efficiency and accuracy is the same as in SPAGxE<sub>CCT</sub>.

## Introduction of SPAGxEmix<sub>CCT-local</sub>
SPAGxEmix<sub>CCT-local</sub> is a G×E analysis framework which is applicable to include local ancestry for analyses of multi-way admixed populations. 

SPAGxEmix<sub>CCT</sub> can further be extended to SPAGxEmix<sub>CCT-local</sub>, which can efficiently identify ancestry-specific GxE effects by incorporating local ancestry. Moreover, we propose SPAGxEmix<sub>CCT-local-global</sub>, which combines the p values from SPAGxEmix<sub>CCT<sub> and SPAGxEmix<sub>CCT-local</sub> and can serve as an optimal unified approach across various cross-ancestry genetic architectures. Similar as SPAGxEmix<sub>CCT-local</sub>, SPAGxEmix<sub>CCT-local<sub> contains two main steps. 

- In step 1, SPAGxEmix<sub>CCT</sub> fits a covariates-only model and then calculates model residuals. Details can be found in step 1 of SPAGxE<sub>CCT</sub>.
  
- In step 2, SPAGxEmix<sub>CCT</sub> identifies genetic variants with marginal G×E effect on the trait of interest. First, SPAGxEmix<sub>CCT</sub> estimates the ancestry-specific allele frequencies of variants being test. Then, SPAGxEmix<sub>CCT-local</sub> tests for ancestry-specific marginal genetic effects via ancestry-specific score statistics. If the ancestry-specific marginal genetic effect is not significant, we use ancestry-specific S<sub>G×E(mix)</sub> as the test statistics to characterize marginal G×E effect. Otherwise, statistics S<sub>G×E(mix)</sub> is updated to genotype-adjusted test statistics. The hybrid strategy to balance the computational efficiency and accuracy is the same as in SPAGxE<sub>CCT</sub>.
  

  
## UK Biobank data analysis results

In the paper **A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits (to be updated)**, we have applied SPAGxE<sub>CCT</sub> to analyze time-to-event traits in UB Biobank. For the SPAGxE<sub>CCT</sub> analyses, 281,299 White British individuals were included. For the SPAGxEmix<sub>CCT</sub> analyses, 338,044 individuals from all ancestries were included. As a universal analysis framework, we also evaluated the performance of SPAGxE<sub>CCT</sub> in ordinal categorical, binary, and quantitative trait analysis.  

**Summary statistics of time-to-event traits in UK Biobank is available [here](https://zenodo.org/records/11571404).**

## Reference
See **A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits** (to be updated) for more details about SPAGxE<sub>CCT</sub> and SPAGxEmix<sub>CCT</sub>.







