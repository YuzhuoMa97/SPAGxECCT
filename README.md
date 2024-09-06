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
Current version is 0.2.3. For older version and version update information, plesase refer to OldVersions.  

Please do not hesitate to contact me (yuzhuoma@stu.pku.edu.cn) if you meet any problem. Suggestions or comments are also welcome.

## Introduction of SPAGxE<sub>CCT</sub> 

**SPAGxE<sub>CCT</sub> is a G×E analysis framework designed to handle a wide range of complex traits with intricate structures, including time-to-event, ordinal categorical, binary, quantitative, longitudinal, and other complex traits.** The framework involves two main steps:

Step 1: SPAGxE<sub>CCT</sub> fits a covariates-only model to calculate model residuals. These covariates include, but are not limited to, confounding factors such as age, gender, SNP-derived principal components (PCs), and environmental factors. The specifics of the model and residuals vary depending on the trait type. Since the covariates-only model is genotype-independent, it only needs to be fitted once across a genome-wide analysis.

Step 2: SPAGxE<sub>CCT</sub> identifies genetic variants with marginal G×E effects on the trait of interest. First, marginal genetic effects are tested using score statistics. If the marginal genetic effect is not significant, S<sub>G×E</sub> is used as the test statistic to characterize the marginal G×E effect. If significant, S<sub>G×E</sub> is updated to genotype-adjusted test statistics. To balance computational efficiency and accuracy, SPAGxE<sub>CCT</sub> employs a hybrid strategy combining normal distribution approximation and SPA to calculate p-values, as used in previous studies such as [SAIGE](https://saigegit.github.io/SAIGE-doc/) and [SPAGE](https://github.com/WenjianBI/SPAGE). For variants with significant marginal genetic effects, SPAGxE<sub>CCT</sub> additionally calculates p value through Wald test and uses Cauchy combination (CCT) to combine p values from Wald test and the proposed genotype-adjusted test statistics.



![plot](https://github.com/YuzhuoMa97/SPAGxECCT/blob/main/workflow/workflow_SPAGxECCT_MYZ.png)

## Introduction of SPAGxEmix<sub>CCT</sub>

As an extension of SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub> is a G×E analysis framework which is applicable to include individuals from multiple ancestries or multi-way admixed populations. 

SPAGxE<sub>CCT</sub> assumes that genotypes across different individuals follow an identical binomial distribution, a valid assumption in homogeneous populations. However, this assumption may not hold in study cohorts comprising individuals from diverse ancestries. To address this limitation, SPAGxEmix<sub>CCT</sub> modifies the approach to allow for different allele frequencies while maintaining binomial distributions for genotypes. SPAGxEmix<sub>CCT</sub> estimates individual-level allele frequencies using SNP-derived principal components (PCs) and raw genotypes. Similar to SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub> involves two main steps:

Step 1: SPAGxEmix<sub>CCT</sub> fits a genotype-independent (covariates-only) model and calculates the model residuals. Detailed information is provided in Step 1 of SPAGxE<sub>CCT</sub>.

Step 2: SPAGxEmix<sub>CCT</sub> identifies genetic variants with marginal G×E effects on the trait of interest. It first estimates the individual-level allele frequencies of the tested variants using SNP-derived PCs and raw genotypes. Next, SPAGxEmix<sub>CCT</sub> evaluates marginal genetic effects using score statistics. If the marginal genetic effect is not significant, S<sub>G×E(mix)</sub> is used as the test statistic to characterize the marginal G×E effect. If significant, S<sub>G×E(mix)</sub> is updated to genotype-adjusted test statistics. The hybrid strategy to balance computational efficiency and accuracy follows the approach used in SPAGxE<sub>CCT</sub>.

## Introduction of SPAGxEmix<sub>CCT-local</sub>
SPAGxEmix<sub>CCT-local</sub> is a G×E analysis framework designed to efficiently and accurately identify ancestry-specific G×E effects by incorporating local ancestry in multi-way admixed populations.

SPAGxEmix<sub>CCT-local</sub> extends SPAGxEmix<sub>CCT</sub> by integrating local ancestry information to enhance the precision of ancestry-specific G×E effect detection. Furthermore, we introduce SPAGxEmix<sub>CCT-local-global</sub>, which combines p-values from both SPAGxEmix<sub>CCT</sub> and SPAGxEmix<sub>CCT-local</sub>, offering an optimal unified approach for various cross-ancestry genetic architectures. As with SPAGxEmix<sub>CCT</sub>, SPAGxEmix<sub>CCT-local</sub> involves two main steps:

- In Step 1, SPAGxEmix<sub>CCT-local</sub> fits a genotype-independent (covariates-only) model and calculates the model residuals. Details are described in Step 1 of SPAGxE<sub>CCT</sub>.
  
- In Step 2, SPAGxEmix<sub>CCT-local</sub> identifies genetic variants with marginal G×E effectS on the trait of interest. First, SPAGxEmix<sub>CCT</sub> estimates the ancestry-specific allele frequencies for the variants being tested. Next, SPAGxEmix<sub>CCT-local</sub> evaluates ancestry-specific marginal genetic effects using ancestry-specific score statistics. If an ancestry-specific marginal genetic effect is not significant, we use ancestry-specific S<sub>G×E(mix)</sub> as the test statistic to characterize the ancestry-specific marginal G×E effect. If significant, S<sub>G×E(mix)</sub> is updated to ancestry-specific genotype-adjusted test statistics. The hybrid strategy to balance computational efficiency and accuracy follows the approach used in SPAGxE<sub>CCT</sub>.

**Compared to traditional statistical testing methods that account for local ancestry, SPAGxEmix<sub>CCT-local</sub> offers much greater computational efficiency.** SPAGxEmix<sub>CCT-local</sub> leverages local ancestry information to estimate the distribution of ancestry-specific genotypes and the null distribution of test statistics. For most tests (ancestry-specific genetic main effect p-values > 0.001) in a genome-wide G×E analysis, SPAGxEmix<sub>CCT-local</sub> utilizes residuals from a genotype-independent model (fitted only once across the genome-wide analysis) to construct test statistics. In contrast, traditional methods must incorporate local ancestry as covariates and fit a separate model for each variant, which is computationally burdensome for genome-wide analyses.

## The computational efficiency can be greatly enhanced through incorporating polygenic scores (PGSs) as covariates with fixed effects.

For SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub>, and SPAGxEmix<sub>CCT-local</sub>, incorporating polygenic scores (PGSs) as covariates with fixed effects can significantly enhance computational efficiency. When PGSs are available, a genotype-independent model can be fitted across genome-wide analyses, allowing us to use regular score statistics as test statistics, followed by a hybrid test employing normal approximation and SPA. By including PGSs as covariates, we eliminate the need for constructing statistics through matrix projection or linear regression using genotype data, thereby improving computational efficiency. Additionally, incorporating PGSs can further boost the statistical power of SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub>, and SPAGxEmix<sub>CCT-local</sub> by accounting for polygenic effects, as recent studies have demonstrated that adjusting for PGSs enhances statistical power.

In our approach, we implement a two-stage strategy, denoted as SPAGxE<sub>CCT</sub>-PGS, SPAGxEmix<sub>CCT</sub>-PGS, and SPAGxEmix<sub>CCT-local</sub>-PGS. In stage 1, we perform an initial round of GWAS using tools such as SAIGE and calculate the PGS based on summary statistics, or alternatively, use PGSs obtained from other databases. In stage 2, PGSs are included as additional covariates in a genome-wide GxE analysis using SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub>, or SPAGxEmix<sub>CCT-local</sub>. When applying this strategy to test for ancestry-specific GxE effects with SPAGxEmix<sub>CCT-local</sub>, it is important to use methods that improve PGS construction by first disentangling such mosaics through local ancestry inference and then taking inferred local ancestry into PGS construction, such as the pPRS method.

This strategy improves both computational efficiency and statistical power.

## UK Biobank data analysis results

In the paper **A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits (to be updated)**, we have applied SPAGxE<sub>CCT</sub> to analyze time-to-event traits in UB Biobank. For the SPAGxE<sub>CCT</sub> analyses, 281,299 White British individuals were included. For the SPAGxEmix<sub>CCT</sub> analyses, 338,044 individuals from all ancestries were included. As a universal analysis framework, we also evaluated the performance of SPAGxE<sub>CCT</sub> in ordinal categorical, binary, and quantitative trait analysis.  

**Summary statistics of time-to-event and binary traits in UK Biobank is available [here](https://zenodo.org/records/11571404).**

## Codes for simulation studies

The primary codes used in the SPAGxECCT project's simulation studies are available in the **simulation studies** folder and have also been uploaded to **Code Ocean**.

These codes cover GxE analyses for time-to-event, binary, quantitative, and ordinal categorical traits. For analyses of homogeneous populations, please assess the **simulation studies/SPAGxECCT** folder. For analyses of heterogeneous or admixed populations, please refer to the **simulation studies/SPAGxEmixCCT** folder. If you are using SPAGxEmix<sub>CCT-local</sub> to identify ancestry-specific G×E effects by incorporating local ancestry in analyses of admixed populations, please refer to the **simulation studies/SPAGxEmixCCT-local** folder. 

With the provided codes, users can reproducibly conduct simulation studies of GxE analyses, including generating genotype data, phenotype data, calculating genetic principal components (PCs), and performing p-value calculations using our proposed methods of SPAGxE<sub>CCT</sub> (for testing GxE effects in homogeneous populations), SPAGxEmix<sub>CCT</sub> (for testing GxE effects in heterogeneous or admixed populations), and SPAGxEmix<sub>CCT-local</sub> (for testing ancestry-specific GxE effects by incorporating local ancestry into analyses of admixed populations).


## Reference
See **A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits** (to be updated) for more details about SPAGxE<sub>CCT</sub>, SPAGxEmix<sub>CCT</sub>, and SPAGxEmix<sub>CCT-local</sub>.







