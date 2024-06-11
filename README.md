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
Current version is 0.0.1. For older version and version update information, plesase refer to OldVersions.  

Please do not hesitate to contact me (yuzhuoma@stu.pku.edu.cn) if you meet any problem. Suggestions or comments are also welcome.

## SPAGxECCT introduction
SPAGxECCT is a a G×E analysis framework which is applicable to a wide variety of complex traits with intricate structures (e.g. time-to-event and ordinal categorical, binary, quantitative, and longitudinal traits). SPAGxECCT contains two main steps (Figure 1). 

- In step 1, SPAGxECCT fits a covariates-only model and then calculates model residuals. The covariates include, but are not limited to, confounding factors such as age, gender, SNP-derived principal components (PCs), and environmental factors. The model specification and the corresponding model residuals vary depending on the type of trait. In Material and methods section and Appendix, we demonstrated regression models to fit time-to-event traits, binary traits, and ordinal categorical traits, along with the corresponding model residuals. As the covariates-only model is genotype-independent, the model fitting and residuals calculation are only required once across a genome-wide analysis.
  
- In step 2, SPAGxECCT identifies genetic variants with marginal G×E effect on the trait of interest. First, SPAGxECCT tests for marginal genetic effect via score statistic. If the marginal genetic effect is not significant, we use S_(G×E) as the test statistics to characterize marginal G×E effect, where E_i  i≤n denote the environmental factors. Otherwise, statistics S_(G×E) is updated to S ̃_(G×E).

To balance the computational efficiency and accuracy, SPAGxECCT employs a hybrid strategy to combine normal distribution approximation and SPA to calculate p values, as in previous studies. For variants with significant marginal genetic effect, SPAGxECCT additionally calculates p value through Wald test and then uses Cauchy combination (CCT) to combine p values from Wald test and the proposed test statistics.




In the paper **A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits (to be updated)**, we applied SPAGxECCT to analyze time-to-event traits in UB Biobank. As a universal analysis framework, we also evaluated the performance of SPAGxECCT in ordinal categorical, binary, and quantitative trait analysis.  

**Summary statistics of time-to-event traits in UK Biobank is available here.**



## Reference
See **A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits** (to be updated) for more details about SPAGxECCT.







