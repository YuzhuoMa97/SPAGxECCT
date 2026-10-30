# Statistical policy and comparison scope

SPAGxE 1.2.0 provides one implementation of SPAGxE, SPAGxE+, SPAGxEmix and
SPAGxEmix_local, following the [published paper](https://www.nature.com/articles/s41467-025-57887-3)
with documented input and numerical corrections and guarded partial-normal
CGF acceleration for all four methods.

| Component | Release policy | Reason/limitation |
|---|---|---|
| Branch B | Residualize R on intercept and G, then E*Rtilde | Preserves the published definition; further weight projection changes the test |
| CGF | Guarded partial-normal SPA; exact binomial tails and moment-matched central Gaussian block | Adds an explicit approximation; fixed 1.5-IQR partition and 1e-3 full-CGF diagnostic, with automatic full evaluation |
| Near-mean threshold | Normal approximation at default absolute Z threshold 2 | Paper-defined normal/tail-SPA strategy |
| Individual AF | Minor-oriented pooled/OLS/carrier-logistic, mapped to ALT | Two-sided bounds and strict 90% criterion follow the supplied supplementary method |
| Wald | Native full-interaction model for eligible unrelated raw traits in B | No independent-sample Wald for related or residual-only input |
| Cox ties | Breslow for native null and full models | Explicit single convention; historical Efron comparisons require matching refits |
| Root solving | Bracketing, log tails, support checks | Numerical safeguards, no second method mode |
| CCT | Stable log-domain combination, explicit weights | No blanket clipping of near-one P; exact P=1 has its documented convention |
| Local G/H | Bin(H,q), H-G, H*q, shared complete-H mask | Respects known ancestry copies rather than diploidizing each ancestry |

Individual AF uses pooled frequency when completed-cohort minor MAC is at most
20. Otherwise OLS is accepted only when more than 90% of fitted AF is in [0,1],
then clipped to that interval. Otherwise eligible PC slopes select a carrier
logistic fit. The documented low-count/constant-outcome routing uses pooled AF.
A logistic fit that exhausts its iteration budget, fails its linear solve, or
has non-finite coefficients is rejected rather than replaced by pooled AF:
association P fields are NA with status 9 (AF_FIT_FAILED), and AF exports record
status -2. This makes an unavailable reference model explicit. There is no separate
greater-than-20 carrier rule. Observed QC MAC serves a different purpose.
See [AF exports](AF_MODEL_FORMAT.md).

Default global/local CCT weights each group by one half.
--cct-per-ancestry-equal instead equally weights available components.
Local-only cross-ancestry CCT is a different output.

## Environmental coding

E is used as supplied. In the paper's Branch-B working approximation, E+c can
change its plug-in variance/SPA leg even where the observed homogeneous score
and prospective interaction Wald are unchanged. Choose and report environmental
origin/units before analysis, not after inspecting P.

A further covariance-metric projection along Rtilde can reduce a fixed working
variance and remove this origin dependence. These are useful mathematical
properties, but alone do not establish better calibration of the complete
data-dependent SPA/CCT: Rtilde itself depends on observed G. This release keeps
the published weight definition. The declared partial-CGF acceleration is
documented separately in [partial SPA](PARTIAL_SPA.md); it is not an additional
weight projection. Changing the working weights would require a separate
method revision with calibration, power and sensitivity studies.

## Agreement and scope

Historical comparisons must match masks, alleles, environmental coding, PCs,
AF routing, residual/null models, ties, thresholds and CCT weights. A difference
is not automatically a solver error; smaller P does not prove greater accuracy.

The full-CGF baseline closely reproduced original R in stated common domains.
The current guarded partial approximation may differ even with matched inputs;
compare it against an independent implementation of the declared partial policy
and report full-reference differences separately. The 1e-3 diagnostic is not a
P-value error theorem. Historical defects
such as 2-G where H-G is required, or 2*q where H*q is required, are not copied.
Rank-deficient and missing-H cases require separate checks. Universal equality
to every historical default is not claimed.

SPAGxE+ is the published homogeneous related-sample score method. Mix and local
require unrelated subjects. The unpublished individual-AF plus cross-individual
GRM GxE extension is excluded. Local input does not infer ancestry or integrate
arbitrary uncertain-H/joint ancestry-genotype posteriors.

Example agreement and source-archive checks are not universal GWAS-tail
calibration or a measured large-biobank speed comparison.
See [release scope](RELEASE_STATUS.md), [ownership](../OWNERSHIP.md) and
[third-party notices](../THIRD_PARTY_NOTICES.md).
