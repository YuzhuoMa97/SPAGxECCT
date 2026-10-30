# Implementation choices for the official C++ release

The release preserves the published SPAGxE weight definitions and uses one
statistical policy with declared guarded partial-normal CGF acceleration. A different implementation is not a reason, by
itself, to replace an equation or declare either complete test superior.
The incorporated-source provenance remains in ../THIRD_PARTY_NOTICES.md.

| Item | Official release decision | Evidence and tradeoff |
|---|---|---|
| Branch B | Use E*R0 after residualizing R on intercept and G | Reproduces the published definition |
| Further weight projection | Do not add silently | Can reduce a fixed working variance and remove E-origin dependence, but changes the working test |
| CGF | Guarded partial-normal CGF beyond the normal threshold, in all four methods | Reduces repeated central CGF work; adds an approximation with full-CGF candidate-root diagnostics and automatic full evaluation |
| AF policy | Minor-oriented pooled/OLS/carrier-logistic with two-sided bounds and strict 90% criterion | Matches the stated method; carrier count and minor MAC are not interchangeable |
| Failed AF fits | NA/status 9; model export -2 | Finite but unconverged coefficients are not a valid fitted reference law |
| Wald | Eligible unrelated raw-trait Branch B only | Residual-only and GRM paths do not acquire a fictitious independent-sample Wald leg |
| Cox ties | Breslow for the native fits | Historical R default Efron comparisons must explicitly match ties |
| Root and tails | Bracket/support checks and log-domain arithmetic | Numerical safeguards preserve valid tails and expose failed computations |
| CCT | Explicit weights, stable log-domain combination | Legal P below 1 is not blanket-clipped; exact P=1 uses the documented endpoint policy |
| Local G/H | Known H, Bin(H,q), H-G, H*q, rank-aware Wald | Respect ancestry copy counts and stated input/missingness contracts |

The extra covariance-metric projection has useful mathematical properties.
They do not alone prove that a complete data-dependent SPA/CCT has better
type-I error or power in all settings: R0 itself depends on observed G.
The paper's Branch-B working approximation can depend on environmental origin.
Choose and report E coding before inspecting association P values.

The partial-CGF acceleration is an explicit implementation revision, described
in [partial SPA](PARTIAL_SPA.md). It preserves the first two working cumulants
but not every higher cumulant or P value. Fixed 1.5-IQR partitioning and 1e-3
full-CGF diagnostics define one policy; automatic full evaluation is its
protection path, not a second user-selectable method. Diagnostic tolerance is
not a theorem about P-value error. Performance claims require measured workloads.
A further weight revision would require separate calibration, power and coding
sensitivity studies; it is not silently included with this acceleration.
Current code is not known to be universally worse or better than an alternative
complete implementation.

The four methods have distinct scopes: homogeneous unrelated SPAGxE,
homogeneous GRM-adjusted SPAGxE+, individual-AF unrelated SPAGxEmix, and
known-H local-ancestry unrelated SPAGxEmix_local.
The combined individual-AF plus cross-subject-GRM interaction extension is not
provided. Local input does not infer ancestry or integrate arbitrary uncertain H.

R agreement requires matched samples, masks, alleles, E coding, PCs/AF routing,
residuals/null fits, ties, thresholds and CCT weights. Boundary repairs are not
required to reproduce historical defects. Smaller P does not establish accuracy.
See [method policy](METHOD_POLICY.md), [equations](METHOD_EQUATIONS.md),
[release scope](RELEASE_STATUS.md) and [validation](validation/README.md).


Coverage of the three global interaction workflows does not mean that every
alternative approximation, selector or historical R extension is preserved.
The single documented statistical policy is deliberate; see
[feature coverage](FEATURE_COVERAGE.md). Historical reference comparisons and
current executable regressions are different evidence layers. Matched-domain
agreement is not an assertion of all-default, all-input equality.
