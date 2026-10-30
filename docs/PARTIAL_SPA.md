# Guarded partial-normal saddlepoint calculation

All four C++ methods use the `guarded_partial_normal` policy beyond the
normal/SPA threshold. It accelerates repeated cumulant generating function
(CGF) evaluations by replacing central contributions with a Gaussian block
whose mean and variance match the binomial working model. Tail contributions
retain their binomial form. This adds an approximation; it is not a promise
that every P value is identical to a full-binomial SPA calculation.

This is one official policy. Full-CGF evaluation is an automatic protection
path within it, not a second user-selectable method. The paper-defined Branch B
weight remains E times the genotype-adjusted residual. No additional weight
projection or individual-AF plus GRM extension is introduced.

## Working moments and partition

For the actual interaction score weights w, let h=2 for the three global
methods and h=H for each local ancestry. The working genotype law has allele
frequency q, which may vary between individuals in SPAGxEmix. Partition the
random contributions with h>0, w!=0 and 0<q<1 using type-7 quartiles of w.
Weights outside [Q1-1.5 IQR, Q3+1.5 IQR] form the exact binomial tail set O;
remaining eligible weights form the Gaussian central set C. The multiplier
is fixed at 1.5 and is not iteratively shrunk to force outliers.

Zero-copy and zero-weight terms do not affect the partition. A deterministic
q=1 contribution shifts the score mean; q=0 contributes zero. In particular,
H=0 subjects cannot change the ancestry-specific quartiles. If IQR is zero
but valid outliers exist, the partial candidate is still checked normally.
No outliers or no central random variance selects full evaluation.

The central moments and partial CGF are

$$
\mu_C=\sum_{i\in C}h_iq_iw_i,\qquad
V_C=\sum_{i\in C}h_iq_i(1-q_i)w_i^2,
$$

$$
K_P(t)=d t+\sum_{i\in O}h_i\log(1-q_i+q_i e^{t w_i})
       +\mu_Ct+\tfrac12 V_Ct^2,
$$

where d is the sum of deterministic h*w contributions with q=1. Its derivatives
use the exact outlier binomial derivatives plus mu_C+t*V_C and V_C. Thus the
full and partial laws have the same first two cumulants at t=0, including
individual q and local H. Their higher cumulants generally differ.

SPAGxE+ retains its existing relatedness score-variance adjustment before the
SPA calculation. Moment matching preserves the independent variance used by
that adjustment; it does not create a joint individual-AF/GRM model.

## Full-CGF diagnostics and protection

For each partial candidate saddlepoint, evaluate the full working CGF once.
The checks examine the absolute CGF difference, standardized full-root residual
and relative curvature difference. All three must satisfy the fixed 1e-3
diagnostic threshold on both tails. A failed diagnostic or unsuitable numerical
candidate selects full-CGF evaluation for both tails. The decision is made
together so that the two-sided result does not mix independently accepted
partial and full tails.

Precisely, at a candidate root t for target s, the diagnostics are
`abs(K_F(t)-K_P(t))`, `abs(K_F'(t)-s)/sqrt(K_F''(t))`, and
`abs(K_F''(t)-K_P''(t))/K_F''(t)`, where F denotes the full reference.
Finite positive full curvature is required. Small machine-roundoff allowances
are added to the three thresholds: with `e=64*epsilon`, these are
`e*(1+abs(K_F)+abs(K_P)+abs(t*s))`,
`e*(1+abs(s)+abs(K_F'))/sqrt(K_F'')`, and `e`, respectively.

Original finite-support boundary or exterior targets also select the full
calculation. A partial CGF with V_C>0 has unbounded support; finite-binomial
endpoint masses cannot be inserted into that Gaussian mixture. Full support
handling then preserves valid opposite tails and the documented endpoint atoms.

The 1e-3 threshold is a numerical diagnostic criterion, **not a bound on
absolute or relative P-value error** and not a proof of tail calibration. It
limits the checked discrepancy at the candidate roots. The approximation and
its protection path require independent numerical checks and empirical
calibration appropriate to the scientific application.

Selecting full evaluation is not the same as a numerical failure followed by
a normal approximation. Existing SPA_STATUS failure codes retain that separate
meaning. Settings receipts record the statistical policy; they should be kept
with analysis results and executable checksums.

## Interpretation and validation

The normal branch near the mean, marginal genetic screen, native null model,
Wald eligibility and CCT weights are unchanged. Individual-AF model selection
also remains minor-oriented pooled/OLS/carrier-logistic with explicit failure
reporting. Stable conversion from logistic linear predictors avoids premature
saturation but does not change the intended AF model.

An independent partial-policy oracle must reproduce the specified partition,
moments, candidate-root diagnostics and automatic protection, not merely allow
a wider tolerance against an old full-CGF reference. Full-reference comparisons
measure the effect of the new approximation and are reported separately from
same-definition correctness checks. Unchanged AF/null/score quantities and
format/thread equivalence retain their strict checks.

Replacing many central terms in every root iteration can reduce CGF work.
The full-CGF diagnostic adds evaluations, and protected cases also pay for a
full solve. Wall-clock gains depend on sample size, tail count, effect strength,
trait, AF, local H, and fallback frequency. A universal speed advantage or
biobank-scale benchmark is not implied by the implementation.

See [equations](METHOD_EQUATIONS.md), [statistical policy](METHOD_POLICY.md),
[implementation choices](IMPLEMENTATION_CHOICE.md) and
[validation evidence](validation/README.md).
