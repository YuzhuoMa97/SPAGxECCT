# SPAGxE method equations

The native executable fits a genotype-independent covariate model once and
then tests variants. The environmental main effect and relevant ancestry PCs
must enter that model. Let its residual vector be $R$, the environment be
$E$, and the tested genotype dosage be $G$.

The executable preserves the paper's Branch-B weight and uses one guarded
partial-normal cumulant generating function (CGF) policy for all four methods.
The normal approximation is used when $|Z|\le2$, and guarded partial SPA
beyond that cutoff. Automatic full-CGF evaluation protects unsuitable partial
candidates; it is not a second user-selected profile. The normal/SPA cutoff
and the central-Gaussian CGF acceleration are distinct approximations.

## SPAGxE for a homogeneous unrelated population

With allele frequency $q$, the retrospective reference law is
$G_i\sim\operatorname{Bin}(2,q)$. For a fixed weight vector $w$:

$$
S_w=\sum_iG_iw_i,\quad
\mu_w=2q\sum_iw_i,\quad
V_w=2q(1-q)\sum_iw_i^2.
$$

The marginal genetic screen uses $w=R$. Its normal P value selects one of
two branches using the default threshold $\epsilon=0.001$.

For Branch A, when the marginal P value exceeds $\epsilon$:

$$
\lambda=\frac{\sum_iE_iR_i^2}{\sum_iR_i^2},\qquad
w_i=(E_i-\lambda)R_i.
$$

For Branch B, when the marginal P value is at most $\epsilon$:

$$
W=(\mathbf1,G),\qquad
\widetilde R=R-W(W^TW)^{-1}W^TR,\qquad
w_i=E_i\widetilde R_i.
$$

The implementation uses a numerically stable projection rather than forming
a dense $n\times n$ projection matrix. For raw phenotype analyses, Branch B
also fits the full prospective model with $G$ and $G\times E$. Its interaction
Wald P value is combined with the SPA P value by CCT. Residual-only input
does not provide the phenotype/model information needed for that Wald refit.

The original homogeneous marginal screen assumes that fitted residuals sum
to zero. The native implementation explicitly centers scores, which also
supports externally supplied residuals with a nonzero mean. Such residuals
are outside the original zero-mean assumption, so this correction can change
branch decisions relative to the old R code.

## Working CGF and guarded partial SPA

The full working CGF, used as the reference and protection calculation, is

$$
K(t)=\sum_i2\log\{1-q_i+q_i\exp(tw_i)\}.
$$

For SPAGxE, every $q_i=q$. For SPAGxEmix they may differ. The default partial
calculation replaces central terms by their Gaussian contribution:

$$
K_P(t)=\sum_{i\in O}h_i\log\{1-q_i+q_i\exp(tw_i)\}
       +\mu_Ct+\tfrac12V_Ct^2,\qquad
\mu_C=\sum_{i\in C}h_iq_iw_i,\quad
V_C=\sum_{i\in C}h_iq_i(1-q_i)w_i^2.
$$

Here $h_i=2$ for global methods and $h_i=H_i$ locally. The display above
assumes random terms; deterministic $q_i=1$ terms add their linear CGF
contribution separately. Type-7 quartiles and fixed 1.5-IQR Tukey fences
partition actual score weights with $h_i>0$, $w_i\ne0$ and $0<q_i<1$ into
exact outlier terms $O$ and Gaussian central terms $C$. The multiplier is
not iteratively shrunk. The first two cumulants are preserved; higher
cumulants generally change. No outliers or no central random variance
selects full evaluation. Zero IQR can still define valid outliers.
A full-CGF diagnostic checks both candidate
roots at tolerance $10^{-3}$; unsuitable candidates or original finite-support
boundary targets trigger full evaluation of both tails. This threshold is
not a mathematical bound on P-value error. See [partial SPA](PARTIAL_SPA.md)
for the complete policy.

A bracketed root solver finds $t$ satisfying the selected $K'(t)=s$.
The signed saddlepoint and curvature
terms are

$$
r=\operatorname{sign}(t)\sqrt{2\{ts-K(t)\}},\qquad
v=t\sqrt{K''(t)},\qquad
r^*=r+\frac{\log(v/r)}r.
$$

The two tails are evaluated around the fitted score mean and added on the
log scale. Near singularities or failed roots produce an explicit status
and the documented fallback; a numerical failure is not changed into P=0.

When the SPA branch is entered, a pure binomial CGF also permits analytic
support checks. Write $h_i=2$ globally or $h_i=H_i$ locally, let
$I=\{i:0<q_i<1\}$ and $C=\sum_{i:q_i=1}h_iw_i$. Its support endpoints are
$L=C+\sum_{i\in I}h_i\min(w_i,0)$ and
$U=C+\sum_{i\in I}h_i\max(w_i,0)$. Deterministic $q_i=0/1$ terms shift the
support but do not constrain an endpoint event. A strict support-exterior
tail is exactly zero under this working law; the valid opposite tail is
retained. At $U$, the inclusive upper-tail atom is the product of $q_i^{h_i}$
for positive weights and $(1-q_i)^{h_i}$ for negative weights; $L$ reverses
those choices. Zero weights impose no constraint. These endpoint products
are accumulated on the log scale, with a scale-dependent floating-point
rounding allowance when comparing scores to bounds.

The same guarded policy is used by the homogeneous, individual-AF and local
kernels. A partial CGF with $V_C>0$ has unbounded support, so finite binomial
endpoint masses are used only after full evaluation has been selected. The
normal branch selected near the mean remains a separate approximation.

Full CGF still means a saddlepoint approximation, not exact enumeration of
the score distribution.

Dosages enter the observed score. The working CGF uses the population
binomial reference law; it does not replace that law by each sample's
posterior genotype-probability distribution.

## SPAGxEmix for population structure

SPAGxEmix estimates an individual allele frequency $q_i$ using the specified
PCs. The model uses an intercept and PCs, with pooled AF for very low MAC,
an OLS AF fit where appropriate, and a carrier-logistic fallback. AF fitting
is performed in the minor direction and mapped back to the reported allele
direction. The PCs must also adjust the phenotype null model.

The AF cascade first uses pooled AF when MAC is at most 20. Otherwise it
fits genotype on an intercept and the PCs by OLS. In this implementation,
OLS is retained only when **more than 90%** of its fitted AF values lie in
[0,1]; retained values are clipped to that interval. If at least 10% lie
outside either bound, the estimator selects PCs with OLS P<.05 and fits a
carrier-logistic model, converting carrier probability $p_i$ to
$q_i=1-\sqrt{1-p_i}$. No selected PC or a constant carrier response uses
pooled AF. A selected fit that is invalid or fails to converge is rejected:
association P values are NA with status 9 (AF_FIT_FAILED), and AF exports
record status -2. A nonconstant response with 20 or fewer carriers is still
eligible for logistic fitting once MAC exceeds 20.

The two-sided bounds and strict 90% rule follow the supplied supplementary
method. The original R code instead counts only negative OLS values and
switches when that proportion is **greater than** .1. It can therefore
select a different AF model on otherwise valid minor-allele data. These
are algorithmic differences, not saddlepoint solver tolerances. See [method policy](METHOD_POLICY.md) for comparison requirements.

Set $v_i=2q_i(1-q_i)$. Then

$$
\mu_w=2\sum_iq_iw_i,\qquad
V_w=\sum_iv_iw_i^2,\qquad
\lambda=\frac{\sum_iv_iE_iR_i^2}{\sum_iv_iR_i^2}.
$$

The same paper-defined Branch-B projection and guarded partial CGF are used.

## SPAGxE+ for relatedness

GRM values use the relationship scale: twice the kinship coefficient.
With $\Phi$ denoting that matrix, homogeneous SPAGxE+ uses

$$
V_w=2q(1-q)w^T\Phi w,\qquad
\lambda=\frac{R^T\Phi(E\circ R)}{R^T\Phi R}.
$$

Its final P is SPA-only; no independent-sample Wald test is attached to
related samples. Let $V_0=2q(1-q)\sum_iw_i^2$. The score adjustment sends
$\mu_w\pm|S_w-\mu_w|\sqrt{V_0/V_w}$ to the guarded partial CGF of the
independent-binomial reference law. Its moment-matched central block preserves
$V_0$; it does not add an individual-AF plus GRM method.
SPAGxEmix and local methods in this release require unrelated samples.

## Local ancestry

For ancestry $a$, $H_{ia}\in\{0,1,2\}$ is its chromosome-copy count and
$0\le G_{ia}\le H_{ia}$ is its ALT dosage. Known copy counts sum to two
across ancestries. The reference law is $G_{ia}\sim\operatorname{Bin}(H_{ia},q_a)$.

Estimate ALT AF from observed G and corresponding known H. Missing G is
imputed by $Hq$; H=0 is a structural zero. Unknown H is not imputed. All
ancestry tests at a variant use the same complete-H subjects after applying
the phenotype mask.

If ALT AF exceeds .5, the tested minor dosage is **H-G**, not 2-G. Output
columns `TESTED_ALLELE` and `MINOR_FLIPPED` record this direction.

For a tested ancestry:

$$
\mu_w=q\sum_iH_iw_i,\qquad
V_w=q(1-q)\sum_iH_iw_i^2,\qquad
K(t)=\sum_iH_i\log\{1-q+q\exp(tw_i)\}.
$$

Branch A uses $\lambda=\sum_iH_iE_iR_i^2/\sum_iH_iR_i^2$; Branch B uses
$E\circ\widetilde R$. The full-model Wald refit includes only H columns
that increase the covariate-design rank, because the intercept already
captures $\sum_aH_{ia}=2$. It then reduces any redundant nuisance columns,
including a redundant G main effect, while preserving G:E as the final
tested target. The target must increase the nuisance rank by one. Thus
G proportional to H can still have an estimable G:E effect; a redundant
main-effect design is not by itself a reason to discard that interaction.
An all-zero-H ancestry has no available trials and produces no test.

The local full-CGF protection calculation has finite support
$[\sum_iH_i\min(w_i,0),\sum_iH_i\max(w_i,0)]$.
A reflected target strictly outside that interval has tail probability
zero, while the valid opposite tail remains a SPA calculation. At a support
endpoint, the inclusive tail includes its atom; it is not replaced by zero.
The partial CGF with positive Gaussian-block variance has unbounded support;
original support-boundary targets select full evaluation for both tails.
Local central mean and variance use H, including structural zero contribution
for H=0, rather than pretending that every subject has two ancestry copies.

Local ancestry CCT combines valid ancestry-specific final P values and
reports the number included. It is distinct from global/local CCT, which
also includes a global association result. The dedicated combination
command makes the grouping weights explicit. Ordinary diploid HWE is a
diagnostic of summed global hardcalls, not a Bin(H,q) local test.

## Compatibility and interpretation

This executable provides one documented statistical policy. Partial-CGF
acceleration, formula changes and input-defect
corrections can change historical outputs; universal equality with every historical
default is not an acceptance claim. See [method policy](METHOD_POLICY.md).

`BETA_Gx*` and `SE_Gx*` are working-score summaries based on the reported
score variance. They are not the prospective interaction coefficient and
standard error from the separate Wald fit. A final CCT P value also need
not correspond to the normal score Z; the output separates these quantities.

The environment is used in its supplied coding. In the published Branch-B
working approximation, changing E to E+c leaves the observed homogeneous
score unchanged but can change its plug-in variance and final P value.
Choose and report environmental coding before analysis; do not select
centering after inspecting P values. The [method policy](METHOD_POLICY.md) explains this limitation and why changes
to the weight require separate method validation.
