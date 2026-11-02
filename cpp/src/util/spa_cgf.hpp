// Binomial CGF reductions for shared AF, individual AF and known trial counts.
//
// A weighted Binomial(h,p) contribution has alpha=(1-p)+p*exp(t*r):
// K=h*log(alpha), K'=h*r*pi and K''=h*r*r*pi*(1-pi).
// The implementation obtains pi and its complement as e/alpha and
// (1-p)/alpha, avoiding subtraction of nearly equal values.
//
// K12 evaluates derivatives for Newton iterations. KFull reuses that exact
// derivative path and evaluates terminal K separately. KFullExact retains
// scalar log1p/expm1 for diagnostics that need relative accuracy in K itself.
//
// Derivative exponent arguments saturate at +/-708. Terminal K restores the
// correct directional linear excess; the exact-relative diagnostic instead
// uses a factored expression above 700. These are distinct numerical policies.
// Scalar, AVX2 and AVX-512 paths keep their existing summation/lane order.
// Frequencies must be in [0,1], trial counts nonnegative, and input products
// finite. Empty reductions return zero. Runtime dispatch is architecture-aware.

#pragma once

#include <cmath>

namespace spa_cgf {

// ══════════════════════════════════════════════════════════════════════
// Return types
// ══════════════════════════════════════════════════════════════════════

// Newton-loop payload.  For the reductions, K1 and K2 are sums over subjects.
struct Cgf12 {
    double K1;   // sum_i K'_i(t)
    double K2;   // sum_i K''_i(t)
};

// Terminal payload.  For the reductions, all three are sums over subjects.
struct Cgf012 {
    double K0;   // sum_i K_i(t)
    double K1;   // sum_i K'_i(t)
    double K2;   // sum_i K''_i(t)
};

// ══════════════════════════════════════════════════════════════════════
// Domain clamp
// ══════════════════════════════════════════════════════════════════════
//
// t·r is clamped to +/- kTrClamp before `exp`, in every tier.
//
// Two reasons, both about agreement rather than about accuracy:
//
//   1. Unclamped, `std::exp(800)` returns +inf, giving alpha = +inf,
//      inva = 0 and K' = h·r·(+inf)·0 = NaN.  The vectorized `avx2_exp_pd` /
//      `avx512_exp_pd` clamp their argument internally to +/-709, so they would
//      return a finite result where the scalar tier returns NaN.  Clamping
//      explicitly, at a bound strictly inside the vectorized kernels' own, makes
//      all three tiers agree by construction and removes the NaN.
//
//   2. Inside the bound the clamp is inert, and at the bound the kernel is still
//      exact: at t·r = 708 with p = 0.3 it returns K'' = 1.5435e-307, which
//      agrees with a long-double evaluation to every printed digit.
//
//      PAST the bound the clamp does change the value, and deliberately so.  It
//      returns the derivatives evaluated AT the bound rather than their true,
//      smaller magnitudes: K'' saturates at ~1.5e-307·r^2 whereas the true K''
//      keeps decaying and underflows to zero near t·r = 745.  The alternative at
//      those arguments is not a more accurate number but NaN — unclamped,
//      exp(800) = +inf makes alpha = +inf, inva = 0 and pi = inf·0 = NaN — and a
//      Newton loop that divides by K'' is far better served by a finite positive
//      saturation than by NaN.  In absolute terms the substitution is below
//      1e-307 in a quantity the solver compares against O(1) variances.
//
//      The saturation also keeps K'' STRICTLY positive in the far tails, which
//      the unclamped algebra does not: at t·r = -800 an unclamped exp underflows
//      to exactly zero, giving e = 0 and K'' = 0, and a zero curvature is a
//      division by zero in the Newton step.
//
// exp(708) = 3.0e307 and exp(-708) = 3.3e-308 are both finite and normal, so no
// intermediate of the clamped kernel is subnormal.
//
// K is NOT clamped: `logAlpha` grows linearly in t·r and is overflow-free by
// construction, so the terminal K stays correct past the bound even though the
// derivatives have saturated. That asymmetry is intentional and correct.
//
// A non-finite t·r is outside the contract of these kernels.
inline constexpr double kTrClamp = 708.0;

inline double clampTr(double x) noexcept {
    return (x > kTrClamp) ? kTrClamp : ((x < -kTrClamp) ? -kTrClamp : x);
}

// ══════════════════════════════════════════════════════════════════════
// log(alpha), accurate at t = 0 and overflow-free at both extremes
// ══════════════════════════════════════════════════════════════════════
//
// Returns log(alpha) = log((1-p) + p·e^x).  Multiply by h to obtain K.
//
// ── N1: why expm1/log1p, and where they must NOT be used ─────────────
//
// Near x = 0, alpha = 1 + delta with delta = p·(lambda - 1) far below 1.
// Rounding alpha to a double discards delta's low bits, so log(alpha) inherits a
// relative error of order eps/delta — measured at up to 3.05e-7, nine orders of
// magnitude worse than the K' and K'' the same routine returns.  Forming delta
// directly as p·expm1(x) and taking log1p(delta) never materializes 1 + delta
// and restores full precision at identical cost.
//
// The critical caveat: expm1 is used ONLY for delta, never to recover lambda.
// Reconstructing lambda as 1 + expm1(x) cancels catastrophically for x very
// negative — expm1 returns a value near -1 with absolute accuracy ~eps, so the
// sum carries absolute accuracy ~eps while lambda itself is exponentially
// small.  K' and K'' consequently keep plain `exp`, which is harmless for them
// because alpha = u + e is a sum of two non-negative quantities and never
// cancels.
//
// ── The two branches ─────────────────────────────────────────────────
//
// The log1p form is used over the whole region the saddlepoint solver visits,
// |x| <= kLogAlphaReflect = 700, on both sides of zero.  It is correct
// throughout: at x -> -inf, expm1 saturates at -1 and the result is
// log1p(-p) = log(u), which is the exact limit.
//
// Beyond x = 700 the form cannot be used because expm1(x) overflows at
// x ~ 709.8.  There the exponential is factored out instead:
//
//     log(u + p·e^x) = x + log(p + u·e^{-x})
//
// The parenthesis is a sum of two non-negative quantities, so it never cancels,
// and e^{-x} merely underflows towards zero, leaving x + log(p) — the exact
// limit.  This branch is deliberately NOT written as x + log1p(u·expm1(-x)):
// that spelling reintroduces 1 - u implicitly and loses p entirely whenever
// u rounds to exactly 1.0, i.e. for every p below 2^-53, returning -infinity
// where the true value is finite.
//
// Splitting at 700 rather than at 0 matters for accuracy, not only for
// overflow.  For small positive x the reflected spelling would evaluate
// x + log(...) as a difference of two nearly equal quantities and lose about
// four digits, defeating N1 on the t > 0 side; the log1p branch keeps full
// precision there.
//
// ── The two endpoints ────────────────────────────────────────────────
//
// p == 1 and p == 0 are handled first, and must be: with p == 1 the x <= 0
// branch degenerates to log1p(expm1(x)), which returns -inf once expm1(x) rounds
// to exactly -1, and with p == 0 the x > 0 branch degenerates to
// x + log1p(expm1(-x)), which returns -inf for the same reason.  Both are
// exactly representable instead:
//
//     p == 1  =>  alpha = e^x,  log(alpha) = x     (G == h almost surely)
//     p == 0  =>  alpha = 1,    log(alpha) = 0     (G == 0 almost surely)
//
// Preconditions: p in [0, 1].  Outside that range alpha may vanish by
// cancellation and the result is unspecified.
// Above this argument expm1 overflows and the factored branch takes over.
// expm1 overflows at x ~ 709.78; p·expm1(700) <= 1.02e304 is comfortably finite.
inline constexpr double kLogAlphaReflect = 700.0;

inline double logAlpha(double x, double p) noexcept {
    const double u = 1.0 - p;
    if (u == 0.0) return x;              // p == 1: alpha = e^x exactly
    if (p == 0.0) return 0.0;            // p == 0: alpha = 1 exactly
    if (x <= kLogAlphaReflect) return std::log1p(p * std::expm1(x));
    return x + std::log(p + u * std::exp(-x));
}

// ══════════════════════════════════════════════════════════════════════
// log(alpha), the vectorizable spelling — the production path
// ══════════════════════════════════════════════════════════════════════
//
// Returns log(alpha) = log((1-p) + p·e^x) formed the way the SIMD tiers form it,
// so the scalar tier differs from them only by the vendored vector exp and log
// and not by the algebra.  See the terminal-K section at the top of this header
// for the measurement that justifies preferring this to `logAlpha`, and
// `binom*KFullExact` for the entry point that keeps `logAlpha`.
//
// Three properties are relied on and are not obvious:
//
//   1. THE CLAMP IS UNDONE FOR K.  `alpha` is built from the clamped argument,
//      because that is what k12 does and what keeps every intermediate normal.
//      Past the clamp the derivatives saturate deliberately (see kTrClamp) but K
//      does not: log(u + p·e^x) -> x + log(p) as x -> +inf, i.e. K keeps growing
//      linearly.  Adding back `x - clamp(x)` recovers that exactly:
//      log(alpha_clamped) + (x - 708) = 708 + log(p + u·e^-708) + x - 708, which
//      differs from the true x + log(p + u·e^-x) by under 1e-307.
//
//   2. THE CORRECTION IS DIRECTIONAL, AND VANISHES AT p == 0.  For x < -708 with
//      u > 0 the true value is log(u + p·e^x) -> log(u), which
//      log(alpha_clamped) already gives to full precision since
//      p·e^-708 <= 3.3e-308 is negligible beside u >= 2^-53.  Adding the deficit
//      there would be wrong.  When u == 0 exactly (p == 1) the true value is x
//      itself, and then the deficit is exactly what recovers it.  And at p == 0
//      alpha is identically 1 and log(alpha) is identically 0, so NO correction
//      applies however large x is.  Hence: the full excess when u == 0, its
//      positive part when 0 < p < 1, and nothing at all when p == 0.
//
//   3. alpha IS NEVER ZERO AND NEVER SUBNORMAL.  alpha = u + p·e^{clamp(x)} with
//      u >= 0 and, for p > 0, p·e^{clamp(x)} >= p·3.3e-308.  u == 0 forces p == 1
//      (for any double p < 1, Sterbenz gives 1 - p >= 2^-53), and then
//      alpha = e^{clamp(x)} >= 3.3e-308, which is normal.  This matters because
//      `avx512_log_pd` decomposes its argument by masking the exponent field and
//      would return a wrong answer, not an infinity, on a subnormal.
namespace detail {

// Endpoint-aware correction is kept separate from exp/log evaluation so the
// terminal-CGF clamp policy has one explicit scalar definition.
inline double terminalClampCorrection(double complement,
                                      double frequency,
                                      double excess) noexcept {
    double correction = 0.0;
    if (complement == 0.0)                    correction = excess;
    else if (frequency > 0.0 && excess > 0.0) correction = excess;
    return correction;
}

}  // namespace detail

inline double logAlphaFast(double x, double p) noexcept {
    const double u  = 1.0 - p;
    const double xc = clampTr(x);
    const double a  = u + p * std::exp(xc);
    const double excess = x - xc;              // 0 strictly inside the clamp
    const double add = detail::terminalClampCorrection(u, p, excess);
    return std::log(a) + add;
}

// ══════════════════════════════════════════════════════════════════════
// Per-subject scalar cores
// ══════════════════════════════════════════════════════════════════════
//
// One core covers all three variants: binomUniform is h == 2 with p shared,
// binomIndiv is h == 2 with p per-subject, binomHapcount is p == q shared with
// h per-subject.
//
// ── The alpha == 0 guard ─────────────────────────────────────────────
//
// alpha = u + e with u >= 0 and e >= 0, so alpha == 0 requires both to vanish.
// u == 0 means p == 1 exactly (for any double p < 1, Sterbenz gives
// 1 - p >= 2^-53 > 0), and e == 0 then means lambda underflowed.  The exact law
// in that corner is G == h almost surely, hence pi = 1 and pi·(1-pi) = 0.
// Substituting (e, alpha) <- (1, 1) while leaving u == 0 reproduces both values
// exactly, and is the only degenerate input this core can receive.
//
// h == 0 needs no special case: K', K'' and K are all proportional to h, and
// every factor they multiply is finite, so all three return exactly zero.
inline Cgf12 subjectK12(double t, double r, double p, double h) noexcept {
    const double u  = 1.0 - p;
    const double tr = clampTr(t * r);
    const double e  = p * std::exp(tr);

    double a  = u + e;                   // sum of non-negatives — never cancels
    double en = e;
    if (a == 0.0) { a = 1.0; en = 1.0; } // p == 1 and lambda underflowed: pi = 1

    const double inva = 1.0 / a;
    const double pi   = en * inva;       // pi   = e/alpha  in [0, 1]
    const double omp  = u  * inva;       // 1-pi = u/alpha  in [0, 1] — a quotient,
                                         //        never the subtraction 1 - pi
    const double hr = h * r;
    return Cgf12{hr * pi, hr * r * pi * omp};
}

// Terminal core: adds K = h·log(alpha).  See `logAlpha` for the N1 treatment.
inline Cgf012 subjectKFull(double t, double r, double p, double h) noexcept {
    const Cgf12 d = subjectK12(t, r, p, h);
    return Cgf012{h * logAlpha(t * r, p), d.K1, d.K2};
}

// ══════════════════════════════════════════════════════════════════════
// Production entry points — sums over the supplied subject set
// ══════════════════════════════════════════════════════════════════════
//
// Each dispatches to the widest SIMD tier the host supports.  `n == 0` returns
// all zeros.  Preconditions: `resid` has n entries; allele frequencies lie in
// [0, 1]; hapcounts are non-negative.

// binomUniform — G_i ~ Binomial(2, p), residual weight resid[i].
//
// The p == 1 and p == 0 endpoints are resolved in closed form before the loop is
// entered, so the vector body may assume 0 < p < 1 and therefore
// alpha >= u > 0.  That removes the per-lane alpha == 0 guard from the hot path
// of the shared-frequency variant.
Cgf12  binomUniformK12(double t, const double *resid, int n, double p) noexcept;
Cgf012 binomUniformKFull(double t, const double *resid, int n, double p) noexcept;

// binomIndiv — G_i ~ Binomial(2, af[i]), residual weight resid[i].
//
// af varies per lane, so the alpha == 0 guard is applied per lane here; it
// cannot be hoisted.
Cgf12  binomIndivK12(double t, const double *resid, const double *af, int n) noexcept;
Cgf012 binomIndivKFull(double t, const double *resid, const double *af, int n) noexcept;

// binomHapcount — G_i ~ Binomial(hap[i], q), residual weight resid[i].
//
// q is shared, so its endpoints are hoisted exactly as in binomUniform; hap
// varies per lane but needs no guard, every cumulant being proportional to it.
Cgf12  binomHapcountK12(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
Cgf012 binomHapcountKFull(
    double t, const double *resid, const double *hap, int n, double q) noexcept;

// ══════════════════════════════════════════════════════════════════════
// The exact-relative-accuracy terminal entry points
// ══════════════════════════════════════════════════════════════════════
//
// Identical to `binom*KFull` in K' and K'' — the same dispatched k12 kernels —
// but accumulating K through `logAlpha`, i.e. N1's log1p(p*expm1(x)) spelling,
// in a scalar loop.  Relatively accurate in K unconditionally, and about six
// times slower.
//
// Production uses `binom*KFull`.  These exist because the relative-accuracy
// property is worth keeping proved and reachable: `tests/spa_cgf_test.cpp` pins
// N1 on them, and any future consumer that needs K itself rather than
// zeta*s - K (a moment-generating-function evaluation, say, or a diagnostic that
// reports K directly) has an entry point that delivers it.
Cgf012 binomUniformKFullExact(
    double t, const double *resid, int n, double p) noexcept;
Cgf012 binomIndivKFullExact(
    double t, const double *resid, const double *af, int n) noexcept;
Cgf012 binomHapcountKFullExact(
    double t, const double *resid, const double *hap, int n, double q) noexcept;

// ══════════════════════════════════════════════════════════════════════
// Tier-specific variants — exposed for cross-tier equality testing
// ══════════════════════════════════════════════════════════════════════
//
// Production code should call the dispatched entry points above.  These are
// declared so `tests/spa_cgf_test.cpp` can compare the tiers against each other
// and against the long-double reference on every tier the host supports.  A
// variant whose ISA the host lacks must never be called: guard with
// `simdLevelValue()`, following `tests/lanc_simd_test.cpp`.
namespace tier {

Cgf12 binomUniformK12_scalar(double t, const double *resid, int n, double p) noexcept;
Cgf12 binomIndivK12_scalar(
    double t, const double *resid, const double *af, int n) noexcept;
Cgf12 binomHapcountK12_scalar(
    double t, const double *resid, const double *hap, int n, double q) noexcept;

#if defined(__x86_64__) || defined(_M_X64)
Cgf12 binomUniformK12_avx2(double t, const double *resid, int n, double p) noexcept;
Cgf12 binomUniformK12_avx512(double t, const double *resid, int n, double p) noexcept;
Cgf12 binomIndivK12_avx2(
    double t, const double *resid, const double *af, int n) noexcept;
Cgf12 binomIndivK12_avx512(
    double t, const double *resid, const double *af, int n) noexcept;
Cgf12 binomHapcountK12_avx2(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
Cgf12 binomHapcountK12_avx512(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
#endif

// The K half of kFull, as a bare reduction.  binomUniform and binomIndiv return
// sum_i log(alpha_i) with the common factor h == 2 left to the caller, which is
// bit-identical to folding it in because multiplication by two is exact;
// binomHapcount returns sum_i hap_i * log(alpha_i), h varying per subject.
//
// Preconditions match the k12 kernels: binomUniform and binomHapcount require
// 0 < p < 1, their endpoints being resolved by the public entry point.
double binomUniformK0_scalar(double t, const double *resid, int n, double p) noexcept;
double binomIndivK0_scalar(
    double t, const double *resid, const double *af, int n) noexcept;
double binomHapcountK0_scalar(
    double t, const double *resid, const double *hap, int n, double q) noexcept;

#if defined(__x86_64__) || defined(_M_X64)
double binomUniformK0_avx2(double t, const double *resid, int n, double p) noexcept;
double binomUniformK0_avx512(double t, const double *resid, int n, double p) noexcept;
double binomIndivK0_avx2(
    double t, const double *resid, const double *af, int n) noexcept;
double binomIndivK0_avx512(
    double t, const double *resid, const double *af, int n) noexcept;
double binomHapcountK0_avx2(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
double binomHapcountK0_avx512(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
#endif

}  // namespace tier

// Runtime SIMD tier actually selected: 0 = scalar, 1 = AVX2, 2 = AVX-512.
// Mirrors `simdLevel()` without exposing the enum to test translation units.
int simdLevelValue();

}  // namespace spa_cgf
