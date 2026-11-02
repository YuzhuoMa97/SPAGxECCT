// Independent-binomial saddlepoint calculation with a central Gaussian block.
// The exact block permits shared AF, individual AF, or known haplotype counts.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

#include "util/spa.hpp"
#include "util/spa_cgf.hpp"
#include "util/binomial_support.hpp"

namespace binomial_cgf {

struct Context {
    const double *resid = nullptr;
    const double *af = nullptr;
    int nOutlier = 0;
    double q = 0.0;
    double mean = 0.0;
    double var = 0.0;
    // Known ancestry copies (0, 1 or 2). Individual AF and hap are exclusive.
    const double *hap = nullptr;
};

// usedGaussian describes the returned result, not merely an attempted solve.
// The three diagnostics compare the full CGF at each partial saddlepoint.
// They control approximation error, not the calibration error of SPA itself.
struct ApproximationInfo {
    bool attemptedGaussian = false;
    bool usedGaussian = false;
    bool fullFallback = false;
    bool fullSupportFallback = false;
    int exactCount = 0;
    int centralCount = 0;
    double maxK0Error = 0.0;
    double maxStandardizedResidual = 0.0;
    double maxRelativeK2Error = 0.0;
};

inline spa::K12 k12(double t, const Context &c, double s) noexcept {
    const spa_cgf::Cgf12 d = c.hap
        ? spa_cgf::binomHapcountK12(t, c.resid, c.hap, c.nOutlier, c.q)
        : c.af ? spa_cgf::binomIndivK12(t, c.resid, c.af, c.nOutlier)
               : spa_cgf::binomUniformK12(t, c.resid, c.nOutlier, c.q);
    return spa::K12{d.K1 + c.mean + c.var * t - s, d.K2 + c.var};
}

inline spa::K012 kFull(double t, const Context &c, double s) noexcept {
    const spa_cgf::Cgf012 d = c.hap
        ? spa_cgf::binomHapcountKFull(t, c.resid, c.hap, c.nOutlier, c.q)
        : c.af ? spa_cgf::binomIndivKFull(t, c.resid, c.af, c.nOutlier)
               : spa_cgf::binomUniformKFull(t, c.resid, c.nOutlier, c.q);
    return spa::K012{d.K0 + t * (c.mean + 0.5 * c.var * t),
                     d.K1 + c.mean + c.var * t - s,
                     d.K2 + c.var};
}

namespace detail {

inline spa::BinomialSupport supportOf(const Context &c, bool computeAtomMasses = true) noexcept {
    spa::BinomialSupport support(c.mean);
    for (int i = 0; i < c.nOutlier; ++i)
        support.add(c.resid[i], c.hap ? c.hap[i] : 2.0,
                    c.af ? c.af[i] : c.q, computeAtomMasses);
    return support;
}

inline bool approximationAccepted(const spa::Saddle &sd,
                                  const Context &full, double target,
                                  ApproximationInfo *info) noexcept {
    const spa::K012 ref = kFull(sd.zeta, full, target);
    if (!std::isfinite(ref.k0) || !std::isfinite(ref.k1) ||
        !std::isfinite(ref.k2) || !(ref.k2 > 0.0) ||
        !std::isfinite(sd.K0) || !std::isfinite(sd.K2)) return false;

    const double err0 = std::abs(ref.k0 - sd.K0);
    const double err1 = std::abs(ref.k1) / std::sqrt(ref.k2);
    const double err2 = std::abs(ref.k2 - sd.K2) / ref.k2;
    if (info) {
        info->maxK0Error = std::max(info->maxK0Error, err0);
        info->maxStandardizedResidual =
            std::max(info->maxStandardizedResidual, err1);
        info->maxRelativeK2Error = std::max(info->maxRelativeK2Error, err2);
    }
    // The absolute K0 check controls the exponential part of the SPA tail;
    // the other checks guard root displacement and the curvature correction.
    // Roundoff floors accommodate cancellation in a nonzero score mean.
    constexpr double tolerance = 1e-3;
    constexpr double eps = 64.0 * std::numeric_limits<double>::epsilon();
    const double floor0 = eps * (1.0 + std::abs(ref.k0) +
                                std::abs(sd.K0) + std::abs(sd.zeta * target));
    const double floor1 = eps * (1.0 + std::abs(target) +
                                std::abs(ref.k1 + target)) / std::sqrt(ref.k2);
    return std::isfinite(err0) && std::isfinite(err1) && std::isfinite(err2) &&
           err0 <= tolerance + floor0 &&
           err1 <= tolerance + floor1 && err2 <= tolerance + eps;
}

} // namespace detail

// Tails are reflected about K'(0). GRM adjustment belongs in the centered
// deviation supplied by the caller; it must not rescale this mean.
// A failed partial solve or diagnostic recomputes BOTH tails with the full
// independent-binomial CGF. Existing full-CGF failure/status semantics remain.
inline spa::Result twoSidedSpa(
    const Context &c,
    double sMean,
    double absDev,
    double indepVar,
    double zNorm,
    double rtol = 1e-6,
    double tailScale = 1.0,
    bool supportAware = false,
    const Context *fullReference = nullptr,
    ApproximationInfo *info = nullptr
) noexcept {
    const bool guardedPartial = c.var > 0.0 && fullReference != nullptr;
    if (info) {
        *info = ApproximationInfo{};
        info->attemptedGaussian = c.var > 0.0;
        info->usedGaussian = c.var > 0.0;
        info->exactCount = c.nOutlier;
        info->centralCount = fullReference
            ? std::max(0, fullReference->nOutlier - c.nOutlier) : 0;
    }
    const auto fullFallback = [&](bool supportBoundary) {
        if (info) {
            info->fullFallback = true;
            info->fullSupportFallback = supportBoundary;
            info->usedGaussian = false;
        }
        return twoSidedSpa(*fullReference, sMean, absDev, indepVar, zNorm,
                           rtol, tailScale, true, nullptr, nullptr);
    };

    // A Gaussian surrogate has unbounded support. At an original support
    // endpoint, select the FULL law before applying its binomial endpoint mass.
    if (guardedPartial && fullReference->var == 0.0) {
        // This check selects a law; only the subsequent full solve needs the
        // endpoint probabilities. Bounds and rounding tolerances are identical.
        const spa::BinomialSupport original = detail::supportOf(*fullReference, false);
        double ignored = 0.0;
        if (original.tail(sMean + absDev, false, ignored, tailScale == 1.0) ||
            original.tail(sMean - absDev, true, ignored, tailScale == 1.0))
            return fullFallback(true);
    }

    double zeta0 = 0.0;
    if (indepVar > 0.0) {
        const double e = absDev / indepVar;
        if (std::isfinite(e) && e > 0.0) zeta0 = std::min(e, 1.2);
    }

    double logp[2];
    spa::Status st[2];
    const bool bounded = supportAware && c.var == 0.0;
    const spa::BinomialSupport support = bounded
        ? detail::supportOf(c) : spa::BinomialSupport(c.mean);
    for (int k = 0; k < 2; ++k) {
        const bool lowerTail = (k == 1);
        const double s = lowerTail ? sMean - absDev : sMean + absDev;
        if (bounded && support.tail(s, lowerTail, logp[k], tailScale == 1.0)) {
            st[k] = spa::Status::SpaOk;
            continue;
        }
        spa::SolveOpts opt;
        opt.init = lowerTail ? -zeta0 : zeta0;
        opt.rtol = rtol;
        const spa::Saddle sd = spa::solveSaddlepoint(
            s, [&](double t) { return k12(t, c, s); },
            [&](double t) { return kFull(t, c, s); }, opt);
        if (guardedPartial &&
            (!spa::statusIsUsable(sd.status) ||
             !detail::approximationAccepted(sd, *fullReference, s, info)))
            return fullFallback(false);

        spa::Status stLog = spa::Status::SpaOk;
        logp[k] = spa::bnTailLog(sd.zeta, s, sd.K0, sd.K2,
                               lowerTail, stLog, tailScale);
        st[k] = spa::worseStatus(sd.status, stLog);
        if (guardedPartial &&
            (!spa::statusIsUsable(st[k]) || !std::isfinite(logp[k])))
            return fullFallback(false);
    }
    return spa::assemble(logp[0], logp[1], st[0], st[1], zNorm);
}

inline double statusCode(spa::Status s) noexcept {
    return static_cast<double>(static_cast<uint8_t>(s));
}

} // namespace binomial_cgf
