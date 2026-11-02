// Independent checks of partial CGFs, partitioning and the full-law safeguard.
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "util/partial_spa.hpp"
#include "tinytest.hpp"

namespace {
using ld = long double;
struct Reference { ld k0 = 0, k1 = 0, k2 = 0; };

Reference scalarReference(const Eigen::VectorXd &w, double q,
                          const Eigen::VectorXd *af,
                          const Eigen::VectorXd *hap,
                          const partial_spa::Partition &part, ld t) {
    std::vector<bool> exact(static_cast<size_t>(w.size()), false);
    for (uint32_t i : part.exactPositions) exact[i] = true;
    Reference out;
    for (Eigen::Index i = 0; i < w.size(); ++i) {
        const ld r = w[i], p = af ? (*af)[i] : q;
        const ld h = hap ? (*hap)[i] : 2;
        if (r == 0 || h == 0 || p == 0) continue;
        if (!exact[static_cast<size_t>(i)]) {
            const ld mu = h * p * r, var = h * p * (1 - p) * r * r;
            out.k0 += mu * t + var * t * t / 2;
            out.k1 += mu + var * t;
            out.k2 += var;
        } else if (p == 1) {
            out.k0 += h * t * r;
            out.k1 += h * r;
        } else {
            const ld e = std::exp(t * r), a = 1 - p + p * e;
            out.k0 += h * std::log(a);
            out.k1 += h * r * p * e / a;
            out.k2 += h * r * r * p * (1 - p) * e / (a * a);
        }
    }
    return out;
}

Eigen::VectorXd regularWeights(int n = 1000, double central = 0.02,
                               double extreme = 4.0) {
    Eigen::VectorXd w(n + 2);
    for (int i = 0; i < n; ++i) w[i] = (i % 2 ? central : -central);
    w[n] = -extreme;
    w[n + 1] = extreme;
    return w;
}

spa::Result guarded(const partial_spa::Workspace &ws, double dev,
                    binomial_cgf::ApproximationInfo &info) {
    return binomial_cgf::twoSidedSpa(
        ws.partial, ws.fullMean, dev, ws.fullVariance,
        dev / std::sqrt(ws.fullVariance), 1e-9, 1.0, true, &ws.full, &info);
}

spa::Result full(const partial_spa::Workspace &ws, double dev) {
    return binomial_cgf::twoSidedSpa(
        ws.full, ws.fullMean, dev, ws.fullVariance,
        dev / std::sqrt(ws.fullVariance), 1e-9, 1.0, true);
}
} // namespace

TEST(partial_cgf_matches_independent_scalar_law) {
    Eigen::VectorXd w = regularWeights(100, 0.2, 4);
    Eigen::VectorXd af(w.size()), hap(w.size());
    for (Eigen::Index i = 0; i < w.size(); ++i) {
        af[i] = 0.05 + 0.003 * (i % 80);
        hap[i] = static_cast<double>(i % 3);
    }
    af[7] = 0.0;
    af[8] = 1.0;
    for (int law = 0; law < 3; ++law) {
        const Eigen::VectorXd *a = law == 1 ? &af : nullptr;
        const Eigen::VectorXd *h = law == 2 ? &hap : nullptr;
        const auto part = partial_spa::split(w, h, a);
        REQUIRE(!part.allExact);
        partial_spa::Workspace ws;
        partial_spa::bind(part, w, 0.21, a, h, ws);
        for (double t : {-0.7, -0.03, 0.0, 0.04, 0.8}) {
            const auto ref = scalarReference(w, 0.21, a, h, part, t);
            const auto actual = binomial_cgf::kFull(t, ws.partial, 0.31);
            const auto shortActual = binomial_cgf::k12(t, ws.partial, 0.31);
            CHECK_CLOSE(actual.k0, static_cast<double>(ref.k0), 3e-12, 2e-13);
            CHECK_CLOSE(actual.k1, static_cast<double>(ref.k1) - 0.31, 3e-12, 2e-13);
            CHECK_CLOSE(actual.k2, static_cast<double>(ref.k2), 3e-12, 2e-13);
            CHECK_CLOSE(shortActual.k1, actual.k1, 2e-14, 2e-13);
            CHECK_CLOSE(shortActual.k2, actual.k2, 2e-14, 2e-13);
        }
        const auto zero = binomial_cgf::kFull(0, ws.partial, 0);
        CHECK_NEAR(zero.k0, 0.0, 1e-14);
        CHECK_CLOSE(zero.k1, ws.fullMean, 2e-14, 2e-13);
        CHECK_CLOSE(zero.k2, ws.fullVariance, 2e-14, 2e-13);
    }
}

TEST(partial_derivatives_match_finite_differences) {
    const auto w = regularWeights(100, 0.2, 4);
    const auto part = partial_spa::split(w);
    partial_spa::Workspace ws;
    partial_spa::bind(part, w, 0.17, nullptr, nullptr, ws);
    const double dt = 1e-5;
    for (double t : {-0.3, 0.1, 0.9}) {
        const auto l = binomial_cgf::kFull(t - dt, ws.partial, 0);
        const auto m = binomial_cgf::kFull(t, ws.partial, 0);
        const auto r = binomial_cgf::kFull(t + dt, ws.partial, 0);
        CHECK_CLOSE((r.k0 - l.k0) / (2 * dt), m.k1, 2e-8, 2e-9);
        CHECK_CLOSE((r.k1 - l.k1) / (2 * dt), m.k2, 2e-8, 2e-9);
    }
}

TEST(type7_tukey_partition_and_degenerate_cases) {
    Eigen::VectorXd w(8);
    w << -10, 1, 2, 3, 4, 5, 6, 20;
    const auto p = partial_spa::split(w);
    // Type-7 Q1=1.75, Q3=5.25: fences -3.5 and 10.5.
    REQUIRE(p.exactPositions.size() == 2);
    CHECK(p.exactPositions[0] == 0);
    CHECK(p.exactPositions[1] == 7);
    CHECK(p.centralPositions.size() == 6);
    CHECK_NEAR(p.centralWeightSum, 21.0, 0.0);
    CHECK_NEAR(p.centralWeightSq, 91.0, 0.0);
    Eigen::VectorXd binary(100);
    binary.setOnes();
    binary[0] = -20;
    const auto b = partial_spa::split(binary);
    CHECK(!b.allExact); // Zero IQR with a true outlier keeps a central block.
    CHECK(b.exactPositions.size() == 1);
    CHECK(partial_spa::split(Eigen::VectorXd::Ones(100)).allExact);
    CHECK(partial_spa::split(Eigen::VectorXd::Zero(100)).allExact);
    CHECK(partial_spa::split(Eigen::VectorXd(0)).allExact);
    Eigen::VectorXd balanced(100);
    for (int i = 0; i < 100; ++i) balanced[i] = i % 2 ? 1 : -1;
    CHECK(partial_spa::split(balanced).allExact);
}

TEST(inactive_rows_do_not_select_outliers) {
    const auto w = regularWeights(100, 0.02, 4);
    Eigen::VectorXd h = Eigen::VectorXd::Constant(w.size(), 2.0);
    const auto p = partial_spa::split(w, &h);
    Eigen::VectorXd extended(w.size() + 50), extH(w.size() + 50);
    extended.head(w.size()) = w;
    extH.head(w.size()) = h;
    for (int i = 0; i < 50; ++i) {
        extended[w.size() + i] = (i % 2 ? 1e100 : -1e100);
        extH[w.size() + i] = 0;
    }
    const auto e = partial_spa::split(extended, &extH);
    CHECK(p.exactPositions == e.exactPositions);
    CHECK(p.nActive == e.nActive);
    partial_spa::Workspace a, b;
    partial_spa::bind(p, w, 0.2, nullptr, &h, a);
    partial_spa::bind(e, extended, 0.2, nullptr, &extH, b);
    for (double t : {-0.2, 0.2}) {
        const auto x = binomial_cgf::kFull(t, a.partial, 0);
        const auto y = binomial_cgf::kFull(t, b.partial, 0);
        CHECK_NEAR(x.k0, y.k0, 1e-13);
        CHECK_NEAR(x.k1, y.k1, 1e-13);
        CHECK_NEAR(x.k2, y.k2, 1e-13);
    }
    // Deterministic AF rows remain in the score mean but not the fences.
    Eigen::VectorXd af = Eigen::VectorXd::Constant(extended.size(), 0.2);
    for (int i = 0; i < 50; ++i) {
        extended[w.size() + i] = 1.0;
        af[w.size() + i] = i % 2 ? 1.0 : 0.0;
    }
    const auto f = partial_spa::split(extended, nullptr, &af);
    CHECK(p.exactPositions == f.exactPositions);
    partial_spa::bind(f, extended, 0, &af, nullptr, b);
    CHECK_NEAR(b.fullMean, 50.0, 1e-12);
    CHECK_NEAR(binomial_cgf::kFull(0, b.partial, 0).k1, 50.0, 1e-12);
}

TEST(central_variance_is_not_total_minus_outlier) {
    Eigen::VectorXd w = Eigen::VectorXd::Ones(101);
    w[100] = 1e12;
    const auto p = partial_spa::split(w);
    partial_spa::Workspace ws;
    partial_spa::bind(p, w, 0.2, nullptr, nullptr, ws);
    CHECK_NEAR(ws.partial.var, 32.0, 1e-13);
    CHECK_NEAR(ws.partial.mean, 40.0, 1e-13);
}

TEST(guard_accepts_accurate_partial_tails) {
    const auto w = regularWeights();
    const auto p = partial_spa::split(w);
    partial_spa::Workspace ws;
    partial_spa::bind(p, w, 0.2, nullptr, nullptr, ws);
    binomial_cgf::ApproximationInfo info;
    const double dev = 2 * std::sqrt(ws.fullVariance);
    const auto a = guarded(ws, dev, info);
    const auto b = full(ws, dev);
    CHECK(spa::statusIsUsable(a.status));
    CHECK(info.attemptedGaussian);
    CHECK(info.usedGaussian);
    CHECK(!info.fullFallback);
    CHECK(info.exactCount == 2);
    CHECK(info.centralCount == 1000);
    CHECK(info.maxK0Error <= 1e-3);
    CHECK(info.maxStandardizedResidual <= 1e-3);
    CHECK(info.maxRelativeK2Error <= 1e-3);
    CHECK_NEAR(a.negLog10p, b.negLog10p, 1e-4);
}

TEST(rare_central_guard_recomputes_both_full_tails) {
    const auto w = regularWeights(10000, 1.0, 10.0);
    const auto p = partial_spa::split(w);
    partial_spa::Workspace ws;
    partial_spa::bind(p, w, 1e-4, nullptr, nullptr, ws);
    binomial_cgf::ApproximationInfo info;
    const auto a = guarded(ws, 8.0, info);
    const auto b = full(ws, 8.0);
    CHECK(info.attemptedGaussian);
    CHECK(info.fullFallback);
    CHECK(!info.usedGaussian);
    CHECK(!info.fullSupportFallback);
    CHECK(a.status == b.status);
    CHECK_NEAR(a.negLog10p, b.negLog10p, 0.0);
}

TEST(all_three_laws_use_gaussian_and_guard_rare_tails) {
    for (int law = 0; law < 3; ++law) {
        for (bool rare : {false, true}) {
            const auto w = rare ? regularWeights(10000, 1.0, 10.0)
                                : regularWeights();
            const double q = rare ? 1e-4 : 0.2;
            Eigen::VectorXd af(w.size()), hap(w.size());
            for (Eigen::Index i = 0; i < w.size(); ++i) {
                // Adjacent +/- weights share AF/H, so the score mean is zero.
                af[i] = q * (0.75 + 0.05 * ((i / 2) % 10));
                hap[i] = 1.0 + ((i / 2) % 2);
            }
            af[w.size() - 2] = af[w.size() - 1] = q;
            hap[w.size() - 2] = hap[w.size() - 1] = 2.0;
            const Eigen::VectorXd *a = law == 1 ? &af : nullptr;
            const Eigen::VectorXd *h = law == 2 ? &hap : nullptr;
            const auto p = partial_spa::split(w, h, a);
            partial_spa::Workspace ws;
            partial_spa::bind(p, w, q, a, h, ws);
            REQUIRE(ws.partial.var > 0.0);
            binomial_cgf::ApproximationInfo info;
            const double dev = rare ? 8.0 : 2 * std::sqrt(ws.fullVariance);
            const auto result = guarded(ws, dev, info);
            const auto reference = full(ws, dev);
            CHECK(spa::statusIsUsable(result.status));
            CHECK(std::isfinite(result.negLog10p));
            CHECK(result.negLog10p > 0.0);
            CHECK(info.attemptedGaussian);
            CHECK(info.usedGaussian == !rare);
            CHECK(info.fullFallback == rare);
            CHECK_NEAR(result.negLog10p, reference.negLog10p, rare ? 0.0 : 1e-4);
            std::printf("    law=%s rare=%d gaussian=%d fallback=%d dLog10P=%.9g "
                        "dK0=%.9g stdResidual=%.9g relK2=%.9g\n",
                        law == 0 ? "shared" : law == 1 ? "individual" : "known-H",
                        int(rare), int(info.usedGaussian), int(info.fullFallback),
                        std::abs(result.negLog10p - reference.negLog10p),
                        info.maxK0Error, info.maxStandardizedResidual,
                        info.maxRelativeK2Error);
        }
    }
}

TEST(original_support_selects_full_law_before_endpoint_atoms) {
    Eigen::VectorXd w = Eigen::VectorXd::Constant(101, 0.1);
    w[100] = 1;
    Eigen::VectorXd h = Eigen::VectorXd::Ones(w.size());
    const auto p = partial_spa::split(w, &h);
    partial_spa::Workspace ws;
    partial_spa::bind(p, w, 0.25, nullptr, &h, ws);
    binomial_cgf::ApproximationInfo info;
    const auto a = guarded(ws, ws.fullMean, info); // Lower target is exactly 0.
    const auto b = full(ws, ws.fullMean);
    CHECK(info.fullSupportFallback);
    CHECK(info.fullFallback);
    CHECK(!info.usedGaussian);
    CHECK_NEAR(a.negLog10p, b.negLog10p, 0.0);
    // A Gaussian surrogate alone has no finite binomial support.
    binomial_cgf::Context gaussian;
    gaussian.mean = 2.0;
    gaussian.var = 3.0;
    const auto g = binomial_cgf::twoSidedSpa(gaussian, 2.0, 4.0, 3.0,
                                           4.0 / std::sqrt(3.0), 1e-9, 1, true);
    CHECK(spa::statusIsUsable(g.status));
    CHECK(std::isfinite(g.negLog10p));
    CHECK(g.negLog10p > 1.0);
}

TEST(no_outliers_and_invalid_or_deterministic_laws) {
    Eigen::VectorXd w = Eigen::VectorXd::Ones(10);
    const auto p = partial_spa::split(w);
    partial_spa::Workspace ws;
    partial_spa::bind(p, w, 0.2, nullptr, nullptr, ws);
    CHECK(ws.partial.var == 0);
    CHECK(ws.partial.nOutlier == w.size());
    binomial_cgf::ApproximationInfo info;
    const auto a = guarded(ws, 1, info);
    const auto b = full(ws, 1);
    CHECK(!info.attemptedGaussian);
    CHECK(!info.fullFallback);
    CHECK_NEAR(a.negLog10p, b.negLog10p, 0.0);
    Eigen::VectorXd h = Eigen::VectorXd::Zero(10);
    const auto zero = partial_spa::split(w, &h);
    partial_spa::bind(zero, w, 0.2, nullptr, &h, ws);
    CHECK(ws.fullVariance == 0.0);
    CHECK(ws.fullMean == 0.0);
    bool rejected = false;
    try { partial_spa::split(w, &h, &h); }
    catch (const std::invalid_argument &) { rejected = true; }
    CHECK(rejected);
    rejected = false;
    h[0] = 0.5;
    try { partial_spa::split(w, &h); }
    catch (const std::invalid_argument &) { rejected = true; }
    CHECK(rejected);
    rejected = false;
    h[0] = 1.1;
    try { partial_spa::split(w, nullptr, &h); }
    catch (const std::invalid_argument &) { rejected = true; }
    CHECK(rejected);
}

TEST(support_precheck_preserves_boundaries_without_endpoint_logs) {
    Eigen::VectorXd w(8), h(8), q(8);
    w << -3, 4, -2, 0, 7, -5, 11, -13;
    h << 2, 1, 0, 2, 1, 2, 0, 1;
    q << .2, .7, .4, 1, 0, 1, .5, .1;
    for (bool local : {false, true}) {
        binomial_cgf::Context c;
        c.resid=w.data(); c.nOutlier=8; c.q=.3;
        if (local) c.hap=h.data(); else c.af=q.data();
        const auto actual=binomial_cgf::detail::supportOf(c);
        const auto bounds=binomial_cgf::detail::supportOf(c,false);
        CHECK(actual.lower==bounds.lower);
        CHECK(actual.upper==bounds.upper);
        CHECK(actual.scale==bounds.scale);
        CHECK(actual.valid==bounds.valid);
        for (double target : {static_cast<double>(actual.lower)-1,
                              static_cast<double>(actual.lower), 0.0,
                              static_cast<double>(actual.upper),
                              static_cast<double>(actual.upper)+1})
            for (bool lower : {false,true}) {
                double pa=0,pb=0;
                CHECK(actual.tail(target,lower,pa)==bounds.tail(target,lower,pb));
            }
    }
}

TEST(inactive_large_weights_do_not_overflow_unused_uniform_caches) {
    Eigen::VectorXd w=Eigen::VectorXd::Ones(40);
    w[0]=1e200; w[1]=20;
    Eigen::VectorXd h=Eigen::VectorXd::Ones(40);h[0]=0;
    Eigen::VectorXd q=Eigen::VectorXd::Constant(40,.2);q[0]=0;
    partial_spa::Workspace ws;
    auto p=partial_spa::split(w,&h);
    partial_spa::bind(p,w,.2,nullptr,&h,ws);
    CHECK(std::isfinite(ws.fullVariance));
    CHECK(std::isfinite(ws.partial.var));
    p=partial_spa::split(w,nullptr,&q);
    partial_spa::bind(p,w,0,&q,nullptr,ws);
    CHECK(std::isfinite(ws.fullVariance));
    CHECK(std::isfinite(ws.partial.var));
}

TINYTEST_MAIN
