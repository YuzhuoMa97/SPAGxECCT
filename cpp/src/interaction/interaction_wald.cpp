// interaction_wald.cpp — Branch-B Wald leg of SPAGxE_CCT (per-marker G×E dispatcher).
//
// SPAGxE_CCT (Ma et al., Nat. Commun. 2025) combines, for the small fraction of
// variants with a significant marginal genetic effect (Branch B, p_marg ≤ ε), a
// retrospective saddlepoint G×E p-value with a prospective Wald p-value of the
// interaction coefficient in the FULL model  trait ~ [covar] + g + g:E , then
// reports  P = CCT(p_spa, p_wald)  (Cauchy combination; the log-domain
// math::cauchyCombineLog10, over the two magnitudes).
//
// This file assembles the full-interaction design  M = [covar | g | g∘E]  per
// marker (the interaction g∘E is the appended last column), drops incomplete-
// case rows, and calls the matching standard-model Wald fitter in
// src/util/wald.hpp to obtain the MAGNITUDE  L = −log10 P  of the two-sided
// Wald p-value of the G:E coefficient (the log-domain implementation; the fitters
// returned a linear p until then, and that p was the last quantity in this
// path with an underflow ceiling).
// The fitters themselves (OLS / logistic IRLS / Breslow-Cox Newton /
// proportional-odds Fisher scoring) live in namespace wald so any method can
// share them without depending on SPAGxE.
//
// A singular information matrix, non-convergence, or a degenerate design returns
// NaN, which the caller folds into the Cauchy combination as a fall-back to the
// SPA p (never as p = 0 — missing-test convention).

#include "interaction/interaction_wald.hpp"

#include "util/wald.hpp" // wald::lastCoef*Log10P, wald::OrdinalInfo

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace interaction_wald {

namespace {
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Keep an intercept in the rank calculation even for Cox/ordinal models:
// their baseline hazard/thresholds absorb constants. Select original nuisance
// columns in order, so the tested interaction keeps its original coefficient
// and units. Twice-reorthogonalized, normalized columns make the independence
// decision relative to each column rather than its arbitrary measurement unit.
bool reduceNuisanceColumns(Eigen::MatrixXd &design) {
    constexpr double tolerance = 1e-10;
    const Eigen::Index n = design.rows();
    const Eigen::Index target = design.cols() - 1;
    Eigen::MatrixXd basis(n, std::min(n, design.cols()));
    basis.col(0).setConstant(1.0 / std::sqrt(static_cast<double>(n)));
    Eigen::Index rank = 1;
    std::vector<Eigen::Index> kept;

    auto independentResidual = [&](Eigen::Index column, Eigen::VectorXd &residual) {
        const double scale = design.col(column).cwiseAbs().maxCoeff();
        if (!(scale > 0.0)) return false;
        residual = design.col(column) / scale;
        residual /= residual.stableNorm();
        for (int pass = 0; pass < 2; ++pass)
            residual -= basis.leftCols(rank) *
                        (basis.leftCols(rank).transpose() * residual).eval();
        return residual.stableNorm() > tolerance;
    };

    for (Eigen::Index j = 0; j < target; ++j) {
        Eigen::VectorXd residual;
        if (rank == n || !independentResidual(j, residual)) continue;
        basis.col(rank++) = residual / residual.stableNorm();
        kept.push_back(j);
    }
    Eigen::VectorXd targetResidual;
    if (rank == n || !independentResidual(target, targetResidual)) return false;

    Eigen::MatrixXd reduced(n, static_cast<Eigen::Index>(kept.size()) + 1);
    for (size_t j = 0; j < kept.size(); ++j)
        reduced.col(static_cast<Eigen::Index>(j)) = design.col(kept[j]);
    reduced.col(reduced.cols() - 1) = design.col(target);
    design.swap(reduced);
    return true;
}
} // namespace

// ══════════════════════════════════════════════════════════════════════════
// Dispatcher — assemble  M = [covar | g | g∘E], drop NaN rows, fit per trait.
// ══════════════════════════════════════════════════════════════════════════

double waldInteractionLog10P(
    const WaldData &wd,
    const Eigen::Ref<const Eigen::VectorXd> &g,
    const Eigen::Ref<const Eigen::VectorXd> &E,
    wald::OrdinalInfo ordInfo,
    bool reduceNuisance
) {
    if (wd.trait == TraitType::None) return kNaN;
    const Eigen::Index n = wd.covar.rows();
    const Eigen::Index c = wd.covar.cols();
    if (g.size() != n || E.size() != n) return kNaN;

    // M = [covar | g | g∘E]; the interaction g∘E is the last column.
    Eigen::MatrixXd M(n, c + 2);
    if (c > 0) M.leftCols(c) = wd.covar;
    M.col(c) = g;
    M.col(c + 1) = g.array() * E.array();

    const bool isCox = (wd.trait == TraitType::Cox);

    // Complete-case rows: no NaN in M or in the trait's response(s).
    std::vector<Eigen::Index> keep;
    keep.reserve(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i) {
        bool ok = M.row(i).allFinite();
        if (ok) {
            if (isCox)
                ok = std::isfinite(wd.time[i]) && std::isfinite(wd.event[i]);
            else
                ok = std::isfinite(wd.y[i]);
        }
        if (ok) keep.push_back(i);
    }
    const Eigen::Index m = static_cast<Eigen::Index>(keep.size());
    if (m == 0) return kNaN;

    Eigen::MatrixXd Mf(m, M.cols());
    for (Eigen::Index i = 0; i < m; ++i) Mf.row(i) = M.row(keep[i]);
    if (reduceNuisance && !reduceNuisanceColumns(Mf)) return kNaN;

    switch (wd.trait) {
    case TraitType::Linear:
    case TraitType::Logistic: {
        // Prepend the intercept: Z = [1 | covar | g | g∘E].
        Eigen::MatrixXd Z(m, Mf.cols() + 1);
        Z.col(0).setOnes();
        Z.rightCols(Mf.cols()) = Mf;
        Eigen::VectorXd yf(m);
        for (Eigen::Index i = 0; i < m; ++i) yf[i] = wd.y[keep[i]];
        return (wd.trait == TraitType::Linear) ? wald::lastCoefLinearLog10P(Z, yf)
                                               : wald::lastCoefLogisticLog10P(Z, yf);
    }
    case TraitType::Cox: {
        Eigen::VectorXd tf(m), ef(m);
        for (Eigen::Index i = 0; i < m; ++i) {
            tf[i] = wd.time[keep[i]];
            ef[i] = wd.event[keep[i]];
        }
        return wald::lastCoefCoxLog10P(tf, ef, Mf);
    }
    case TraitType::Ordinal: {
        Eigen::VectorXi yf(m);
        for (Eigen::Index i = 0; i < m; ++i)
            yf[i] = static_cast<int>(std::lround(wd.y[keep[i]]));
        return wald::lastCoefOrdinalLog10P(yf, Mf, ordInfo);
    }
    default:
        return kNaN;
    }
}

} // namespace interaction_wald
