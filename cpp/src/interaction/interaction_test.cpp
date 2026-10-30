// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE GxE tests with guarded partial-normal saddlepoint acceleration.
// Copyright (c) 2026 Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
// Adapted in part from GRAB infrastructure; see THIRD_PARTY_NOTICES.
//
// One genotype-independent null model supplies residuals for each phenotype.
// Branch A uses (E-lambda)*R; Branch B projects R on [1,G] then uses E*R0.
// The marginal genetic screen and small interaction scores use normal tails.
// Beyond the configured Z threshold, the central weights use a Gaussian CGF
// and outliers retain their binomial law, with full-CGF diagnostic safeguards.
// A sparse GRM changes the score variance and target
// rescaling; raw-trait Wald+CCT is enabled only for unrelated Branch B.

#include "util/version.hpp"
#include "interaction/interaction_test.hpp"

#include "genotype/genotype_source.hpp" // readGenotypeSampleIds, openGenotypeCatalog, GenotypeCatalog
#include "io/sparse_grm.hpp"          // SparseGRM
#include "io/subject_data.hpp"        // SubjectData, extractPhenoVec/Mat, PerPhenoInfo
#include "interaction/interaction_wald.hpp"     // interaction_wald::{WaldData, waldInteractionLog10P}
#include "allele_model/binomial_cgf.hpp"      // binomial_cgf::{Context, twoSidedSpa, statusCode}
#include "allele_model/frequency_model.hpp"        // AFContext, computeAFVec (SPAGxEmix)
#include "util/logging.hpp"           // infoMsg, warnMsg
#include "util/probability_math.hpp"       // math::pnorm, zFromNegLog10P, cauchyCombineLog10
#include "util/spa.hpp"               // spa::normalTwoSidedLog, spa::Status
#include "util/text_stream.hpp"
#include "util/null_model.hpp"        // nullmodel::fitAll and friends
#include "util/partial_spa.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Relatedness rescales the centered score before the independence-CGF SPA.
// Branch-A partitions are cached because shared-AF weights do not change
// between markers. Branch-B weights are partitioned only when SPA is needed.
spa::Result evaluateSharedFrequencyScore(
    const Eigen::VectorXd &weight, double wSum, double wSq, double grmQuad,
    double rawScore, double q, double spaCutoff, double &zScore, double &outVar,
    const partial_spa::Partition *cachedPartition = nullptr
) {
    const double f = 2.0 * q * (1.0 - q);
    outVar = f * grmQuad;
    if (!(outVar > 0.0) || !std::isfinite(outVar)) {
        zScore = 0.0;
        outVar = 0.0;
        return {kNaN, spa::Status::NaNoTest};
    }
    const double mean = 2.0 * q * wSum;
    zScore = (rawScore - mean) / std::sqrt(outVar);
    if (!std::isfinite(zScore)) return {kNaN, spa::Status::NaNoTest};
    if (std::abs(zScore) <= spaCutoff) return spa::normalBranch(zScore);
    const double independentVariance = f * wSq;
    const double ratio = std::sqrt(wSq / grmQuad);
    if (!(ratio > 0.0) || !std::isfinite(ratio)) return {kNaN, spa::Status::NaNoTest};
    partial_spa::Partition partition;
    if (!cachedPartition) {
        partition = partial_spa::split(weight);
        cachedPartition = &partition;
    }
    partial_spa::Workspace workspace;
    partial_spa::bind(*cachedPartition, weight, q, nullptr, nullptr, workspace);
    return binomial_cgf::twoSidedSpa(workspace.partial, mean,
        std::abs(rawScore - mean) * ratio, independentVariance, zScore,
        1e-6, 1.0, true, &workspace.full);
}

spa::Result evaluateIndividualFrequencyScore(
    const Eigen::VectorXd &weight, const Eigen::VectorXd &af,
    const Eigen::VectorXd &genoVariance, double rawScore, double spaCutoff,
    double &zScore, double &outVar
) {
    const double mean = 2.0 * weight.dot(af);
    outVar = (weight.array().square() * genoVariance.array()).sum();
    if (!(outVar > 0.0) || !std::isfinite(outVar)) {
        zScore = 0.0;
        outVar = 0.0;
        return {kNaN, spa::Status::NaNoTest};
    }
    zScore = (rawScore - mean) / std::sqrt(outVar);
    if (!std::isfinite(zScore)) return {kNaN, spa::Status::NaNoTest};
    if (std::abs(zScore) <= spaCutoff) return spa::normalBranch(zScore);
    const auto partition = partial_spa::split(weight, nullptr, &af);
    partial_spa::Workspace workspace;
    partial_spa::bind(partition, weight, 0.0, &af, nullptr, workspace);
    return binomial_cgf::twoSidedSpa(workspace.partial, mean,
        std::abs(rawScore - mean), outVar, zScore,
        1e-6, 1.0, true, &workspace.full);
}

// Resolve the trait family of one phenotype, matching nullmodel::fitAll's own
// dispatch: a Cox spec (TIME:EVENT) → Cox; otherwise the explicit
// --regression-model, or per-column inference in Auto mode.  Drives the
// Branch-B Wald refit; `None` disables it.
interaction_wald::TraitType resolveTraitType(
    const SubjectData &sd,
    const nullmodel::PhenoSpec &spec,
    nullmodel::RegressionModel regModel
) {
    using TT = interaction_wald::TraitType;
    using RM = nullmodel::RegressionModel;
    if (nullmodel::isCoxSpec(spec)) return TT::Cox;
    RM m = regModel;
    if (m == RM::Auto)
        m = nullmodel::inferModelFromColumn(sd.getColumn(spec.yColumn), spec.yColumn,
                                            sd.usedIIDs())
                .model;
    switch (m) {
    case RM::Linear:   return TT::Linear;
    case RM::Logistic: return TT::Logistic;
    case RM::Cox:      return TT::Cox;
    case RM::Ordinal:  return TT::Ordinal;
    default:           return TT::None;
    }
}

// Recode a raw phenotype column into the form the Wald fitter expects.  The two-
// sided interaction Wald p is invariant to the recode direction, so only the
// family matters: Logistic → {0,1} (below the midpoint → 0), Ordinal → {0..J−1}
// (shift by the minimum), Linear/Cox → unchanged (OLS is affine-invariant; the
// Cox response is time/event, handled separately).  NaN entries are preserved
// (the fitter drops those rows).
Eigen::VectorXd recodeResponse(const Eigen::VectorXd &col, interaction_wald::TraitType t) {
    using TT = interaction_wald::TraitType;
    if (t == TT::Logistic) {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (Eigen::Index i = 0; i < col.size(); ++i)
            if (std::isfinite(col[i])) { lo = std::min(lo, col[i]); hi = std::max(hi, col[i]); }
        const double mid = 0.5 * (lo + hi);
        Eigen::VectorXd out(col.size());
        for (Eigen::Index i = 0; i < col.size(); ++i)
            out[i] = std::isnan(col[i]) ? col[i] : (col[i] > mid ? 1.0 : 0.0);
        return out;
    }
    if (t == TT::Ordinal) {
        double lo = std::numeric_limits<double>::infinity();
        for (Eigen::Index i = 0; i < col.size(); ++i)
            if (std::isfinite(col[i])) lo = std::min(lo, col[i]);
        Eigen::VectorXd out(col.size());
        for (Eigen::Index i = 0; i < col.size(); ++i)
            out[i] = std::isnan(col[i]) ? col[i] : std::round(col[i] - lo);
        return out;
    }
    return col; // Linear (affine-invariant); Cox uses time/event, not y
}

} // namespace

// ══════════════════════════════════════════════════════════════════════
// InteractionScan
// ══════════════════════════════════════════════════════════════════════

double InteractionScan::phiQuad(const Eigen::VectorXd &x) const {
    if (m_grm)
        return m_grm->quadForm(x.data(), static_cast<uint32_t>(x.size()));
    return x.squaredNorm();
}

InteractionScan::InteractionScan(
    Eigen::VectorXd resid,
    std::vector<std::string> envNames,
    std::vector<Eigen::VectorXd> envVecs,
    std::shared_ptr<const SparseGRM> grm,
    std::shared_ptr<const interaction_wald::WaldData> wald,
    double marginalCutoff,
    double spaCutoff,
    std::shared_ptr<const AFData> afData
)
    : m_resid(std::move(resid)),
      m_envNames(std::move(envNames)),
      m_grm(std::move(grm)),
      m_wald(std::move(wald)),
      m_afData(std::move(afData)),
      m_marginalCutoff(marginalCutoff),
      m_spaCutoff(spaCutoff)
{
    if (m_afData && m_grm)
        throw std::runtime_error("Individual-AF analysis does not accept a GRM.");
    const Eigen::Index N = m_resid.size();
    if (N < 2 || !m_resid.allFinite() || m_envNames.size() != envVecs.size())
        throw std::runtime_error("SPAGxE: inconsistent or non-finite residual/environment data.");
    for (const auto &env : envVecs)
        if (env.size() != N || !env.allFinite())
            throw std::runtime_error("SPAGxE: environment must align with finite residuals.");
    m_resid2 = m_resid.array().square();
    m_residSum = m_resid.sum();
    m_R_phiQuad = phiQuad(m_resid); // RᵀΦR (= ΣR² in the base case)
    const int nEnv = static_cast<int>(m_envNames.size());
    const bool mix = static_cast<bool>(m_afData);

    if (mix) {
        const auto &X = m_afData->onePlusPCs;
        if (m_afData->nPC < 1 || X.rows() != N ||
            X.cols() != 1 + m_afData->nPC || !X.allFinite() ||
            X.rows() <= X.cols())
            throw std::runtime_error("SPAGxEmix: invalid per-individual AF design.");
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
        if (qr.rank() != X.cols() || !m_afData->XtX_inv_Xt.allFinite() ||
            m_afData->XtX_inv_Xt.rows() != X.cols() ||
            m_afData->XtX_inv_Xt.cols() != N ||
            m_afData->sqrt_XtX_inv_diag.size() != X.cols() ||
            !m_afData->sqrt_XtX_inv_diag.allFinite())
            throw std::runtime_error("SPAGxEmix: rank-deficient or inconsistent PC design.");
        m_afVec.resize(N);
        m_afGeno.resize(N);
        m_wVec.resize(N);
        m_wScratch.resize(N);
    }

    // m_W is the GEMM/GEMV operand: one product against a genotype batch yields
    // the marginal score (col 0) and every per-env score column at once.
    //   base:  m_W = [ R | w_1 | … | w_nEnv ],  w_e = (E_e − λ_e)∘R  (λ_e fixed)
    //   mix:   m_W = [ R | E_1∘R | … ],         λ_e is per-marker so only the
    //          naive weight E_e∘R is precomputable (S_GxE formed per marker).
    m_W.resize(N, 1 + nEnv);
    m_W.col(0) = m_resid;

    m_envs.reserve(nEnv);
    for (int e = 0; e < nEnv; ++e) {
        EnvData ed;
        ed.name = m_envNames[e];
        ed.E = std::move(envVecs[e]);

        if (mix) {
            // Per-individual AF: λ is variance-weighted and per-marker, so the
            // Branch-A weight cannot be precomputed.  Only the marker-independent
            // naive weight E∘R feeds the GEMV (rawScores[1+e] = Σ G_i E_i R_i).
            ed.lambda = 0.0;
            ed.wSum = ed.wSq = ed.wPhiQuad = 0.0;
            m_W.col(1 + e) = (ed.E.array() * m_resid.array()).matrix();
            m_envs.push_back(std::move(ed));
            continue;
        }

        // λ = RᵀΦR_E / RᵀΦR, R_E = R∘E.  RᵀΦR_E is the bilinear Φ-form obtained
        // by polarization (== Σ E R² in the base case; kept as the direct sum
        // there so the base numbers match the pure-Σ path exactly).
        double R_GRM_RE;
        if (m_grm) {
            const Eigen::VectorXd RE = (m_resid.array() * ed.E.array()).matrix();
            const Eigen::VectorXd RpRE = m_resid + RE;
            R_GRM_RE = 0.5 * (phiQuad(RpRE) - m_R_phiQuad - phiQuad(RE));
        } else {
            R_GRM_RE = (m_resid.array().square() * ed.E.array()).sum();
        }
        ed.lambda = (m_R_phiQuad > 0.0) ? R_GRM_RE / m_R_phiQuad : 0.0;
        ed.w = (ed.E.array() - ed.lambda) * m_resid.array();
        ed.wSum = ed.w.sum();
        ed.wSq = ed.w.squaredNorm();
        ed.partition = partial_spa::split(ed.w);
        ed.wPhiQuad = phiQuad(ed.w); // wᵀΦw (= Σw² in the base case)
        m_W.col(1 + e) = ed.w;
        m_envs.push_back(std::move(ed));
    }
}

std::unique_ptr<MethodBase> InteractionScan::clone() const {
    return std::make_unique<InteractionScan>(*this);
}

std::string InteractionScan::getHeaderColumns() const {
    // Leading marginal block LOG10P_G Z_G BETA_G SE_G SPA_STATUS_G
    // (normal-approx, so its Z is already p-consistent and needs no Z_Norm, and
    // its status is the constant 1 NORMAL wherever a statistic exists), then a
    // 9-wide G×E block per environment: the final p as the magnitude LOG10P_Gx
    // (CCT in Branch B, SPA otherwise); the Branch-B Wald magnitude
    // LOG10P_Wald_Gx (NaN when no Wald ran); the p-consistent Z and the
    // raw-score Z_Norm (which differ in the tails); BETA/SE; and the
    // saddlepoint outcome SPA_STATUS_Gx.  The linear P_G / P_Gx columns left in
    // the log-domain implementation.
    std::string h = "\tLOG10P_G\tZ_G\tBETA_G\tSE_G\tSPA_STATUS_G";
    for (const auto &n : m_envNames)
        h += "\tLOG10P_Gx" + n + "\tLOG10P_Wald_Gx" + n + "\tZ_Gx" + n +
             "\tZ_Norm_Gx" + n + "\tBETA_Gx" + n + "\tSE_Gx" + n +
             "\tSPA_STATUS_Gx" + n + "\tLOG10P_SPA_Gx" + n + "\tBRANCH_Gx" + n;
    return h;
}

namespace {

// Nine output cells for one environment. The final result combines the
// score and prospective Wald legs in the log domain when a valid raw-trait
// Wald refit is available. A missing score leg is never hidden by Wald.
// SPA_STATUS always describes the score leg, and LOG10P_SPA preserves it
// separately from the final CCT result.
void writeEnvBlock(
    std::vector<double> &out,
    int base,
    const spa::Result &ts,
    double lWald,
    bool waldEnabled,
    double z,
    double var
) {
    out[base + 6] = binomial_cgf::statusCode(ts.status);
    if (!(var > 0.0)) return;

    double negLog = ts.negLog10p;
    if (waldEnabled && std::isfinite(ts.negLog10p)) {
        // A successful Wald leg must not conceal a failed score test.
        // The CCT in the log domain, over two magnitudes.
        double Ls[2] = {ts.negLog10p, lWald};
        negLog = math::cauchyCombineLog10(Ls, 2);
    }
    out[base + 0] = negLog;                      // LOG10P_Gx<E> (final)
    out[base + 1] = lWald;                       // LOG10P_Wald_Gx<E> (NaN if none)
    // Z from LOG10P_Gx, not from P_Gx: the linear inversion
    // saturated at |Z| = 37.0470962993612 for every L >= 299.698970.  Since
    // Stage 5, `negLog` is the CCT's own magnitude, so the saturation is gone
    // from the combined branch as well rather than merely displaced onto it.
    out[base + 2] = math::zFromNegLog10P(negLog, z);  // Z_Gx<E> (p-consistent)
    out[base + 3] = z;                           // Z_Norm_Gx<E> (raw score z)
    out[base + 4] = z / std::sqrt(var);           // centered-score BETA_Gx<E>
    out[base + 7] = ts.negLog10p; // score-test p, before CCT
    out[base + 5] = 1.0 / std::sqrt(var);        // SE_Gx<E>
}

}  // namespace

void InteractionScan::evalMarker(
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    const Eigen::VectorXd &rawScores,
    std::vector<double> &out
) {
    const int nEnv = static_cast<int>(m_envNames.size());
    out.assign(static_cast<size_t>(5 + 9 * nEnv), kNaN);
    // Pre-set so that a marker leaving the marginal block empty still says why
    //: no statistic exists there, it is not a failed one.
    out[4] = binomial_cgf::statusCode(spa::Status::NaNoTest);   // SPA_STATUS_G

    const double q = altFreq;

    // ── Marginal genetic block (always the normal approximation) ────────
    // S_G = Σ (G_i - 2q) R_i; centering also supports external residuals,
    // variance 2q(1−q)·RᵀΦR.
    const double varSG = m_R_phiQuad * 2.0 * q * (1.0 - q);
    double pMarg = 1.0; // degenerate variance ⇒ take Branch A (no projection)
    if (varSG > 0.0) {
        const double sG = rawScores[0] - 2.0 * q * m_residSum;
        const double sdG = std::sqrt(varSG);
        const double zG = sG / sdG;
        pMarg = 2.0 * math::pnorm(-std::abs(zG));
        // LOG10P_G in the log domain, not −log10(pMarg): the linear tail is
        // exactly zero past |zG| = 38.6 and its logarithm would be the +Inf
        // invariant C1 forbids.  pMarg itself stays linear because its only
        // other consumer is the Branch A / Branch B routing against
        // --spagxe-marginal-cutoff, an input threshold.
        out[0] = -spa::normalTwoSidedLog(zG) / math::kLn10;   // LOG10P_G
        out[1] = zG;          // Z_G (normal ⇒ p-consistent)
        out[2] = sG / varSG;  // BETA_G
        out[3] = 1.0 / sdG;   // SE_G
        // 1 NORMAL: this block never attempts a saddlepoint, and the output contract
        // covers exactly that case.  Var(S_G) <= 0 leaves the pre-set
        // 8 NA_NO_TEST — there is no statistic to report, not a failed one.
        out[4] = binomial_cgf::statusCode(spa::Status::Normal);   // SPA_STATUS_G
    }

    const bool branchB = (pMarg <= m_marginalCutoff);

    // Branch B: genotype-adjusted residual R̃ = R − α − β·G, projecting the
    // marginal genetic effect out.  [α, β]ᵀ = ([1,G]ᵀ[1,G])⁻¹ [1,G]ᵀ R.
    // R̃ is env-independent, so it is formed once per marker.  (Unweighted even
    // for +; the GRM re-enters only through the variance — paper eq. 3.)
    Eigen::VectorXd R0;
    if (branchB) {
        const double n = static_cast<double>(m_resid.size());
        const double meanG = GVec.mean();
        const Eigen::VectorXd centeredG = GVec.array() - meanG;
        const double ssG = centeredG.squaredNorm();
        if (ssG > 0.0 && std::isfinite(ssG)) {
            const double beta = centeredG.dot(m_resid) / ssG;
            R0 = m_resid.array() - m_residSum / n - beta * centeredG.array();
        

        }
        // zero genotype variance (monomorphic G): leave R0 empty ⇒ per-env NaN below.
    }

    // Branch B fires the prospective Wald leg only on the base (non-GRM) path
    // (SPAGxE+ keeps no Wald — paper) and only when raw phenotype data is
    // present (fit mode; residual mode leaves m_wald null → SPA-only).
    const bool waldEnabled = branchB && m_wald && !m_grm &&
                             m_wald->trait != interaction_wald::TraitType::None;

    for (int e = 0; e < nEnv; ++e) {
        const int base = 5 + 9 * e;
        out[base + 8] = branchB ? 2.0 : 1.0;
        // A marker that never reaches the score test still gets a status: a
        // bare row of NA says nothing about why, and this family's NA are not
        // hypothetical (see the SPAGxEmix note in evalMarkerMix).
        out[base + 6] = binomial_cgf::statusCode(spa::Status::NaNoTest);

        double z = 0.0, var = 0.0, score;
        double lWald = kNaN;
        spa::Result ts{kNaN, spa::Status::NaNoTest};
        if (!branchB) {
            // Branch A: precomputed λ-orthogonalised weight and its Φ-form.
            const EnvData &ed = m_envs[e];
            score = rawScores[1 + e]; // Σ G_i w_{e,i}
            ts = evaluateSharedFrequencyScore(ed.w, ed.wSum, ed.wSq, ed.wPhiQuad, score, q,
                              m_spaCutoff, z, var, &ed.partition);
        } else {
            if (R0.size() == 0) continue; // degenerate G ⇒ leave NaN
            // Branch B per-marker weight u = E ∘ R̃; Ŝ_GxE = Σ G_i u_i.
            Eigen::VectorXd u;
            u = m_envs[e].E.array() * R0.array(); // paper / original R Branch B
            score = GVec.dot(u);
            ts = evaluateSharedFrequencyScore(u, u.sum(), u.squaredNorm(), phiQuad(u), score, q,
                              m_spaCutoff, z, var);
            // Prospective Wald p of the G:E coefficient in trait ~ X+E+g+g:E,
            // over the same environment used for the score weight.
            if (waldEnabled)
                lWald = interaction_wald::waldInteractionLog10P(*m_wald, GVec, m_envs[e].E);
        }
        writeEnvBlock(out, base, ts, lWald, waldEnabled, z, var);
    }
}

// ══════════════════════════════════════════════════════════════════════
// evalMarkerMix — SPAGxEmix per-individual-AF marker evaluation
// ══════════════════════════════════════════════════════════════════════
//
// Mirrors evalMarker but with the per-individual genotype law G_i ~ Bin(2, q̂_i)
// (Ma et al. 2025, SPAGxEmix_CCT_one_SNP).  q̂_i is estimated per marker via the
// SPAmix AF cascade (computeAFVec); the marginal screen is re-centred by
// 2·Σ R_i q̂_i, λ is variance-weighted and recomputed per marker, and every SPA
// uses the per-individual afVec (evaluateIndividualFrequencyScore).  Branch B reuses the exact
// same genotype-adjusted-residual projection and the Phase-3 Wald + CCT leg.
void InteractionScan::evalMarkerMix(
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    const Eigen::VectorXd &rawScores,
    std::vector<double> &out
) {
    const int nEnv = static_cast<int>(m_envNames.size());
    out.assign(static_cast<size_t>(5 + 9 * nEnv), kNaN);
    // Pre-set so that a marker leaving the marginal block empty still says why
    //: no statistic exists there, it is not a failed one.
    out[4] = binomial_cgf::statusCode(spa::Status::NaNoTest);   // SPA_STATUS_G

    const double q = altFreq;

    // ── Per-individual ALT frequency q̂_i (SPAmix cascade) ──────────────
    // [1|PCs] design + OLS matrices from m_afData; identical to SPAmix.
    AFContext ctx{
        m_afData->onePlusPCs,
        m_afData->XtX_inv_Xt,
        m_afData->sqrt_XtX_inv_diag,
        m_afData->onePlusPCs.rightCols(m_afData->nPC),
        static_cast<int>(m_resid.size()),
        m_afData->nPC
    };
    const bool flipAF = q > 0.5;
    const double minorFreq = flipAF ? 1.0 - q : q;
    bool afFitted = false;
    if (flipAF) {
        m_afGeno = 2.0 - GVec.array();
        afFitted = computeAFVec(m_afGeno, minorFreq, ctx, m_afVec, true);
    } else {
        afFitted = computeAFVec(GVec, minorFreq, ctx, m_afVec);
    }
    if (!afFitted || !m_afVec.allFinite() ||
        (m_afVec.array() < 0.0).any() || (m_afVec.array() > 1.0).any()) {
        out[4] = binomial_cgf::statusCode(spa::Status::AfFitFailed);
        for (int e = 0; e < nEnv; ++e)
            out[5 + 9 * e + 6] = binomial_cgf::statusCode(spa::Status::AfFitFailed);
        return;
    }
    m_wVec = 2.0 * m_afVec.array() * (1.0 - m_afVec.array());  // 2q̂(1−q̂) → m_wVec

    // ── Marginal genetic block (always the normal approximation) ────────
    // S_G = Σ G_i R_i, re-centred: E[S_G] = 2·Σ R_i q̂_i, Var = Σ R_i²·2q̂(1−q̂).
    // AF fitting is minor-oriented; scores and output effects remain ALT-oriented.
    const double sG = rawScores[0];
    const double meanSG = 2.0 * m_resid.dot(m_afVec);
    const double varSG = m_resid2.dot(m_wVec);
    double pMarg = 1.0; // degenerate variance ⇒ take Branch A (no projection)
    if (varSG > 0.0) {
        const double sdG = std::sqrt(varSG);
        const double zG = (sG - meanSG) / sdG;
        pMarg = 2.0 * math::pnorm(-std::abs(zG));
        // The magnitude in the log domain; see the base path above.
        out[0] = -spa::normalTwoSidedLog(zG) / math::kLn10;   // LOG10P_G
        out[1] = zG;                    // Z_G (normal ⇒ p-consistent)
        out[2] = (sG - meanSG) / varSG; // BETA_G
        out[3] = 1.0 / sdG;             // SE_G
        out[4] = binomial_cgf::statusCode(spa::Status::Normal);   // SPA_STATUS_G
    }

    const bool branchB = (pMarg <= m_marginalCutoff);

    // Branch B: genotype-adjusted residual R̃ = R − α − β·G (env-independent).
    Eigen::VectorXd R0;
    if (branchB) {
        const double n = static_cast<double>(m_resid.size());
        const double meanG = GVec.mean();
        const Eigen::VectorXd centeredG = GVec.array() - meanG;
        const double ssG = centeredG.squaredNorm();
        if (ssG > 0.0 && std::isfinite(ssG)) {
            const double beta = centeredG.dot(m_resid) / ssG;
            R0 = m_resid.array() - m_residSum / n - beta * centeredG.array();
        

        }
    }

    // Raw-trait fit mode provides the prospective Wald leg.
    const bool waldEnabled = branchB && m_wald &&
                             m_wald->trait != interaction_wald::TraitType::None;

    for (int e = 0; e < nEnv; ++e) {
        const int base = 5 + 9 * e;
        out[base + 8] = branchB ? 2.0 : 1.0;
        out[base + 6] = binomial_cgf::statusCode(spa::Status::NaNoTest);

        const EnvData &ed = m_envs[e];
        double z = 0.0, var = 0.0, score;
        double lWald = kNaN;
        spa::Result ts{kNaN, spa::Status::NaNoTest};
        if (!branchB) {
            // Branch A: variance-weighted λ (per marker) on the per-individual q̂.
            //   λ = Σ 2q̂(1−q̂) E R² / Σ 2q̂(1−q̂) R²
            const double denom = (m_wVec.array() * m_resid2.array()).sum();
            const double numer = (m_wVec.array() * ed.E.array() * m_resid2.array()).sum();
            const double lambda = (denom > 0.0) ? numer / denom : 0.0;
            // S_GxE = S2 − λ·S_G  (S2 = Σ G_i E_i R_i from the GEMV col 1+e).
            score = rawScores[1 + e] - lambda * sG;
            m_wScratch = (ed.E.array() - lambda) * m_resid.array(); // weight (E−λ)∘R
            ts = evaluateIndividualFrequencyScore(m_wScratch, m_afVec, m_wVec, score,
                                 m_spaCutoff, z, var);
        } else {
            if (R0.size() == 0) continue; // degenerate G ⇒ leave NaN
            // Branch B per-marker weight u = E ∘ R̃; Ŝ_GxE = Σ G_i u_i.
            Eigen::VectorXd u;
            u = ed.E.array() * R0.array(); // paper / original R Branch B
            score = GVec.dot(u);
            ts = evaluateIndividualFrequencyScore(u, m_afVec, m_wVec, score, m_spaCutoff,
                                 z, var);
            if (waldEnabled)
                lWald = interaction_wald::waldInteractionLog10P(*m_wald, GVec, ed.E);
        }
        writeEnvBlock(out, base, ts, lWald, waldEnabled, z, var);
    }
}

void InteractionScan::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int /*markerInChunkIdx*/,
    std::vector<double> &result
) {
    const Eigen::VectorXd rawScores = m_W.transpose() * GVec;
    if (m_afData)
        evalMarkerMix(GVec, altFreq, rawScores, result);
    else
        evalMarker(GVec, altFreq, rawScores, result);
}

void InteractionScan::getResultBatch(
    const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
    const std::vector<double> &altFreqs,
    const std::vector<int> & /*chunkIdxs*/,
    std::vector<std::vector<double> > &results
) {
    const int B = static_cast<int>(GBatch.cols());
    results.resize(B);
    // One GEMM supplies the marginal (col 0) and every per-env score column.
    //   base:  scores(b, 1+e) = Σ G_b·w_e         (Branch-A λ-orthogonalised)
    //   mix:   scores(b, 1+e) = Σ G_b·(E_e∘R)     (naive S2; λ applied per marker)
    const Eigen::MatrixXd scores = GBatch.transpose() * m_W; // B × (1+nEnv)
    const bool mix = static_cast<bool>(m_afData);
    for (int b = 0; b < B; ++b) {
        const Eigen::VectorXd rs = scores.row(b).transpose();
        if (mix)
            evalMarkerMix(GBatch.col(b), altFreqs[b], rs, results[b]);
        else
            evalMarker(GBatch.col(b), altFreqs[b], rs, results[b]);
    }
}

// ══════════════════════════════════════════════════════════════════════
// executeInteractionWorkflow — full workflow (shares the scan driver); serves both
//   --method spagxe    (base / SPAGxE+; pcColNames empty)
//   --method spagxemix (per-individual AF; pcColNames = the --pc-cols columns)
// ══════════════════════════════════════════════════════════════════════

namespace {

void executeInteractionWorkflow(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::vector<std::string> &pcColNames, // empty → base/+, non-empty → mix
    const std::string &grmTableFile,
    const std::string &spgrmGctaFile,
    const GenotypeInput &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double marginalCutoff,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    bool saveResid
) {
    const bool mix = !pcColNames.empty();
    const bool hasGrm = !grmTableFile.empty() || !spgrmGctaFile.empty();
    if (mix && hasGrm)
        throw std::runtime_error("Individual-AF analysis does not accept a GRM.");
    // Display label (messages) and output-file suffix.  Base and SPAGxE+ share
    // the "SPAGxE" file suffix; the "+" is display only.  Mix writes ".SPAGxEmix".
    const std::string dispLabel = mix ? "SPAGxEmix" : (hasGrm ? "SPAGxE+" : "SPAGxE");
    const char *fileSuffix = mix ? "SPAGxEmix" : "SPAGxE";
    if (hasGrm) infoMsg("GRM SPA adjustment: score variance-ratio rescaling with a guarded partial-normal independence CGF");
    const int nPC = static_cast<int>(pcColNames.size());

    const bool fitPath = !phenoNameSpec.empty();
    nullmodel::RegressionModel regModel{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        regModel = nullmodel::parseRegressionModel(regressionModelStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
        infoMsg("%s: fitting %s null model for %zu phenotype(s)",
                dispLabel.c_str(), nullmodel::regressionModelName(regModel),
                phenoSpecs.size());
    }

    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = readGenotypeSampleIds(geno);
    SubjectData sd(std::move(famIIDs));
    if (fitPath) {
        sd.loadPhenoFile(phenoFile, nullmodel::columnsNeeded(phenoSpecs));
    } else {
        sd.loadResidOne(phenoFile, residNames);
        if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile);
    }
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    // (+) variant: restrict to subjects present in the sparse GRM (the GRM
    // supplies the retrospective genotype covariance among relatives).
    if (hasGrm)
        sd.setGrmSubjects(
            SparseGRM::parseSubjectIDs(grmTableFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();

    // Environment columns must enter the null model in fit mode (so R ⟂ E),
    // exactly as SAGELD requires --envir-name ⊆ --covar-name.
    if (fitPath) {
        for (const auto &en : envNames)
            if (std::find(covarNames.begin(), covarNames.end(), en) == covarNames.end())
                throw std::runtime_error(
                    dispLabel + ": environment '" + en +
                    "' must also appear in --covar-name (it enters the "
                    "genotype-independent null model).");
        // SPAGxEmix: the PC columns feeding the per-individual AF model must
        // also adjust the null model (matching SPAGxEmix_CCT's topPCs ⊆ Cova).
        if (mix)
            for (const auto &pc : pcColNames)
                if (std::find(covarNames.begin(), covarNames.end(), pc) == covarNames.end())
                    throw std::runtime_error(
                        dispLabel + ": PC column '" + pc +
                        "' must also appear in --covar-name (the per-individual "
                        "allele-frequency PCs adjust the null model).");
    }
    // Drop subjects with a missing environment value (a no-op in fit mode,
    // where E is a covariate; necessary in residual mode, where E is only
    // used to form λ and is not otherwise NA-filtered).
    sd.dropNaInColumns(envNames);
    if (mix) sd.dropNaInColumns(pcColNames); // AF design must be complete
    const uint32_t N = sd.nUsed();
    infoMsg("  %u subjects in union mask", N);

    // Covariate union (PCs + env main effects).  In fit mode it is both the
    // null-model design and — reused per marker — the Branch-B Wald design; the
    // per-phenotype trait family (below) selects the Wald fitter.
    Eigen::MatrixXd covarUnion;
    std::vector<interaction_wald::TraitType> traitTypes;
    if (fitPath) {
        if (!covarNames.empty()) covarUnion = sd.getColumns(covarNames);
        else covarUnion.resize(N, 0);
        traitTypes.reserve(phenoSpecs.size());
        for (const auto &spec : phenoSpecs)
            traitTypes.push_back(resolveTraitType(sd, spec, regModel));
    }

    // SPAGxEmix: union-level [1 | PCs] design for the per-individual AF cascade.
    Eigen::MatrixXd unionOnePlusPCs;
    if (mix) {
        unionOnePlusPCs.resize(N, 1 + nPC);
        unionOnePlusPCs.col(0).setOnes();
        unionOnePlusPCs.rightCols(nPC) = sd.getColumns(pcColNames);
        infoMsg("SPAGxEmix: per-individual AF from %d PC(s)", nPC);
    }

    if (fitPath) {
        nullmodel::EngineOptions eo;
        eo.nthreads = nthreads;
        auto fits = nullmodel::fitAll(sd, phenoSpecs, regModel, covarUnion, eo);
        std::vector<Eigen::VectorXd> rs;
        std::vector<std::string> ns;
        rs.reserve(fits.size());
        ns.reserve(fits.size());
        for (auto &f : fits) {
            infoMsg("  Fitted '%s': %d subjects after NaN removal", f.name.c_str(),
                    f.nUsedRows);
            rs.push_back(std::move(f.residuals));
            ns.push_back(f.name);
        }
        if (saveResid) {
            std::vector<nullmodel::NullModelFit> dumpFits(rs.size());
            for (size_t i = 0; i < rs.size(); ++i) {
                dumpFits[i].name = ns[i];
                dumpFits[i].residuals = rs[i];
                dumpFits[i].nUsedRows = static_cast<int>(rs[i].size());
            }
            nullmodel::writeResidualsFile(outPrefix + ".null.resid", sd, dumpFits,
                                          compression, compressionLevel);
        }
        sd.setResidualsFromFit(std::move(rs), std::move(ns));
    }

    // Environment vectors in union-dense order (aligned to the residuals).
    std::vector<Eigen::VectorXd> envUnion;
    envUnion.reserve(envNames.size());
    for (const auto &en : envNames)
        envUnion.push_back(sd.getColumn(en));

    // (+) variant: load the union-level sparse GRM once; re-index per phenotype.
    std::shared_ptr<const SparseGRM> unionGrm;
    if (hasGrm) {
        auto subjIDs = sd.usedIIDs();
        unionGrm = std::make_shared<const SparseGRM>(
            SparseGRM::load(grmTableFile, spgrmGctaFile, subjIDs, sd.famIIDs()));
        infoMsg("Sparse GRM: %zu entries (SPAGxE+ variance correction)",
                unionGrm->nnz());
    }

    auto genoData =
        openGenotypeCatalog(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);
    infoMsg("Genotype data: %u markers, %u subjects", genoData->nMarkers(),
            genoData->nSubjUsed());

    auto phenoInfos = sd.buildPerColumnMasks();
    const int K = sd.residOneCols();
    if (K > 1) infoMsg("Multi-column residual file: %d phenotypes", K);

    std::vector<PhenoTask> tasks(K);
    for (int rc = 0; rc < K; ++rc) {
        const auto &pi = phenoInfos[rc];
        const Eigen::VectorXd residualUnion = K == 1 ? sd.residuals() : Eigen::VectorXd(sd.residMatrix().col(rc));
        Eigen::VectorXd phenoResid = extractPhenoVec(residualUnion, pi);
        std::vector<Eigen::VectorXd> envVecs;
        envVecs.reserve(envUnion.size());
        for (const auto &eu : envUnion)
            envVecs.push_back((pi.nUsed != N) ? extractPhenoVec(eu, pi) : eu);

        // Per-phenotype GRM: for a single phenotype the union GRM is used
        // directly; for K > 1 the union entries are re-indexed to the
        // phenotype's dense subject order (in the scan driver).
        std::shared_ptr<const SparseGRM> phenoGrm;
        if (hasGrm) {
            if (pi.nUsed == N) {
                phenoGrm = unionGrm;
            } else {
                const auto &u2l = pi.unionToLocal;
                std::vector<SparseGRM::Entry> pEntries;
                for (const auto &en : unionGrm->entries()) {
                    const uint32_t li = u2l[en.row], lj = u2l[en.col];
                    if (li != UINT32_MAX && lj != UINT32_MAX)
                        pEntries.push_back({li, lj, en.value});
                }
                phenoGrm = std::make_shared<const SparseGRM>(
                    SparseGRM::fromEntries(pi.nUsed, std::move(pEntries)));
            }
        }

        // Per-phenotype Branch-B Wald data: only in fit mode on the base (non-
        // GRM) path — SPAGxE+ keeps no Wald, residual mode has no raw phenotype.
        // Extract the raw phenotype / covariate design into the same dense
        // subject order as phenoResid (identical (K>1)?extract:union pattern),
        // so the Wald design aligns with the engine's per-phenotype genotype.
        std::shared_ptr<const interaction_wald::WaldData> phenoWald;
        if (fitPath && !hasGrm) {
            const interaction_wald::TraitType tt = traitTypes[rc];
            if (tt != interaction_wald::TraitType::None) {
                auto wd = std::make_shared<interaction_wald::WaldData>();
                wd->trait = tt;
                wd->covar = (pi.nUsed != N) ? extractPhenoMat(covarUnion, pi) : covarUnion;
                const auto &spec = phenoSpecs[rc];
                if (tt == interaction_wald::TraitType::Cox) {
                    Eigen::VectorXd tU = sd.getColumn(spec.timeColumn);
                    Eigen::VectorXd eU = sd.getColumn(spec.eventColumn);
                    wd->time = (pi.nUsed != N) ? extractPhenoVec(tU, pi) : tU;
                    wd->event = (pi.nUsed != N) ? extractPhenoVec(eU, pi) : eU;
                } else {
                    Eigen::VectorXd yU = recodeResponse(sd.getColumn(spec.yColumn), tt);
                    wd->y = (pi.nUsed != N) ? extractPhenoVec(yU, pi) : yU;
                }
                phenoWald = std::move(wd);
            }
        }

        // SPAGxEmix: per-phenotype [1|PCs] design + OLS matrices for computeAFVec.
        // Extracted into the phenotype's dense subject order (as covarUnion is),
        // so the AF cascade aligns with the residual / environment / genotype.
        std::shared_ptr<const AFData> phenoAF;
        if (mix) {
            auto af = std::make_shared<AFData>();
            af->nPC = nPC;
            af->onePlusPCs =
                (pi.nUsed != N) ? extractPhenoMat(unionOnePlusPCs, pi) : unionOnePlusPCs;
            const Eigen::MatrixXd &X = af->onePlusPCs;
            prepareAFLinearDesign(X, af->XtX_inv_Xt, af->sqrt_XtX_inv_diag);
            phenoAF = std::move(af);
        }

        tasks[rc].phenoName = pi.name;
        tasks[rc].method = std::make_unique<InteractionScan>(
            std::move(phenoResid), envNames, std::move(envVecs), phenoGrm,
            std::move(phenoWald), marginalCutoff, spaCutoff, std::move(phenoAF));
        tasks[rc].unionToLocal = pi.unionToLocal;
        tasks[rc].nUsed = pi.nUsed;
        infoMsg("  Phenotype '%s': %u subjects", pi.name.c_str(), pi.nUsed);
    }

    infoMsg("Running %s marker tests (%d thread(s), %d phenotype(s))...",
            dispLabel.c_str(), nthreads, K);
    multiPhenoEngine(*genoData, tasks, outPrefix, fileSuffix, compression,
                     compressionLevel, nthreads, missingCutoff, minMafCutoff,
                     minMacCutoff, hweCutoff);
    TextWriter receipt(outPrefix + ".analysis.settings.tsv");
    receipt.write("key\tvalue\nsoftware\tSPAGxE\nversion\t" SPAGXE_VERSION "\ncgf\tguarded_partial_normal\ncentral_IQR_multiplier\t1.5\npartial_guard_tolerance\t0.001\npartial_fallback\tboth_tails_full_binomial\ngrm_adjustment\tscore\nowner_contact\tyuzhuoma@amss.ac.cn\n");
    receipt.close();
}

} // namespace

// ══════════════════════════════════════════════════════════════════════
// scanHomogeneousInteractions / scanAdmixedInteractions — thin public entry points over executeInteractionWorkflow
// ══════════════════════════════════════════════════════════════════════

void scanHomogeneousInteractions(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::string &grmTableFile,
    const std::string &spgrmGctaFile,
    const GenotypeInput &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double marginalCutoff,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    bool saveResid
) {
    executeInteractionWorkflow(phenoFile, residNames, envNames, /*pcColNames=*/{}, grmTableFile,
                  spgrmGctaFile, geno, outPrefix, compression, compressionLevel,
                  marginalCutoff, spaCutoff, nthreads, nSnpPerChunk,
                  missingCutoff, minMafCutoff, minMacCutoff, hweCutoff, keepFile,
                  removeFile, regressionModelStr, phenoNameSpec, covarFile, covarNames,
                  saveResid);
}

void scanAdmixedInteractions(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::vector<std::string> &pcColNames,
    const GenotypeInput &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double marginalCutoff,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    bool saveResid
) {
    executeInteractionWorkflow(phenoFile, residNames, envNames, pcColNames, {},
                  {}, geno, outPrefix, compression, compressionLevel,
                  marginalCutoff, spaCutoff, nthreads, nSnpPerChunk,
                  missingCutoff, minMafCutoff, minMacCutoff, hweCutoff, keepFile,
                  removeFile, regressionModelStr, phenoNameSpec, covarFile, covarNames,
                  saveResid);
}



