// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE interaction scans for shared or individual allele frequencies.
// Published method: Ma, Zhao, Zhang and Bi, Nature Communications (2025).
// Maintainer: Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
// Incorporated source and license scopes: THIRD_PARTY_NOTICES.md.
//
// A genotype-independent null fit supplies R. The marginal genetic screen
// selects Branch A: (E-lambda)*R, or Branch B: E*R0, where R0 removes intercept
// and G from R. The SPA uses guarded central-Gaussian/binomial acceleration.
//
// SPAGxE+ changes score variance through a sparse GRM and remains SPA-only.
// SPAGxEmix fits individual allele frequencies and does not accept a GRM.
// Eligible unrelated raw-trait Branch B combines SPA and prospective Wald.
// Results use -log10(P), with explicit branch and computation-status fields.
#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "engine/marker.hpp"          // MethodBase
#include "genotype/genotype_source.hpp" // GenotypeInput
#include "interaction/interaction_wald.hpp"     // interaction_wald::WaldData (Branch-B Wald leg)
#include "util/partial_spa.hpp"

class SparseGRM; // held by shared_ptr; full type only needed in the .cpp

// ══════════════════════════════════════════════════════════════════════
// SPAGxEmix — per-individual allele frequency (admixture)
// ══════════════════════════════════════════════════════════════════════
//
// SPAGxEmix (Ma et al. 2025, `SPAGxEmix_CCT`) replaces the single scalar MAF q
// of the base test with a per-individual ALT frequency q̂_i estimated from the
// SNP-derived PCs. The shared estimator fits in the minor direction,
// maps fitted AF back to ALT, and follows the published pooled/OLS/logistic
// cascade (see src/allele_model/frequency_model.hpp and docs/AF_MODEL_FORMAT.md).  The retrospective genotype law becomes
// G_i ~ Bin(2, q̂_i), so the marginal screen is re-centred by 2·Σ R_i q̂_i, the
// λ-orthogonalisation weight is variance-weighted and recomputed PER MARKER
// (λ = Σ 2q̂_i(1−q̂_i) E_i R_i² / Σ 2q̂_i(1−q̂_i) R_i²), and the SPA uses the
// per-individual afVec rather than the scalar-q path. This is the unrelated
// individual-AF method; the separate SPAGxE+ method uses a shared scalar AF.
//
// AFData holds the [1|PCs] design and the OLS matrices needed by computeAFVec
// (indiv_af.hpp).  One per phenotype, shared read-only across worker clones via
// shared_ptr<const>; null → the base / SPAGxE+ scalar-q path.
struct AFData {
    Eigen::MatrixXd onePlusPCs;        // N × (1+nPC) = [1 | PCs]
    Eigen::MatrixXd XtX_inv_Xt;        // (1+nPC) × N = (XᵀX)⁻¹Xᵀ  (OLS solve)
    Eigen::VectorXd sqrt_XtX_inv_diag; // (1+nPC) = sqrt(diag((XᵀX)⁻¹))  (PC SE)
    int nPC = 0;
};

// ══════════════════════════════════════════════════════════════════════
// InteractionScan — MethodBase adapter (one clone per worker thread)
// ══════════════════════════════════════════════════════════════════════
//
// Output schema, with one block per environment:
//   col 0..4        :  LOG10P_G  Z_G  BETA_G  SE_G  SPA_STATUS_G
//                                                    marginal genetic block
//                                                    (always normal-approx)
//   col 5 + 9·e+0   :  LOG10P_Gx<Ee>  final G×E p-value, as the magnitude
//                                     −log10(P).  Branch A (and the SPAGxE+
//                                     GRM / residual-mode paths): the
//                                     SPA/normal score p.  Branch B, base path:
//                                     CCT(p_spa, p_wald).
//   col 5 + 9·e+1   :  LOG10P_Wald_Gx<Ee>  magnitude of the Branch-B Wald p of
//                                     the G:E coefficient (NaN in Branch A, the
//                                     GRM path, and residual mode — no Wald there).
//   col 5 + 9·e+2   :  Z_Gx<Ee>       p-consistent z: 2·pnorm(−|Z|)==10^(−LOG10P)
//   col 5 + 9·e+3   :  Z_Norm_Gx<Ee>  raw score z (Score/√Var, SPA-uncalibrated)
//   col 5 + 9·e+4   :  BETA_Gx<Ee>    score-derived interaction effect
//   col 5 + 9·e+5   :  SE_Gx<Ee>      score-derived interaction SE
//   col 5 + 9·e+6   :  SPA_STATUS_Gx<Ee>  saddlepoint outcome, as
//                                     static_cast<uint8_t>(spa::Status)
// Additional cells +7/+8 are LOG10P_SPA_Gx (before CCT) / BRANCH_Gx (1 or 2).
// resultSize() = 5 + 9·nEnv.
//
// LOG10P_Gx is a magnitude on every path, never −log10 of a linear p.  Where
// the reported p IS the saddlepoint p it comes from the log-domain tail
// assembly, so it survives past the point at which the
// linear-scale sum of the two tails underflows to zero.  Where Branch B adds a
// Wald leg and the reported p is CCT(p_spa, p_wald), the combination is taken
// by math::cauchyCombineLog10 over the two magnitudes,
// so the saddlepoint leg enters exactly and the Cauchy statistic — whose terms
// are 1/(π p) and which therefore overflows for p ≤ 1e-308 — is never formed.
// The linear P_Gx column was deleted in Stage 8; a consumer that
// needs it takes 10^(−LOG10P_Gx).  The Wald leg arrives as a magnitude too
// in log space (wald::lastCoef*Log10P), so no leg of the combination has an
// underflow ceiling any more.
//
// The marginal block carries SPA_STATUS_G and, LOG10P_G.  It is
// always the normal approximation, never a saddlepoint, so the status is the
// constant `1 NORMAL` wherever the block exists at all and `8 NA_NO_TEST`
// where Var(S_G) ≤ 0 — which is exactly the distinction a bare row of NA does
// not make, and the reason the output contract gives every reported test a status.
// LOG10P_G is formed by spa::normalTwoSidedLog rather than as −log10 of the
// linear tail, so it keeps rising past |Z_G| = 38.6 where the linear p flushes
// to exactly zero and its logarithm would be the +Inf invariant C1 forbids.
class InteractionScan : public MethodBase {
  public:
    // resid    — the null-model residual R (per-phenotype dense, length N).
    // envNames — environment column names, one 9-wide output block each.
    // envVecs  — the environment vectors E_e (aligned to resid), one per env.
    // grm      — per-phenotype sparse GRM for the (+) variant; null → base.
    // wald     — per-phenotype Branch-B Wald data (raw phenotype + covariates).
    //            null, or trait==None, or a non-null grm ⇒ no Wald leg (the
    //            SPAGxE+ GRM path and residual mode stay SPA-only).
    // marginalCutoff — ε, the Branch A / Branch B routing threshold.
    // spaCutoff      — r (--spa-z-threshold), the |z| normal↔SPA switch.
    // afData   — per-individual-AF design (SPAGxEmix).  Non-null ⇒ mix mode:
    //            q̂_i is estimated per marker via computeAFVec, the marginal
    //            screen / λ / SPA all use afVec. A GRM cannot be supplied with
    //            afData. Null ⇒ the base / SPAGxE+ scalar-q path.
    InteractionScan(
        Eigen::VectorXd resid,
        std::vector<std::string> envNames,
        std::vector<Eigen::VectorXd> envVecs,
        std::shared_ptr<const SparseGRM> grm,
        std::shared_ptr<const interaction_wald::WaldData> wald,
        double marginalCutoff,
        double spaCutoff,
        std::shared_ptr<const AFData> afData = nullptr
    );

    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 5 + 9 * static_cast<int>(m_envNames.size());
    }

    std::string getHeaderColumns() const override;

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result
    ) override;

    // Batched: one GEMM  scores = GBatchᵀ · [R | w_1 | … | w_nEnv]  supplies the
    // marginal and every Branch-A G×E score; per-marker Branch B (rare) uses the
    // genotype column directly.
    void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        const std::vector<int> &chunkIdxs,
        std::vector<std::vector<double> > &results
    ) override;

  private:
    // Retrospective quadratic form  xᵀΦx: SparseGRM::quadForm when a GRM is
    // present, Σx² otherwise (Φ = identity → base path).
    double phiQuad(const Eigen::VectorXd &x) const;

    // Per-environment precomputed Branch-A machinery.  In the base / SPAGxE+
    // scalar-q path λ and the weight w = (E−λ)∘R are marker-independent and
    // precomputed here.  In mix mode λ is per-marker (variance-weighted on the
    // per-individual q̂_i), so only E and the naive weight E∘R (the m_W column)
    // are meaningful; lambda/w/wSum/wSq/wPhiQuad are left unset.
    struct EnvData {
        std::string name;
        Eigen::VectorXd E;         // environment vector (Branch B needs it raw)
        double lambda;             // λ = RᵀΦR_E / RᵀΦR       (base only)
        Eigen::VectorXd w;         // Branch-A weight (E − λ)∘R (base only)
        double wSum;               // Σ w  (≈ 0 by residual orthogonality; base)
        double wSq;                // Σ w²         (independence quadratic form; base)
        double wPhiQuad;           // wᵀΦw         (GRM quadratic form; == wSq base)
        partial_spa::Partition partition; // cached shared-AF Branch-A split
    };

    // Evaluate one marker given the pre-multiplied raw scores.  The m_W columns
    // differ by mode, so the raw-score interpretation follows suit:
    //   base:  rawScores[0]=Σ G_i R_i,  rawScores[1+e]=Σ G_i w_{e,i}  (Branch-A score)
    //   mix:   rawScores[0]=Σ G_i R_i,  rawScores[1+e]=Σ G_i E_{e,i} R_i (naive S2;
    //          the Branch-A score S_GxE = S2 − λ_e·S_G is formed per marker).
    void evalMarker(
        const Eigen::Ref<const Eigen::VectorXd> &GVec,
        double altFreq,
        const Eigen::VectorXd &rawScores,
        std::vector<double> &out
    );

    // Per-individual-AF marker evaluation (SPAGxEmix; m_afData != null).
    void evalMarkerMix(
        const Eigen::Ref<const Eigen::VectorXd> &GVec,
        double altFreq,
        const Eigen::VectorXd &rawScores,
        std::vector<double> &out
    );

    Eigen::VectorXd m_resid;                    // R (null-model residual)
    Eigen::VectorXd m_resid2;                   // R∘R (marginal / λ variance factor)
    double m_residSum;                          // Σ R
    double m_R_phiQuad;                         // RᵀΦR (marginal variance factor; base)
    Eigen::MatrixXd m_W;                        // N × (1+nEnv): base [R|w_e], mix [R|E_e∘R]
    std::vector<std::string> m_envNames;        // env column names (header order)
    std::vector<EnvData> m_envs;                // per-env Branch-A state
    std::shared_ptr<const SparseGRM> m_grm;     // (+) variant; null → base / mix
    std::shared_ptr<const interaction_wald::WaldData> m_wald; // Branch-B Wald; null → SPA-only
    std::shared_ptr<const AFData> m_afData;     // SPAGxEmix per-individual AF; null → scalar-q
    double m_marginalCutoff;                    // ε
    double m_spaCutoff;                         // r

    // Per-thread scratch for the mix path (freshly sized per clone; reused
    // across markers).  m_afVec = q̂_i, m_wVec = 2 q̂_i(1−q̂_i).
    Eigen::VectorXd m_afVec;
    Eigen::VectorXd m_afGeno; // minor-allele genotype for stable AF estimation
    Eigen::VectorXd m_wVec;
    Eigen::VectorXd m_wScratch; // per-marker Branch-A weight (E_e − λ_e)∘R
};

// ══════════════════════════════════════════════════════════════════════
// scanHomogeneousInteractions — full workflow entry point
// ══════════════════════════════════════════════════════════════════════
//
// Mirrors the shared scan workflow: load SubjectData, fit or load the
// genotype-independent residual, extract each --envir-name column, build one
// PhenoTask per phenotype, and drive multiPhenoEngine.  When a sparse GRM is
// supplied (--sp-grm-*, optional) the (+) variance correction engages; absent,
// the base unrelated path runs.
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
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    const std::string &regressionModelStr = {},
    const std::string &phenoNameSpec = {},
    const std::string &covarFile = {},
    const std::vector<std::string> &covarNames = {},
    bool saveResid = false
);

// ══════════════════════════════════════════════════════════════════════
// scanAdmixedInteractions — SPAGxEmix workflow entry point (per-individual AF)
// ══════════════════════════════════════════════════════════════════════
//
// A separate --method spagxemix.  Identical to scanHomogeneousInteractions except that it builds a
// per-phenotype [1|PCs] AF design (from the --pc-cols columns, which must also
// appear in --covar-name so they adjust the null model) and threads it into each
// InteractionScan, switching the marginal screen / λ / SPA to the per-individual
// q̂_i (SPAGxEmix_CCT). This entry point does not accept a sparse GRM.
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
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    const std::string &regressionModelStr = {},
    const std::string &phenoNameSpec = {},
    const std::string &covarFile = {},
    const std::vector<std::string> &covarNames = {},
    bool saveResid = false
);
