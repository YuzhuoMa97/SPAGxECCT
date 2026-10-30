// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE local ancestry: H-weighted scores and binomial-hapcount SPA.
// Copyright (c) 2026 Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
#include "local/local_kernel.hpp"
#include "util/probability_math.hpp"
#include "util/spa.hpp"
#include "util/spa_cgf.hpp"
#include "util/binomial_support.hpp"
#include "util/partial_spa.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace spagxe_local {
namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();
double statusCode(spa::Status s) { return static_cast<double>(static_cast<uint8_t>(s)); }

// Preserve the raw-trait design and append only local-ancestry columns that
// increase its rank. The intercept handles sum_k H_ik = 2, so at most K-1
// ancestry columns survive; constant/repeated H columns do not poison Wald.
interaction_wald::WaldData ancestryAdjustedWald(
    const interaction_wald::WaldData &raw, const Eigen::MatrixXd &allH) {
    interaction_wald::WaldData wd = raw;
    const Eigen::Index n = allH.rows();
    if (raw.covar.rows() != n)
        throw std::invalid_argument("Local SPAGxE: Wald covariates do not align with G/H.");
    Eigen::MatrixXd design(n, 1 + raw.covar.cols());
    design.col(0).setOnes();
    if (raw.covar.cols()) design.rightCols(raw.covar.cols()) = raw.covar;
    auto rank = [](const Eigen::MatrixXd &a) {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a);
        qr.setThreshold(1e-10);
        return qr.rank();
    };
    Eigen::Index previousRank = rank(design);
    for (Eigen::Index k = 0; k < allH.cols(); ++k) {
        Eigen::MatrixXd candidate(n, design.cols() + 1);
        candidate.leftCols(design.cols()) = design;
        candidate.col(design.cols()) = allH.col(k);
        const Eigen::Index candidateRank = rank(candidate);
        if (candidateRank > previousRank) {
            design.swap(candidate);
            previousRank = candidateRank;
        }
    }
    wd.covar = design.rightCols(design.cols() - 1);
    return wd;
}

spa::Result localScore(const Eigen::VectorXd &g, const Eigen::VectorXd &h,
                       const Eigen::VectorXd &weight, double q,
                       double spaCutoff,
                       double &z, double &variance) {
    const double mean = q * h.dot(weight);
    variance = q * (1.0 - q) * h.dot(weight.array().square().matrix());
    if (!(variance > 0.0) || !std::isfinite(variance)) {
        z = kNaN;
        return {kNaN, spa::Status::NaNoTest};
    }
    const double centered = g.dot(weight) - mean;
    z = centered / std::sqrt(variance);
    if (!std::isfinite(z)) return {kNaN, spa::Status::NaNoTest};
    if (std::abs(z) <= spaCutoff) return spa::normalBranch(z);

    const auto partition = partial_spa::split(weight, &h);
    partial_spa::Workspace workspace;
    partial_spa::bind(partition, weight, q, nullptr, &h, workspace);
    return binomial_cgf::twoSidedSpa(workspace.partial, mean,
        std::abs(centered), variance, z, 1e-6, 1.0, true, &workspace.full);
}
} // namespace

LocalKernel::LocalKernel(Eigen::VectorXd resid,
                         std::vector<std::string> envNames,
                         std::vector<Eigen::VectorXd> envVecs,
                         std::shared_ptr<const interaction_wald::WaldData> wald,
                         double marginalCutoff, double spaCutoff)
    : m_resid(std::move(resid)), m_envNames(std::move(envNames)),
      m_envs(std::move(envVecs)), m_wald(std::move(wald)),
      m_marginalCutoff(marginalCutoff), m_spaCutoff(spaCutoff) {
    if (m_resid.size() < 2 || !m_resid.allFinite() || m_envNames.empty() ||
        m_envNames.size() != m_envs.size())
        throw std::invalid_argument("Local SPAGxE: invalid residual/environment dimensions.");
    for (const auto &e : m_envs)
        if (e.size() != m_resid.size() || !e.allFinite())
            throw std::invalid_argument("Local SPAGxE: environments must be finite and aligned.");
    if (!std::isfinite(marginalCutoff) || marginalCutoff < 0 || marginalCutoff > 1 ||
        !std::isfinite(spaCutoff) || spaCutoff < 0)
        throw std::invalid_argument("Local SPAGxE: invalid test thresholds.");
}
int LocalKernel::resultSize() const { return 5 + 9*static_cast<int>(m_envs.size()); }
std::string LocalKernel::getHeaderColumns() const {
    std::string h = "\tLOG10P_G\tZ_G\tBETA_G\tSE_G\tSPA_STATUS_G";
    for (const auto &e : m_envNames) {
        h += "\tLOG10P_Gx"+e+"\tLOG10P_Wald_Gx"+e+"\tZ_Gx"+e+
             "\tZ_Norm_Gx"+e+"\tBETA_Gx"+e+"\tSE_Gx"+e+
             "\tSPA_STATUS_Gx"+e+"\tLOG10P_SPA_Gx"+e+"\tBRANCH_Gx"+e;
    }
    return h;
}
LocalResult LocalKernel::evaluate(const Eigen::VectorXd &G, const Eigen::VectorXd &H,
                                  const Eigen::MatrixXd &allH) const {
    const Eigen::Index n = m_resid.size();
    if (G.size()!=n || H.size()!=n || allH.rows()!=n || !H.allFinite() || !allH.allFinite())
        throw std::invalid_argument("Local SPAGxE: finite G/H ancestry counts must align with residuals.");
    LocalResult out;
    out.values.assign(static_cast<size_t>(resultSize()), kNaN);
    out.altFreq=kNaN; out.mac=0; out.missingRate=0;
    out.values[4] = statusCode(spa::Status::NaNoTest);
    for (size_t e=0;e<m_envs.size();++e) {
        out.values[5+9*e+6] = statusCode(spa::Status::NaNoTest);
        out.values[5+9*e+8] = 1;
    }
    Eigen::VectorXd g = G;
    double observedG=0, observedH=0;
    for (Eigen::Index i=0;i<n;++i) {
        if (H[i]<0 || H[i]>2 || H[i]!=std::round(H[i]))
            throw std::invalid_argument("Local SPAGxE: H must contain integer haplotype counts 0, 1, or 2.");
        if (H[i]==0 && std::isnan(g[i])) g[i]=0;
        if (std::isnan(g[i])) { out.missingRate += 1.0/static_cast<double>(n); continue; }
        if (!std::isfinite(g[i]) || g[i]<-1e-10 || g[i]>H[i]+1e-10)
            throw std::invalid_argument("Local SPAGxE: ancestry-specific dosage must satisfy 0 <= G <= H.");
        observedG += g[i]; observedH += H[i];
    }
    if (!(observedH>0)) return out;
    out.altFreq = std::max(0.0, std::min(1.0, observedG/observedH));
    out.mac = H.sum()*std::min(out.altFreq,1.0-out.altFreq);
    for (Eigen::Index i=0;i<n;++i) if (std::isnan(g[i])) g[i]=H[i]*out.altFreq;
    out.flipped = out.altFreq>.5;
    const double q=out.flipped ? 1-out.altFreq : out.altFreq;
    if (out.flipped) g = H-g;
    if (!(q>0 && q<1)) return out;
    const double marginalVariance=q*(1-q)*H.dot(m_resid.array().square().matrix());
    if (!(marginalVariance>0)) return out;
    const double marginalScore=(g-q*H).dot(m_resid);
    const double zG=marginalScore/std::sqrt(marginalVariance);
    const auto marginal=spa::normalBranch(zG);
    out.values[0]=marginal.negLog10p;
    out.values[1]=zG; out.values[2]=marginalScore/marginalVariance;
    out.values[3]=1/std::sqrt(marginalVariance); out.values[4]=statusCode(marginal.status);
    const bool branchB = 2*math::pnorm(-std::abs(zG))<=m_marginalCutoff;
    Eigen::VectorXd r0;
    if (branchB) {
        const Eigen::VectorXd cg=g.array()-g.mean();
        const double ss=cg.squaredNorm();
        if (ss>0) r0=m_resid.array()-m_resid.mean()-cg.array()*(cg.dot(m_resid)/ss);
    }
    const bool waldEnabled=branchB && m_wald && m_wald->trait!=interaction_wald::TraitType::None;
    interaction_wald::WaldData wd;
    if (waldEnabled) wd=ancestryAdjustedWald(*m_wald,allH);
    for (size_t e=0;e<m_envs.size();++e) {
        const size_t b=5+9*e;
        out.values[b+8]=branchB?2:1;
        Eigen::VectorXd weight;
        if (branchB) {
            if (r0.size()==0) continue;
            weight=m_envs[e].array()*r0.array();
        } else {
            const double denominator=H.dot(m_resid.array().square().matrix());
            const double lambda=(H.array()*m_envs[e].array()*m_resid.array().square()).sum()/denominator;
            weight=(m_envs[e].array()-lambda)*m_resid.array();
        }
        double z=kNaN,variance=kNaN;
        const auto score=localScore(g,H,weight,q,m_spaCutoff,z,variance);
        out.values[b+6]=statusCode(score.status);
        if (!std::isfinite(score.negLog10p) || !(variance>0)) continue;
        double wald=kNaN,final=score.negLog10p;
        if (waldEnabled) {
            wald=interaction_wald::waldInteractionLog10P(
                wd,g,m_envs[e],wald::OrdinalInfo::Observed,true);
            if (std::isfinite(wald)) {
                const double legs[2]={score.negLog10p,wald};
                final=math::cauchyCombineLog10(legs,2);
            }
        }
        out.values[b]=final; out.values[b+1]=wald;
        out.values[b+2]=math::zFromNegLog10P(final,z); out.values[b+3]=z;
        out.values[b+4]=z/std::sqrt(variance); out.values[b+5]=1/std::sqrt(variance);
        out.values[b+7]=score.negLog10p;
    }
    return out;
}
} // namespace spagxe_local
