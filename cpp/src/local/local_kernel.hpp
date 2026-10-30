// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE local-ancestry GxE kernel.
// Copyright (c) 2026 Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
#pragma once

#include "interaction/interaction_wald.hpp"
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace spagxe_local {
struct LocalResult {
    std::vector<double> values; // five marginal cells, nine per environment
    double altFreq;
    double mac;
    double missingRate;
    bool flipped = false; // statistics test the ancestry-specific minor allele
};

// A pure C++ conditional-binomial interaction test, G_i | H_i ~ Bin(H_i,q).
// Residuals/environments/raw traits are in the phenotype's complete-case order.
// evaluate() is const: all marker scratch is local and worker-safe.
class LocalKernel {
  public:
    LocalKernel(Eigen::VectorXd resid,
                std::vector<std::string> envNames,
                std::vector<Eigen::VectorXd> envVecs,
                std::shared_ptr<const interaction_wald::WaldData> wald,
                double marginalCutoff, double spaCutoff);
    int resultSize() const;
    std::string getHeaderColumns() const;
    LocalResult evaluate(const Eigen::VectorXd &G,
                         const Eigen::VectorXd &H,
                         const Eigen::MatrixXd &allH) const;
  private:
    Eigen::VectorXd m_resid;
    std::vector<std::string> m_envNames;
    std::vector<Eigen::VectorXd> m_envs;
    std::shared_ptr<const interaction_wald::WaldData> m_wald;
    double m_marginalCutoff, m_spaCutoff;
};
} // namespace spagxe_local
