// Public method scope: individual-AF and relatedness inputs are separate.
#include "interaction/interaction_test.hpp"
#include "io/sparse_grm.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
int checks = 0;
void require(bool condition, const std::string &label) {
    ++checks;
    if (!condition) throw std::runtime_error(label);
}
}

int main() {
    try {
        constexpr int n = 32;
        Eigen::VectorXd r(n), e(n), g(n);
        auto af = std::make_shared<AFData>();
        af->nPC = 1;
        af->onePlusPCs.resize(n, 2);
        af->onePlusPCs.col(0).setOnes();
        af->onePlusPCs.col(1) = Eigen::VectorXd::LinSpaced(n, -1, 1);
        const Eigen::MatrixXd inverse =
            (af->onePlusPCs.transpose() * af->onePlusPCs).inverse();
        af->XtX_inv_Xt = inverse * af->onePlusPCs.transpose();
        af->sqrt_XtX_inv_diag = inverse.diagonal().cwiseSqrt();
        std::vector<SparseGRM::Entry> entries;
        for (int i = 0; i < n; ++i) {
            r[i] = std::sin(.7 * i) + .1 * std::cos(.13 * i);
            e[i] = std::cos(.37 * i) + .02 * i;
            g[i] = i % 3;
            entries.push_back({static_cast<uint32_t>(i), static_cast<uint32_t>(i), 1});
        }
        entries.push_back({1, 0, .125});
        const auto grm = std::make_shared<const SparseGRM>(
            SparseGRM::fromEntries(n, std::move(entries)));

        bool rejected = false;
        std::string reason;
        try {
            InteractionScan invalid(r, {"E"}, {e}, grm, nullptr, 0, 100, af);
        } catch (const std::runtime_error &error) {
            rejected = true;
            reason = error.what();
        }
        require(rejected, "individual-AF with GRM must be rejected by the API");
        require(reason.find("GRM") != std::string::npos,
                "the API rejection identifies the incompatible input");

        for (int mode = 0; mode < 3; ++mode) {
            InteractionScan allowed(r, {"E"}, {e}, mode == 1 ? grm : nullptr,
                                nullptr, 0, 100, mode == 2 ? af : nullptr);
            require(allowed.resultSize() == 14, "allowed method output schema");
            std::vector<double> result;
            Eigen::VectorXd copy = g;
            allowed.getResultVec(copy, g.mean() / 2, 0, result);
            require(result.size() == 14, "allowed method evaluates a marker");
            require(std::isfinite(result[5]), "allowed method has a finite interaction P");
            require(allowed.clone()->resultSize() == 14, "allowed worker clone retains schema");
        }
        std::cout << "PASS " << checks << " method scope checks\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAIL after " << checks << " checks: " << error.what() << '\n';
        return 1;
    }
}
