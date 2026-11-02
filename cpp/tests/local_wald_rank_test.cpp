// Local Wald: nuisance redundancy must not discard an identifiable G:E target.
// The reference fits use manually specified full-rank designs, not the
// production rank-selection helper. No reference runtime or data is needed.
#include "interaction/interaction_wald.hpp"
#include "util/wald.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
int checks = 0;
void require(bool ok, const std::string &name) {
    ++checks;
    if (!ok) throw std::runtime_error(name);
}
void close(double got, double expected, const std::string &name) {
    require(std::isfinite(got) && std::isfinite(expected) &&
            std::abs(got - expected) <= 2e-7 * (1 + std::abs(expected)), name);
}

struct Fixture {
    Eigen::VectorXd x, e, h1, h2, g, target, linear, binary, ordinal, time, event;
    explicit Fixture(int n = 480) {
        x.resize(n); e.resize(n); h1.resize(n); h2.resize(n); g.resize(n);
        target.resize(n); linear.resize(n); binary.resize(n); ordinal.resize(n);
        time.resize(n); event.resize(n);
        uint32_t state = 20261026;
        auto uniform = [&]() {
            state = 1664525u * state + 1013904223u;
            return (double(state) + 0.5) / 4294967296.0;
        };
        for (int i = 0; i < n; ++i) {
            x[i] = std::sin(i * .31) + .2 * std::cos(i * .17);
            e[i] = std::cos(i * .23) + .3 * std::sin(i * .13);
            const int a = i % 3, b = (i / 3) % 3;
            h1[i] = (a == 0) + (b == 0);
            h2[i] = (a == 1) + (b == 1);
            g[i] = .23 * h1[i];
            target[i] = g[i] * e[i];
            const double eta = .15 * x[i] + .22 * e[i] +
                               .18 * h1[i] - .1 * h2[i] + .6 * target[i];
            linear[i] = eta + 1.8 * (uniform() - .5);
            binary[i] = uniform() < 1. / (1. + std::exp(-eta));
            const double u = uniform();
            const double latent = eta + std::log(u / (1 - u));
            ordinal[i] = latent < -.65 ? 0 : (latent < .7 ? 1 : 2);
            const double failure = -std::log(uniform()) * std::exp(-eta);
            const double censoring = 1.5 + 2 * uniform();
            time[i] = std::min(failure, censoring);
            event[i] = failure <= censoring;
        }
    }
    interaction_wald::WaldData data(interaction_wald::TraitType trait, bool redundant) const {
        interaction_wald::WaldData out;
        out.trait = trait;
        out.covar.resize(x.size(), redundant ? 7 : 3);
        out.covar.col(0) = x;
        out.covar.col(1) = e;
        if (redundant) {
            out.covar.col(2) = h1;
            out.covar.col(3) = h2;
            out.covar.col(4) = 2 * h1 + .5 * h2;
            out.covar.col(5).setConstant(3);
            out.covar.col(6) = 4 * x - 2 * e + h2;
        } else {
            // [x,e,h2,g] has the same nuisance span as [x,e,h1,h2].
            out.covar.col(2) = h2;
        }
        if (trait == interaction_wald::TraitType::Linear) out.y = linear;
        if (trait == interaction_wald::TraitType::Logistic) out.y = binary;
        if (trait == interaction_wald::TraitType::Ordinal) out.y = ordinal;
        out.time = time;
        out.event = event;
        return out;
    }
    double reference(interaction_wald::TraitType trait, int start = 0) const {
        const int n = int(x.size()) - start;
        Eigen::MatrixXd design(n, 5);
        design.col(0) = x.tail(n); design.col(1) = e.tail(n);
        design.col(2) = h1.tail(n); design.col(3) = h2.tail(n);
        design.col(4) = target.tail(n);
        if (trait == interaction_wald::TraitType::Cox)
            return wald::lastCoefCoxLog10P(time.tail(n), event.tail(n), design);
        if (trait == interaction_wald::TraitType::Ordinal)
            return wald::lastCoefOrdinalLog10P(ordinal.tail(n).cast<int>(), design);
        Eigen::MatrixXd withIntercept(n, 6);
        withIntercept.col(0).setOnes(); withIntercept.rightCols(5) = design;
        if (trait == interaction_wald::TraitType::Linear)
            return wald::lastCoefLinearLog10P(withIntercept, linear.tail(n));
        return wald::lastCoefLogisticLog10P(withIntercept, binary.tail(n));
    }
};

void allFamilies() {
    using TT = interaction_wald::TraitType;
    Fixture f;
    const TT traits[] = {TT::Linear, TT::Logistic, TT::Cox, TT::Ordinal};
    for (TT trait : traits) {
        const auto full = f.data(trait, false);
        auto repeated = f.data(trait, true);
        const double expected = f.reference(trait);
        const double actual = interaction_wald::waldInteractionLog10P(
            repeated, f.g, f.e, wald::OrdinalInfo::Observed, true);
        close(actual, expected, "redundant nuisance matches manual full-rank fit");
        close(interaction_wald::waldInteractionLog10P(full, f.g, f.e), expected,
              "default dispatcher matches full-rank equivalent design");
        close(interaction_wald::waldInteractionLog10P(full, f.g, f.e),
              interaction_wald::waldInteractionLog10P(
                  full, f.g, f.e, wald::OrdinalInfo::Observed, false),
              "default false is unchanged");

        // A column only becomes redundant after dropping an incomplete row.
        repeated.covar(0, 4) += 10;
        if (trait == TT::Cox)
            repeated.time[0] = std::numeric_limits<double>::quiet_NaN();
        else
            repeated.y[0] = std::numeric_limits<double>::quiet_NaN();
        close(interaction_wald::waldInteractionLog10P(
                  repeated, f.g, f.e, wald::OrdinalInfo::Observed, true),
              f.reference(trait, 1), "rank reduction follows complete-case filtering");

        auto unidentifiable = f.data(trait, true);
        unidentifiable.covar.conservativeResize(Eigen::NoChange, 8);
        unidentifiable.covar.col(7) = f.target;
        require(std::isnan(interaction_wald::waldInteractionLog10P(
                    unidentifiable, f.g, f.e, wald::OrdinalInfo::Observed, true)),
                "target in nuisance span is rejected");
        const Eigen::VectorXd constant = Eigen::VectorXd::Ones(f.e.size());
        require(std::isnan(interaction_wald::waldInteractionLog10P(
                    full, f.g, constant, wald::OrdinalInfo::Observed, true)),
                "constant environment makes target equal genotype nuisance");
    }
}

void relativeRankTolerance() {
    using TT = interaction_wald::TraitType;
    Fixture f;
    auto wd = f.data(TT::Linear, true);
    Eigen::VectorXd independent(f.x.size());
    for (int i = 0; i < independent.size(); ++i) independent[i] = std::sin(i * .73);
    wd.covar.col(4) = f.h1 + 1e-12 * independent;
    close(interaction_wald::waldInteractionLog10P(
              wd, f.g, f.e, wald::OrdinalInfo::Observed, true),
          f.reference(TT::Linear), "relative 1e-10 drops numerical nuisance redundancy");
    wd.covar.col(4) *= 1e5;
    close(interaction_wald::waldInteractionLog10P(
              wd, f.g, f.e, wald::OrdinalInfo::Observed, true),
          f.reference(TT::Linear), "rank decision is invariant to nuisance units");
    wd.covar.col(4) = f.target + 1e-12 * independent;
    require(std::isnan(interaction_wald::waldInteractionLog10P(
                wd, f.g, f.e, wald::OrdinalInfo::Observed, true)),
            "target below relative rank threshold is rejected");
}
} // namespace

int main() {
    try {
        allFamilies();
        relativeRankTolerance();
        std::cout << checks << " local Wald rank checks passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "FAIL " << e.what() << "\n";
        return 1;
    }
}
