// SPDX-License-Identifier: GPL-3.0-or-later
// Exact support and endpoint atoms for sums of independent weighted binomials.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>

namespace spa {
struct BinomialSupport {
    long double lower = 0, upper = 0;
    long double logLowerMass = 0, logUpperMass = 0, scale = 0;
    bool valid = true;

    explicit BinomialSupport(double shift = 0) noexcept
        : lower(shift), upper(shift), scale(std::abs(shift)), valid(std::isfinite(shift)) {}

    void add(double weight, double count, double p, bool computeAtomMasses = true) noexcept {
        if (!std::isfinite(weight) || !std::isfinite(count) || !std::isfinite(p) ||
            count < 0 || p < 0 || p > 1) { valid = false; return; }
        if (count == 0 || weight == 0) return;
        const long double amount = static_cast<long double>(weight) * count;
        scale += std::abs(amount);
        // Degenerate q=0/1 variables shift the support, without constraining
        // an endpoint event or adding log(0) to its probability.
        if (p == 0) return;
        if (p == 1) { lower += amount; upper += amount; return; }
        lower += std::min(0.L, amount);
        upper += std::max(0.L, amount);
        if (!computeAtomMasses) return;
        const long double yes = count * std::log(static_cast<long double>(p));
        const long double no = count * std::log1p(-static_cast<long double>(p));
        logLowerMass += weight > 0 ? no : yes;
        logUpperMass += weight > 0 ? yes : no;
    }

    // A finite-support tail is not a failed root solve. Strictly beyond the
    // support its probability is zero; at an endpoint use the inclusive atom.
    // Interior points still require SPA. The rounding allowance scales with
    // the statistic (no fixed absolute epsilon or lattice correction).
    // If a later BN variance-ratio transformation is requested, endpoint atoms
    // are not silently transformed: leave them to the existing solver policy.
    bool tail(double target, bool lowerTail, double &logp,
              bool endpointAtoms = true) const noexcept {
        if (!valid || !std::isfinite(target)) return false;
        const long double x = target;
        const long double tol = 8 * std::numeric_limits<double>::epsilon() *
            std::max(scale, std::abs(x));
        const double zero = -std::numeric_limits<double>::infinity();
        if (x < lower - tol) { logp = lowerTail ? zero : 0.; return true; }
        if (x > upper + tol) { logp = lowerTail ? 0. : zero; return true; }
        if (!endpointAtoms) return false;
        if (std::abs(x - lower) <= tol) {
            logp = lowerTail ? static_cast<double>(logLowerMass) : 0.; return true;
        }
        if (std::abs(x - upper) <= tol) {
            logp = lowerTail ? 0. : static_cast<double>(logUpperMass); return true;
        }
        return false;
    }
};
} // namespace spa
