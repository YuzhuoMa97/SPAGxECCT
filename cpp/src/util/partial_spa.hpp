// A shared, deterministic central-Gaussian partition for interaction SPA.
#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "allele_model/binomial_cgf.hpp"

namespace partial_spa {

struct Partition {
    std::vector<uint32_t> exactPositions;
    std::vector<uint32_t> centralPositions;
    Eigen::VectorXd exactWeights;
    bool allExact = true;
    int nActive = 0;
    int sourceSize = 0;
    double centralWeightSum = 0.0;
    double centralWeightSq = 0.0;
    double totalWeightSum = 0.0;
    double totalWeightSq = 0.0;
};

struct Workspace {
    std::vector<double> afGather;
    std::vector<double> hapGather;
    binomial_cgf::Context partial;
    binomial_cgf::Context full;
    double fullMean = 0.0;
    double fullVariance = 0.0;
};

namespace detail {

inline void validate(const Eigen::VectorXd &weight,
                     const Eigen::VectorXd *hap,
                     const Eigen::VectorXd *af) {
    if (hap && af)
        throw std::invalid_argument("Partial SPA cannot combine hap counts and individual AF");
    if (weight.size() > std::numeric_limits<int>::max())
        throw std::length_error("Partial SPA sample count exceeds the kernel limit");
    if ((hap && hap->size() != weight.size()) ||
        (af && af->size() != weight.size()))
        throw std::invalid_argument("Partial SPA vectors have inconsistent sizes");
    for (Eigen::Index i = 0; i < weight.size(); ++i) {
        if (!std::isfinite(weight[i]))
            throw std::invalid_argument("Partial SPA weights must be finite");
        if (hap && (!std::isfinite((*hap)[i]) || (*hap)[i] < 0.0 ||
                    (*hap)[i] > 2.0 || std::floor((*hap)[i]) != (*hap)[i]))
            throw std::invalid_argument("Partial SPA hap counts must be 0, 1 or 2");
        if (af && (!std::isfinite((*af)[i]) || (*af)[i] < 0.0 || (*af)[i] > 1.0))
            throw std::invalid_argument("Partial SPA allele frequencies must be in [0,1]");
    }
}

inline bool active(Eigen::Index i, const Eigen::VectorXd &weight,
                   const Eigen::VectorXd *hap, const Eigen::VectorXd *af) {
    return weight[i] != 0.0 && (!hap || (*hap)[i] > 0.0) &&
           (!af || ((*af)[i] > 0.0 && (*af)[i] < 1.0));
}

// Type-7 sample quantiles use linear interpolation at (n-1)*p. Selection
// avoids sorting all subjects on each variant-specific local/AF partition.
inline long double quantile(std::vector<double> &x, long double p) {
    const long double pos = static_cast<long double>(x.size() - 1) * p;
    const size_t lo = static_cast<size_t>(pos);
    std::nth_element(x.begin(), x.begin() + lo, x.end());
    const long double a = x[lo];
    const long double fraction = pos - lo;
    if (fraction == 0.0L) return a;
    const long double b = *std::min_element(x.begin() + lo + 1, x.end());
    return a + fraction * (b - a);
}

inline double finiteMoment(long double x) {
    const double d = static_cast<double>(x);
    if (!std::isfinite(d))
        throw std::overflow_error("Partial SPA moment exceeds floating-point range");
    return d;
}

} // namespace detail

// Partition the actual score weights, not the unmultiplied phenotype residual.
// Only random, nonzero contributions determine the Tukey 1.5-IQR fences.
// An empty outlier block retains the full law: a large N alone cannot justify
// normal approximation when the allele is rare or weights are concentrated.
inline Partition split(const Eigen::VectorXd &weight,
                       const Eigen::VectorXd *hap = nullptr,
                       const Eigen::VectorXd *af = nullptr) {
    detail::validate(weight, hap, af);
    Partition out;
    out.sourceSize = static_cast<int>(weight.size());
    std::vector<double> eligible;
    eligible.reserve(static_cast<size_t>(weight.size()));
    long double totalSum = 0.0L, totalSq = 0.0L;
    for (Eigen::Index i = 0; i < weight.size(); ++i) {
        const long double w = weight[i];
        if (!hap && !af) {
            totalSum += w;
            totalSq += w * w;
        }
        if (detail::active(i, weight, hap, af)) eligible.push_back(weight[i]);
    }
    out.nActive = static_cast<int>(eligible.size());
    out.totalWeightSum = detail::finiteMoment(totalSum);
    out.totalWeightSq = detail::finiteMoment(totalSq);
    if (!eligible.empty()) {
        const long double q1 = detail::quantile(eligible, 0.25L);
        const long double q3 = detail::quantile(eligible, 0.75L);
        const long double width = 1.5L * (q3 - q1);
        const long double lower = q1 - width;
        const long double upper = q3 + width;
        int centralActive = 0;
        for (Eigen::Index i = 0; i < weight.size(); ++i) {
            if (detail::active(i, weight, hap, af) &&
                (weight[i] < lower || weight[i] > upper)) {
                out.exactPositions.push_back(static_cast<uint32_t>(i));
            } else {
                out.centralPositions.push_back(static_cast<uint32_t>(i));
                if (detail::active(i, weight, hap, af)) ++centralActive;
            }
        }
        out.allExact = out.exactPositions.empty() || centralActive == 0;
    }
    if (out.allExact) {
        out.exactPositions.resize(static_cast<size_t>(weight.size()));
        for (Eigen::Index i = 0; i < weight.size(); ++i)
            out.exactPositions[static_cast<size_t>(i)] = static_cast<uint32_t>(i);
        out.centralPositions.clear();
    }
    out.exactWeights.resize(static_cast<Eigen::Index>(out.exactPositions.size()));
    for (size_t j = 0; j < out.exactPositions.size(); ++j)
        out.exactWeights[static_cast<Eigen::Index>(j)] = weight[out.exactPositions[j]];
    long double sum = 0.0L, sq = 0.0L;
    if (!hap && !af)
        for (uint32_t i : out.centralPositions) {
            const long double w = weight[i];
            sum += w;
            sq += w * w;
        }
    out.centralWeightSum = detail::finiteMoment(sum);
    out.centralWeightSq = detail::finiteMoment(sq);
    return out;
}

// Bind a partition to the current marker. The referenced weight/AF/H vectors
// and partition must remain alive and unchanged until the solve completes.
// Central moments are accumulated directly, never as total minus outlier:
// subtraction would lose a small central variance beside a very large outlier.
inline void bind(const Partition &partition, const Eigen::VectorXd &weight,
                 double q, const Eigen::VectorXd *af,
                 const Eigen::VectorXd *hap, Workspace &ws) {
    if (hap && af)
        throw std::invalid_argument("Partial SPA cannot combine hap counts and individual AF");
    if (partition.sourceSize != weight.size() ||
        (hap && hap->size() != weight.size()) ||
        (af && af->size() != weight.size()))
        throw std::invalid_argument("Partial SPA partition does not match the input vectors");
    if (!af && (!std::isfinite(q) || q < 0.0 || q > 1.0))
        throw std::invalid_argument("Partial SPA shared allele frequency must be in [0,1]");
    ws.full = binomial_cgf::Context{};
    ws.full.resid = weight.data();
    ws.full.af = af ? af->data() : nullptr;
    ws.full.hap = hap ? hap->data() : nullptr;
    ws.full.q = q;
    ws.full.nOutlier = static_cast<int>(weight.size());

    long double fullMean = 0.0L, fullVariance = 0.0L;
    long double mean = 0.0L, variance = 0.0L;
    if (!af && !hap) {
        fullMean = 2.0L * q * partition.totalWeightSum;
        fullVariance = 2.0L * q * (1.0L - q) * partition.totalWeightSq;
        mean = 2.0L * q * partition.centralWeightSum;
        variance = 2.0L * q * (1.0L - q) * partition.centralWeightSq;
    } else {
        for (Eigen::Index i = 0; i < weight.size(); ++i) {
            const long double w = weight[i];
            const long double p = af ? (*af)[i] : q;
            const long double h = hap ? (*hap)[i] : 2.0L;
            fullMean += h * p * w;
            fullVariance += h * p * (1.0L - p) * w * w;
        }
        for (uint32_t i : partition.centralPositions) {
            const long double w = weight[i];
            const long double p = af ? (*af)[i] : q;
            const long double h = hap ? (*hap)[i] : 2.0L;
            mean += h * p * w;
            variance += h * p * (1.0L - p) * w * w;
        }
    }
    ws.fullMean = detail::finiteMoment(fullMean);
    ws.fullVariance = detail::finiteMoment(fullVariance);
    ws.partial = ws.full;
    if (partition.allExact || !(variance > 0.0L)) return;

    ws.partial.resid = partition.exactWeights.data();
    ws.partial.nOutlier = static_cast<int>(partition.exactPositions.size());
    ws.partial.mean = detail::finiteMoment(mean);
    ws.partial.var = detail::finiteMoment(variance);
    if (af) {
        ws.afGather.resize(partition.exactPositions.size());
        for (size_t j = 0; j < partition.exactPositions.size(); ++j)
            ws.afGather[j] = (*af)[partition.exactPositions[j]];
        ws.partial.af = ws.afGather.data();
    }
    if (hap) {
        ws.hapGather.resize(partition.exactPositions.size());
        for (size_t j = 0; j < partition.exactPositions.size(); ++j)
            ws.hapGather[j] = (*hap)[partition.exactPositions[j]];
        ws.partial.hap = ws.hapGather.data();
    }
}

} // namespace partial_spa
