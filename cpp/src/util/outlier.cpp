// IQR partition selection and residual gathering.

#include "util/outlier.hpp"
#include "util/logging.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace {

// R-style type-7 quantile on a pre-sorted vector.
double type7QuantileOnSorted(
    const std::vector<double> &sorted,
    double p
) {
    const int n = static_cast<int>(sorted.size());
    double h = (n - 1) * p;
    int j = static_cast<int>(h);
    double g = h - j;
    if (j + 1 >= n) return sorted[n - 1];
    return sorted[j] + g * (sorted[j + 1] - sorted[j]);
}

std::vector<uint32_t> selectOutlierPositions(
    const Eigen::VectorXd &residuals,
    double lowerQuartile,
    double upperQuartile,
    double iqr,
    double multiplier
) {
    const double lo = lowerQuartile - multiplier * iqr;
    const double hi = upperQuartile + multiplier * iqr;
    const int count = static_cast<int>(residuals.size());
    std::vector<uint32_t> positions;
    for (int i = 0; i < count; ++i)
        if (residuals[i] < lo || residuals[i] > hi)
            positions.push_back(static_cast<uint32_t>(i));
    return positions;
}

OutlierData gatherPartition(
    const Eigen::VectorXd &residuals,
    std::vector<uint32_t> outlierPositions
) {
    const int count = static_cast<int>(residuals.size());
    std::vector<bool> isOutlier(count, false);
    for (uint32_t row : outlierPositions)
        isOutlier[row] = true;
    std::vector<uint32_t> centralPositions;
    centralPositions.reserve(static_cast<size_t>(count) - outlierPositions.size());
    for (int i = 0; i < count; ++i)
        if (!isOutlier[i])
            centralPositions.push_back(static_cast<uint32_t>(i));

    OutlierData partition;
    partition.posOutlier = std::move(outlierPositions);
    partition.posNonOutlier = std::move(centralPositions);
    const int nOutlier = static_cast<int>(partition.outlierCount());
    const int nCentral = static_cast<int>(partition.centralCount());
    partition.residOutlier.resize(nOutlier);
    partition.residNonOutlier.resize(nCentral);
    partition.resid2NonOutlier.resize(nCentral);
    for (int i = 0; i < nOutlier; ++i)
        partition.residOutlier[i] = residuals[partition.posOutlier[i]];
    for (int i = 0; i < nCentral; ++i) {
        const double value = residuals[partition.posNonOutlier[i]];
        partition.residNonOutlier[i] = value;
        partition.resid2NonOutlier[i] = value * value;
    }
    return partition;
}

} // anonymous namespace

OutlierData detectOutliers(
    const Eigen::VectorXd &resid,
    double outlierRatio
) {
    const int N = static_cast<int>(resid.size());

    std::vector<double> sorted(resid.data(), resid.data() + N);
    std::sort(sorted.begin(), sorted.end());

    double q25 = type7QuantileOnSorted(sorted, 0.25);
    double q75 = type7QuantileOnSorted(sorted, 0.75);
    double iqr = q75 - q25;

    auto findOutliers = [&](double r) {
        return selectOutlierPositions(resid, q25, q75, iqr, r);
    };

    auto posOutlier = findOutliers(outlierRatio);
    // For discrete-valued residuals (e.g. binary phenotypes whose logistic
    // residuals take only two values around ±p), every shrink of the cutoff
    // band still leaves the band strictly enclosing all observed residual
    // values, so the loop would never terminate.  Two safeguards:
    //   1. IQR <= 0 → cutoffs cannot adapt; bail out with an empty outlier
    //      set and emit [WARN].
    //   2. Otherwise cap the shrink loop at 200 iterations and warn if no
    //      outlier was ever found.
    constexpr int kMaxShrinkIter = 200;
    if (iqr <= 0.0 && posOutlier.empty()) {
        warnMsg(
            "  detectOutliers: residual IQR is zero (degenerate distribution; "
            "typically a binary phenotype with extreme imbalance); proceeding "
            "with no outliers"
        );
    } else {
        int shrinkIter = 0;
        while (posOutlier.empty() && shrinkIter < kMaxShrinkIter) {
            outlierRatio *= 0.8;
            posOutlier = findOutliers(outlierRatio);
            ++shrinkIter;
        }
        if (posOutlier.empty()) {
            warnMsg(
                "  detectOutliers: no outliers found after %d shrink iterations "
                "(IQR ratio %.3g); residual distribution lacks tails "
                "(e.g. binary phenotype); proceeding with no outliers",
                shrinkIter, outlierRatio
            );
        }
    }

    return gatherPartition(resid, std::move(posOutlier));
}
