// Residual partitions for partial-normal saddlepoint evaluation.
//
// Splits a residual vector into outlier and non-outlier subsets using the
// 1.5×IQR rule (or any other multiplier).  The non-outlier sub-vector lets
// downstream SPA implementations replace an O(N) empirical CGF sum with a
// closed-form Gaussian approximation of equal mean and variance.
#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <vector>

// ======================================================================
// Outlier data — precomputed from residuals
// ======================================================================

struct OutlierData {
    std::vector<uint32_t> posOutlier;
    std::vector<uint32_t> posNonOutlier;
    Eigen::VectorXd residOutlier;     // resid[posOutlier]
    Eigen::VectorXd residNonOutlier;  // resid[posNonOutlier]
    Eigen::VectorXd resid2NonOutlier; // resid[posNonOutlier]^2

    // Position vectors determine the partition sizes and preserve input order.
    std::size_t outlierCount() const noexcept {
        return posOutlier.size();
    }

    std::size_t centralCount() const noexcept {
        return posNonOutlier.size();
    }
};

// Detect outlier residuals using IQR method.
// The ratio shrinks by 0.8 while searching for outliers, up to 200 steps.
// Degenerate IQR or a search with no tail values returns an empty outlier set.
OutlierData detectOutliers(
    const Eigen::VectorXd &resid,
    double outlierRatio
);
