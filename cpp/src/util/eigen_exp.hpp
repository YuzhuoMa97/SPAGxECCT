// SPAGxE v1: safe Eigen exponentials for legacy MinGW toolchains.
// Copyright (c) 2026 Yuzhuo Ma, yuzhuoma@amss.ac.cn.
// Keep matrix products and other Eigen operations vectorized.
#pragma once
#include <Eigen/Core>
#include <cmath>

namespace spagxe_numeric {

// An owning result also prevents a lazy expression from outliving its input.
// Eigen's AVX Packet4d pexp crashes with the legacy Windows GCC toolchain used
// for v1 validation (fault traced to Eigen AVX/MathFunctions.h). A lambda has
// no packetOp, so unaryExpr evaluates scalar std::exp on those builds only.
template <typename Derived>
Eigen::ArrayXd eigenSafeExp(const Eigen::ArrayBase<Derived> &values) {
#if defined(_WIN32) && defined(__GNUC__) && !defined(__clang__) && (__GNUC__ <= 10)
    return values.unaryExpr([](double value) { return std::exp(value); }).eval();
#else
    return values.exp().eval();
#endif
}
} // namespace spagxe_numeric
