// simd_dispatch.hpp — Runtime SIMD tier detection
//
// On x86_64: detects AVX-512 / AVX2 / scalar at process startup via
// __builtin_cpu_supports().  Hot kernels dispatch to the widest available
// SIMD implementation, enabling a single binary that exploits AVX-512 when
// present and falls back gracefully.
//
// On non-x86 (e.g. arm64): always returns Scalar.  AVX kernels are not
// compiled in; callers fall through to the scalar variant.
#pragma once

enum class SimdLevel : int { Scalar = 0, AVX2 = 1, AVX512 = 2 };

namespace simd_dispatch_detail {

inline SimdLevel detectLevel() {
#if defined(__x86_64__) || defined(_M_X64)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512vl"))
        return SimdLevel::AVX512;
    if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"))
        return SimdLevel::AVX2;
    return SimdLevel::Scalar;
#else
    return SimdLevel::Scalar;
#endif
}

} // namespace simd_dispatch_detail

inline SimdLevel simdLevel() {
    static const SimdLevel level = simd_dispatch_detail::detectLevel();
    return level;
}
