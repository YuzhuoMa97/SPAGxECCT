// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE local ancestry I/O, Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
// Native count and posterior-expectation input for the SPAGxE streaming API.
#pragma once
#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace spagxe_local {
struct LocalInputOptions {
    // Exactly one: SPAGxE G/H manifest, SQRLOC01 prefix, or BGEN manifest.
    // Manifest header: ancestry genotype hapcount counted_allele (TAB).
    std::string manifest;
    std::string binaryPrefix;
    std::string bgenManifest;
    uint64_t maxBlockBytes = 268435456;
};
struct LocalVariant {
    // CHROM POS ID REF ALT. G is always canonical ALT count, including
    // explicit H-G recoding when the source declares REF counted allele.
    std::array<std::string, 5> metadata;
    std::vector<Eigen::VectorXd> genotype, hapcount;
    uint64_t index = 0;
    bool hasDosage = false; // fractional G or non-point-mass posterior source
};
class LocalSource {
public:
    explicit LocalSource(const LocalInputOptions &options);
    ~LocalSource();
    LocalSource(LocalSource &&) noexcept;
    LocalSource &operator=(LocalSource &&) noexcept;
    LocalSource(const LocalSource &) = delete;
    LocalSource &operator=(const LocalSource &) = delete;
    const std::vector<std::string> &sampleIds() const;
    const std::vector<std::string> &ancestries() const;
    const std::string &sourceCountedAllele() const;
    uint64_t nVariants() const; // zero means not known without scanning (text).
    void align(const std::vector<std::string> &iids);
    bool next(LocalVariant &variant);
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
} // namespace spagxe_local
