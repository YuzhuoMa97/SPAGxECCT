// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
// Native known-haplotype-count BGEN collection input for SPAGxE.
#pragma once

#include "local/local_io.hpp"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace spagxe_local {

// A collection consists of one biallelic, unphased Layout-2 file per ancestry.
// Its per-sample ploidy is the known ancestry haplotype count H (0, 1 or 2).
// Only posterior expectations enter LocalVariant; the association kernel does
// not integrate a genotype/ancestry posterior distribution.
class LocalBgenCollection {
public:
    LocalBgenCollection(const std::string &manifest, uint64_t maxBlockBytes);
    ~LocalBgenCollection();
    LocalBgenCollection(const LocalBgenCollection &) = delete;
    LocalBgenCollection &operator=(const LocalBgenCollection &) = delete;
    const std::vector<std::string> &sampleIds() const;
    const std::vector<std::string> &ancestries() const;
    const std::string &sourceCountedAllele() const;
    uint64_t nVariants() const;
    // BGEN probability categories always count the second listed allele.
    // Return canonical ALT directly, avoiding a lossy H-(H-G) round trip when
    // a manifest requests REF-count interchange. The declaration is retained
    // separately by sourceCountedAllele().
    bool next(LocalVariant &variant);
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

} // namespace spagxe_local
