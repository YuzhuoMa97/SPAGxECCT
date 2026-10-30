// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
#pragma once
#include "local/local_io.hpp"
#include <cstdint>
#include <string>

namespace spagxe_local {
struct PackOptions {
    LocalInputOptions input;
    std::string outPrefix;
    std::string bimFile; // optional six-column CHROM ID CM POS ALT REF
    std::string famFile; // optional six-column FID IID PAT MAT SEX PHENO
    uint32_t blockVariants = 64;
    int compressionLevel = 3;
};
// Stream validated integer G/H into the SQRLOC01 count container. Output is
// canonical ALT. Non-point posterior input is rejected, including integer
// posterior means. Existing output files are never overwritten.
void packLocalCounts(const PackOptions &options);
} // namespace spagxe_local
