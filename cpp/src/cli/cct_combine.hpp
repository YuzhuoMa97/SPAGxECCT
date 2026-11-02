// SPDX-License-Identifier: GPL-3.0-or-later
// Native global/local CCT utility, Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
#pragma once
#include <string>
#include <vector>
namespace spagxe_cct {
struct Options {
    std::string globalFile, localFile, phenotype, outPrefix, compression;
    std::vector<std::string> envNames;
    double globalWeight=.5, localWeight=.5;
    bool perAncestryEqual=false;
    int compressionLevel=0;
};
void run(const Options &options);
} // namespace spagxe_cct
