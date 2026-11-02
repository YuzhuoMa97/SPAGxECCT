// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE native local ancestry workflow, Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
#pragma once
#include "local/local_io.hpp"
#include <string>
#include <unordered_set>
#include <vector>
namespace spagxe_local {
struct LocalWorkflowOptions {
    LocalInputOptions input;
    std::string phenoFile,covarFile,phenoNameSpec,regressionModel;
    std::vector<std::string> residNames,envNames,covarNames;
    std::string keepFile,removeFile,extractFile,excludeFile,outPrefix,compression;
    std::unordered_set<std::string> chromosomes;
    int compressionLevel=0,nthreads=1,nSnpPerChunk=64;
    double marginalCutoff=0.001,spaCutoff=2.0;
    double missingCutoff=0.1,minMafCutoff=1e-5,minMacCutoff=10,hweCutoff=0;
    bool saveResid=false;
};
void scanHomogeneousInteractionsLocal(const LocalWorkflowOptions &options);
} // namespace spagxe_local
