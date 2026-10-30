// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE standalone command-line interface.
// SPAGxE algorithm: Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
// Derived in part from the GRAB command-line infrastructure; see THIRD_PARTY_NOTICES.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>

namespace cli {
struct Args {
    std::string method; // Empty until dispatch: SPAGxE is the analysis default.
    std::string helpTopic;
    bool showVersion = false;
    std::string phenoFile, covarFile, covarName, phenoName, residName;
    std::string regressionModel;
    bool saveResid = false;
    std::string pcCols = "PC1,PC2,PC3,PC4";
    bool pcColsExplicit = false;
    std::string envName;
    double spagxeMarginalCutoff = 0.001;
    std::string bfilePrefix, pfilePrefix, vcfFile, bcfFile, bgenFile, bgenRefMode;
    std::string grmTableFile, spGrmPlink2File;
    std::string localManifest, localPrefix, localBgenManifest;
    uint64_t localMaxBlockBytes = 268435456;
    bool localMaxBlockBytesExplicit = false;
    std::string cctGlobal, cctLocal, cctPhenotype, cctWeights = "0.5,0.5";
    bool cctPerAncestryEqual = false, cctWeightsExplicit = false;
    std::unordered_set<std::string> explicitFlags;
    std::string outPrefix, compression;
    std::string extractFile, excludeFile, chrSpec, keepFile, removeFile;
    bool calAfCoef = false;
    bool intPheno = false;
    bool packLocal = false;
    std::string localBimFile, localFamFile;
    uint32_t localBlockVariants = 64;
    int localZstdLevel = 3;
    bool localBlockVariantsExplicit = false, localZstdLevelExplicit = false;
    double spaCutoff = 2.0;
    double missingCutoff = 0.1;
    double minMafCutoff = 1e-5;
    double minMacCutoff = 10.0;
    double hweCutoff = 0.0;
    double hardCallThreshold = 0.1;
    bool hardCallThresholdExplicit = false;
    int nthread = 1;
    int nSnpPerChunk = 8192;
    int compressionLevel = 0;
    bool compressionLevelExplicit = false;
    uint64_t seed = 0; // Accepted for upstream CLI compatibility; analysis is deterministic.
};
int run(int argc, char *argv[]);
Args parseArgs(int argc, char *argv[]);
void printHelp(const std::string &topic);
void printVersion();
} // namespace cli
