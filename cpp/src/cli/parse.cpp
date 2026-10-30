// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE argument parser. Derived in part from GRAB; see THIRD_PARTY_NOTICES.
#include "cli/cli.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace cli {
namespace {
[[noreturn]] void invalid(const std::string &message) {
    throw std::invalid_argument(message);
}
double parseDouble(const std::string &value, const std::string &flag) {
    try {
        size_t pos = 0;
        double parsed = std::stod(value, &pos);
        if (pos != value.size() || !std::isfinite(parsed))
            invalid(flag + " requires a finite numeric value, got '" + value + "'.");
        return parsed;
    } catch (const std::invalid_argument &) {
        invalid(flag + " requires a finite numeric value, got '" + value + "'.");
    } catch (const std::out_of_range &) {
        invalid(flag + " value is out of range: '" + value + "'.");
    }
}
int parseInt(const std::string &value, const std::string &flag) {
    try {
        size_t pos = 0;
        int parsed = std::stoi(value, &pos);
        if (pos != value.size()) invalid(flag + " requires an integer, got '" + value + "'.");
        return parsed;
    } catch (const std::invalid_argument &) {
        invalid(flag + " requires an integer, got '" + value + "'.");
    } catch (const std::out_of_range &) {
        invalid(flag + " value is out of range: '" + value + "'.");
    }
}
uint64_t parseSeed(const std::string &value, const std::string &flag) {
    try {
        if (value.empty() || value.front() == '-') invalid(flag + " requires a non-negative integer.");
        size_t pos = 0;
        auto parsed = std::stoull(value, &pos);
        if (pos != value.size()) invalid(flag + " requires a non-negative integer.");
        return static_cast<uint64_t>(parsed);
    } catch (const std::invalid_argument &) {
        invalid(flag + " requires a non-negative integer, got '" + value + "'.");
    } catch (const std::out_of_range &) {
        invalid(flag + " value is out of range: '" + value + "'.");
    }
}
} // namespace

Args parseArgs(int argc, char *argv[]) {
    Args a;
    std::unordered_set<std::string> seen;
    auto markSeen = [&seen](const std::string &flag) {
        if (!seen.insert(flag).second) invalid(flag + " was specified more than once.");
    };
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            a.helpTopic = i + 1 < argc && argv[i + 1][0] != '-' ? argv[++i] : "__short__";
            return a;
        }
        if (arg == "--version" || arg == "-V") { a.showVersion = true; return a; }
        auto next = [&]() -> std::string {
            markSeen(arg);
            if (i + 1 >= argc || std::string(argv[i + 1]).rfind("--", 0) == 0 ||
                std::string(argv[i + 1]) == "-h" || std::string(argv[i + 1]) == "-V")
                invalid("Missing value for " + arg + ".");
            const std::string value = argv[++i];
            if (value.empty()) invalid("Empty value for " + arg + ".");
            return value;
        };
        if (arg == "--method") a.method = next();
        else if (arg == "--pheno") a.phenoFile = next();
        else if (arg == "--covar") a.covarFile = next();
        else if (arg == "--covar-name") a.covarName = next();
        else if (arg == "--pheno-name") a.phenoName = next();
        else if (arg == "--resid-name") a.residName = next();
        else if (arg == "--regression-model") a.regressionModel = next();
        else if (arg == "--save-resid") { markSeen(arg); a.saveResid = true; }
        else if (arg == "--algorithm-profile" || arg == "--spa-mode" ||
                 arg == "--spa-exact" || arg == "--outlier-iqr-multiplier" ||
                 arg == "--grm-spa-adjust")
            invalid(arg + " was removed: this release uses one guarded partial-normal SPA implementation and score-based GRM adjustment. Remove this option; --spa-z-threshold remains available (default 2).");
        else if (arg == "--pc-cols") { a.pcCols = next(); a.pcColsExplicit = true; }
        else if (arg == "--envir-name") a.envName = next();
        else if (arg == "--spagxe-marginal-cutoff") a.spagxeMarginalCutoff = parseDouble(next(), arg);
        else if (arg == "--bfile") a.bfilePrefix = next();
        else if (arg == "--pfile") a.pfilePrefix = next();
        else if (arg == "--vcf") a.vcfFile = next();
        else if (arg == "--bcf") a.bcfFile = next();
        else if (arg == "--bgen") {
            a.bgenFile = next();
            if (i + 1 >= argc) invalid("--bgen requires FILE followed by ref-first, ref-last, or ref-unknown.");
            a.bgenRefMode = argv[++i];
            if (a.bgenRefMode != "ref-first" && a.bgenRefMode != "ref-last" && a.bgenRefMode != "ref-unknown")
                invalid("--bgen requires FILE followed by ref-first, ref-last, or ref-unknown.");
        }
        else if (arg == "--grm-file" || arg == "--sp-grm-grab") {
            if (!a.grmTableFile.empty())
                invalid("--grm-file was specified more than once, including a legacy alias.");
            a.grmTableFile = next();
        }
        else if (arg == "--sp-grm-plink2") a.spGrmPlink2File = next();
        else if (arg == "--pack-local") { markSeen(arg); a.packLocal = true; }
        else if (arg == "--local-bim") a.localBimFile = next();
        else if (arg == "--local-fam") a.localFamFile = next();
        else if (arg == "--local-block-variants") {
            const auto value = parseSeed(next(), arg);
            if (value < 1 || value > 65536) invalid("--local-block-variants must be in 1..65536.");
            a.localBlockVariants = static_cast<uint32_t>(value); a.localBlockVariantsExplicit = true;
        }
        else if (arg == "--local-zstd-level") { a.localZstdLevel = parseInt(next(), arg); a.localZstdLevelExplicit = true; }
        else if (arg == "--local-manifest") a.localManifest = next();
        else if (arg == "--local-bgen-manifest") a.localBgenManifest = next();
        else if (arg == "--local-prefix" || arg == "--local-binary") {
            if (!a.localPrefix.empty()) invalid("--local-prefix and --local-binary select the same input and cannot be combined.");
            a.localPrefix = next();
        }
        else if (arg == "--local-max-block-bytes") {
            a.localMaxBlockBytes = parseSeed(next(), arg); a.localMaxBlockBytesExplicit = true;
        }
        else if (arg == "--cct-global") a.cctGlobal = next();
        else if (arg == "--cct-local") a.cctLocal = next();
        else if (arg == "--cct-phenotype") a.cctPhenotype = next();
        else if (arg == "--cct-weights") { a.cctWeights = next(); a.cctWeightsExplicit = true; }
        else if (arg == "--cct-per-ancestry-equal") { markSeen(arg); a.cctPerAncestryEqual = true; }
        else if (arg == "--out") a.outPrefix = next();
        else if (arg == "--spa-z-threshold") a.spaCutoff = parseDouble(next(), arg);
        else if (arg == "--geno") a.missingCutoff = parseDouble(next(), arg);
        else if (arg == "--maf") a.minMafCutoff = parseDouble(next(), arg);
        else if (arg == "--mac") a.minMacCutoff = parseDouble(next(), arg);
        else if (arg == "--hwe") a.hweCutoff = parseDouble(next(), arg);
        else if (arg == "--hard-call-threshold") { a.hardCallThreshold = parseDouble(next(), arg); a.hardCallThresholdExplicit = true; }
        else if (arg == "--threads") a.nthread = parseInt(next(), arg);
        else if (arg == "--chunk-ksnp") {
            double ksnp = parseDouble(next(), arg);
            double halves = ksnp * 2.0;
            if (ksnp < 0.5 || ksnp > std::numeric_limits<int>::max() / 1024.0 ||
                std::fabs(halves - std::round(halves)) > 1e-9)
                invalid("--chunk-ksnp must be a positive multiple of 0.5 within the supported integer range.");
            a.nSnpPerChunk = static_cast<int>(std::llround(halves) * 512);
        }
        else if (arg == "--seed") a.seed = parseSeed(next(), arg);
        else if (arg == "--extract") a.extractFile = next();
        else if (arg == "--exclude") a.excludeFile = next();
        else if (arg == "--chr") a.chrSpec = next();
        else if (arg == "--keep") a.keepFile = next();
        else if (arg == "--remove") a.removeFile = next();
        else if (arg == "--compression") a.compression = next();
        else if (arg == "--compression-level") { a.compressionLevel = parseInt(next(), arg); a.compressionLevelExplicit = true; }
        else if (arg == "--cal-af-coef") { markSeen(arg); a.calAfCoef = true; }
        else if (arg == "--int-pheno") { markSeen(arg); a.intPheno = true; }
        else invalid("Unknown option '" + arg + "'. Run 'SPAGxE --help options'.");
    }
    a.explicitFlags = std::move(seen);
    return a;
}
} // namespace cli
