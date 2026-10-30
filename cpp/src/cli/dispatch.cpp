// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE standalone validation and dispatch.
// Adapted from GRAB CLI infrastructure; see THIRD_PARTY_NOTICES.
#include "util/version.hpp"
#include "cli/cli.hpp"
#include "cli/cct_combine.hpp"
#include "genotype/genotype_source.hpp"
#include "io/subject_data.hpp"
#include "local/local_workflow.hpp"
#include "local/local_pack.hpp"
#include "interaction/interaction_test.hpp"
#include "allele_model/frequency_model.hpp"
#include "util/int_pheno.hpp"
#include "util/logging.hpp"
#include "util/null_model.hpp"
#include "util/text_stream.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace cli {
namespace {
void need(bool condition, const std::string &message) {
    if (!condition) throw std::invalid_argument(message);
}
void require(const std::string &value, const std::string &flag, const std::string &context) {
    need(!value.empty(), flag + " is required for " + context + ".");
}
std::string lower(std::string value) {
    for (auto &ch : value) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return value;
}
std::string trim(const std::string &value) {
    auto first = value.find_first_not_of(" \t\r\n");
    return first == std::string::npos ? "" : value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}
std::vector<std::string> splitComma(const std::string &value, const std::string &flag) {
    std::vector<std::string> result;
    std::unordered_set<std::string> seen;
    size_t begin = 0;
    do {
        auto end = value.find(',', begin);
        std::string token = trim(value.substr(begin, end == std::string::npos ? end : end - begin));
        need(!token.empty(), flag + " contains an empty list entry.");
        need(seen.insert(token).second, flag + " contains duplicate entry '" + token + "'.");
        result.push_back(token);
        if (end == std::string::npos) break;
        begin = end + 1;
    } while (true);
    return result;
}
std::string joinComma(const std::vector<std::string> &values) {
    std::string result;
    for (const auto &value : values) { if (!result.empty()) result += ','; result += value; }
    return result;
}
std::unordered_set<std::string> parseChrSpec(const std::string &value) {
    std::unordered_set<std::string> result;
    for (const auto &token : splitComma(value, "--chr")) {
        auto dash = token.find('-');
        if (dash == std::string::npos) { result.insert(token); continue; }
        auto allDigits = [](const std::string &s) {
            return !s.empty() && std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isdigit(c); });
        };
        auto left = token.substr(0, dash), right = token.substr(dash + 1);
        need(allDigits(left) && allDigits(right), "Invalid --chr numeric range '" + token + "'.");
        const long long lo = std::stoll(left), hi = std::stoll(right);
        need(lo <= hi && hi <= std::numeric_limits<int>::max() && hi - lo <= 10000,
             "Invalid or excessive --chr range '" + token + "'.");
        for (long long chr = lo; chr <= hi; ++chr) result.insert(std::to_string(chr));
    }
    return result;
}
GenotypeInput resolveGeno(const Args &args, const std::string &context) {
    int count = !args.bfilePrefix.empty() + !args.pfilePrefix.empty() + !args.vcfFile.empty() +
                !args.bcfFile.empty() + !args.bgenFile.empty();
    need(count == 1, "Exactly one of --bfile, --pfile, --vcf, --bcf, or --bgen is required for " + context + ".");
    GenotypeInput result;
    result.extractFile = args.extractFile;
    result.excludeFile = args.excludeFile;
    result.hardCallThreshold = args.hardCallThreshold;
    if (!args.chrSpec.empty()) result.chrFilter = parseChrSpec(args.chrSpec);
    if (!args.bfilePrefix.empty()) { result.format = GenotypeFormat::Plink; result.path = args.bfilePrefix; }
    else if (!args.pfilePrefix.empty()) { result.format = GenotypeFormat::Pgen; result.path = args.pfilePrefix; }
    else if (!args.vcfFile.empty()) { result.format = GenotypeFormat::Vcf; result.path = args.vcfFile; }
    else if (!args.bcfFile.empty()) { result.format = GenotypeFormat::Bcf; result.path = args.bcfFile; }
    else {
        result.format = GenotypeFormat::Bgen;
        result.path = args.bgenFile;
        result.bgenAltFirst = args.bgenRefMode != "ref-first";
        if (args.bgenRefMode == "ref-unknown")
            warnMsg("BGEN ref-unknown: the last allele is provisional REF; the first allele is the tested ALT.");
    }
    return result;
}
void validateCommon(Args &args) {
    need(args.nthread >= 1, "--threads must be at least 1.");
    need(args.nSnpPerChunk >= 512 && args.nSnpPerChunk % 512 == 0,
         "--chunk-ksnp must be a positive multiple of 0.5.");
    need(args.missingCutoff >= 0 && args.missingCutoff <= 1, "--geno must be in [0,1].");
    need(args.minMafCutoff >= 0 && args.minMafCutoff <= 0.5, "--maf must be in [0,0.5].");
    need(args.minMacCutoff >= 0, "--mac must be non-negative.");
    need(args.hweCutoff >= 0 && args.hweCutoff <= 1, "--hwe must be in [0,1].");
    need(args.hardCallThreshold >= 0 && args.hardCallThreshold < 0.5, "--hard-call-threshold must be in [0,0.5).");
    need(args.spagxeMarginalCutoff >= 0 && args.spagxeMarginalCutoff <= 1,
         "--spagxe-marginal-cutoff must be in [0,1].");
    need(args.spaCutoff >= 0.01, "--spa-z-threshold must be at least 0.01.");
    need(args.grmTableFile.empty() || args.spGrmPlink2File.empty(),
         "--grm-file and --sp-grm-plink2 are mutually exclusive.");
    need(args.localMaxBlockBytes >= 4096, "--local-max-block-bytes must be at least 4096.");
    const bool local = lower(args.method) == "spagxelocal" || lower(args.method) == "spagxemix_local" || args.packLocal;
    const bool cct = lower(args.method) == "spagxecct";
    need(args.packLocal || (args.localBimFile.empty() && args.localFamFile.empty() &&
                          !args.localBlockVariantsExplicit && !args.localZstdLevelExplicit),
         "--local-bim, --local-fam, --local-block-variants and --local-zstd-level require --pack-local.");
    need(args.localZstdLevel >= -5 && args.localZstdLevel <= 22, "--local-zstd-level must be in -5..22.");
    need(cct || (args.cctGlobal.empty() && args.cctLocal.empty() && args.cctPhenotype.empty() &&
                 !args.cctWeightsExplicit && !args.cctPerAncestryEqual),
         "--cct-* options require --method SPAGxEcct.");
    need(local || (args.localManifest.empty() && args.localPrefix.empty() && args.localBgenManifest.empty() && !args.localMaxBlockBytesExplicit),
         "--local-manifest, --local-prefix and --local-max-block-bytes require --method SPAGxElocal.");
    if (local) {
        const int sources = !args.localManifest.empty() + !args.localPrefix.empty() + !args.localBgenManifest.empty();
        need(sources == 1, "SPAGxEmix_local requires one of --local-manifest, --local-prefix or --local-bgen-manifest.");
        need(!args.localMaxBlockBytesExplicit || !args.localPrefix.empty() || !args.localBgenManifest.empty() || !args.localManifest.empty(),
             "--local-max-block-bytes applies to local binary/BGEN decoded-block memory.");
        need(args.bfilePrefix.empty() && args.pfilePrefix.empty() && args.vcfFile.empty() && args.bcfFile.empty() && args.bgenFile.empty(),
             "SPAGxElocal uses paired ancestry G/H input and does not accept --bfile/--pfile/--vcf/--bcf/--bgen.");
        need(args.grmTableFile.empty() && args.spGrmPlink2File.empty(),
             "SPAGxElocal does not implement local-ancestry relatedness correction; GRM input is not supported.");
        need(!args.pcColsExplicit, "SPAGxElocal does not fit individual AF from --pc-cols; use --covar-name for null-model PCs.");
        need(!args.hardCallThresholdExplicit, "SPAGxElocal does not use --hard-call-threshold; its HWE diagnostic uses summed ancestry genotypes.");
    }
    need(args.compression.empty() || args.compression == "gz" || args.compression == "zst",
         "--compression must be gz or zst.");
    if (args.compressionLevelExplicit) {
        need(!args.compression.empty(), "--compression-level requires --compression gz|zst.");
        need(args.compressionLevel >= 1 && args.compressionLevel <= (args.compression == "gz" ? 9 : 22),
             "--compression-level must be in 1..9 for gz or 1..22 for zst.");
    } else if (!args.compression.empty()) args.compressionLevel = args.compression == "gz" ? 6 : 3;
}
std::string codecSuffix(const Args &args) {
    return args.compression.empty() ? "" : "." + args.compression;
}
void logCommand(int argc, char *argv[], const Args &args) {
    infoMsg("SPAGxE " SPAGXE_VERSION " | Yuzhuo Ma | yuzhuoma@amss.ac.cn");
    infoMsg("SPAGxE algorithms: Yuzhuo Ma. Incorporated components: THIRD_PARTY_NOTICES.md.");
    std::cerr << "Command:";
    for (int i = 0; i < argc; ++i) std::cerr << ' ' << std::quoted(std::string(argv[i]));
    std::cerr << '\n';
    if (!args.method.empty()) infoMsg("Analysis method: %s", args.method.c_str());
    if (!args.grmTableFile.empty() || !args.spGrmPlink2File.empty())
        infoMsg("GRM-corrected mode: %s; SPA variance adjustment: score",
                "SPAGxE+");
    infoMsg("Threads: %d; markers/chunk: %d; QC: missing<=%g, MAF>=%g, MAC>=%g, HWE P>=%g",
            args.nthread, args.nSnpPerChunk, args.missingCutoff, args.minMafCutoff, args.minMacCutoff, args.hweCutoff);
    infoMsg("CGF: guarded partial normal (IQR 1.5, guard 0.001); normal/SPA threshold: %g", args.spaCutoff);
    if (args.seed) warnMsg("--seed is retained for CLI compatibility; SPAGxE analysis is deterministic and does not use random numbers.");
}
} // namespace

int run(int argc, char *argv[]) {
    if (argc < 2) { printHelp("__short__"); return 1; }
    const auto wallStart = std::chrono::steady_clock::now();
    const auto cpuStart = std::clock();
    try {
        Args args = parseArgs(argc, argv);
        if (args.showVersion) { printVersion(); return 0; }
        if (!args.helpTopic.empty()) { printHelp(args.helpTopic); return 0; }
        validateCommon(args);
        need(int(args.calAfCoef) + int(args.intPheno) + int(args.packLocal) <= 1,
             "--cal-af-coef, --int-pheno and --pack-local are mutually exclusive.");
        need(!(args.calAfCoef || args.intPheno || args.packLocal) || args.method.empty(), "Utility modes cannot be combined with --method.");
        auto columns = [](const std::string &s, const std::string &flag) {
            return s.empty() ? std::vector<std::string>{} : splitComma(s, flag);
        };
        auto pcNames = columns(args.pcCols, "--pc-cols");
        auto covarNames = columns(args.covarName, "--covar-name");
        auto phenoNames = columns(args.phenoName, "--pheno-name");
        auto residNames = columns(args.residName, "--resid-name");
        require(args.outPrefix, "--out", "all modes");
        if (args.packLocal) {
            const std::unordered_set<std::string> allowed = {"--pack-local", "--local-manifest", "--local-prefix", "--local-binary", "--local-bgen-manifest", "--local-max-block-bytes", "--local-bim", "--local-fam", "--local-block-variants", "--local-zstd-level", "--out"};
            for (const auto &flag : args.explicitFlags)
                need(allowed.count(flag) != 0, "--pack-local does not accept analysis option " + flag + ".");
            spagxe_local::PackOptions options;
            options.input.manifest = args.localManifest;
            options.input.binaryPrefix = args.localPrefix;
            options.input.bgenManifest = args.localBgenManifest;
            options.input.maxBlockBytes = args.localMaxBlockBytes;
            options.outPrefix = args.outPrefix;
            options.bimFile = args.localBimFile; options.famFile = args.localFamFile;
            options.blockVariants = args.localBlockVariants;
            options.compressionLevel = args.localZstdLevel;
            args.method = "pack-local"; logCommand(argc, argv, args);
            spagxe_local::packLocalCounts(options);
            infoMsg("Packed hardcall counts: %s.localbed/.localbim/.localfam", args.outPrefix.c_str());
        } else if (lower(args.method) == "spagxecct") {
            const std::unordered_set<std::string> allowed = {"--method", "--cct-global", "--cct-local", "--cct-phenotype",
                "--cct-weights", "--cct-per-ancestry-equal", "--envir-name", "--out", "--compression", "--compression-level"};
            for (const auto &flag : args.explicitFlags)
                need(allowed.count(flag) != 0, "SPAGxEcct does not accept analysis option " + flag + ".");
            require(args.cctGlobal, "--cct-global", "SPAGxEcct");
            require(args.cctLocal, "--cct-local", "SPAGxEcct");
            require(args.cctPhenotype, "--cct-phenotype", "SPAGxEcct");
            require(args.envName, "--envir-name", "SPAGxEcct");
            need(!(args.cctWeightsExplicit && args.cctPerAncestryEqual),
                 "--cct-weights and --cct-per-ancestry-equal are mutually exclusive.");
            const auto comma = args.cctWeights.find(',');
            need(comma != std::string::npos && args.cctWeights.find(',', comma+1) == std::string::npos,
                 "--cct-weights requires two positive group weights: GLOBAL,LOCAL.");
            const std::array<std::string,2> weights = {trim(args.cctWeights.substr(0,comma)), trim(args.cctWeights.substr(comma+1))};
            auto weight = [](const std::string &s) {
                size_t consumed = 0; double value;
                try { value = std::stod(s, &consumed); }
                catch (...) { throw std::invalid_argument("--cct-weights requires finite positive numeric weights."); }
                need(consumed == s.size() && std::isfinite(value) && value > 0,
                     "--cct-weights requires finite positive numeric weights.");
                return value;
            };
            spagxe_cct::Options options;
            options.globalFile = args.cctGlobal; options.localFile = args.cctLocal;
            options.phenotype = args.cctPhenotype; options.outPrefix = args.outPrefix;
            options.envNames = splitComma(args.envName, "--envir-name");
            options.globalWeight = weight(weights[0]); options.localWeight = weight(weights[1]);
            options.perAncestryEqual = args.cctPerAncestryEqual;
            options.compression = args.compression; options.compressionLevel = args.compressionLevel;
            args.method = "SPAGxEcct"; logCommand(argc, argv, args);
            spagxe_cct::run(options);
        } else if (args.intPheno) {
            require(args.phenoFile, "--pheno", "--int-pheno");
            need(args.bfilePrefix.empty() && args.pfilePrefix.empty() && args.vcfFile.empty() && args.bcfFile.empty() && args.bgenFile.empty(),
                 "--int-pheno does not accept genotype input.");
            need(args.covarFile.empty() && args.covarName.empty() && args.residName.empty() && args.envName.empty() &&
                 args.regressionModel.empty() && !args.saveResid && !args.pcColsExplicit &&
                 args.grmTableFile.empty() && args.spGrmPlink2File.empty() && args.compression.empty(),
                 "--int-pheno accepts --pheno, --pheno-name and --out; incompatible analysis or compression options were supplied.");
            logCommand(argc, argv, args);
            TextWriter::assertWritable(args.outPrefix + ".txt");
            runIntPheno(args.phenoFile, args.outPrefix + ".txt", phenoNames);
        } else if (args.calAfCoef) {
            need(!pcNames.empty(), "--cal-af-coef requires --pc-cols.");
            need(!args.phenoFile.empty() || !args.covarFile.empty(), "--cal-af-coef requires --pheno or --covar containing PC columns.");
            need(args.residName.empty() && args.phenoName.empty() && args.envName.empty() && args.regressionModel.empty() &&
                 !args.saveResid && args.grmTableFile.empty() && args.spGrmPlink2File.empty(),
                 "--cal-af-coef does not accept phenotype fitting, environment, or GRM options.");
            need(!args.compressionLevelExplicit, "--cal-af-coef uses the codec default level; --compression-level applies to association output.");
            auto geno = resolveGeno(args, "--cal-af-coef");
            logCommand(argc, argv, args);
            auto path = args.outPrefix + ".afc" + codecSuffix(args);
            TextWriter::assertWritable(path);
            exportFrequencyModels(pcNames, args.phenoFile, args.covarFile, geno, path, args.nthread,
                        args.nSnpPerChunk, args.missingCutoff, args.minMafCutoff, args.minMacCutoff,
                        args.hweCutoff, args.keepFile, args.removeFile);
        } else {
            const auto method = lower(args.method.empty() ? "SPAGxE" : args.method);
            need(method == "spagxe" || method == "spagxe+" || method == "spagxemix" || method == "spagxelocal" || method == "spagxemix_local",
                 "--method must be SPAGxE, SPAGxE+, SPAGxEmix or SPAGxEmix_local.");
            if (method == "spagxe+") need(!args.grmTableFile.empty() || !args.spGrmPlink2File.empty(), "SPAGxE+ requires a sparse GRM.");
            args.method = (method == "spagxe" || method == "spagxe+") ? "SPAGxE" : (method == "spagxemix" ? "SPAGxEmix" : "SPAGxEmix_local");
            require(args.phenoFile, "--pheno", args.method);
            require(args.envName, "--envir-name", args.method);
            auto envNames = splitComma(args.envName, "--envir-name");
            const bool local = args.method == "SPAGxEmix_local";
            GenotypeInput geno;
            if (!local) geno = resolveGeno(args, args.method);
            need(args.phenoName.empty() || args.residName.empty(), "--pheno-name and --resid-name are mutually exclusive.");
            need(args.method == "SPAGxEmix" || !args.pcColsExplicit,
                 "--pc-cols is for SPAGxEmix / --cal-af-coef. For null-model PC adjustment use --covar-name.");
            if (args.method == "SPAGxEmix") {
                need(!pcNames.empty(), "SPAGxEmix requires --pc-cols.");
                need(args.grmTableFile.empty() && args.spGrmPlink2File.empty(),
                     "SPAGxEmix supports unrelated samples; GRM input is not supported.");
            }
            // Materialize the documented all-covariates default, so the engine
            // can verify that each environment (and mix PC) enters the null model.
            if (covarNames.empty() && !args.covarFile.empty()) covarNames = SubjectData::readColumnNames(args.covarFile);
            std::string effectiveCovar = args.covarFile.empty() && !covarNames.empty() ? args.phenoFile : args.covarFile;
            if (args.phenoName.empty() && args.residName.empty()) {
                auto names = SubjectData::readColumnNames(args.phenoFile);
                std::unordered_set<std::string> excluded(covarNames.begin(), covarNames.end());
                excluded.insert(envNames.begin(), envNames.end());
                if (args.method == "SPAGxEmix") excluded.insert(pcNames.begin(), pcNames.end());
                for (const auto &name : names) if (!excluded.count(name)) phenoNames.push_back(name);
                need(!phenoNames.empty(), "No phenotype columns remain after excluding covariates, environments and PCs; specify --pheno-name or --resid-name.");
                args.phenoName = joinComma(phenoNames);
                infoMsg("Automatically selected phenotype columns: %s", args.phenoName.c_str());
            }
            auto model = nullmodel::parseRegressionModel(args.regressionModel);
            if (!args.phenoName.empty()) {
                for (const auto &environment : envNames)
                    need(std::find(covarNames.begin(), covarNames.end(), environment) != covarNames.end(),
                         "Environment '" + environment + "' must appear in --covar-name for the null model.");
                if (args.method == "SPAGxEmix") for (const auto &pc : pcNames)
                    need(std::find(covarNames.begin(), covarNames.end(), pc) != covarNames.end(),
                         "PC '" + pc + "' must appear in --covar-name for SPAGxEmix.");
                nullmodel::parsePhenoSpecList(model, args.phenoName);
            } else {
                need(model == nullmodel::RegressionModel::Auto, "A non-auto --regression-model requires --pheno-name.");
                need(!args.saveResid, "--save-resid requires --pheno-name.");
            }
            logCommand(argc, argv, args);
            TextWriter::assertWritable(args.outPrefix + ".write_test");
            if (local) {
                spagxe_local::LocalWorkflowOptions options;
                options.input.manifest = args.localManifest;
                options.input.binaryPrefix = args.localPrefix;
                options.input.bgenManifest = args.localBgenManifest;
                options.input.maxBlockBytes = args.localMaxBlockBytes;
                options.phenoFile = args.phenoFile; options.covarFile = effectiveCovar;
                options.phenoNameSpec = args.phenoName; options.regressionModel = args.regressionModel;
                options.residNames = residNames; options.envNames = envNames; options.covarNames = covarNames;
                options.keepFile = args.keepFile; options.removeFile = args.removeFile;
                options.extractFile = args.extractFile; options.excludeFile = args.excludeFile;
                if (!args.chrSpec.empty()) options.chromosomes = parseChrSpec(args.chrSpec);
                options.outPrefix = args.outPrefix; options.compression = args.compression;
                options.compressionLevel = args.compressionLevel; options.nthreads = args.nthread;
                options.nSnpPerChunk = args.nSnpPerChunk; options.marginalCutoff = args.spagxeMarginalCutoff;
                options.spaCutoff = args.spaCutoff;
                options.missingCutoff = args.missingCutoff; options.minMafCutoff = args.minMafCutoff;
                options.minMacCutoff = args.minMacCutoff; options.hweCutoff = args.hweCutoff;
                options.saveResid = args.saveResid;
                spagxe_local::scanHomogeneousInteractionsLocal(options);
            } else if (args.method == "SPAGxE")
                scanHomogeneousInteractions(args.phenoFile, residNames, envNames, args.grmTableFile,
                          args.spGrmPlink2File, geno, args.outPrefix, args.compression,
                          args.compressionLevel, args.spagxeMarginalCutoff, args.spaCutoff,
                          args.nthread, args.nSnpPerChunk, args.missingCutoff,
                          args.minMafCutoff, args.minMacCutoff, args.hweCutoff, args.keepFile,
                          args.removeFile, args.regressionModel, args.phenoName, effectiveCovar,
                          covarNames, args.saveResid);
            else
                scanAdmixedInteractions(args.phenoFile, residNames, envNames, pcNames, geno,
                             args.outPrefix, args.compression, args.compressionLevel,
                             args.spagxeMarginalCutoff, args.spaCutoff,
                             args.nthread, args.nSnpPerChunk, args.missingCutoff, args.minMafCutoff,
                             args.minMacCutoff, args.hweCutoff, args.keepFile, args.removeFile,
                             args.regressionModel, args.phenoName, effectiveCovar, covarNames,
                             args.saveResid);
        }
        const auto wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wallStart).count();
        const auto cpu = static_cast<double>(std::clock() - cpuStart) / CLOCKS_PER_SEC;
        infoMsg("Wall time: %.3f seconds; CPU time: %.3f seconds", wall, cpu);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "[ERROR] " << error.what() << '\n';
        return 1;
    }
}
} // namespace cli
