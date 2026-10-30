// Per-individual allele-frequency fitting, storage and reconstruction.
// One minor-oriented paper cascade is shared by online and exported models.
#include "allele_model/frequency_model.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "util/probability_math.hpp"
#include "util/text_scanner.hpp"
#include "util/text_stream.hpp"

namespace {
bool validAF(const Eigen::Ref<const Eigen::VectorXd> &af) {
    return af.allFinite() && !(af.array() < 0.0).any() &&
           !(af.array() > 1.0).any() &&
           (af.array() * (1.0 - af.array())).sum() > 0.0;
}

// g and frequency are already oriented to the cohort's minor allele.
// Keep fitted OLS values in out to avoid a second matrix multiply online.
AFModel fitMinorAF(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    double frequency,
    const AFContext &ctx,
    Eigen::Ref<Eigen::VectorXd> out,
    bool complementOutput = false
) {
    constexpr double pcPValueCutoff = 0.05;
    constexpr double outsideRatioCutoff = 0.1;
    constexpr double macCutoff = 20.0;
    const int n = ctx.N, p = ctx.nPC;
    if (n <= 0 || p < 0 || g.size() != n || out.size() != n ||
        ctx.onePlusPCs.rows() != n || ctx.onePlusPCs.cols() != p + 1 ||
        ctx.XtX_inv_Xt.rows() != p + 1 || ctx.XtX_inv_Xt.cols() != n ||
        ctx.sqrt_XtX_inv_diag.size() != p + 1 ||
        ctx.PCs.rows() != n || ctx.PCs.cols() != p ||
        !g.allFinite() || (g.array() < 0).any() || (g.array() > 2).any() ||
        !std::isfinite(frequency) || frequency < 0 || frequency > 0.5)
        throw std::runtime_error("AF fitting requires finite minor-oriented diploid G, frequency in [0,0.5], and matching design dimensions");

    AFModel model;
    model.status = 0;
    model.betas = Eigen::VectorXd::Zero(p + 1);
    auto pooled = [&]() {
        model.status = 0;
        model.betas.setZero();
        out.setConstant(complementOutput ? 1.0 - frequency : frequency);
    };
    auto failed = [&]() {
        model.status = -2;
        model.betas.setZero();
        out.setConstant(std::numeric_limits<double>::quiet_NaN());
    };
    // Count the actual minor-oriented, mean-imputed vector. This is
    // mathematically 2*N*frequency, without cancellation in 1-ALT_FREQ
    // spuriously turning a high-ALT hardcall MAC of 20 into 20+roundoff.
    if (g.sum() <= macCutoff || n <= p + 1) {
        pooled();
        return model;
    }

    const Eigen::VectorXd coef = ctx.XtX_inv_Xt * g;
    out.noalias() = ctx.onePlusPCs * coef;
    if (!out.allFinite() || !coef.allFinite()) {
        failed();
        return model;
    }
    int outside = 0;
    for (int i = 0; i < n; ++i)
        if (out[i] < 0.0 || out[i] > 2.0) ++outside;

    if (static_cast<double>(outside) / n < outsideRatioCutoff) {
        model.status = 1;
        model.betas = coef;
        out *= 0.5;
        for (int i = 0; i < n; ++i) out[i] = std::max(0.0, std::min(1.0, out[i]));
        if (complementOutput) out = 1.0 - out.array();
    } else {
        const double df = static_cast<double>(n - p - 1);
        const double sqrtS2 = std::sqrt((g - out).squaredNorm() / df);
        std::vector<int> selected;
        for (int j = 0; j < p; ++j) {
            const double se = ctx.sqrt_XtX_inv_diag[j + 1] * sqrtS2;
            const double t = coef[j + 1] / se;
            if (2.0 * math::pt(-std::abs(t), df) < pcPValueCutoff)
                selected.push_back(j);
        }
        if (selected.empty()) {
            pooled();
            return model;
        }

        Eigen::VectorXd carriers = Eigen::VectorXd::Zero(n);
        int count = 0;
        for (int i = 0; i < n; ++i)
            if (g[i] > 0.5) { carriers[i] = 1.0; ++count; }
        // The minor MAC gate was already applied. Carrier count <=20 is not
        // a second MAC definition; only constant binary outcomes cannot fit.
        if (count == 0 || count == n) {
            pooled();
            return model;
        }
        Eigen::MatrixXd selectedPCs(n, selected.size());
        for (size_t j = 0; j < selected.size(); ++j)
            selectedPCs.col(j) = ctx.PCs.col(selected[j]);
        bool converged = false;
        const Eigen::VectorXd beta = math::logisticRegressionBeta(
            selectedPCs, carriers, 1e-6, 100, &converged);
        if (!converged || !beta.allFinite()) {
            failed();
            return model;
        }
        model.status = 2;
        model.betas[0] = beta[0];
        for (size_t j = 0; j < selected.size(); ++j)
            model.betas[selected[j] + 1] = beta[j + 1];
        if (!model.betas.allFinite()) {
            failed();
            return model;
        }
        if (complementOutput) model.status = 3;
        getAFVecFromModel(model, frequency, ctx.onePlusPCs, n, out);
        model.status = 2;
    }
    if (!validAF(out)) failed();
    return model;
}
} // namespace

void prepareAFLinearDesign(
    const Eigen::MatrixXd &x,
    Eigen::MatrixXd &inverseTimesX,
    Eigen::VectorXd &sqrtDiag
) {
    if (!x.allFinite() || x.cols() < 1 || x.rows() <= x.cols())
        throw std::runtime_error("AF PC design must be finite with more subjects than columns");
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(x);
    if (qr.rank() != x.cols())
        throw std::runtime_error("AF PC design is rank deficient; remove duplicate or constant PCs");
    const Eigen::Index p = x.cols();
    const Eigen::MatrixXd r = qr.matrixR().topLeftCorner(p, p);
    const Eigen::MatrixXd inverseR = r.triangularView<Eigen::Upper>().solve(
        Eigen::MatrixXd::Identity(p, p));
    const Eigen::MatrixXd left = qr.colsPermutation() * inverseR;
    const Eigen::MatrixXd thinQ = qr.householderQ() *
        Eigen::MatrixXd::Identity(x.rows(), p);
    inverseTimesX = left * thinQ.transpose();
    const Eigen::MatrixXd inverse = left * left.transpose();
    sqrtDiag = inverse.diagonal().cwiseMax(0.0).cwiseSqrt();
    if (!inverseTimesX.allFinite() || !sqrtDiag.allFinite())
        throw std::runtime_error("AF PC regression design produced non-finite values");
}

AFModel computeAFModel(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    double altFreq,
    const AFContext &ctx
) {
    if (!std::isfinite(altFreq) || altFreq < 0 || altFreq > 1)
        throw std::runtime_error("AF model ALT frequency must be in [0,1]");
    if (ctx.N <= 0)
        throw std::runtime_error("AF model requires at least one sample");
    const bool flip = altFreq > 0.5;
    Eigen::VectorXd fitted(ctx.N);
    Eigen::VectorXd minor;
    if (flip) minor = 2.0 - g.array();
    AFModel model = flip ? fitMinorAF(minor, 1.0 - altFreq, ctx, fitted, true)
                         : fitMinorAF(g, altFreq, ctx, fitted);
    if (flip && model.status == 1) {
        // Store ordinary ALT-scale OLS coefficients: 2 - X*beta_minor.
        model.betas = -model.betas;
        model.betas[0] += 2.0;
    } else if (flip && model.status == 2) {
        // Complementing this nonlinear AF transform cannot be encoded merely
        // by negating logistic coefficients. Preserve minor betas explicitly.
        model.status = 3;
    }
    return model;
}

bool computeAFVec(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    double minorFreq,
    const AFContext &ctx,
    Eigen::Ref<Eigen::VectorXd> out,
    bool complementOutput
) {
    return fitMinorAF(g, minorFreq, ctx, out, complementOutput).status != -2;
}

void getAFVecFromModel(
    const AFModel &model,
    double altFreq,
    const Eigen::MatrixXd &onePlusPCs,
    int n,
    Eigen::Ref<Eigen::VectorXd> out
) {
    if (model.status == -2)
        throw std::runtime_error("AF model fitting failed (status -2); no individual AF is available");
    if (model.status == -1)
        throw std::runtime_error("AF model is not fitted or was filtered by QC (status -1)");
    if (model.status < 0 || model.status > 3)
        throw std::runtime_error("Unknown AF model status: " + std::to_string(int(model.status)));
    if (n <= 0 || out.size() != n || !std::isfinite(altFreq) || altFreq < 0 || altFreq > 1)
        throw std::runtime_error("AF reconstruction requires matching sample count and ALT frequency in [0,1]");
    if (model.status == 0) {
        out.setConstant(altFreq);
        return;
    }
    if (onePlusPCs.rows() != n || onePlusPCs.cols() != model.betas.size() ||
        !onePlusPCs.allFinite() || !model.betas.allFinite())
        throw std::runtime_error("AF reconstruction requires finite coefficients and a matching finite PC design");
    out.noalias() = onePlusPCs * model.betas;
    if (model.status == 1) {
        out *= 0.5;
        for (int i = 0; i < n; ++i) out[i] = std::max(0.0, std::min(1.0, out[i]));
    } else {
        for (int i = 0; i < n; ++i) {
            const double eta = out[i];
            // Under HWE, 1-AF = sqrt(1-sigmoid(eta)). Work on the
            // log scale so neither sigmoid rounding to one nor subtracting
            // a tiny minor AF destroys a representable allele frequency.
            const double softplus = eta >= 0.0
                ? eta + std::log1p(std::exp(-eta))
                : std::log1p(std::exp(eta));
            out[i] = model.status == 3 ? std::exp(-0.5 * softplus)
                                       : -std::expm1(-0.5 * softplus);
        }
    }
}

// Storage layout: status byte plus (1+nPC) doubles; see header for statuses.
static IndivAFWriter::Mode inferMode(const std::string &path) {
    const auto n = path.size();
    if (n > 4 && path.compare(n - 4, 4, ".bin") == 0) return IndivAFWriter::Mode::Binary;
    return IndivAFWriter::Mode::Text;
}

// ======================================================================
// IndivAFWriter — constructor
// ======================================================================

IndivAFWriter::IndivAFWriter(
    const std::string &outputFile,
    uint64_t nBimMarkers,
    int nPC
)
    : m_mode(inferMode(outputFile)),
      m_nPC(nPC),
      m_recordSize(static_cast<long long>(sizeof(int8_t)) +
                   static_cast<long long>(1 + nPC) * static_cast<long long>(sizeof(double)))
{
    if (m_mode == Mode::Binary) {
        // Each raw slot starts unfitted; skipped/filtered slots must never
        // masquerade as a completed pooled-frequency fit (status 0).
        const long long totalSize = static_cast<long long>(nBimMarkers) * m_recordSize;

        {
            std::ofstream pre(outputFile, std::ios::binary | std::ios::trunc);
            if (!pre.is_open()) throw std::runtime_error("Cannot create binary AF file: " + outputFile);
            const long long recordsPerChunk = std::max(1LL, (1LL << 20) / m_recordSize);
            std::vector<char> zeros(static_cast<size_t>(recordsPerChunk * m_recordSize), 0);
            for (long long i = 0; i < recordsPerChunk; ++i)
                zeros[static_cast<size_t>(i * m_recordSize)] = static_cast<char>(-1);
            for (long long rem = totalSize; rem > 0;) {
                const std::streamsize n =
                    static_cast<std::streamsize>(std::min<long long>(rem, static_cast<long long>(zeros.size())));
                pre.write(zeros.data(), n);
                rem -= n;
            }
        }

        m_binOut.open(outputFile, std::ios::binary | std::ios::in | std::ios::out);
        if (!m_binOut.is_open()) throw std::runtime_error("Cannot open binary AF file for random write: " + outputFile);

    } else {
        m_writer = std::make_unique<TextWriter>(outputFile);

        // Write TSV header
        std::string hdr = "#STATUS";
        for (int j = 0; j <= nPC; ++j)
            hdr += "\tBETA" + std::to_string(j);
        hdr += "\n";
        m_writer->write(hdr);
        m_writer->write("#AFMODEL_SCHEMA\t2\t-2=AF_fit_failed\t-1=unfitted_or_QC_filtered\t0=pooled_ALT\t1=ALT_OLS\t2=logistic\t3=complement_logistic\n");
    }
}

// ======================================================================
// IndivAFWriter — write one record
// ======================================================================

void IndivAFWriter::write(
    uint64_t genoIndex,
    int8_t status,
    const Eigen::VectorXd &betas
) {
    if (status < -2 || status > 3 || betas.size() != 1 + m_nPC || !betas.allFinite())
        throw std::runtime_error("Cannot write AF model: invalid status or coefficient vector");
    if (m_mode == Mode::Binary) {
        const long long filePos = static_cast<long long>(genoIndex) * m_recordSize;
        m_binOut.seekp(filePos);
        m_binOut.write(reinterpret_cast<const char *>(&status), sizeof(int8_t));
        m_binOut.write(
            reinterpret_cast<const char *>(betas.data()),
            static_cast<std::streamsize>((1 + m_nPC) * sizeof(double))
        );

    } else {
        char buf[64];
        std::string line;
        line.reserve(64 + 26 * (1 + m_nPC));

        std::snprintf(buf, sizeof(buf), "%d", static_cast<int>(status));
        line += buf;
        for (int j = 0; j <= m_nPC; ++j) {
            line += '\t';
            std::snprintf(buf, sizeof(buf), "%.17g", betas[j]);
            line += buf;
        }
        line += '\n';

        m_writer->write(line);
    }
}

// ======================================================================
// IndivAFWriter — close
// ======================================================================

void IndivAFWriter::close() {
    if (m_closed) return;
    m_closed = true;

    if (m_mode == Mode::Binary && m_binOut.is_open())m_binOut.close();
    else if (m_writer)m_writer->close();
}

IndivAFWriter::~IndivAFWriter()
{
    close();
}

// ======================================================================
// IndivAFReader
// ======================================================================

IndivAFReader::IndivAFReader(
    const std::string &binFile,
    int nPC
)
    : m_nPC(nPC),
      m_recordSize(static_cast<long long>(sizeof(int8_t)) +
                   static_cast<long long>(1 + nPC) * static_cast<long long>(sizeof(double)))
{
    m_in.open(binFile, std::ios::binary);
    if (!m_in.is_open()) throw std::runtime_error("Cannot open binary AF file: " + binFile);
}

IndivAFReader::~IndivAFReader()
{
    m_in.close();
}

bool IndivAFReader::read(
    uint64_t genoIndex,
    AFModel &model
) {
    const long long filePos = static_cast<long long>(genoIndex) * m_recordSize;
    m_in.seekg(filePos);
    if (!m_in.good()) throw std::runtime_error("Binary AF seek failed at genoIndex=" + std::to_string(genoIndex));

    int8_t st = 0;
    m_in.read(reinterpret_cast<char *>(&st), sizeof(int8_t));
    model.status = st;
    model.betas.resize(1 + m_nPC);
    m_in.read(reinterpret_cast<char *>(model.betas.data()), static_cast<std::streamsize>((1 + m_nPC) * sizeof(double)));
    if (!m_in.good()) throw std::runtime_error("Binary AF read failed at genoIndex=" + std::to_string(genoIndex));
    if (st < -2 || st > 3 || !model.betas.allFinite())
        throw std::runtime_error("Invalid binary AF model status/coefficients at genoIndex=" + std::to_string(genoIndex));
    return st > 0;
}

// ======================================================================
// ======================================================================
// loadAFModels — helpers (binary + text/gz) and dispatcher
// ======================================================================

static std::vector<AFModel>loadAFModelsBinary(
    const std::string &path,
    int nPC,
    uint32_t nMarkers,
    const std::vector<
        uint64_t> &genoIndices
) {
    std::vector<AFModel> models(nMarkers);
    IndivAFReader reader(path, nPC);
    for (uint32_t fi = 0; fi < nMarkers; ++fi)
        reader.read(genoIndices[fi], models[fi]);
    return models;
}

static std::vector<AFModel> loadAFModelsText(
    const std::string &path,
    int nPC,
    uint32_t nExpected
) {
    std::vector<AFModel> models;
    models.reserve(nExpected);

    uint32_t lineNo = 0;
    auto parseLine = [&](std::string &line) {
        ++lineNo;
        if (text::skipLine(line)) return;
        text::TokenScanner ts(line);
        auto sv = ts.nextView();
        if (sv.empty()) return;
        char *end;
        long status = std::strtol(sv.data(), &end, 10);
        if (end != sv.data() + sv.size() || status < -2 || status > 3) throw std::runtime_error(path + " line " + std::to_string(lineNo) +
                                                       ": expected STATUS BETA0 ...");
        AFModel m;
        m.status = static_cast<int8_t>(status);
        m.betas.resize(1 + nPC);
        for (int j = 0; j <= nPC; ++j) {
            ts.skipWS();
            if (ts.atEnd()) throw std::runtime_error(
                          path + " line " + std::to_string(lineNo) + ": expected " +
                          std::to_string(1 + nPC) + " beta values"
            );
            const char *cur = ts.pos();
            m.betas[j] = std::strtod(cur, &end);
            if (end == cur) throw std::runtime_error(
                          path + " line " + std::to_string(lineNo) + ": expected " +
                          std::to_string(1 + nPC) + " beta values"
            );
            ts.p = end;
        }
        ts.skipWS();
        if (!ts.atEnd() || !m.betas.allFinite())
            throw std::runtime_error(path + " line " + std::to_string(lineNo) + ": invalid or extra AF coefficients");
        models.push_back(std::move(m));
    };

    TextReader reader(path);
    std::string line;
    while (reader.getline(line))
        parseLine(line);

    return models;
}

std::vector<AFModel>loadAFModels(
    const std::string &path,
    int nPC,
    uint32_t nMarkers,
    const std::vector<uint64_t> &
    genoIndices
) {
    const auto len = path.size();
    std::vector<AFModel> models;
    if (len > 4 && path.compare(len - 4, 4, ".bin") == 0)models = loadAFModelsBinary(path, nPC, nMarkers, genoIndices);
    else models = loadAFModelsText(path, nPC, nMarkers);

    if (models.size() != nMarkers)throw std::runtime_error(
                  "AF model count (" + std::to_string(models.size()) +
                  ") does not match marker count (" + std::to_string(nMarkers) + ")"
    );
    return models;
}

// ======================================================================
// exportFrequencyModels — pre-compute per-marker AF models (top-level orchestration)
// ======================================================================

#include <atomic>
#include <chrono>
#include <ctime>
#include <exception>
#include <mutex>
#include <sstream>
#include <thread>

#include "genotype/genotype_source.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"

namespace {

void imputeMissing(
    Eigen::Ref<Eigen::VectorXd> g,
    int N,
    double altFreq
) {
    const double mean = 2.0 * altFreq;
    for (int i = 0; i < N; ++i)
        if (std::isnan(g[i])) g[i] = mean;
}

} // anonymous namespace

std::vector<AFModel> computeAFModelsInMemory(
    const GenotypeCatalog &plinkData,
    const AFContext &afCtx,
    const std::vector<uint32_t> &genoToFlat,
    int nthread,
    const AFQCThresholds &qc
) {
    if (nthread < 1 || qc.missing < 0 || qc.missing > 1 || qc.maf < 0 || qc.maf > 0.5 ||
        qc.mac < 0 || qc.hwe < 0 || qc.hwe > 1 ||
        !std::isfinite(qc.missing) || !std::isfinite(qc.maf) ||
        !std::isfinite(qc.mac) || !std::isfinite(qc.hwe))
        throw std::runtime_error("Invalid AF export thread count or QC thresholds");
    const auto &markerInfo = plinkData.markerInfo();
    const auto &chunkIndices = plinkData.chunkIndices();
    const size_t nTotalChunks = chunkIndices.size();
    const uint32_t nMarkers = static_cast<uint32_t>(markerInfo.size());
    const int N = afCtx.N;

    std::vector<AFModel> models(nMarkers);
    for (auto &m : models) {
        m.status = -1;
        m.betas = Eigen::VectorXd::Zero(1 + afCtx.nPC);
    }

    std::atomic<size_t> nextChunk(0);
    std::exception_ptr workerError;
    std::mutex errorMutex;

    const double hweLog10Cutoff = qc.hwe > 0 ? -std::log10(qc.hwe) : 0.0;
    auto workerFn = [&]() { try {
        auto cursor = plinkData.makeCursor();
        Eigen::VectorXd GVec(N);
        std::vector<uint32_t> missing;

        while (true) {
            const size_t ci = nextChunk.fetch_add(1);
            if (ci >= nTotalChunks) break;

            const auto &chunk = chunkIndices[ci];
            if (!chunk.empty()) cursor->beginSequentialBlock(chunk.front());

            for (uint64_t gIdx : chunk) {
                if (gIdx >= genoToFlat.size())
                    throw std::runtime_error("AF genotype index exceeds marker mapping");
                const uint32_t flatIdx = genoToFlat[gIdx];
                if (flatIdx == UINT32_MAX) continue;
                if (flatIdx >= models.size())
                    throw std::runtime_error("AF marker mapping exceeds result count");
                const GenoStats stats = cursor->getGenotypes(gIdx, GVec, missing);
                const bool pass = missing.size() < static_cast<size_t>(N) &&
                    std::isfinite(stats.altFreq) && std::isfinite(stats.missingRate) &&
                    std::isfinite(stats.maf) && std::isfinite(stats.mac) &&
                    !(stats.missingRate > qc.missing || stats.maf < qc.maf || stats.mac < qc.mac ||
                      (qc.hwe > 0 && (!std::isfinite(stats.log10pHwe) || stats.log10pHwe > hweLog10Cutoff)));
                if (!pass) continue;
                imputeMissing(GVec, N, stats.altFreq);
                models[flatIdx] = computeAFModel(GVec, stats.altFreq, afCtx);
            }
            infoMsg("  AF chunk %zu/%zu done", ci + 1, nTotalChunks);
        }
    } catch (...) {
        std::lock_guard<std::mutex> lk(errorMutex);
        if (!workerError) workerError = std::current_exception();
    }};

    const int nt = std::min(nthread, static_cast<int>(nTotalChunks));
    std::vector<std::thread> workers;
    workers.reserve(nt);
    for (int i = 0; i < nt; ++i)
        workers.emplace_back(workerFn);
    for (auto &w : workers)
        w.join();

    if (workerError) std::rethrow_exception(workerError);

    int cntUnfitted = 0, cntFailed = 0, cnt0 = 0, cnt1 = 0, cnt2 = 0, cnt3 = 0;
    for (uint32_t i = 0; i < nMarkers; ++i) {
        if (models[i].status == 1)++cnt1;
        else if (models[i].status == 2)++cnt2;
        else if (models[i].status == 3)++cnt3;
        else if (models[i].status == -1)++cntUnfitted;
        else if (models[i].status == -2)++cntFailed;
        else ++cnt0;
    }
    infoMsg("AF models: failed=%d, unfitted/QC-filtered=%d, pooled=%d, OLS=%d, logistic=%d, complement-logistic=%d",
            cntFailed, cntUnfitted, cnt0, cnt1, cnt2, cnt3);

    return models;
}

void exportFrequencyModels(
    const std::vector<std::string> &pcColNames,
    const std::string &phenoFile,
    const std::string &covarFile,
    const GenotypeInput &geno,
    const std::string &outputFile,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const auto wallStart = std::chrono::steady_clock::now();
    const std::clock_t cpuStart = std::clock();

    // Load pheno / covar files and extract PC columns.
    infoMsg("Loading PC data for AF model (%zu columns)", pcColNames.size());
    auto famIIDs = readGenotypeSampleIds(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile);
    if (!covarFile.empty()) sd.loadCovar(covarFile, pcColNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();
    sd.dropNaInColumns(pcColNames);

    const int N = static_cast<int>(sd.nUsed());
    const int nPC = static_cast<int>(pcColNames.size());
    Eigen::MatrixXd PCs = sd.getColumns(pcColNames);
    infoMsg("  %u subjects, %d PCs", sd.nUsed(), nPC);

    Eigen::MatrixXd onePlusPCs(N, 1 + nPC);
    onePlusPCs.col(0).setOnes();
    onePlusPCs.rightCols(nPC) = PCs;
    Eigen::MatrixXd XtX_inv_Xt;
    Eigen::VectorXd sqrt_XtX_inv_diag;
    prepareAFLinearDesign(onePlusPCs, XtX_inv_Xt, sqrt_XtX_inv_diag);

    const AFContext afCtx{onePlusPCs, XtX_inv_Xt, sqrt_XtX_inv_diag, PCs, N, nPC};

    auto genoData = openGenotypeCatalog(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);
    infoMsg("  %u subjects, %u markers", genoData->nSubjUsed(), genoData->nMarkers());

    const auto &markerInfo = genoData->markerInfo();
    const auto &chunkIndices = genoData->chunkIndices();
    const size_t nTotalChunks = chunkIndices.size();
    const uint32_t nMarkers = static_cast<uint32_t>(markerInfo.size());
    const uint32_t nBimMarkers = genoData->nMarkers();

    // Build genoIndex → flat index map.
    std::vector<uint32_t> genoToFlat(genoData->nMarkers(), UINT32_MAX);
    for (uint32_t fi = 0; fi < nMarkers; ++fi)
        genoToFlat[markerInfo[fi].genoIndex] = fi;

    infoMsg("Computing AF models (%d thread(s), %zu chunks)...", nthread, nTotalChunks);
    const AFQCThresholds qc{missingCutoff, minMafCutoff, minMacCutoff, hweCutoff};
    const std::vector<AFModel> allModels = computeAFModelsInMemory(*genoData, afCtx, genoToFlat, nthread, qc);

    infoMsg("Writing output: %s", outputFile.c_str());
    IndivAFWriter writer(outputFile, nBimMarkers, nPC);
    for (uint32_t fi = 0; fi < nMarkers; ++fi) {
        const auto &mi = markerInfo[fi];
        writer.write(mi.genoIndex, allModels[fi].status, allModels[fi].betas);
    }
    writer.close();

    const auto wallEnd = std::chrono::steady_clock::now();
    const double wallSec = std::chrono::duration<double>(wallEnd - wallStart).count();
    const double cpuSec = static_cast<double>(std::clock() - cpuStart) / CLOCKS_PER_SEC;
    infoMsg("Wall time: %.1f seconds, CPU time: %.1f seconds", wallSec, cpuSec);
}
