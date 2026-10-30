// Exported AF models and the online minor-oriented paper cascade must agree.
#include "allele_model/frequency_model.hpp"
#include "util/probability_math.hpp"
#include "util/spa.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
int checks = 0;
void require(bool ok, const std::string &name) {
    ++checks;
    if (!ok) throw std::runtime_error(name);
}
struct Design {
    Eigen::MatrixXd pcs, x, inverseTimesX;
    Eigen::VectorXd diag;
    explicit Design(const Eigen::VectorXd &pc) : pcs(pc), x(pc.size(), 2) {
        x.col(0).setOnes(); x.col(1) = pc;
        prepareAFLinearDesign(x, inverseTimesX, diag);
    }
    AFContext context() const { return {x, inverseTimesX, diag, pcs, int(x.rows()), 1}; }
};
void close(const Eigen::VectorXd &a, const Eigen::VectorXd &b, const std::string &name) {
    require(a.size() == b.size() && a.allFinite() && b.allFinite() &&
            (a - b).cwiseAbs().maxCoeff() < 3e-12, name);
}
template<class F> void throws(F fn, const std::string &name) {
    bool error = false;
    try { fn(); } catch (const std::runtime_error &) { error = true; }
    require(error, name);
}
AFModel compare(const Eigen::VectorXd &g, const Design &design, int status,
                const std::string &name) {
    const double q = g.mean() / 2;
    const AFContext ctx = design.context();
    const AFModel model = computeAFModel(g, q, ctx);
    require(model.status == status, name + " model status");
    Eigen::VectorXd reconstructed(g.size()), online(g.size());
    getAFVecFromModel(model, q, design.x, int(g.size()), reconstructed);
    if (q > .5) {
        const Eigen::VectorXd minor = 2 - g.array();
        require(computeAFVec(minor, 1 - q, ctx, online, true), name + " online complement fit succeeds");
    } else computeAFVec(g, q, ctx, online);
    close(reconstructed, online, name + " reconstructed equals online AF");
    require(model.betas.size() == 2 && model.betas.allFinite(), name + " coefficients");
    return model;
}

void stableLogisticMappings(std::vector<AFModel> &saved) {
    const std::vector<double> eta = {-700, -100, -40, -10, 0, 10, 40, 100, 700};
    Eigen::MatrixXd x(eta.size(), 2);
    x.col(0).setOnes();
    for (size_t i = 0; i < eta.size(); ++i) x(i, 1) = eta[i];
    AFModel model{2, Eigen::VectorXd(2)};
    model.betas << 0, 1;
    Eigen::VectorXd minor(eta.size()), complement(eta.size());
    getAFVecFromModel(model, .25, x, int(eta.size()), minor);
    saved.push_back(model);
    model.status = 3;
    getAFVecFromModel(model, .75, x, int(eta.size()), complement);
    saved.push_back(model);
    for (size_t i = 0; i < eta.size(); ++i) {
        const long double t = eta[i];
        const long double logDenom = t >= 0 ? t + std::log1p(std::exp(-t))
                                            : std::log1p(std::exp(t));
        const long double minorReference = -std::expm1(-.5L * logDenom);
        const long double complementReference = std::exp(-.5L * logDenom);
        require(minor[i] > 0 && complement[i] > 0,
                "representable logistic AF tails do not collapse to zero");
        require(std::abs(static_cast<long double>(minor[i]) / minorReference - 1) < 3e-14L,
                "minor AF agrees with extended-precision tail reference");
        require(std::abs(static_cast<long double>(complement[i]) / complementReference - 1) < 3e-14L,
                "complement AF agrees with extended-precision tail reference");
        require(std::abs(minor[i] + complement[i] - 1) < 3e-16,
                "minor and complement AF sum to one");
        if (i) {
            require(minor[i] >= minor[i - 1], "minor AF is monotone in logistic predictor");
            require(complement[i] <= complement[i - 1], "complement AF is monotone in logistic predictor");
        }
    }
    require(std::abs(minor[4] - (1 - std::sqrt(.5))) < 2e-16 &&
            std::abs(complement[4] - std::sqrt(.5)) < 2e-16,
            "zero-predictor HWE AF has the closed-form value");
    require(minor[2] < 1e-17 && complement[6] > 2e-9,
            "small minor and small complement frequencies survive previous cancellation cases");
}

void linearDesignStability() {
    const int n = 200;
    Eigen::MatrixXd x(n, 3);
    x.col(0).setOnes();
    for (int i = 0; i < n; ++i) {
        const double value = -1.0 + 2.0 * i / (n - 1);
        x(i, 1) = value;
        x(i, 2) = value + 1e-7 * std::sin(1.37 * i);
    }
    Eigen::MatrixXd inverseTimesX;
    Eigen::VectorXd diag;
    prepareAFLinearDesign(x, inverseTimesX, diag);
    const Eigen::Vector3d coefficients(.6, .2, -.1);
    const Eigen::VectorXd observed = x * coefficients;
    const Eigen::VectorXd recovered = inverseTimesX * observed;
    require((recovered - coefficients).cwiseAbs().maxCoeff() < 2e-7,
            "nearly collinear PC design recovers known coefficients");
    require((x * recovered - observed).cwiseAbs().maxCoeff() < 2e-8,
            "nearly collinear PC fitted dosage remains accurate");
    require(diag.allFinite() && (diag.array() > 0).all(),
            "nearly collinear PC coefficient uncertainty remains finite");
    Eigen::MatrixXd bad = x;
    bad.col(2) = bad.col(1);
    throws([&]() { prepareAFLinearDesign(bad, inverseTimesX, diag); },
           "duplicate PCs are rejected by shared online/export design preparation");
    bad = x;
    bad(1, 2) = std::numeric_limits<double>::quiet_NaN();
    throws([&]() { prepareAFLinearDesign(bad, inverseTimesX, diag); },
           "non-finite PC design is rejected");
    bad = Eigen::MatrixXd::Identity(3, 3);
    throws([&]() { prepareAFLinearDesign(bad, inverseTimesX, diag); },
           "PC design without residual degrees of freedom is rejected");
}

std::vector<AFModel> boundaries() {
    std::vector<AFModel> saved;
    const int n = 400;
    for (int test = 0; test < 3; ++test) {
        Eigen::VectorXd pc(n), g = Eigen::VectorXd::Zero(n);
        if (test == 0) {
            for (int i = 0; i < n; ++i) pc[i] = i < 60 ? 3 : (i < 100 ? 1 : 0);
            for (int i = 0; i < 58; ++i) g[i] = 2;
            for (int i = 60; i < 98; ++i) g[i] = 2;
            for (int i = 100; i < 110; ++i) g[i] = 1;
        } else if (test == 1) {
            for (int i = 0; i < n; ++i) pc[i] = i < 40 ? -1 : (i < 360 ? 0 : 1);
            for (int i = 0; i < 3; ++i) g[i] = 1;
            for (int i = 40; i < 50; ++i) g[i] = 1;
            for (int i = 360; i < 395; ++i) g[i] = 2;
        } else {
            pc = Eigen::VectorXd::LinSpaced(n, -1, 1);
            g[0] = 2;
            for (int i = 372; i <= 398; i += 2) g[i] = 2;
        }
        pc.array() -= pc.mean();
        pc /= std::sqrt(pc.squaredNorm() / (n - 1));
        Design design(pc);
        saved.push_back(compare(g, design, 2, "AF boundary " + std::to_string(test)));
        saved.push_back(compare((2 - g.array()).matrix(), design, 3,
                                "high ALT logistic " + std::to_string(test)));
        if (test == 2) {
            require((g.array() > .5).count() <= 20 && g.sum() > 20,
                    "carrier count differs from minor MAC");
            compare((.95 * g).eval(), design, 2, "fractional low-carrier dosage");
        }
    }
    Eigen::VectorXd pc = Eigen::VectorXd::LinSpaced(n, -1, 1);
    Design design(pc);
    const Eigen::VectorXd dosage = .6 + .3 * pc.array();
    const AFModel low = compare(dosage, design, 1, "OLS dosage");
    const AFModel high = compare((2 - dosage.array()).matrix(), design, 1, "high ALT OLS dosage");
    require(std::abs(high.betas[0] - (2 - low.betas[0])) < 1e-13 &&
            std::abs(high.betas[1] + low.betas[1]) < 1e-13, "ALT OLS coefficient complement");
    Eigen::VectorXd rare = Eigen::VectorXd::Zero(n);
    rare.head(10).setConstant(2);
    compare(rare, design, 0, "minor MAC exactly 20");
    compare((2 - rare.array()).matrix(), design, 0, "high ALT minor MAC exactly 20");
    for (double displacement : {-1e-8, 1e-8}) {
        const Eigen::VectorXd near = Eigen::VectorXd::Constant(n, (20 + displacement) / n);
        const int status = displacement < 0 ? 0 : 1;
        compare(near, design, status, "minor dosage MAC just below/above 20");
        compare((2 - near.array()).matrix(), design, status,
                "high ALT dosage MAC just below/above 20");
    }
    return saved;
}

void convergenceFailures(std::vector<AFModel> &saved) {
    const int n = 400;
    const Eigen::VectorXd pc = Eigen::VectorXd::LinSpaced(n, -2, 2);
    Design design(pc);
    Eigen::VectorXd separated(n), ordinary(n);
    for (int i = 0; i < n; ++i) {
        separated[i] = pc[i] > 0 ? 1 : 0;
        ordinary[i] = std::fmod(i * .618033988749895, 1.0) <
                          1.0 / (1.0 + std::exp(-3.0 * pc[i])) ? 1 : 0;
    }
    bool converged = false;
    const Eigen::VectorXd normal = math::logisticRegressionBeta(
        design.pcs, ordinary, 1e-6, 100, &converged);
    require(converged && normal.allFinite(), "ordinary logistic fit converges");
    close(normal, math::logisticRegressionBeta(design.pcs, ordinary),
          "diagnostic leaves converged logistic coefficients unchanged");
    const Eigen::VectorXd exhausted = math::logisticRegressionBeta(
        design.pcs, ordinary, 1e-6, 1, &converged);
    require(!converged && exhausted.allFinite(), "iteration exhaustion is explicit despite finite coefficients");
    const Eigen::VectorXd unbounded = math::logisticRegressionBeta(
        design.pcs, separated, 1e-6, 100, &converged);
    require(!converged && unbounded.allFinite(), "complete separation is not reported as converged");
    compare(ordinary, design, 2, "converged logistic after diagnostics");
    compare((2 - ordinary.array()).matrix(), design, 3, "converged complement-logistic after diagnostics");
    for (bool flip : {false, true}) {
        const Eigen::VectorXd g = flip ? (2 - separated.array()).matrix().eval() : separated;
        const double q = g.mean() / 2;
        const AFModel failed = computeAFModel(g, q, design.context());
        require(failed.status == -2, "failed fit remains distinct from pooled AF in either allele direction");
        require(failed.betas.isZero(), "failed AF record has no usable coefficients");
        saved.push_back(failed);
        Eigen::VectorXd out(n);
        const Eigen::VectorXd minor = flip ? (2 - g.array()).matrix().eval() : g;
        require(!computeAFVec(minor, std::min(q, 1 - q), design.context(), out),
                "online AF reports failed fitting");
        require(out.array().isNaN().all(), "failed online AF values are all NaN");
        throws([&]() { getAFVecFromModel(failed, q, design.x, n, out); },
               "failed model cannot be reconstructed as a usable AF");
    }
    require(spa::statusIsNA(spa::Status::AfFitFailed) &&
            !spa::statusIsUsable(spa::Status::AfFitFailed) &&
            !spa::statusIsFallback(spa::Status::AfFitFailed),
            "AF failure is NA rather than a usable or normal-fallback result");
    require(std::string(spa::statusName(spa::Status::AfFitFailed)) == "AF_FIT_FAILED",
            "AF failure has an explicit diagnostic name");
    require(spa::worseStatus(spa::Status::Normal, spa::Status::AfFitFailed) == spa::Status::AfFitFailed,
            "AF failure is preserved when statuses combine");
}

struct MockCatalog : GenotypeCatalog {
    std::vector<Eigen::VectorXd> columns;
    std::vector<GenoStats> stats;
    std::vector<MarkerInfo> markers;
    std::vector<std::vector<uint64_t>> chunks;
    explicit MockCatalog(int n) {
        Eigen::VectorXd common(n);
        for (int i = 0; i < n; ++i) common[i] = i % 3;
        columns.assign(8, common);
        for (int i = 0; i < n / 10; ++i) columns[1][i] = std::numeric_limits<double>::quiet_NaN();
        columns[2].setZero(); columns[2].head(2).setOnes();
        columns[3].setZero(); columns[3].head(8).setOnes();
        columns[4].setOnes();
        columns[5].setConstant(std::numeric_limits<double>::quiet_NaN());
        columns[7][0] = std::numeric_limits<double>::quiet_NaN();
        chunks = {{0, 1, 2, 3}, {4, 5, 6, 7}};
        for (int j = 0; j < 8; ++j) {
            int observed = 0; double sum = 0;
            for (int i = 0; i < n; ++i) if (std::isfinite(columns[j][i])) { ++observed; sum += columns[j][i]; }
            const double q = observed ? sum / (2 * observed) : 0;
            const double maf = std::min(q, 1 - q);
            stats.push_back({q, sum, 1 - double(observed) / n, j == 4 ? 12.0 : 0.0,
                             maf, 2 * observed * maf, false});
            markers.push_back({"1", uint32_t(j + 1), "m" + std::to_string(j), "A", "G", uint64_t(j)});
        }
    }
    uint32_t nMarkers() const override { return uint32_t(markers.size()); }
    uint32_t nSubjUsed() const override { return uint32_t(columns[0].size()); }
    uint32_t nSubjInFile() const override { return nSubjUsed(); }
    std::string_view chr(uint64_t i) const override { return markers[i].chrom; }
    std::string_view markerId(uint64_t i) const override { return markers[i].id; }
    uint32_t pos(uint64_t i) const override { return markers[i].pos; }
    std::string_view ref(uint64_t i) const override { return markers[i].ref; }
    std::string_view alt(uint64_t i) const override { return markers[i].alt; }
    const std::vector<MarkerInfo> &markerInfo() const override { return markers; }
    const std::vector<std::vector<uint64_t>> &chunkIndices() const override { return chunks; }
    struct Cursor : VariantCursor {
        const MockCatalog &owner;
        explicit Cursor(const MockCatalog &catalog) : owner(catalog) {}
        void beginSequentialBlock(uint64_t) override {}
        GenoStats getGenotypes(uint64_t i, Eigen::Ref<Eigen::VectorXd> out,
                                std::vector<uint32_t> &missing) override {
            out = owner.columns[i]; missing.clear();
            for (int j = 0; j < out.size(); ++j) if (std::isnan(out[j])) missing.push_back(j);
            return owner.stats[i];
        }
        void getGenotypesSimple(uint64_t, Eigen::Ref<Eigen::VectorXd>) override {
            throw std::runtime_error("QC export must use the stats-producing cursor");
        }
    };
    std::unique_ptr<VariantCursor> makeCursor() const override { return std::make_unique<Cursor>(*this); }
};

void qcAndStorage(const std::vector<AFModel> &models, const std::string &directory) {
    const int n = 120;
    Design design(Eigen::VectorXd::LinSpaced(n, -1, 1));
    MockCatalog catalog(n);
    std::vector<uint32_t> mapping = {0, 1, 2, 3, 4, 5, UINT32_MAX, 7};
    const AFQCThresholds qc{.05, .02, 10, 1e-6};
    const auto fitted = computeAFModelsInMemory(catalog, design.context(), mapping, 2, qc);
    require(fitted[0].status >= 0 && fitted[7].status >= 0, "QC accepts eligible markers");
    for (int i = 1; i <= 6; ++i) require(fitted[i].status == -1, "QC/skip leaves explicit unfitted status");
    Eigen::VectorXd imputed = catalog.columns[7];
    imputed[0] = 2 * catalog.stats[7].altFreq;
    const AFModel direct = computeAFModel(imputed, catalog.stats[7].altFreq, design.context());
    require(fitted[7].status == direct.status, "QC mean-imputation model status");
    close(fitted[7].betas, direct.betas, "QC fits the mean-imputed genotype");

    std::filesystem::create_directories(directory);
    for (const std::string suffix : {".bin", ".txt", ".txt.gz"}) {
        const std::string path = directory + "/models" + suffix;
        { IndivAFWriter writer(path, models.size() + 1, 1);
          for (size_t i = 0; i < models.size(); ++i) writer.write(i, models[i].status, models[i].betas);
          if (suffix != ".bin") writer.write(models.size(), -1, Eigen::VectorXd::Zero(2)); }
        std::vector<uint64_t> indices(models.size() + 1);
        for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;
        const auto loaded = loadAFModels(path, 1, uint32_t(indices.size()), indices);
        for (size_t i = 0; i < models.size(); ++i) {
            require(loaded[i].status == models[i].status, "status survives storage");
            close(loaded[i].betas, models[i].betas, "betas survive storage");
            if (loaded[i].status == -2) {
                Eigen::VectorXd rejected(n);
                throws([&]() { getAFVecFromModel(loaded[i], .1, design.x, n, rejected); },
                       "stored failed model is rejected after loading");
            } else {
                Eigen::MatrixXd extreme(3, 2);
                extreme << 1, -700, 1, 40, 1, 700;
                Eigen::VectorXd before(3), after(3);
                getAFVecFromModel(models[i], .1, extreme, 3, before);
                getAFVecFromModel(loaded[i], .1, extreme, 3, after);
                require((before.array() == after.array()).all(),
                        "stored logistic model preserves extreme AF reconstruction exactly");
            }
        }
        require(loaded.back().status == -1, "unwritten binary/QC text slot is unfitted");
    }
    Eigen::VectorXd out(n);
    AFModel invalid{42, Eigen::VectorXd::Zero(2)};
    throws([&]() { getAFVecFromModel(invalid, .1, design.x, n, out); }, "unknown status is rejected");
    invalid.status = -1;
    throws([&]() { getAFVecFromModel(invalid, .1, design.x, n, out); }, "unfitted status is rejected");
    invalid.status = 2; invalid.betas[0] = std::numeric_limits<double>::quiet_NaN();
    throws([&]() { getAFVecFromModel(invalid, .1, design.x, n, out); }, "invalid coefficients are rejected");
    const std::string bad = directory + "/bad.txt";
    { std::ofstream f(bad); f << "256\t0\t0\n"; }
    throws([&]() { loadAFModels(bad, 1, 1, {0}); }, "status narrowing overflow is rejected");
}
} // namespace

int main(int argc, char **argv) {
    try {
        auto models = boundaries();
        stableLogisticMappings(models);
        linearDesignStability();
        convergenceFailures(models);
        qcAndStorage(models, argc > 1 ? argv[1] : "tmp/af_model_test");
        std::cout << checks << " AF model/online/QC/storage checks passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "FAIL " << e.what() << "\n"; return 1;
    }
}
