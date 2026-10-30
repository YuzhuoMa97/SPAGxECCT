// indiv_af.hpp — Per-individual allele-frequency model:
//   compute, write (binary / text / gz), and reconstruct AF vectors.
//   Also: exportFrequencyModels — pre-compute AF models for all markers.
//
// Binary format (.bin):
//   No header. Records indexed by raw .bim line (genoIndex).
//   filePos = genoIndex * recordSize
//   recordSize = sizeof(int8_t) + (1+nPC)*sizeof(double)
//   status int8_t: -2 = failed AF fit, -1 = not fitted/QC filtered,
//                  0 = pooled ALT frequency,
//                  1 = ALT OLS, 2 = logistic, 3 = complement-logistic
//   betas  (1+nPC) doubles: intercept first, then one per PC
//   Unprocessed markers: status=-1 and zero betas (not a fitted pooled model).
//
// Text / gz format:
//   Header: #STATUS\tBETA0\tBETA1\t...\n
//   One row per filtered marker, tab-separated, doubles as %.17g.
//
// Shared between SPAmixAF (pre-compute step) and SPAmix (on-the-fly).
#pragma once

#include "genotype/genotype_source.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

// ======================================================================
// AFModel — result of fitting one marker
// ======================================================================

struct AFModel {
    int8_t status;         // -2=failed fit, -1=unfitted, 0=pooled ALT, 1=ALT OLS, 2=logistic,
                           // 3=1-logisticAF(minor beta); same record size
    Eigen::VectorXd betas; // length (1+nPC); all zeros for status 0
};

// ======================================================================
// AFContext — pre-computed OLS matrices for AF estimation
// ======================================================================

struct AFContext {
    const Eigen::MatrixXd &onePlusPCs;        // (N × 1+nPC) = [1 | PCs]
    const Eigen::MatrixXd &XtX_inv_Xt;        // (1+nPC × N) = (X'X)^{-1} X'
    const Eigen::VectorXd &sqrt_XtX_inv_diag; // (1+nPC) = sqrt(diag((X'X)^{-1}))
    Eigen::Ref<const Eigen::MatrixXd> PCs;    // (N × nPC) — accepts MatrixXd or block
    int N;
    int nPC;
};

// Shared rank-checked QR preparation for online and exported AF models.
// Avoid forming X'X, which would square the PC design's condition number.
void prepareAFLinearDesign(
    const Eigen::MatrixXd &x,
    Eigen::MatrixXd &inverseTimesX,
    Eigen::VectorXd &sqrtDiag
);

// ======================================================================
// computeAFModel — fit model, return status + betas (for pre-compute storage)
//
// Cascade:
// Fit internally to min(ALT, REF), including the MAC <=20 gate.
//   minor MAC <= 20         → status 0 (uniform altFreq)
//   OLS fit good            → status 1 (full ALT-scale OLS coefficients)
//   OLS bad, no sig PCs     → status 0
//   logistic on sig PCs     → status 2 if ALT is minor, status 3 otherwise
// Status 3 stores minor-carrier betas: ALT q = 1-[1-sqrt(1-sigmoid(X beta))].
// Non-convergent/invalid fitted models have status -2, not pooled frequency.
// The paper's pooled branches remain unchanged. No second <=20-carrier gate
// is imposed; all-zero/all-one logistic outcomes use the existing pooled branch.
// ======================================================================

AFModel computeAFModel(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    double altFreq,
    const AFContext &ctx
);

// ======================================================================
// computeAFVec — fused compute: directly writes per-individual AF into out
//
// Same shared paper cascade as computeAFModel. This low-level entry point
// requires g and minorFreq already minor-oriented. Set complementOutput for
// ALT frequency >0.5; logistic AF is then reconstructed directly in the ALT
// direction without subtracting a rounded minor AF. Reuses fitted OLS values
// rather than multiplying the design again.
// ======================================================================

// Returns false with an all-NaN output when the selected AF fit fails.
bool computeAFVec(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    double minorFreq,
    const AFContext &ctx,
    Eigen::Ref<Eigen::VectorXd> out,
    bool complementOutput = false
);

// ======================================================================
// getAFVecFromModel — reconstruct per-individual AF from a stored model
//
//  model    — AFModel as returned by computeAFModel / loaded from file
//  altFreq  — overall alt allele frequency (used for status 0)
//  onePlusPCs — (N × 1+nPC)
//  N        — sample size
//  out      — output vector (size N), caller-allocated
// Rejects status -2 (failed fit), -1 (unfitted/QC filtered), unknown statuses and invalid model
// dimensions/coefficients. Status 3 explicitly reconstructs complement-logistic.
// ======================================================================

void getAFVecFromModel(
    const AFModel &model,
    double altFreq,
    const Eigen::MatrixXd &onePlusPCs,
    int N,
    Eigen::Ref<Eigen::VectorXd> out
);

// ======================================================================
// IndivAFWriter — write all marker AF models to disk
//
// Usage:
//   IndivAFWriter w(outputFile, nBimMarkers, nPC, markerInfo);
//   w.write(genoIndex, chr, bp, status, betas);  // called per marker
//   w.close();
// ======================================================================

// Forward-declared from plink.hpp — avoid circular include.
// Caller passes CHR / BP as strings / uint32_t directly.

class IndivAFWriter {
  public:
    enum class Mode { Binary, Text };

// Infer mode from outputFile extension: .bin → Binary, else Text.
// Text mode supports .gz, .zst, or plain via TextWriter.
    IndivAFWriter(
        const std::string &outputFile,
        uint64_t nBimMarkers,           // total .bim lines (binary pre-alloc)
        int nPC
    );

    ~IndivAFWriter();

// Write one record.  For Binary, genoIndex determines the file offset.
    void write(
        uint64_t genoIndex,
        int8_t status,
        const Eigen::VectorXd &betas
    );

    void close();

    Mode mode() const {
        return m_mode;
    }

  private:
    Mode m_mode;
    int m_nPC;
    long long m_recordSize; // valid only for Binary

// Binary output
    std::fstream m_binOut;

// Text output (TextWriter handles gz/zst/plain)
    std::unique_ptr<class TextWriter> m_writer;

    bool m_closed = false;
};

// ======================================================================
// IndivAFReader — load one record from a binary AF file by genoIndex
//
//   model.betas will have length (1+nPC).
// ======================================================================

class IndivAFReader {
  public:
    IndivAFReader(
        const std::string &binFile,
        int nPC
    );
    ~IndivAFReader();

// Seek to genoIndex and fill model. Returns true for a fitted non-pooled model.
// Always inspect status: 0 is a usable pooled estimate, -1 is unfitted and -2
// is a failed fit that must not be used for association.
    bool read(
        uint64_t genoIndex,
        AFModel &model
    );

  private:
    std::ifstream m_in;
    int m_nPC;
    long long m_recordSize;
};

// ======================================================================
// loadAFModels — load pre-computed AF models from file
//
// Dispatches by extension:
//   .bin  → binary random-access (seeks by genoIndices[fi])
//   .gz / other → tab-separated text, read sequentially
//
// genoIndices[fi] = markerInfo[fi].genoIndex (used only for .bin).
// Always returns exactly nMarkers models in flat marker order.
// ======================================================================

std::vector<AFModel>loadAFModels(
    const std::string &path,
    int nPC,
    uint32_t nMarkers,
    const std::vector<uint64_t> &
    genoIndices
);

// ======================================================================
// exportFrequencyModels — pre-compute per-marker AF models and write to disk.
//
// The eigenvector file provides both subject IDs (first column) and PCs.
//
// Utility: SPAGxE --cal-af-coef --bfile PREFIX --covar PCS
//          --pc-cols PC1,PC2 --out OUTPUT
// ======================================================================

// Compute AF models for all markers in parallel, storing results in memory.
// genoToFlat maps raw genoIndex → flat marker index (UINT32_MAX = skip).
// Callers must also include "genotype/bed_source.hpp".
class GenotypeCatalog;
struct AFQCThresholds {
    double missing = 1.0;
    double maf = 0.0;
    double mac = 0.0;
    double hwe = 0.0;
};
std::vector<AFModel> computeAFModelsInMemory(
    const GenotypeCatalog &plinkData,
    const AFContext &afCtx,
    const std::vector<uint32_t> &genoToFlat,
    int nthread,
    const AFQCThresholds &qc = AFQCThresholds{}
);

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
    const std::string &keepFile = {},
    const std::string &removeFile = {}

);
