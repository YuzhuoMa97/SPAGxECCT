// geno_data.hpp — Abstract interfaces for genotype data backends
//
// GenotypeCatalog:   shared metadata + cursor factory (constructed once, main thread)
// VariantCursor: per-thread lightweight genotype decoder
//
// Concrete backends (BedCatalog, PgenCatalog, VcfCatalog, BgenCatalog)
// derive from these interfaces.  The marker engine and helper functions
// operate through the abstract interface so they work with any backend.
#pragma once

#include "genotype/equilibrium.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

class VariantCursor;

// ======================================================================
// GenotypeInput — lightweight descriptor for the genotype data source
// ======================================================================

// Vcf and Bcf share the same htslib-based backend but are kept as distinct
// values so that error messages, help output, and "Options in effect"
// logging echo the flag the user actually typed.  Content-based validation
// in VcfCatalog enforces that --vcf rejects BCF2 input and --bcf rejects VCF
// text, matching plink2's --vcf / --bcf behaviour.
enum class GenotypeFormat { Plink, Pgen, Vcf, Bcf, Bgen };

struct GenotypeInput {
    GenotypeFormat format;
    std::string path;        // bfilePrefix for Plink/Pgen, full file for Vcf/Bcf/Bgen
    std::string extractFile; // --extract: SNP include list
    std::string excludeFile; // --exclude: SNP exclude list
    std::unordered_set<std::string> chrFilter; // --chr: chromosome include set
    // BGEN-only: alleles[0]=ALT (set by --bgen FILE ref-last or ref-unknown;
    // false ⇒ --bgen FILE ref-first).
    bool bgenAltFirst = false;

    // plink2 --hard-call-threshold: dosage→hard-call cutoff for HWE counts
    // (default 0.1).  Propagated into the GenotypeCatalog and read by the cursors.
    double hardCallThreshold = 0.1;

    // Return a descriptive flag label, e.g. "--bfile d_bed".
    std::string flagLabel() const {
        const char *flag = nullptr;
        switch (format) {
        case GenotypeFormat::Plink: flag = "--bfile"; break;
        case GenotypeFormat::Pgen:  flag = "--pfile"; break;
        case GenotypeFormat::Vcf:   flag = "--vcf";   break;
        case GenotypeFormat::Bcf:   flag = "--bcf";   break;
        case GenotypeFormat::Bgen:  flag = "--bgen";  break;
        }
        return std::string(flag) + " " + path;
    }

};

// Return a descriptive label for the GRM source, e.g. "--grm-file e_grm.grab".
inline std::string grmFlagLabel(
    const std::string &tripletFile,
    const std::string &gctaFile
) {
    if (!tripletFile.empty()) return "--grm-file " + tripletFile;
    if (!gctaFile.empty()) return "--sp-grm-plink2 " + gctaFile;
    return "";
}

// Extract all sample IIDs from the genotype file metadata
// (plink .fam, pgen .psam, VCF header, BGEN sample block).
std::vector<std::string> readGenotypeSampleIds(const GenotypeInput &spec);

// Construct the appropriate GenotypeCatalog backend from a spec + subject mask.
std::unique_ptr<class GenotypeCatalog> openGenotypeCatalog(
    const GenotypeInput &spec,
    const std::vector<uint64_t> &usedMask,
    uint32_t nSamplesInFile,
    uint32_t nUsed,
    int nMarkersEachChunk = 1024
);

// ======================================================================
// GenotypeCatalog — shared, read-only genotype metadata + cursor factory
// ======================================================================

class GenotypeCatalog {
  public:
    // plink2 --hard-call-threshold, propagated from GenotypeInput by openGenotypeCatalog.
    // Read by the per-format cursors when classifying dosages into hard-calls
    // for the HWE genotype counts (default 0.1).
    double hardCallThreshold = 0.1;

    struct MarkerInfo {
        std::string chrom;
        uint32_t pos;
        std::string id;
        std::string ref;
        std::string alt;
        uint64_t genoIndex; // backend-specific row index
    };

    virtual ~GenotypeCatalog() = default;

    virtual uint32_t nMarkers() const = 0;

    virtual uint32_t nSubjUsed() const = 0;

    virtual uint32_t nSubjInFile() const = 0;

    virtual std::string_view chr(uint64_t i) const = 0;

    virtual std::string_view markerId(uint64_t i) const = 0;

    virtual uint32_t pos(uint64_t i) const = 0;

    virtual std::string_view ref(uint64_t i) const = 0;

    virtual std::string_view alt(uint64_t i) const = 0;

    virtual const std::vector<MarkerInfo> &markerInfo() const = 0;

    virtual const std::vector<std::vector<uint64_t> > &chunkIndices() const = 0;

    // Factory: create a per-thread cursor for genotype decoding.
    virtual std::unique_ptr<VariantCursor> makeCursor() const = 0;

};

// ======================================================================
// VariantCursor — per-thread, lightweight genotype decoder
// ======================================================================

class VariantCursor {
  public:
    virtual ~VariantCursor() = default;

    // Prepare sequential reading starting from a given marker index.
    virtual void beginSequentialBlock(uint64_t firstMarker) = 0;

    // Decode genotype for marker gIndex into caller-owned Eigen vector, and
    // return its QC statistics.  `indexForMissing` receives the dense indices
    // that were written NaN; it is cleared first.
    virtual GenoStats getGenotypes(
        uint64_t gIndex,
        Eigen::Ref<Eigen::VectorXd> out,
        std::vector<uint32_t> &indexForMissing
    ) = 0;

    // Lightweight variant: genotype column only, missing → NaN, no QC stats.
    // `out` may alias a destination matrix column, so a backend that can decode
    // a variant cheaply (pgen's difflist: set the modal genotype, then patch the
    // few that differ) writes the result in place rather than handing the caller
    // a sparse representation to expand.
    virtual void getGenotypesSimple(
        uint64_t gIndex,
        Eigen::Ref<Eigen::VectorXd> out
    ) = 0;

};
