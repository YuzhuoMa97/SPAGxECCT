// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
// Independently implemented BGEN Layout-2 known-H collection decoder.
#include "util/path_io.hpp"
#include "local/local_bgen.hpp"
#include "util/text_stream.hpp"
#include <zlib.h>
#include <zstd.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace spagxe_local {
namespace {
constexpr double missing = std::numeric_limits<double>::quiet_NaN();
using Names = std::vector<std::string>;
using Metadata = std::array<std::string, 5>;

void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error("Local BGEN: " + message);
}

Names splitTabs(std::string text) {
    if (!text.empty() && text.back() == '\r') text.pop_back();
    Names result;
    size_t begin = 0;
    for (;;) {
        const size_t end = text.find('\t', begin);
        result.push_back(text.substr(begin, end == std::string::npos ? end : end - begin));
        if (end == std::string::npos) return result;
        begin = end + 1;
    }
}

void requireUnique(const Names &names, const std::string &kind) {
    std::unordered_set<std::string> seen;
    for (const auto &name : names) {
        require(!name.empty() && name.find_first_of("\t\r\n") == std::string::npos &&
                    seen.insert(name).second,
                "empty, invalid or duplicate " + kind + ": " + name);
    }
}

std::string relativeTo(const std::string &manifest, const std::string &file) {
    require(!file.empty(), "empty file path in manifest");
    return pathio::resolveManifestPath(manifest, file);
}

struct PosteriorRow {
    Metadata metadata;
    Eigen::VectorXd dosage, hapcount;
    bool uncertain = false;
};

// The file reader is independent of the diploid global-genotype reader. In
// particular, haploid and zero-copy ancestry records have different numbers
// of stored probability values and must not be decoded as diploid genotypes.
class AncestryBgen {
public:
    Names ids;
    uint32_t markerCount = 0;

    AncestryBgen(const std::string &path, uint64_t blockBudget)
        : input(pathio::nativePath(path), std::ios::binary), file(path), budget(blockBudget) {
        require(bool(input), "cannot open " + file);
        input.seekg(0, std::ios::end);
        const auto end = input.tellg();
        require(end >= std::streampos(0), "cannot determine file size: " + file);
        fileSize = uint64_t(end);
        input.seekg(0);
        const uint64_t offset = integer(4), headerLength = integer(4);
        require(headerLength >= 20 && headerLength <= (1u << 20) && offset >= headerLength &&
                    uint64_t(4) + offset <= fileSize,
                "invalid header/variant offset: " + file);
        markerCount = uint32_t(integer(4));
        samples = uint32_t(integer(4));
        require(markerCount > 0 && samples > 0 && uint64_t(samples) + 10 <= budget,
                "empty dimensions or sample count exceeds decoded-block budget: " + file);
        char magic[4];
        bytes(magic, 4);
        require(std::string(magic, 4) == "bgen" || std::string(magic, 4) == std::string(4, '\0'),
                "invalid file magic: " + file);
        input.seekg(std::streamoff(4 + headerLength - 4));
        const uint32_t flags = uint32_t(integer(4));
        compression = flags & 3u;
        require(compression <= 2 && ((flags >> 2) & 15u) == 2 &&
                    (flags & 0x7fffffc0u) == 0,
                "requires Layout 2 with raw/zlib/zstd and zero reserved flags: " + file);
        require((flags & 0x80000000u) != 0,
                "known-H collection requires embedded sample IDs: " + file);
        const uint64_t sampleStart = position();
        const uint64_t sampleBytes = integer(4);
        require(sampleBytes >= 8 + 3 * uint64_t(samples) && sampleBytes <= (uint64_t(1) << 30) &&
                    sampleBytes <= uint64_t(4) + offset - sampleStart && integer(4) == samples,
                "invalid embedded sample block: " + file);
        ids.reserve(samples);
        for (uint32_t i = 0; i < samples; ++i) ids.push_back(string(2, 65535));
        require(position() == sampleStart + sampleBytes, "sample block length mismatch: " + file);
        requireUnique(ids, "sample IID");
        input.seekg(std::streamoff(4 + offset));
    }

    bool next(PosteriorRow &row) {
        if (current == markerCount) {
            require(position() == fileSize, "trailing bytes after declared variants: " + file);
            return false;
        }
        const std::string variantId = string(2, 65535), rsid = string(2, 65535);
        const std::string chromosome = string(2, 65535);
        const uint32_t pos = uint32_t(integer(4));
        require(integer(2) == 2, "multiallelic variants are outside the known-H profile: " + file);
        const std::string ref = string(4, budget), alt = string(4, budget);
        row.metadata = {chromosome, std::to_string(pos), rsid.empty() ? variantId : rsid, ref, alt};
        require(pos > 0 && !chromosome.empty() && !row.metadata[2].empty() &&
                    !ref.empty() && !alt.empty() && ref != alt &&
                    ref.find(',') == std::string::npos && alt.find(',') == std::string::npos,
                "invalid biallelic variant metadata: " + file);
        require(variants.insert(row.metadata[2]).second, "duplicate variant ID: " + row.metadata[2]);

        uint64_t compressedBytes = integer(4), decodedBytes = compressedBytes;
        if (compression != 0) {
            require(compressedBytes >= 4, "compressed block omits decoded length: " + file);
            decodedBytes = integer(4);
            compressedBytes -= 4;
        }
        require(decodedBytes >= uint64_t(samples) + 10 && decodedBytes <= budget,
                "decoded block exceeds --local-max-block-bytes or is too short: " + file);
        const uint64_t packedBound = std::max<uint64_t>(
            uint64_t(compressBound(uLong(decodedBytes))), ZSTD_compressBound(size_t(decodedBytes)));
        require(compressedBytes > 0 && compressedBytes <= packedBound &&
                    compressedBytes <= fileSize - position(),
                "compressed block length exceeds bounds: " + file);
        std::vector<unsigned char> packed(size_t(compressedBytes), 0), raw;
        bytes(reinterpret_cast<char *>(packed.data()), packed.size());
        if (compression == 0) raw = std::move(packed);
        else {
            raw.resize(size_t(decodedBytes));
            if (compression == 1) {
                z_stream stream{};
                require(inflateInit(&stream) == Z_OK, "zlib initialization failed");
                stream.next_in = packed.data();
                stream.avail_in = uInt(packed.size());
                stream.next_out = raw.data();
                stream.avail_out = uInt(raw.size());
                const int result = inflate(&stream, Z_FINISH);
                const bool exact = result == Z_STREAM_END && stream.total_in == packed.size() &&
                                   stream.total_out == raw.size();
                inflateEnd(&stream);
                require(exact, "zlib stream integrity or decoded length mismatch: " + file);
            } else {
                const size_t frame = ZSTD_findFrameCompressedSize(packed.data(), packed.size());
                require(!ZSTD_isError(frame) && frame == packed.size(),
                        "zstd frame is truncated or has trailing bytes: " + file);
                const size_t result = ZSTD_decompress(raw.data(), raw.size(), packed.data(), packed.size());
                require(!ZSTD_isError(result) && result == raw.size(),
                        "zstd decoded length/integrity mismatch: " + file);
            }
        }
        decode(raw, row);
        ++current;
        return true;
    }

private:
    std::ifstream input;
    std::string file;
    uint64_t budget, fileSize = 0;
    uint32_t samples = 0, current = 0, compression = 0;
    std::unordered_set<std::string> variants;

    uint64_t position() {
        const auto pos = input.tellg();
        require(pos >= std::streampos(0), "file read/seek failed: " + file);
        return uint64_t(pos);
    }
    void bytes(char *target, size_t count) {
        require(count <= fileSize - position(), "truncated block/string: " + file);
        input.read(target, std::streamsize(count));
        require(bool(input), "input read failed: " + file);
    }
    uint64_t integer(unsigned count) {
        unsigned char data[8]{};
        bytes(reinterpret_cast<char *>(data), count);
        uint64_t value = 0;
        for (unsigned i = 0; i < count; ++i) value |= uint64_t(data[i]) << (8 * i);
        return value;
    }
    std::string string(unsigned lengthBytes, uint64_t maximum) {
        const uint64_t size = integer(lengthBytes);
        require(size <= maximum && size <= fileSize - position(), "invalid string length: " + file);
        std::string result(size_t(size), '\0');
        if (size) bytes(result.data(), size_t(size));
        return result;
    }
    void decode(const std::vector<unsigned char> &raw, PosteriorRow &row) {
        size_t cursor = 0;
        auto take = [&](unsigned count) {
            require(count <= raw.size() - cursor, "truncated probability header: " + file);
            uint64_t result = 0;
            for (unsigned i = 0; i < count; ++i) result |= uint64_t(raw[cursor++]) << (8 * i);
            return result;
        };
        require(take(4) == samples && take(2) == 2, "probability dimensions mismatch: " + file);
        const unsigned minimumH = unsigned(take(1)), maximumH = unsigned(take(1));
        require(minimumH <= maximumH && maximumH <= 2, "ploidy H must be known 0/1/2: " + file);
        std::vector<unsigned char> hapcount(samples), absent(samples);
        uint64_t probabilityCount = 0;
        unsigned actualMinimum = 2, actualMaximum = 0;
        for (uint32_t i = 0; i < samples; ++i) {
            const unsigned value = unsigned(take(1)), h = value & 63u;
            require((value & 64u) == 0 && h >= minimumH && h <= maximumH,
                    "invalid ploidy/reserved bit: " + file);
            hapcount[i] = static_cast<unsigned char>(h);
            absent[i] = static_cast<unsigned char>((value & 128u) != 0);
            probabilityCount += h;
            actualMinimum = std::min(actualMinimum, h);
            actualMaximum = std::max(actualMaximum, h);
        }
        require(actualMinimum == minimumH && actualMaximum == maximumH,
                "ploidy extrema do not match sample values: " + file);
        require(take(1) == 0, "phased probabilities are outside the known-H profile: " + file);
        const unsigned bits = unsigned(take(1));
        require(bits >= 1 && bits <= 32, "probability precision must be 1..32 bits: " + file);
        const uint64_t totalBits = probabilityCount * bits;
        require(raw.size() - cursor == (totalBits + 7) / 8,
                "probability payload length mismatch: " + file);
        if (totalBits % 8) {
            const unsigned mask = (1u << unsigned(totalBits % 8)) - 1;
            require((raw.back() & ~mask) == 0, "nonzero probability padding bits: " + file);
        }
        const uint64_t denominator = (uint64_t(1) << bits) - 1;
        uint64_t accumulator = 0;
        unsigned available = 0;
        auto probability = [&]() {
            while (available < bits) {
                accumulator |= uint64_t(raw[cursor++]) << available;
                available += 8;
            }
            const uint64_t value = accumulator & denominator;
            accumulator >>= bits;
            available -= bits;
            return value;
        };
        row.dosage.resize(samples);
        row.hapcount.resize(samples);
        row.uncertain = false;
        for (uint32_t i = 0; i < samples; ++i) {
            const unsigned h = hapcount[i];
            uint64_t sum = 0;
            std::array<uint64_t, 3> p = {0, 0, 0};
            for (unsigned g = 0; g < h; ++g) {
                p[g] = probability();
                require(!absent[i] || p[g] == 0, "missing posterior has nonzero stored values: " + file);
                sum += p[g];
            }
            require(sum <= denominator, "stored posterior probabilities sum above one: " + file);
            p[h] = denominator - sum;
            row.hapcount[i] = h;
            if (absent[i] && h != 0) row.dosage[i] = missing;
            else {
                row.dosage[i] = (double(p[1]) + 2 * double(p[2])) / double(denominator);
                unsigned support = 0;
                for (unsigned g = 0; g <= h; ++g) support += p[g] != 0;
                row.uncertain |= support > 1;
            }
        }
    }
};

} // namespace

struct LocalBgenCollection::Impl {
    Names ids, labels;
    std::string counted;
    std::vector<std::unique_ptr<AncestryBgen>> files;
    std::vector<std::vector<size_t>> order;
    uint64_t current = 0;

    Impl(const std::string &manifest, uint64_t budget) {
        TextReader input(manifest);
        std::string text;
        require(input.getline(text) && splitTabs(text) == Names({"ancestry", "bgen", "counted_allele"}),
                "manifest requires TAB header ancestry bgen counted_allele");
        while (input.getline(text)) {
            const auto fields = splitTabs(text);
            require(fields.size() == 3, "manifest requires exactly three TAB fields per ancestry");
            require(!fields[0].empty() && fields[0] != "global" &&
                        fields[0].find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") == std::string::npos,
                    "unsafe/reserved ancestry label");
            labels.push_back(fields[0]);
            require(fields[2] == "ALT" || fields[2] == "REF", "counted allele must be ALT or REF");
            if (counted.empty()) counted = fields[2];
            require(counted == fields[2], "all ancestries must declare the same counted direction");
            files.push_back(std::make_unique<AncestryBgen>(relativeTo(manifest, fields[1]), budget));
        }
        requireUnique(labels, "ancestry label");
        require(!files.empty() && files.size() <= 65535, "ancestry count must be 1..65535");
        ids = files[0]->ids;
        for (const auto &file : files) {
            require(file->markerCount == files[0]->markerCount, "ancestry variant counts differ");
            require(file->ids.size() == ids.size(), "ancestry sample sets differ");
            std::unordered_map<std::string, size_t> lookup;
            for (size_t i = 0; i < file->ids.size(); ++i) lookup.emplace(file->ids[i], i);
            std::vector<size_t> align;
            align.reserve(ids.size());
            for (const auto &id : ids) {
                const auto it = lookup.find(id);
                require(it != lookup.end(), "ancestry file omits IID " + id);
                align.push_back(it->second);
            }
            order.push_back(std::move(align));
        }
    }

    bool next(LocalVariant &out) {
        std::vector<PosteriorRow> rows(files.size());
        size_t available = 0;
        for (size_t a = 0; a < files.size(); ++a) available += files[a]->next(rows[a]);
        if (available == 0) return false;
        require(available == files.size(), "ancestry files have different row counts");
        out.metadata = rows[0].metadata;
        out.index = current++;
        out.hasDosage = false;
        out.genotype.resize(files.size());
        out.hapcount.resize(files.size());
        for (size_t a = 0; a < files.size(); ++a) {
            require(rows[a].metadata == out.metadata, "ancestry variant metadata/order differ");
            out.genotype[a].resize(ids.size());
            out.hapcount[a].resize(ids.size());
            out.hasDosage |= rows[a].uncertain;
            for (size_t i = 0; i < ids.size(); ++i) {
                const size_t j = order[a][i];
                const double h = rows[a].hapcount[j], alt = rows[a].dosage[j];
                out.hapcount[a][i] = h;
                out.genotype[a][i] = alt;
            }
        }
        for (size_t i = 0; i < ids.size(); ++i) {
            double total = 0;
            for (const auto &h : out.hapcount) total += h[i];
            require(total == 2, "complete ancestry haplotype counts must sum to two");
        }
        return true;
    }
};

LocalBgenCollection::LocalBgenCollection(const std::string &path, uint64_t budget)
    : impl(std::make_unique<Impl>(path, budget)) {}
LocalBgenCollection::~LocalBgenCollection() = default;
const Names &LocalBgenCollection::sampleIds() const { return impl->ids; }
const Names &LocalBgenCollection::ancestries() const { return impl->labels; }
const std::string &LocalBgenCollection::sourceCountedAllele() const { return impl->counted; }
uint64_t LocalBgenCollection::nVariants() const { return impl->files[0]->markerCount; }
bool LocalBgenCollection::next(LocalVariant &row) { return impl->next(row); }

} // namespace spagxe_local
