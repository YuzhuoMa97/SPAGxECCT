// variant_filter.cpp — Shared --extract / --exclude logic for all readers

#include "genotype/variant_selection.hpp"

#include "util/text_scanner.hpp"
#include "util/text_stream.hpp"

#include <fstream>
#include <stdexcept>

namespace spagxe_genotype {

std::unordered_set<std::string> readVariantSelection(const std::string &path) {
    TextReader in(path);

    std::unordered_set<std::string> ids;
    std::string line;
    while (in.getline(line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        text::TokenScanner ts(line);
        auto tok = ts.nextView();
        if (tok.empty()) continue;
        ids.emplace(tok);
    }
    return ids;
}

void applyVariantSelection(
    std::vector<GenotypeCatalog::MarkerInfo> &markers,
    const std::string &extractFile,
    const std::string &excludeFile
) {
    const bool anyExtract = !extractFile.empty();
    const bool anyExclude = !excludeFile.empty();
    if (!anyExtract && !anyExclude) return;

    std::unordered_set<std::string> includeSet, excludeSet;
    if (anyExtract) includeSet = readVariantSelection(extractFile);
    if (anyExclude) excludeSet = readVariantSelection(excludeFile);

    std::vector<GenotypeCatalog::MarkerInfo> filtered;
    filtered.reserve(markers.size());
    for (auto &m : markers) {
        if (anyExtract && includeSet.find(m.id) == includeSet.end()) continue;
        if (anyExclude && excludeSet.find(m.id) != excludeSet.end()) continue;
        filtered.push_back(std::move(m));
    }
    markers = std::move(filtered);

    if (markers.empty())
        throw std::runtime_error("No markers remain after --extract/--exclude");
}

} // namespace spagxe_genotype
