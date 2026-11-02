// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE standalone entry point. See the retained provenance notices.
#include "cli/cli.hpp"
#include "util/path_io.hpp"
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <shellapi.h>
#include <cwchar>
#include <string>
#include <vector>
namespace {
bool unicodeLocalMode(int argc, wchar_t **argv) {
    for (int i=1;i<argc;++i) {
        if (std::wcscmp(argv[i],L"--pack-local")==0) return true;
        if (std::wcscmp(argv[i],L"--method")==0 && i+1<argc &&
            (_wcsicmp(argv[i+1],L"SPAGxEmix_local")==0 ||
             _wcsicmp(argv[i+1],L"SPAGxElocal")==0)) return true;
    }
    return false;
}
std::string utf8Argument(const wchar_t *value) {
    const int bytes=WideCharToMultiByte(CP_UTF8,0,value,-1,nullptr,0,nullptr,nullptr);
    if (bytes<=0) return {};
    std::string result(static_cast<size_t>(bytes),'\0');
    WideCharToMultiByte(CP_UTF8,0,value,-1,result.data(),bytes,nullptr,nullptr);
    result.pop_back(); return result;
}
}
#endif
int main(int argc, char *argv[]) {
#ifdef _WIN32
    // Local streams use the Unicode-aware I/O layer. Preserve wide command-line
    // paths here instead of losing non-ACP characters before the parser runs.
    // Other genotype workflows retain their existing native argv contract.
    int wideCount=0;
    wchar_t **wideArgs=CommandLineToArgvW(GetCommandLineW(),&wideCount);
    if (wideArgs && unicodeLocalMode(wideCount,wideArgs)) {
        pathio::setInputEncodingUtf8(true);
        std::vector<std::string> values; values.reserve(static_cast<size_t>(wideCount));
        for (int i=0;i<wideCount;++i) values.push_back(utf8Argument(wideArgs[i]));
        LocalFree(wideArgs);
        std::vector<char *> pointers; pointers.reserve(values.size());
        for (auto &value:values) pointers.push_back(value.data());
        return cli::run(static_cast<int>(pointers.size()),pointers.data());
    }
    if (wideArgs) LocalFree(wideArgs);
    pathio::setInputEncodingUtf8(false);
#endif
    return cli::run(argc,argv);
}
