// SPDX-License-Identifier: GPL-3.0-or-later
// Native paths at the text/local I/O boundary. Copyright (c) 2026 Yuzhuo Ma.
#pragma once

#include <atomic>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <cwchar>
#endif

namespace pathio {
namespace detail {
// -1: API auto detection; 0: legacy native Windows argv; 1: Unicode-aware
// UTF-8 argv. Shared by all threads and translation units in this process.
inline std::atomic<int> inputEncoding{-1};
}

inline void setInputEncodingUtf8(bool utf8) noexcept {
    detail::inputEncoding.store(utf8 ? 1 : 0, std::memory_order_relaxed);
}
inline void setInputEncodingAuto() noexcept {
    detail::inputEncoding.store(-1, std::memory_order_relaxed);
}

// A CLI sets its declared byte encoding once before worker threads start.
// API callers retain auto detection (UTF-8 first, ACP fallback). Manifest
// entry decoding below is always explicitly UTF-8, independent of this mode.
inline std::filesystem::path nativePath(const std::string &value) {
    if (value.find('\0') != std::string::npos)
        throw std::invalid_argument("A file path cannot contain an embedded NUL.");
#ifdef _WIN32
    if (value.empty()) return {};
    if (value.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("File path is too long to decode.");
    const int bytes = static_cast<int>(value.size());
    const int encoding = detail::inputEncoding.load(std::memory_order_relaxed);
    unsigned int page = encoding == 0 ? CP_ACP : CP_UTF8;
    unsigned long flags = encoding == 0 ? 0 : MB_ERR_INVALID_CHARS;
    int length = MultiByteToWideChar(page, flags, value.data(), bytes, nullptr, 0);
    if (!length && encoding < 0) {
        page = CP_ACP;
        flags = 0;
        length = MultiByteToWideChar(page, flags, value.data(), bytes, nullptr, 0);
    }
    if (!length) throw std::invalid_argument("Cannot decode the file path.");
    std::wstring wide(static_cast<size_t>(length), L'\0');
    if (MultiByteToWideChar(page, flags, value.data(), bytes, wide.data(), length) != length)
        throw std::invalid_argument("Cannot decode the file path.");
    return std::filesystem::path(wide);
#else
    return std::filesystem::path(value);
#endif
}

// Join native parent and UTF-8 child as native paths before converting back.
// This avoids concatenating ACP parent bytes with UTF-8 manifest entry bytes.
inline std::string resolveManifestPath(const std::string &manifest,
                                       const std::string &entry) {
    if (entry.empty()) throw std::invalid_argument("Empty file path in manifest.");
    if (entry.find('\0') != std::string::npos)
        throw std::invalid_argument("A file path cannot contain an embedded NUL.");
    std::filesystem::path child = std::filesystem::u8path(entry);
    if (!child.is_absolute()) child = nativePath(manifest).parent_path() / child;
    return child.lexically_normal().u8string();
}

inline FILE *openFile(const std::string &path, const char *mode) {
#ifdef _WIN32
    const auto native = nativePath(path);
    const std::wstring wideMode(mode, mode + std::strlen(mode));
    return ::_wfopen(native.c_str(), wideMode.c_str());
#else
    return std::fopen(path.c_str(), mode);
#endif
}

inline int removeFile(const std::string &path) {
#ifdef _WIN32
    return ::_wremove(nativePath(path).c_str());
#else
    return std::remove(path.c_str());
#endif
}

} // namespace pathio
