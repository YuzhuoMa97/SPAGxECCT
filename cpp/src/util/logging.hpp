// logging.hpp — Thread-safe timestamped stderr logger
//
// Usage:  infoMsg("processed %d markers", n);
// Output: [INFO] 2026-03-28 20:05:56 processed 236 markers
#pragma once

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <mutex>

namespace logging_detail {
inline std::mutex &mutex() {
    static std::mutex m;
    return m;
}

inline void writeStamped(const char *level, const char *msg) {
    std::lock_guard<std::mutex> lk(mutex());
    auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    char ts[20];
    std::strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    std::fprintf(stderr, "[%s] %s %s\n", level, ts, msg);
}

inline void writeFormatted(const char *level, const char *fmt, va_list args) {
    char message[512];
    std::vsnprintf(message, sizeof(message), fmt, args);
    writeStamped(level, message);
}
} // namespace logging_detail

// Function-local static ensures the mutex is initialized exactly once
// even across translation units (C++17 magic statics are thread-safe).
inline void infoMsg(
    const char *fmt,
    ...
) {
    va_list args;
    va_start(args, fmt);
    logging_detail::writeFormatted("INFO", fmt, args);
    va_end(args);
}

// Same shape as infoMsg, but tagged [WARN] for non-fatal anomalies the user
// should notice (e.g. solver failed to converge, fell back to a default).
inline void warnMsg(
    const char *fmt,
    ...
) {
    va_list args;
    va_start(args, fmt);
    logging_detail::writeFormatted("WARN", fmt, args);
    va_end(args);
}
