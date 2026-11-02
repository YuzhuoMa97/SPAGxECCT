// mmap_ro.cpp — read-only memory-mapped file, POSIX / Windows backends.

#include "util/mmap_ro.hpp"

#include <stdexcept>
#include <utility>

#if defined(_WIN32)
#  include <windows.h>
#else
#  include <fcntl.h>
#  include <sys/mman.h>
#  include <sys/stat.h>
#  include <unistd.h>
   // MAP_POPULATE pre-faults pages on Linux; macOS / BSD have no
   // equivalent so we fall back to plain MAP_PRIVATE.
#  ifndef MAP_POPULATE
#    define MAP_POPULATE 0
#  endif
#endif

namespace util {

namespace {

#if defined(_WIN32)
class ScopedHandle {
public:
    explicit ScopedHandle(HANDLE value) noexcept : value_(value) {}
    ~ScopedHandle() {
        if (value_) ::CloseHandle(value_);
    }
    ScopedHandle(const ScopedHandle &) = delete;
    ScopedHandle &operator=(const ScopedHandle &) = delete;

    HANDLE get() const noexcept { return value_; }
    HANDLE release() noexcept { return std::exchange(value_, nullptr); }

private:
    HANDLE value_;
};
#else
class ScopedDescriptor {
public:
    explicit ScopedDescriptor(int value) noexcept : value_(value) {}
    ~ScopedDescriptor() { close(); }
    ScopedDescriptor(const ScopedDescriptor &) = delete;
    ScopedDescriptor &operator=(const ScopedDescriptor &) = delete;

    int get() const noexcept { return value_; }
    void close() noexcept {
        if (value_ >= 0) ::close(std::exchange(value_, -1));
    }

private:
    int value_;
};
#endif

} // namespace

MmapReadOnly::MmapReadOnly(const std::string &path) {
#if defined(_WIN32)
    HANDLE opened = ::CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ,
                                  nullptr, OPEN_EXISTING,
                                  FILE_ATTRIBUTE_NORMAL, nullptr);
    if (opened == INVALID_HANDLE_VALUE)
        throw std::runtime_error("Cannot open file: " + path);
    ScopedHandle file(opened);

    LARGE_INTEGER li;
    if (!::GetFileSizeEx(file.get(), &li))
        throw std::runtime_error("Cannot stat file: " + path);
    const std::size_t fileSize = static_cast<std::size_t>(li.QuadPart);
    if (fileSize == 0)
        throw std::runtime_error("Empty file: " + path);

    HANDLE created = ::CreateFileMappingA(file.get(), nullptr, PAGE_READONLY, 0, 0,
                                          nullptr);
    if (!created)
        throw std::runtime_error("Cannot create file mapping: " + path);
    ScopedHandle mapping(created);

    void *mapped = ::MapViewOfFile(mapping.get(), FILE_MAP_READ, 0, 0, 0);
    if (!mapped)
        throw std::runtime_error("Cannot map view of file: " + path);

    addr_  = mapped;
    size_  = fileSize;
    hFile_ = file.release();
    hMap_  = mapping.release();
#else
    const int opened = ::open(path.c_str(), O_RDONLY);
    if (opened < 0)
        throw std::runtime_error("Cannot open file: " + path);
    ScopedDescriptor file(opened);

    struct stat st;
    if (::fstat(file.get(), &st) != 0)
        throw std::runtime_error("Cannot stat file: " + path);
    const std::size_t fileSize = static_cast<std::size_t>(st.st_size);
    if (fileSize == 0)
        throw std::runtime_error("Empty file: " + path);

    void *mapped = ::mmap(nullptr, fileSize, PROT_READ,
                          MAP_PRIVATE | MAP_POPULATE, file.get(), 0);
    file.close();
    if (mapped == MAP_FAILED)
        throw std::runtime_error("Cannot mmap file: " + path);

    ::madvise(mapped, fileSize, MADV_SEQUENTIAL);

    addr_ = mapped;
    size_ = fileSize;
#endif
}

MmapReadOnly::~MmapReadOnly() noexcept {
#if defined(_WIN32)
    if (addr_)  ::UnmapViewOfFile(addr_);
    if (hMap_)  ::CloseHandle(static_cast<HANDLE>(hMap_));
    if (hFile_) ::CloseHandle(static_cast<HANDLE>(hFile_));
#else
    if (addr_) ::munmap(addr_, size_);
#endif
}

} // namespace util
