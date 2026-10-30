// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
// Independent SQRLOC01 integer-count writer. This is a data conversion tool;
// it does not alter the association statistic or infer local ancestry.
#include "local/local_pack.hpp"
#include "util/path_io.hpp"
#include <zlib.h>
#include <zstd.h>
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace spagxe_local {
namespace {
namespace fs = std::filesystem;
void require(bool ok, const std::string &message) {
    if (!ok) throw std::runtime_error("Local count packing: " + message);
}
std::vector<std::string> fields(const std::string &line) {
    std::istringstream in(line); std::vector<std::string> out; std::string x;
    while (in >> x) out.push_back(x);
    return out;
}
void put(std::ostream &out, uint64_t value, unsigned bytes) {
    for (unsigned i=0; i<bytes; ++i) out.put(char((value >> (8*i)) & 255u));
    require(bool(out), "binary write failed");
}
uint64_t product(uint64_t a, uint64_t b) {
    require(!b || a <= std::numeric_limits<uint64_t>::max()/b, "size overflow");
    return a*b;
}
uint32_t position(const std::string &text) {
    require(!text.empty() && text.find_first_not_of("0123456789")==std::string::npos,
            "BIM POS must be a positive integer");
    uint64_t x=0;
    for (char c:text) {
        require(x <= (uint64_t(UINT32_MAX)-unsigned(c-'0'))/10, "BIM POS exceeds uint32");
        x=x*10+unsigned(c-'0');
    }
    require(x>0, "BIM POS must be positive");
    return static_cast<uint32_t>(x);
}
uint32_t crcFile(const fs::path &file) {
    std::ifstream in(file,std::ios::binary); require(bool(in), "cannot checksum sidecar");
    char buffer[65536]; uLong result=crc32(0L,Z_NULL,0);
    while (in) { in.read(buffer,sizeof(buffer)); const auto n=in.gcount();
        if (n) result=crc32(result,reinterpret_cast<const Bytef*>(buffer),static_cast<uInt>(n)); }
    require(in.eof(), "sidecar checksum read failed");
    return static_cast<uint32_t>(result);
}
void copyBytes(const fs::path &from, const fs::path &to) {
    std::ifstream in(from,std::ios::binary); require(bool(in), "cannot read supplied sidecar");
    std::ofstream out(to,std::ios::binary); require(bool(out), "cannot create temporary sidecar");
    char buffer[65536];
    while (in) { in.read(buffer,sizeof(buffer)); if (in.gcount()) out.write(buffer,in.gcount()); }
    require(in.eof(), "sidecar snapshot read failed"); out.flush();
    require(bool(out), "sidecar snapshot write failed"); out.close();
    require(bool(out), "sidecar snapshot close failed");
}
// The same-filesystem temporary directory is exclusively owned. Remove only
// the known files it contains; never recursively remove a caller's directory.
struct Transaction {
    fs::path dir;
    std::vector<fs::path> published;
    bool complete=false;
    explicit Transaction(const fs::path &prefix) {
        auto parent=prefix.parent_path(); if (parent.empty()) parent=fs::path(".");
        require(fs::is_directory(parent), "output parent directory must exist");
        static std::atomic<uint64_t> sequence{0};
        for (unsigned tries=0; tries<100; ++tries) {
            const auto tick=std::chrono::high_resolution_clock::now().time_since_epoch().count();
            dir=parent/(prefix.filename().native()+fs::path(".pack-tmp-"+std::to_string(tick)+"-"+std::to_string(sequence++)).native());
            std::error_code ec;
            if (fs::create_directory(dir,ec)) return;
            require(!ec || ec==std::errc::file_exists, "cannot reserve temporary output directory");
        }
        throw std::runtime_error("Local count packing: could not reserve unique temporary files");
    }
    ~Transaction() {
        std::error_code ec;
        if (!complete) for (auto i=published.rbegin(); i!=published.rend(); ++i) fs::remove(*i,ec);
        for (const char *name:{"counts","bim","fam"}) fs::remove(dir/name,ec);
        fs::remove(dir,ec);
    }
    // Exclusive creation is used rather than rename-overwrite. If any write
    // fails, the destructor removes only destinations created by this call.
    void publish(const fs::path &from, const fs::path &to) {
#ifdef _WIN32
        int fd=_wopen(to.c_str(),_O_CREAT|_O_EXCL|_O_WRONLY|_O_BINARY,_S_IREAD|_S_IWRITE);
#else
        int fd=::open(to.c_str(),O_CREAT|O_EXCL|O_WRONLY,0666);
#endif
        require(fd>=0, "cannot exclusively create output (it may already exist): "+std::string(std::strerror(errno)));
        published.push_back(to);
#ifdef _WIN32
        FILE *file=_fdopen(fd,"wb");
#else
        FILE *file=fdopen(fd,"wb");
#endif
        if (!file) {
#ifdef _WIN32
            _close(fd);
#else
            ::close(fd);
#endif
            throw std::runtime_error("Local count packing: cannot open reserved output stream");
        }
        try {
            std::ifstream in(from,std::ios::binary); require(bool(in), "cannot reopen temporary output");
            char buffer[65536];
            while (in) { in.read(buffer,sizeof(buffer)); const auto n=in.gcount();
                if (n) require(std::fwrite(buffer,1,static_cast<size_t>(n),file)==static_cast<size_t>(n),"output copy failed"); }
            require(in.eof(), "temporary output read failed");
            require(std::fflush(file)==0, "output flush failed");
            const int status=std::fclose(file); file=nullptr; require(status==0,"output close failed");
        } catch (...) { if(file) std::fclose(file); throw; }
    }
};
unsigned countCode(double value, const char *kind) {
    if (std::isnan(value)) return 3;
    require(std::isfinite(value)&&value>=0&&value<=2&&value==std::floor(value),
            std::string(kind)+" must be integer 0/1/2 or missing; fractional dosage cannot be packed");
    return static_cast<unsigned>(value);
}
void closeChecked(std::ofstream &out,const char *name) {
    out.flush(); require(bool(out),std::string(name)+" flush failed");
    out.close(); require(bool(out),std::string(name)+" close failed");
}
} // namespace

void packLocalCounts(const PackOptions &options) {
    require(!options.outPrefix.empty(), "--out is required");
    require(options.blockVariants>=1&&options.blockVariants<=65536, "block variants must be 1..65536");
    require(options.compressionLevel>=-5&&options.compressionLevel<=22, "compression level must be -5..22");
    const fs::path destination[3]={pathio::nativePath(options.outPrefix+".localbed"),
                                 pathio::nativePath(options.outPrefix+".localbim"),
                                 pathio::nativePath(options.outPrefix+".localfam")};
    for (const auto &file:destination) require(!fs::exists(file), "output already exists; choose a new prefix");
    LocalSource source(options.input);
    Transaction transaction(pathio::nativePath(options.outPrefix));
    const auto tempBed=transaction.dir/"counts",tempBim=transaction.dir/"bim",tempFam=transaction.dir/"fam";
    std::vector<std::string> ids=source.sampleIds();
    if (!options.famFile.empty()) {
        copyBytes(pathio::nativePath(options.famFile),tempFam);
        std::ifstream in(tempFam,std::ios::binary); require(bool(in), "cannot read FAM snapshot");
        std::string line; ids.clear(); std::unordered_set<std::string> seen;
        while (std::getline(in,line)) { const auto row=fields(line);
            require(row.size()==6,"FAM requires six headerless FID IID PAT MAT SEX PHENO fields");
            require(seen.insert(row[1]).second,"duplicate IID in supplied FAM"); ids.push_back(row[1]); }
        require(in.eof()&&!ids.empty(),"empty or unreadable supplied FAM");
        require(ids.size()==source.sampleIds().size(),"supplied FAM must contain the same complete sample set as source"); source.align(ids);
    } else {
        std::ofstream fam(tempFam,std::ios::binary); require(bool(fam),"cannot write FAM");
        for (const auto &id:ids) {
            require(!id.empty()&&id.find_first_of(" \t\r\n\v\f")==std::string::npos,
                    "count-container IID cannot contain whitespace");
            fam<<"0\t"<<id<<"\t0\t0\t0\t-9\n";
        }
        closeChecked(fam,"FAM");
    }
    const uint64_t n=ids.size(),k=source.ancestries().size();
    require(n>0&&n<=UINT32_MAX&&k>0&&k<=65535,"invalid sample/ancestry count");
    const uint64_t plane=(n+3)/4,stride=product(product(2,k),plane);
    const uint64_t capacity=product(stride,options.blockVariants);
    require(capacity<=options.input.maxBlockBytes&&capacity<=std::numeric_limits<size_t>::max(),
            "requested decoded block exceeds memory guard; reduce block variants or raise the guard");
    std::ifstream suppliedBim;
    std::ofstream generatedBim;
    if (!options.bimFile.empty()) {
        copyBytes(pathio::nativePath(options.bimFile),tempBim);
        suppliedBim.open(tempBim,std::ios::binary); require(bool(suppliedBim),"cannot read BIM snapshot");
    } else { generatedBim.open(tempBim,std::ios::binary); require(bool(generatedBim),"cannot write BIM"); }
    std::ofstream bed(tempBed,std::ios::binary); require(bool(bed),"cannot write temporary count container");
    bed.write("SQRLOC01",8); put(bed,n,8); put(bed,0,8); put(bed,k,4);
    put(bed,options.blockVariants,4); put(bed,0,8); put(bed,0,1); // canonical ALT
    for (const auto &label:source.ancestries()) {
        require(!label.empty()&&label.size()<=65535,"invalid ancestry label length");
        put(bed,label.size(),2);bed.write(label.data(),std::streamsize(label.size()));
    }
    std::vector<unsigned char> raw; raw.reserve(static_cast<size_t>(capacity));
    std::vector<uint64_t> offsets; uint64_t markers=0;
    auto flush=[&]() {
        if (raw.empty()) return;
        require(raw.size()%stride==0,"internal count-plane size mismatch");
        offsets.push_back(static_cast<uint64_t>(bed.tellp()));
        std::vector<unsigned char> packed(ZSTD_compressBound(raw.size()));
        const size_t length=ZSTD_compress(packed.data(),packed.size(),raw.data(),raw.size(),options.compressionLevel);
        require(!ZSTD_isError(length),"zstd count-block compression failed");
        put(bed,raw.size()/stride,4);put(bed,raw.size(),8);put(bed,length,8);
        put(bed,crc32(0L,raw.data(),static_cast<uInt>(raw.size())),4);
        bed.write(reinterpret_cast<const char*>(packed.data()),std::streamsize(length));
        require(bool(bed),"compressed block write failed");raw.clear();
    };
    LocalVariant row;
    while (source.next(row)) {
        require(!row.hasDosage,"packing requires integer calls/point posteriors; uncertain or fractional genotypes are not hard-called");
        require(row.genotype.size()==k&&row.hapcount.size()==k,"ancestry dimensions changed");
        if (suppliedBim.is_open()) {
            std::string line; require(bool(std::getline(suppliedBim,line)),"supplied BIM has too few variants");
            const auto f=fields(line);require(f.size()==6,"BIM requires CHROM ID CM POS ALT REF");
            require(f[0]==row.metadata[0]&&f[1]==row.metadata[2]&&position(f[3])==position(row.metadata[1])&&
                    f[4]==row.metadata[4]&&f[5]==row.metadata[3],
                    "supplied BIM variant order/coordinates/ALT/REF differ from source");
        } else generatedBim<<row.metadata[0]<<'\t'<<row.metadata[2]<<"\t0\t"<<row.metadata[1]<<'\t'<<row.metadata[4]<<'\t'<<row.metadata[3]<<'\n';
        const size_t begin=raw.size();raw.resize(begin+static_cast<size_t>(stride),0);
        for (size_t i=0;i<n;++i) {
            unsigned sumH=0;bool missingH=false;
            for (size_t a=0;a<k;++a) {
                require(row.genotype[a].size()==static_cast<Eigen::Index>(n)&&row.hapcount[a].size()==static_cast<Eigen::Index>(n),"sample dimensions changed");
                unsigned gc=countCode(row.genotype[a][i],"G"),hc=countCode(row.hapcount[a][i],"H");
                if (hc==0&&gc==3) gc=0;
                require(gc==3||hc==3||gc<=hc,"G exceeds H");missingH|=hc==3;if(hc!=3)sumH+=hc;
                const unsigned shift=2*unsigned(i%4);
                raw[begin+(2*a)*plane+i/4]|=static_cast<unsigned char>(gc<<shift);
                raw[begin+(2*a+1)*plane+i/4]|=static_cast<unsigned char>(hc<<shift);
            }
            require(missingH||sumH==2,"complete ancestry H must sum to two");
        }
        require(markers<std::numeric_limits<uint64_t>::max(),"too many variants");++markers;
        if (markers%options.blockVariants==0) flush();
    }
    require(markers>0,"source has no variants");flush();
    if (suppliedBim.is_open()) {
        std::string extra; require(!std::getline(suppliedBim,extra)&&suppliedBim.eof(),"supplied BIM has extra variants or cannot be read");
        suppliedBim.close();
    } else closeChecked(generatedBim,"BIM");
    const uint64_t footer=static_cast<uint64_t>(bed.tellp());
    put(bed,offsets.size(),8);for (auto offset:offsets) put(bed,offset,8);
    put(bed,crcFile(tempBim),4);put(bed,crcFile(tempFam),4);
    bed.seekp(16);require(bool(bed),"cannot patch variant count");put(bed,markers,8);
    bed.seekp(32);require(bool(bed),"cannot patch footer offset");put(bed,footer,8);closeChecked(bed,"count container");
    // Publish metadata before the completed count file. No existing path can
    // be replaced even if another caller creates it after initial checks.
    transaction.publish(tempBim,destination[1]);transaction.publish(tempFam,destination[2]);
    transaction.publish(tempBed,destination[0]);transaction.complete=true;
}
} // namespace spagxe_local
