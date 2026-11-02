// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE native local ancestry I/O. Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
// Independently implemented native G/H count and posterior-expectation I/O.
#include "util/path_io.hpp"
#include "local/local_io.hpp"
#include "local/local_bgen.hpp"
#include "util/text_stream.hpp"
#include <zlib.h>
#include <zstd.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace spagxe_local {
namespace {
constexpr double NA = std::numeric_limits<double>::quiet_NaN();
void need(bool ok, const std::string &msg) { if (!ok) throw std::runtime_error("Local ancestry: " + msg); }
std::vector<std::string> tabs(const std::string &s) {
    std::vector<std::string> out; size_t start=0;
    while (true) { size_t end=s.find('\t',start); out.push_back(s.substr(start,end==std::string::npos?end:end-start)); if(end==std::string::npos) break; start=end+1; }
    if(!out.empty()&&!out.back().empty()&&out.back().back()=='\r') out.back().pop_back();
    return out;
}
std::vector<std::string> spaces(const std::string &s) { std::istringstream in(s);std::vector<std::string> out;std::string v;while(in>>v)out.push_back(v);return out; }
void unique(const std::vector<std::string> &v,const std::string &what) {
    std::unordered_set<std::string> seen;
    for(const auto &s:v) need(!s.empty()&&seen.insert(s).second,"empty or duplicate "+what+": "+s);
}
bool safeLabel(const std::string &s) { return !s.empty()&&s!="global"&&s.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")==std::string::npos; }
size_t product(uint64_t a,uint64_t b) {need(!b||a<=std::numeric_limits<size_t>::max()/b,"dimension overflow");return size_t(a*b);}
uint64_t get(std::istream &in,int bytes) {uint64_t x=0;for(int j=0;j<bytes;++j){int c=in.get();need(c!=EOF,"truncated localbed");x|=uint64_t(uint8_t(c))<<(8*j);}return x;}
std::string label(std::istream &in) {size_t n=get(in,2);std::string s(n,'\0');in.read(s.data(),n);need(bool(in),"truncated ancestry label");return s;}
uint32_t crcFile(const std::string &p) {
    std::ifstream in(pathio::nativePath(p),std::ios::binary);need(bool(in),"cannot checksum "+p);
    char buffer[65536];uLong c=crc32(0,nullptr,0);while(in){in.read(buffer,sizeof(buffer));if(in.gcount())c=crc32(c,reinterpret_cast<const Bytef*>(buffer),uInt(in.gcount()));}need(in.eof(),"sidecar checksum read failed");return uint32_t(c);
}
void validMeta(const std::array<std::string,5> &v,const std::string &where) {
    for(const auto &s:v)need(!s.empty()&&s.find_first_of(" \t\r\n")==std::string::npos,"empty/whitespace metadata at "+where);
    need(v[1].find_first_not_of("0123456789")==std::string::npos,"noninteger POS at "+where);
    size_t end=0;unsigned long long pos=0;try {pos=std::stoull(v[1],&end);}catch(...){need(false,"invalid POS at "+where);}
    need(end==v[1].size()&&pos>0&&pos<=UINT32_MAX,"POS out of range at "+where);
    need(v[3]!=v[4]&&v[3].find(',')==std::string::npos&&v[4].find(',')==std::string::npos,"distinct biallelic REF/ALT required at "+where);
}
double count(const std::string &s,const std::string &where,bool dosage) {
    if(s.empty()||s=="NA"||s=="."||s=="-9"||s=="NaN"||s=="nan"||s=="NAN")return NA;
    size_t end=0;double d=0;try{d=std::stod(s,&end);}catch(...){need(false,"invalid count at "+where);}
    need(end==s.size()&&std::isfinite(d)&&d>=0&&d<=2&&(dosage||d==std::floor(d)),"G must be finite in [0,H] and H integer 0/1/2 (or missing) at "+where);return d;
}
struct TableStream {
    std::string path;TextReader in;std::vector<std::string> ids;size_t line=1;std::unordered_set<std::string> seen;bool dosage;
    explicit TableStream(const std::string &p,bool allowDosage):path(p),in(p),dosage(allowDosage) {
        std::string s;need(in.getline(s),"empty G/H stream: "+p);auto f=tabs(s);
        need(f.size()>5&&(f[0]=="CHROM"||f[0]=="CHR")&&f[1]=="POS"&&f[2]=="ID"&&f[3]=="REF"&&f[4]=="ALT","G/H stream requires TAB header CHROM POS ID REF ALT followed by IID: "+p);
        ids.assign(f.begin()+5,f.end());unique(ids,"G/H IID");
    }
    bool next(std::array<std::string,5> &meta,Eigen::VectorXd &values) {
        std::string s;if(!in.getline(s))return false;++line;auto f=tabs(s);std::string where=path+":"+std::to_string(line);
        need(f.size()==ids.size()+5,"G/H column count differs from header at "+where);
        std::copy_n(f.begin(),5,meta.begin());validMeta(meta,where);need(seen.insert(meta[2]).second,"duplicate variant "+meta[2]+" in "+path);
        values.resize(ids.size());for(size_t i=0;i<ids.size();++i)values[i]=count(f[i+5],where,dosage);return true;
    }
};
}

struct LocalSource::Impl {
    LocalInputOptions opt;std::vector<std::string> ids,labels;std::string counted;
    std::vector<size_t> selected;
    std::vector<std::unique_ptr<TableStream>> tables;
    std::vector<std::vector<size_t>> tableOrder;
    std::unique_ptr<LocalBgenCollection> posterior;
    std::ifstream binary;
    uint64_t n=0,m=0,current=0,footer=0,loaded=UINT64_MAX;
    size_t block=0,plane=0,stride=0;
    std::vector<uint64_t> offsets;
    std::vector<unsigned char> raw;
    std::vector<std::array<std::string,5>> metadata;
    explicit Impl(const LocalInputOptions &o):opt(o) {
        const int inputs = int(!o.manifest.empty()) + int(!o.binaryPrefix.empty()) + int(!o.bgenManifest.empty());
        need(inputs == 1, "supply exactly one local G/H, binary-count or BGEN-collection input");
        need(o.maxBlockBytes>=1024&&o.maxBlockBytes<=uint64_t(1)<<30,"max-block-bytes must be 1024..1073741824");
        std::string bgen = o.bgenManifest;
        if (!o.manifest.empty()) {
            TextReader probe(o.manifest);
            std::string first;
            need(probe.getline(first), "empty local manifest");
            if (tabs(first) == std::vector<std::string>({"ancestry", "bgen", "counted_allele"})) bgen = o.manifest;
        }
        if (!bgen.empty()) {
            posterior = std::make_unique<LocalBgenCollection>(bgen, o.maxBlockBytes);
            ids = posterior->sampleIds();
            labels = posterior->ancestries();
            counted = posterior->sourceCountedAllele();
            n = ids.size();
            m = posterior->nVariants();
        } else if (!o.manifest.empty()) openText();
        else openBinary();
        align(ids);
    }
    void align(const std::vector<std::string> &wanted) {
        unique(wanted,"selected IID");need(!wanted.empty(),"no selected samples");std::unordered_map<std::string,size_t> lookup;
        for (size_t i = 0; i < ids.size(); ++i) {
            lookup.emplace(ids[i], i);
        }
        selected.clear();
        selected.reserve(wanted.size());
        for(const auto &id:wanted){auto it=lookup.find(id);need(it!=lookup.end(),"missing IID "+id);selected.push_back(it->second);}
    }
    void openText() {
        TextReader manifest(opt.manifest);std::string s;need(manifest.getline(s),"empty local manifest");
        need(tabs(s)==std::vector<std::string>({"ancestry","genotype","hapcount","counted_allele"}),"manifest header must be ancestry<TAB>genotype<TAB>hapcount<TAB>counted_allele");
        size_t l=1;
        while(manifest.getline(s)) {
            ++l;auto f=tabs(s);need(f.size()==4,"manifest needs four TAB columns at row "+std::to_string(l));
            need(safeLabel(f[0]),"unsafe/reserved ancestry label: "+f[0]);labels.push_back(f[0]);
            need(f[3]=="ALT"||f[3]=="REF","counted_allele must be ALT or REF");
            if (counted.empty()) counted = f[3];
            need(counted == f[3], "mixed counted-allele directions across ancestry files");
            for(int j=1;j<=2;++j){need(!f[j].empty(),"empty G/H path");const std::string p=pathio::resolveManifestPath(opt.manifest,f[j]);tables.push_back(std::make_unique<TableStream>(p,j==1));}
        }
        unique(labels,"ancestry");need(!labels.empty()&&labels.size()<=65535,"ancestry count must be 1..65535");
        ids=tables[0]->ids;n=ids.size();need(n<=UINT32_MAX,"too many local samples");
        std::unordered_map<std::string,size_t> base;for(size_t i=0;i<ids.size();++i)base.emplace(ids[i],i);
        for(const auto &t:tables) {need(t->ids.size()==ids.size(),"G/H sample sets differ");std::vector<size_t> order;order.reserve(n);for(const auto &id:t->ids){auto it=base.find(id);need(it!=base.end(),"G/H sample sets differ at IID "+id);order.push_back(it->second);}tableOrder.push_back(std::move(order));}
    }
    void openBinary() {
        const auto &p=opt.binaryPrefix;binary.open(pathio::nativePath(p+".localbed"),std::ios::binary);need(bool(binary),"cannot open "+p+".localbed");
        char magic[8];binary.read(magic,8);need(binary.gcount()==8&&std::string(magic,8)=="SQRLOC01","invalid SQRLOC01 localbed magic/version");
        n=get(binary,8);m=get(binary,8);uint64_t k=get(binary,4);block=size_t(get(binary,4));footer=get(binary,8);auto direction=get(binary,1);
        need(n>0&&n<=UINT32_MAX&&m>0&&m<=std::numeric_limits<size_t>::max()&&k>0&&k<=65535&&block>0&&direction<=1,"invalid binary dimensions or counted allele");counted=direction?"REF":"ALT";
        plane=size_t((n+3)/4);stride=product(2*k,plane);
        for(uint64_t a=0;a<k;++a){labels.push_back(label(binary));need(safeLabel(labels.back()),"unsafe ancestry label");}unique(labels,"ancestry");
        uint64_t start=uint64_t(binary.tellg());binary.seekg(0,std::ios::end);need(bool(binary),"cannot seek binary input");uint64_t end=uint64_t(binary.tellg());
        need(footer>=start&&footer<end&&end-footer>=16,"invalid footer location");binary.seekg(std::streamoff(footer));auto nb=get(binary,8);
        need(nb==1+(m-1)/block&&nb<=(end-footer-16)/8&&footer+16+nb*8==end,"invalid footer index length");
        for(uint64_t b=0;b<nb;++b){auto off=get(binary,8);need(off>=start&&off<footer&&(offsets.empty()||off>offsets.back()),"invalid block offset");offsets.push_back(off);}offsets.push_back(footer);
        need(get(binary,4)==crcFile(p+".localbim"),"localbim checksum mismatch");need(get(binary,4)==crcFile(p+".localfam"),"localfam checksum mismatch");
        TextReader fam(p+".localfam");std::string line;while(fam.getline(line)){auto f=spaces(line);need(f.size()==6,"localfam requires six headerless columns");ids.push_back(f[1]);}unique(ids,"localfam IID");need(ids.size()==n,"localfam sample count mismatch");
        TextReader bim(p+".localbim");std::unordered_set<std::string> seen;while(bim.getline(line)){auto f=spaces(line);need(f.size()==6,"localbim requires CHROM ID CM POS ALT REF");std::array<std::string,5> v={f[0],f[3],f[1],f[5],f[4]};validMeta(v,p+".localbim");need(seen.insert(v[2]).second,"duplicate localbim variant "+v[2]);metadata.push_back(std::move(v));}need(metadata.size()==m,"localbim variant count mismatch");
    }
    void loadBlock(size_t b) {
        binary.clear();binary.seekg(std::streamoff(offsets[b]));need(bool(binary),"binary seek failed");auto nr=get(binary,4),bytes=get(binary,8),compressed=get(binary,8),crc=get(binary,4);
        need(nr==std::min<uint64_t>(block,m-uint64_t(b)*block)&&bytes==product(nr,stride),"invalid binary raw block dimensions");
        need(bytes<=opt.maxBlockBytes,"decoded block exceeds --local-max-block-bytes; repack smaller blocks or explicitly raise the budget");
        need(offsets[b+1]>=offsets[b]+24&&compressed>0&&compressed==offsets[b+1]-offsets[b]-24&&compressed<=ZSTD_compressBound(size_t(bytes)),"invalid compressed block extent");
        std::vector<unsigned char> packed(static_cast<size_t>(compressed));binary.read(reinterpret_cast<char*>(packed.data()),std::streamsize(compressed));need(bool(binary),"truncated compressed local block");
        auto frame=ZSTD_findFrameCompressedSize(packed.data(),packed.size());need(!ZSTD_isError(frame)&&frame==packed.size(),"invalid/multiple/trailing zstd frame");raw.resize(size_t(bytes));
        auto decoded=ZSTD_decompress(raw.data(),raw.size(),packed.data(),packed.size());need(!ZSTD_isError(decoded)&&decoded==bytes,"local zstd decompression failed");
        need(crc32(0,raw.data(),uInt(bytes))==crc,"local block CRC32 mismatch");
        if(n%4){const unsigned keep=(1u<<(2*(n%4)))-1u;for(size_t v=0;v<size_t(nr);++v)for(size_t a=0;a<2*labels.size();++a)need((raw[v*stride+a*plane+plane-1]&~keep)==0,"nonzero local padding bits");}
        loaded=b;
    }
    bool next(LocalVariant &out) {
        LocalVariant full;full.index=current;
        const size_t k=labels.size();full.genotype.resize(k);full.hapcount.resize(k);
        if (posterior) {
            if (!posterior->next(full)) return false;
        } else if(!tables.empty()) {
            std::vector<Eigen::VectorXd> rows(2*k);size_t available=0;
            for(size_t t=0;t<tables.size();++t){std::array<std::string,5> v;if(tables[t]->next(v,rows[t])){++available;if(!t)full.metadata=v;else need(v==full.metadata,"G/H variant metadata/order mismatch");}}
            if (!available) return false;
            need(available == tables.size(), "G/H row counts differ");
            for(size_t a=0;a<k;++a){full.genotype[a].resize(n);full.hapcount[a].resize(n);for(size_t i=0;i<n;++i){full.genotype[a][tableOrder[2*a][i]]=rows[2*a][i];full.hapcount[a][tableOrder[2*a+1][i]]=rows[2*a+1][i];}}
        } else {
            if (current == m) return false;
            size_t b = size_t(current / block);
            if (loaded != b) loadBlock(b);
            full.metadata = metadata[size_t(current)];
            size_t base = size_t(current % block) * stride;
            for(size_t a=0;a<k;++a){full.genotype[a].resize(n);full.hapcount[a].resize(n);for(size_t i=0;i<n;++i){unsigned g=(raw[base+2*a*plane+i/4]>>(2*(i%4)))&3,h=(raw[base+(2*a+1)*plane+i/4]>>(2*(i%4)))&3;full.genotype[a][i]=g==3?NA:double(g);full.hapcount[a][i]=h==3?NA:double(h);}}
        }
        for(size_t i=0;i<n;++i){double sum=0;bool missingH=false;for(size_t a=0;a<k;++a){double &g=full.genotype[a][i],h=full.hapcount[a][i];need(std::isnan(g)||std::isnan(h)||g<=h,"G>H at variant "+full.metadata[2]+", IID "+ids[i]);full.hasDosage|=std::isfinite(g)&&g!=std::floor(g);if(h==0&&std::isnan(g))g=0;missingH|=std::isnan(h);sum+=h;if(counted=="REF"&&!posterior)g=std::isnan(h)?NA:h-g;}need(missingH||sum==2,"complete ancestry H must sum to two at variant "+full.metadata[2]+", IID "+ids[i]);}
        out.metadata=full.metadata;out.index=current++;out.hasDosage=full.hasDosage;out.genotype.resize(k);out.hapcount.resize(k);
        for(size_t a=0;a<k;++a){out.genotype[a].resize(selected.size());out.hapcount[a].resize(selected.size());for(size_t i=0;i<selected.size();++i){out.genotype[a][i]=full.genotype[a][selected[i]];out.hapcount[a][i]=full.hapcount[a][selected[i]];}}
        return true;
    }
};
LocalSource::LocalSource(const LocalInputOptions &o):impl(std::make_unique<Impl>(o)){}
LocalSource::~LocalSource()=default;
LocalSource::LocalSource(LocalSource &&) noexcept=default;
LocalSource &LocalSource::operator=(LocalSource &&) noexcept=default;
const std::vector<std::string> &LocalSource::sampleIds()const{return impl->ids;}
const std::vector<std::string> &LocalSource::ancestries()const{return impl->labels;}
const std::string &LocalSource::sourceCountedAllele()const{return impl->counted;}
uint64_t LocalSource::nVariants()const{return impl->m;}
void LocalSource::align(const std::vector<std::string> &iids){impl->align(iids);}
bool LocalSource::next(LocalVariant &v){return impl->next(v);}
} // namespace spagxe_local
