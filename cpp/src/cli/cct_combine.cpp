// SPDX-License-Identifier: GPL-3.0-or-later
// Native global/local Cauchy combination with strict marker/trait alignment.
// Copyright (c) 2026 Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
#include "cli/cct_combine.hpp"
#include "util/text_stream.hpp"
#include "util/logging.hpp"
#include "util/probability_math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace spagxe_cct {
namespace {
const double kNaN=std::numeric_limits<double>::quiet_NaN();
void need(bool b,const std::string &s){if(!b)throw std::runtime_error("SPAGxEcct: "+s);}
bool ends(const std::string&s,const std::string &tail){return s.size()>=tail.size()&&s.compare(s.size()-tail.size(),tail.size(),tail)==0;}
std::string numeric(double d){if(std::isnan(d))return "NA";if(std::isinf(d))return "Inf";char b[48];std::snprintf(b,sizeof(b),"%.17g",d);return b;}
std::vector<std::string> tokens(const std::string &line){std::istringstream in(line);std::vector<std::string> out;std::string x;while(in>>x)out.push_back(x);return out;}
double magnitude(const std::string &s,const std::string &context){
    if(s=="NA"||s=="NaN"||s=="nan"||s==".")return kNaN;
    size_t consumed=0;double v;
    try{v=std::stod(s,&consumed);}catch(...){throw std::runtime_error("SPAGxEcct: invalid -log10(P) in "+context+": "+s);}
    need(consumed==s.size()&&!std::isnan(v)&&v>=0,"invalid -log10(P) in "+context+": "+s);
    return v;
}
uint64_t count(const std::string&s,const std::string &context){
    need(!s.empty()&&s.front()!='-',"invalid ancestry count in "+context);
    size_t consumed=0;unsigned long long result;
    try{result=std::stoull(s,&consumed);}catch(...){throw std::runtime_error("SPAGxEcct: invalid ancestry count in "+context);}
    need(consumed==s.size(),"invalid ancestry count in "+context);return result;
}
void checkIdentity(const Options&o){
    need(!o.phenotype.empty()&&o.phenotype.find_first_of("/\\\t\r\n ")==std::string::npos,
         "--cct-phenotype must be a nonempty filename-safe phenotype name.");
    auto stem=[](const std::string&path){std::string s=std::filesystem::path(path).filename().string();
        if(ends(s,".gz"))s.resize(s.size()-3);else if(ends(s,".zst"))s.resize(s.size()-4);return s;};
    const auto global=stem(o.globalFile),local=stem(o.localFile);
    need(ends(global,"."+o.phenotype+".SPAGxE")||ends(global,"."+o.phenotype+".SPAGxEmix"),
         "global filename must end in .PHENOTYPE.SPAGxE, .SPAGxEmix (before compression suffix).");
    need(ends(local,"."+o.phenotype+".SPAGxE_local_CCT"),
         "local filename must end in .PHENOTYPE.SPAGxE_local_CCT (before compression suffix).");
}
struct Table {
    TextReader reader;std::string path;std::vector<std::string> header;
    std::unordered_map<std::string,size_t> positions;
    explicit Table(const std::string&p):reader(p),path(p){std::string line;need(reader.getline(line),"empty input "+p);header=tokens(line);
        for(size_t j=0;j<header.size();++j)need(positions.emplace(header[j],j).second,"duplicate header column "+header[j]+" in "+p);
        for(const auto *key:{"CHROM","POS","ID","REF","ALT"})column(key);
    }
    size_t column(const std::string &name)const{auto it=positions.find(name);need(it!=positions.end(),"missing column "+name+" in "+path);return it->second;}
    bool next(std::vector<std::string>&row){std::string line;if(!reader.getline(line))return false;row=tokens(line);need(row.size()==header.size(),"row has wrong column count in "+path);return true;}
    std::array<std::string,5> key(const std::vector<std::string>&r)const{return {r[column("CHROM")],r[column("POS")],r[column("ID")],r[column("REF")],r[column("ALT")]};}
};
// Same P==1 -> .999 endpoint convention as the original R CCT. A signed
// log-scaled sum supports tails at both P~0 and P~1 without forming 1/P.
double weighted(double a,double b,double wg,double wl){
    if(std::isnan(a)||std::isnan(b))return kNaN;
    need(!(std::isinf(a)&&b==0)&&!(std::isinf(b)&&a==0),"cannot combine exact P=0 and P=1 (singular Cauchy endpoints)");
    const double magnitudes[2]={a,b}, weights[2]={wg,wl};
    return math::cauchyCombineWeightedLog10(magnitudes,weights,2);
}
uint64_t process(const Options&o,TextWriter *writer){
    Table global(o.globalFile),local(o.localFile);
    std::vector<size_t> gc,lc,nc;
    for(const auto&e:o.envNames){gc.push_back(global.column("LOG10P_Gx"+e));
        lc.push_back(local.column("LOG10P_ANCESTRY_CCT_Gx"+e));nc.push_back(local.column("N_ANCESTRIES_Gx"+e));}
    std::unordered_set<std::string> seen;std::vector<std::string> gr,lr;uint64_t rows=0;
    while(true){const bool g=global.next(gr),l=local.next(lr);need(g==l,"global/local files have different numbers of variant rows");if(!g)break;
        const auto key=global.key(gr);need(key==local.key(lr),"variant metadata/order mismatch at row "+std::to_string(rows+1));
        std::string encoded;for(const auto &v:key){encoded+=std::to_string(v.size())+":"+v;}
        need(seen.insert(encoded).second,"duplicate five-field variant key at row "+std::to_string(rows+1));
        std::string output;for(const auto&v:key){if(!output.empty())output+='\t';output+=v;}output+='\t'+o.phenotype;
        for(size_t e=0;e<o.envNames.size();++e){const std::string context="row "+std::to_string(rows+1)+", environment "+o.envNames[e];
            const double gp=magnitude(gr[gc[e]],context),lp=magnitude(lr[lc[e]],context);const auto n=count(lr[nc[e]],context);
            need(n>0||std::isnan(lp),"finite local CCT with zero contributing ancestries in "+context);
            double wg=o.globalWeight,wl=o.localWeight;
            if(o.perAncestryEqual){wg=1.;wl=static_cast<double>(n);}
            const double rawWg=wg,rawWl=wl;
            const double maxWeight=std::max(wg,wl);wg/=maxWeight;wl/=maxWeight;const double sum=wg+wl;wg/=sum;wl/=sum;
            const bool mg=std::isnan(gp),ml=std::isnan(lp);std::string status=mg?(ml?"MISSING_BOTH":"MISSING_GLOBAL"):(ml?"MISSING_LOCAL":"OK");
            const double combined=weighted(gp,lp,rawWg,rawWl);
            output+='\t'+numeric(gp)+'\t'+numeric(lp)+'\t'+numeric(combined)+'\t'+numeric(wg)+'\t'+numeric(wl)+'\t'+std::to_string(n)+'\t'+status;
        }
        if(writer)writer->write(output+"\n");
        ++rows;
    }need(rows>0,"input tables contain no variant rows");return rows;
}
} // namespace
void run(const Options&o){
    need(!o.globalFile.empty()&&!o.localFile.empty()&&!o.outPrefix.empty()&&!o.envNames.empty(),"global/local inputs, phenotype, environments and output are required");
    need(std::isfinite(o.globalWeight)&&std::isfinite(o.localWeight)&&o.globalWeight>0&&o.localWeight>0,"CCT weights must be finite and strictly positive");
    checkIdentity(o);
    const auto path=TextWriter::buildOutputPath(o.outPrefix,o.phenotype,"SPAGxE_global_local_CCT",o.compression);
    for(const auto &input:{o.globalFile,o.localFile})need(std::filesystem::absolute(input).lexically_normal()!=std::filesystem::absolute(path).lexically_normal(),"output path cannot equal an input path");
    // Validate all rows before opening the output, preserving an existing
    // result when mismatched input would otherwise fail partway through.
    const uint64_t verified=process(o,nullptr);TextWriter::assertWritable(path);
    TextWriter writer(path,TextWriter::modeFromString(o.compression),o.compressionLevel);
    std::string header="CHROM\tPOS\tID\tREF\tALT\tPHENOTYPE";
    for(const auto&e:o.envNames)header+="\tLOG10P_GLOBAL_Gx"+e+"\tLOG10P_LOCAL_CCT_Gx"+e+
        "\tLOG10P_GLOBAL_LOCAL_CCT_Gx"+e+"\tWEIGHT_GLOBAL_Gx"+e+"\tWEIGHT_LOCAL_Gx"+e+
        "\tN_ANCESTRIES_Gx"+e+"\tCCT_STATUS_Gx"+e;
    writer.write(header+"\n");need(process(o,&writer)==verified,"inputs changed between validation and output passes");writer.close();
    infoMsg("SPAGxEcct: %llu matched variants; phenotype=%s; weights=%s; output=%s",
        static_cast<unsigned long long>(verified),o.phenotype.c_str(),o.perAncestryEqual?"per-ancestry equal":"explicit group weights",path.c_str());
}
} // namespace spagxe_cct
