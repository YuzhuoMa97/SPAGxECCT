// SPDX-License-Identifier: GPL-3.0-or-later
// SPAGxE native local ancestry workflow, Yuzhuo Ma <yuzhuoma@amss.ac.cn>.
// Null-model orchestration follows SPAGxE's GPL-3.0 inherited infrastructure;
// original source attribution remains in THIRD_PARTY_NOTICES.
#include "util/version.hpp"
#include "local/local_workflow.hpp"
#include "local/local_kernel.hpp"
#include "io/subject_data.hpp"
#include "genotype/variant_selection.hpp"
#include "genotype/equilibrium.hpp"
#include "util/null_model.hpp"
#include "util/logging.hpp"
#include "util/text_stream.hpp"
#include "util/probability_math.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace spagxe_local {
namespace {
constexpr double NA=std::numeric_limits<double>::quiet_NaN();
using TT=interaction_wald::TraitType;
void need(bool b,const std::string &s){if(!b)throw std::runtime_error("SPAGxE_local: "+s);}
std::string numeric(double d) {if(std::isnan(d))return "NA";if(std::isinf(d))return d>0?"Inf":"-Inf";char b[48];std::snprintf(b,sizeof(b),"%.17g",d);return b;}
TT traitType(const SubjectData &sd,const nullmodel::PhenoSpec &spec,nullmodel::RegressionModel model) {
    using RM=nullmodel::RegressionModel;if(nullmodel::isCoxSpec(spec))return TT::Cox;
    if(model==RM::Auto)model=nullmodel::inferModelFromColumn(sd.getColumn(spec.yColumn),spec.yColumn,sd.usedIIDs()).model;
    switch(model){case RM::Linear:return TT::Linear;case RM::Logistic:return TT::Logistic;case RM::Cox:return TT::Cox;case RM::Ordinal:return TT::Ordinal;default:return TT::None;}
}
Eigen::VectorXd recode(Eigen::VectorXd y,TT t) {
    if (t == TT::Linear) return y;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (Eigen::Index i = 0; i < y.size(); ++i) {
        if (std::isfinite(y[i])) {
            lo = std::min(lo, y[i]);
            hi = std::max(hi, y[i]);
        }
    }
    for (Eigen::Index i = 0; i < y.size(); ++i) {
        if (!std::isfinite(y[i])) continue;
        if (t == TT::Logistic) y[i] = y[i] > (lo + hi) * 0.5 ? 1 : 0;
        else if (t == TT::Ordinal) y[i] = std::round(y[i] - lo);
    }
    return y;
}
Eigen::VectorXd subset(const Eigen::VectorXd &v,const std::vector<size_t> &ix) {Eigen::VectorXd out(ix.size());for(size_t i=0;i<ix.size();++i)out[i]=v[ix[i]];return out;}
Eigen::MatrixXd subset(const Eigen::MatrixXd &m,const std::vector<size_t> &ix) {Eigen::MatrixXd out(ix.size(),m.cols());for(size_t i=0;i<ix.size();++i)out.row(i)=m.row(ix[i]);return out;}
std::shared_ptr<const interaction_wald::WaldData> subsetWald(const std::shared_ptr<const interaction_wald::WaldData> &w,const std::vector<size_t> &ix) {
    if (!w) return {};
    auto o = std::make_shared<interaction_wald::WaldData>();
    o->trait = w->trait;
    o->covar = subset(w->covar, ix);
    if(w->trait==TT::Cox){o->time=subset(w->time,ix);o->event=subset(w->event,ix);}else o->y=subset(w->y,ix);return o;
}
struct Task {
    PerPhenoInfo mask;Eigen::VectorXd resid;std::vector<Eigen::VectorXd> env;
    std::shared_ptr<const interaction_wald::WaldData> wald;
    std::unique_ptr<LocalKernel> kernel;
};
std::unique_ptr<LocalKernel> makeKernel(const Task &t,const LocalWorkflowOptions &o,const std::vector<size_t> *ix=nullptr) {
    if(!ix)return std::make_unique<LocalKernel>(t.resid,o.envNames,t.env,t.wald,o.marginalCutoff,o.spaCutoff);
    std::vector<Eigen::VectorXd> env;for(const auto &e:t.env)env.push_back(subset(e,*ix));return std::make_unique<LocalKernel>(subset(t.resid,*ix),o.envNames,std::move(env),subsetWald(t.wald,*ix),o.marginalCutoff,o.spaCutoff);
}
struct VariantLines{std::string local,combined;};
std::string meta(const LocalVariant &v){return v.metadata[0]+"\t"+v.metadata[1]+"\t"+v.metadata[2]+"\t"+v.metadata[3]+"\t"+v.metadata[4];}
VariantLines analyze(const LocalVariant &v,const Task &t,const std::vector<std::string> &labels,const LocalWorkflowOptions &o) {
    const size_t k=labels.size();std::vector<Eigen::VectorXd> g,h;g.reserve(k);h.reserve(k);
    for(size_t a=0;a<k;++a){g.push_back(extractPhenoVec(v.genotype[a],t.mask));h.push_back(extractPhenoVec(v.hapcount[a],t.mask));}
    // Unknown H is not imputed. All ancestry tests for a variant use the
    // same complete-H sample set; G missingness remains per ancestry.
    std::vector<size_t> complete;complete.reserve(t.resid.size());
    for(Eigen::Index i=0;i<t.resid.size();++i){bool ok=true;for(const auto &hh:h)ok&=std::isfinite(hh[i]);if(ok)complete.push_back(size_t(i));}
    const size_t n=complete.size(),hMissing=size_t(t.resid.size())-n;
    std::unique_ptr<LocalKernel> temporary;const LocalKernel *kernel=t.kernel.get();
    if(hMissing&&n>=2){temporary=makeKernel(t,o,&complete);kernel=temporary.get();}
    if(hMissing)for(size_t a=0;a<k;++a){g[a]=subset(g[a],complete);h[a]=subset(h[a],complete);}
    Eigen::MatrixXd allH(n,k);for(size_t a=0;a<k;++a)allH.col(a)=h[a];
    // HWE is a diploid global diagnostic. A Bin(H,q) local branch with
    // variable H has no ordinary diploid HWE test.
    uint32_t homAlt=0,het=0,homRef=0;for(size_t i=0;i<n;++i){double total=0;bool observed=true;for(size_t a=0;a<k;++a){observed&=std::isfinite(g[a][i]);total+=g[a][i];}if(observed){if(total==0)++homRef;else if(total==1)++het;else if(total==2)++homAlt;}}
    double hwe=!v.hasDosage&&(homAlt+het+homRef)>0?equilibriumNegLog10P(het,homAlt,homRef):NA;
    const bool hwePass=o.hweCutoff==0||(!std::isnan(hwe)&&hwe<=-std::log10(o.hweCutoff));
    std::vector<std::vector<double>> cct(o.envNames.size());VariantLines lines;
    for(size_t a=0;a<k;++a){
        double observedG=0,observedH=0,totalH=0;size_t missingG=0;
        for(size_t i=0;i<n;++i){totalH+=h[a][i];if(std::isfinite(g[a][i])){observedG+=g[a][i];observedH+=h[a][i];}else if(h[a][i]>0)++missingG;}
        double q=observedH>0?observedG/observedH:NA,maf=std::isfinite(q)?std::min(q,1-q):NA,mac=maf*totalH,missing=n?double(missingG)/n:NA;
        const bool pass=n>=2&&std::isfinite(q)&&missing<=o.missingCutoff&&maf>=o.minMafCutoff&&mac>=o.minMacCutoff&&hwePass;
        std::vector<double> values(size_t(t.kernel->resultSize()),NA);bool flipped=std::isfinite(q)&&q>0.5;
        if(pass){auto result=kernel->evaluate(g[a],h[a],allH);values=std::move(result.values);q=result.altFreq;mac=result.mac;missing=result.missingRate;flipped=result.flipped;}
        const char *qc=pass?"PASS":(o.hweCutoff>0&&v.hasDosage?"FAIL_HWE_UNAVAILABLE":"FAIL_QC");
        lines.local+=meta(v)+"\t"+labels[a]+"\t"+(flipped?"REF":"ALT")+"\t"+(flipped?"1":"0")+"\t"+std::to_string(n)+"\t"+std::to_string(hMissing)+"\t"+numeric(missing)+"\t"+numeric(q)+"\t"+numeric(mac)+"\t"+numeric(hwe)+"\t"+qc;
        for (double value : values) lines.local += "\t" + numeric(value);
        lines.local += '\n';
        for(size_t e=0;e<cct.size();++e)if(pass&&!std::isnan(values[5+9*e]))cct[e].push_back(values[5+9*e]);
    }
    lines.combined=meta(v)+"\t"+std::to_string(n)+"\t"+std::to_string(hMissing);
    for(const auto &p:cct){double combined=p.empty()?NA:math::cauchyCombineLog10(p.data(),int(p.size()));lines.combined+="\t"+numeric(combined)+"\t"+std::to_string(p.size());}lines.combined+='\n';return lines;
}
} // namespace

void scanHomogeneousInteractionsLocal(const LocalWorkflowOptions &o) {
    need(!o.phenoFile.empty()&&!o.envNames.empty()&&!o.outPrefix.empty(),"phenotype, environment and output are required");
    need(o.nthreads>=1&&o.nSnpPerChunk>=1,"positive threads/chunk size required");
    need(o.missingCutoff>=0&&o.missingCutoff<=1&&o.minMafCutoff>=0&&o.minMafCutoff<=0.5&&o.minMacCutoff>=0&&o.hweCutoff>=0&&o.hweCutoff<=1,"invalid QC threshold");
    LocalSource source(o.input);infoMsg("SPAGxE_local: %zu ancestries; source allele=%s, canonical input=ALT",source.ancestries().size(),source.sourceCountedAllele().c_str());
    const bool fit=!o.phenoNameSpec.empty();auto model=nullmodel::parseRegressionModel(o.regressionModel);std::vector<nullmodel::PhenoSpec> specs;
    if(fit)specs=nullmodel::parsePhenoSpecList(model,o.phenoNameSpec);
    SubjectData sd(source.sampleIds());sd.setGenoLabel("local G/H");
    if(fit)sd.loadPhenoFile(o.phenoFile,nullmodel::columnsNeeded(specs));else {sd.loadResidOne(o.phenoFile,o.residNames);sd.loadPhenoFile(o.phenoFile);}
    if (!o.covarFile.empty()) sd.loadCovar(o.covarFile, o.covarNames);
    sd.setKeepRemove(o.keepFile, o.removeFile);
    sd.finalize();
    sd.dropNaInColumns(o.envNames);
    const uint32_t N=sd.nUsed();need(N>=2,"at least two usable samples required");source.align(sd.usedIIDs());
    Eigen::MatrixXd covar;std::vector<TT> traits;
    if(fit){for(const auto &env:o.envNames)need(std::find(o.covarNames.begin(),o.covarNames.end(),env)!=o.covarNames.end(),"environment must enter the null model: "+env);covar=o.covarNames.empty()?Eigen::MatrixXd(N,0):sd.getColumns(o.covarNames);for(const auto &s:specs)traits.push_back(traitType(sd,s,model));nullmodel::EngineOptions e;e.nthreads=o.nthreads;auto fits=nullmodel::fitAll(sd,specs,model,covar,e);if(o.saveResid)nullmodel::writeResidualsFile(o.outPrefix+".null.resid",sd,fits,o.compression,o.compressionLevel);std::vector<Eigen::VectorXd> residuals;std::vector<std::string> names;for(auto &f:fits){residuals.push_back(std::move(f.residuals));names.push_back(std::move(f.name));}sd.setResidualsFromFit(std::move(residuals),std::move(names));}
    auto masks=sd.buildPerColumnMasks();std::vector<Task> tasks(masks.size());std::vector<std::unique_ptr<TextWriter>> outputs,combined;
    std::vector<Eigen::VectorXd> envUnion;for(const auto &e:o.envNames)envUnion.push_back(sd.getColumn(e));
    for(size_t p=0;p<masks.size();++p){auto &t=tasks[p];t.mask=masks[p];Eigen::VectorXd r=masks.size()==1?sd.residuals():Eigen::VectorXd(sd.residMatrix().col(p));t.resid=extractPhenoVec(r,t.mask);for(const auto &e:envUnion)t.env.push_back(extractPhenoVec(e,t.mask));
        if(fit){auto wd=std::make_shared<interaction_wald::WaldData>();wd->trait=traits[p];wd->covar=extractPhenoMat(covar,t.mask);if(traits[p]==TT::Cox){wd->time=extractPhenoVec(sd.getColumn(specs[p].timeColumn),t.mask);wd->event=extractPhenoVec(sd.getColumn(specs[p].eventColumn),t.mask);}else wd->y=extractPhenoVec(recode(sd.getColumn(specs[p].yColumn),traits[p]),t.mask);t.wald=std::move(wd);}
        t.kernel=makeKernel(t,o);
        auto outPath=TextWriter::buildOutputPath(o.outPrefix,t.mask.name,"SPAGxE_local",o.compression),combinedPath=TextWriter::buildOutputPath(o.outPrefix,t.mask.name,"SPAGxE_local_CCT",o.compression);
        outputs.push_back(std::make_unique<TextWriter>(outPath,TextWriter::modeFromString(o.compression),o.compressionLevel));combined.push_back(std::make_unique<TextWriter>(combinedPath,TextWriter::modeFromString(o.compression),o.compressionLevel));
        outputs.back()->write("CHROM\tPOS\tID\tREF\tALT\tANCESTRY\tTESTED_ALLELE\tMINOR_FLIPPED\tN_TESTED\tN_H_MISSING\tMISS_RATE\tALT_FREQ\tMAC\tLOG10P_HWE_GLOBAL\tQC_STATUS"+t.kernel->getHeaderColumns()+"\n");std::string header="CHROM\tPOS\tID\tREF\tALT\tN_TESTED\tN_H_MISSING";for(const auto &env:o.envNames)header+="\tLOG10P_ANCESTRY_CCT_Gx"+env+"\tN_ANCESTRIES_Gx"+env;combined.back()->write(header+"\n");
        infoMsg("SPAGxE_local phenotype '%s': %u subjects",t.mask.name.c_str(),t.mask.nUsed);
    }
    const auto include=o.extractFile.empty()?std::unordered_set<std::string>{}:spagxe_genotype::readVariantSelection(o.extractFile),exclude=o.excludeFile.empty()?std::unordered_set<std::string>{}:spagxe_genotype::readVariantSelection(o.excludeFile);
    // Retained decoded genotypes per text processing batch are at most
    // approximately 64 MiB (one variant is always permitted). Binary reader
    // decoded-block budget is separate and enforced before allocation.
    const uint64_t perVariant=uint64_t(N)*source.ancestries().size()*2*sizeof(double);
    const size_t batchSize=std::max<size_t>(1,std::min<uint64_t>(o.nSnpPerChunk,(64ULL<<20)/std::max<uint64_t>(perVariant,1)));
    uint64_t read=0,selected=0;bool eof=false;
    while(!eof){std::vector<LocalVariant> batch;batch.reserve(batchSize);while(batch.size()<batchSize){LocalVariant v;if(!source.next(v)){eof=true;break;}++read;if((!o.chromosomes.empty()&&!o.chromosomes.count(v.metadata[0]))||(!o.extractFile.empty()&&!include.count(v.metadata[2]))||exclude.count(v.metadata[2]))continue;batch.push_back(std::move(v));++selected;}
        if (batch.empty()) continue;
        std::vector<std::vector<VariantLines>> lines(tasks.size(), std::vector<VariantLines>(batch.size()));
        std::atomic<size_t> next{0};
        std::exception_ptr failure;
        std::mutex errorMutex;
        auto worker=[&](){try{while(true){size_t j=next.fetch_add(1);if(j>=batch.size()*tasks.size())break;size_t p=j/batch.size(),v=j%batch.size();lines[p][v]=analyze(batch[v],tasks[p],source.ancestries(),o);}}catch(...){std::lock_guard<std::mutex> lock(errorMutex);if(!failure)failure=std::current_exception();next.store(batch.size()*tasks.size());}};
        const size_t workers=std::min<size_t>(o.nthreads,batch.size()*tasks.size());std::vector<std::thread> pool;for(size_t j=1;j<workers;++j)pool.emplace_back(worker);worker();for(auto &w:pool)w.join();if(failure)std::rethrow_exception(failure);
        for (size_t p = 0; p < tasks.size(); ++p) {
            for (const auto &line : lines[p]) {
                outputs[p]->write(line.local);
                combined[p]->write(line.combined);
            }
        }
        infoMsg("SPAGxE_local: processed %llu selected variants (%llu read)",
                static_cast<unsigned long long>(selected),
                static_cast<unsigned long long>(read));
    }
    need(selected>0,"no variants remain after chromosome/extract/exclude filters");for(auto &out:outputs)out->close();for(auto &out:combined)out->close();
    TextWriter receipt(o.outPrefix+".local.settings.tsv");receipt.write("key\tvalue\nsoftware\tSPAGxE\nversion\t" SPAGXE_VERSION "\ncgf\tguarded_partial_normal\ncentral_IQR_multiplier\t1.5\npartial_guard_tolerance\t0.001\npartial_fallback\tboth_tails_full_binomial\nmethod\tSPAGxEmix_local\nsoftware_algorithm_owner\tyuzhuoma@amss.ac.cn\nsource_counted_allele\t"+source.sourceCountedAllele()+"\ncanonical_input_counted_allele\tALT\ntest_direction\tancestry-specific minor allele; see TESTED_ALLELE\nunknown_H_policy\tcomplete-H samples shared across all ancestries per variant\nmissing_G_policy\tancestry mean H*q; H=0 is structural zero\nHWE_policy\tglobal diploid hardcall diagnostic, not ancestry-specific HWE\nancestry_CCT\tequal-weight per-ancestry final GxE P values; failed QC/tests omitted, count reported\nvariants_read\t"+std::to_string(read)+"\nvariants_selected\t"+std::to_string(selected)+"\nthreads\t"+std::to_string(o.nthreads)+"\nprocessing_batch_variants\t"+std::to_string(batchSize)+"\n");receipt.close();
}
} // namespace spagxe_local
