// Independent local-ancestry score/CGF regression checks.
#include "local/local_kernel.hpp"
#include "util/spa_cgf.hpp"
#include "util/probability_math.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace {
int checks=0;
void require(bool condition,const std::string &name) {
    ++checks;
    if(!condition) throw std::runtime_error(name);
}
void close(double got,double expected,double tol,const std::string &name) {
    require(std::isfinite(got)&&std::isfinite(expected)&&
        std::abs(got-expected)<=tol*(1+std::abs(expected)),name);
}
Eigen::MatrixXd ancestry(const Eigen::VectorXd &h) {
    Eigen::MatrixXd all(h.size(),2); all.col(0)=h;all.col(1)=2-h.array();return all;
}
spagxe_local::LocalResult eval(const Eigen::VectorXd&r,const Eigen::VectorXd&e,
                              const Eigen::VectorXd&g,const Eigen::VectorXd&h,
                              bool branchB,double cutoff=100,bool exact=true,
                              std::shared_ptr<const interaction_wald::WaldData> wd=nullptr) {
    spagxe_local::LocalKernel kernel(r,{"E"},{e},wd,branchB?1:0,cutoff);
    require(kernel.resultSize()==14,"output dimensions");
    require(kernel.getHeaderColumns().find("LOG10P_SPA_GxE")!=std::string::npos,"SPA provenance header");
    return kernel.evaluate(g,h,ancestry(h));
}
void normalOracle() {
    Eigen::VectorXd h(12),g(12),r(12),e(12);
    h<<0,1,2,1,2,0,1,2,1,2,1,2;
    g<<0,0,1,1,0,0,0,1,0,1,0,1;
    r<<.3,-1.4,.8,.2,-.5,.1,1.7,.6,-.4,.5,-.9,.4;
    e<<-.8,.1,.4,1.2,.3,-.1,.5,-1.1,.7,.4,-.2,.6;
    const double q=g.sum()/h.sum();
    for(bool b:{false,true}) {
        const auto got=eval(r,e,g,h,b);
        require(!got.flipped,"minor fixture does not flip");
        close(got.altFreq,q,1e-14,"ancestry AF uses total H");
        close(got.mac,g.sum(),1e-14,"minor count");
        const double mg=(g-q*h).dot(r),vg=q*(1-q)*(h.array()*r.array().square()).sum();
        close(got.values[1],mg/std::sqrt(vg),1e-12,"H-weighted marginal score");
        Eigen::VectorXd weight;
        if(b) {
            // Independent explicit 2x2 normal equations, unlike centered production OLS.
            Eigen::MatrixXd x(g.size(),2);x.col(0).setOnes();x.col(1)=g;
            const Eigen::VectorXd r0=r-x*((x.transpose()*x).inverse()*x.transpose()*r);
            weight=e.array()*r0.array();
        } else {
            const double lambda=(h.array()*e.array()*r.array().square()).sum()/
                                (h.array()*r.array().square()).sum();
            weight=(e.array()-lambda)*r.array();
        }
        const double v=q*(1-q)*(h.array()*weight.array().square()).sum();
        const double s=(g-q*h).dot(weight), z=s/std::sqrt(v);
        close(got.values[8],z,1e-12,"local A/B normal z oracle");
        close(got.values[12],-std::log10(std::erfc(std::abs(z)/std::sqrt(2.))),1e-12,"local normal p");
        close(got.values[9],s/v,1e-12,"local score beta");
        close(got.values[9]/got.values[10],got.values[8],1e-12,"beta/SE");
        close(got.values[13],b?2:1,0,"branch trace");
        const auto swapped=eval(r,e,(h-g).eval(),h,b);
        require(swapped.flipped,"major allele flips H-G");
        for(size_t j=0;j<got.values.size();++j)
            if(std::isfinite(got.values[j])) close(swapped.values[j],got.values[j],2e-12,"complement -> same tested minor allele");
    }
    const auto original=eval(r,e,g,h,false);
    Eigen::VectorXd missing=g;missing[5]=std::numeric_limits<double>::quiet_NaN();
    const auto noTrial=eval(r,e,missing,h,false);
    close(noTrial.missingRate,0,0,"missing G with H0 is known zero");
    close(noTrial.values[8],original.values[8],1e-14,"H0 missing contributes zero");
    missing=g;missing[2]=std::numeric_limits<double>::quiet_NaN();
    const auto imputed=eval(r,e,missing,h,false);
    const double qObs=(g.sum()-g[2])/(h.sum()-h[2]);
    close(imputed.altFreq,qObs,1e-14,"AF excludes missing G trial counts");
    close(imputed.missingRate,1./g.size(),1e-14,"missing G rate");
    Eigen::VectorXd filled=g;filled[2]=h[2]*qObs;
    const auto manual=eval(r,e,filled,h,false);
    close(imputed.values[8],manual.values[8],1e-12,"missing dosage imputed H*q");
}
void cgfEnumeration() {
    Eigen::VectorXd h(5),w(5);h<<0,1,2,1,2;w<<4,-.3,.8,-1.1,.5;
    const double q=.31;
    for(double t:{-.7,-.05,0.,.03,.8}) {
        // Enumerate 3^5 genotype states, assigning zero probability to g>H.
        long double m0=0,m1=0,m2=0;
        for(int code=0;code<243;++code) {
            int rem=code;long double prob=1,score=0;
            for(int i=0;i<5;++i) {
                int g=rem%3;rem/=3;
                if(g>h[i]) {prob=0;break;}
                const double comb=(h[i]==2&&g==1)?2:1;
                prob*=comb*std::pow(q,g)*std::pow(1-q,h[i]-g);
                score+=g*w[i];
            }
            const long double mass=prob*std::exp(static_cast<long double>(t)*score);
            m0+=mass;m1+=mass*score;m2+=mass*score*score;
        }
        const auto got=spa_cgf::binomHapcountKFull(t,w.data(),h.data(),5,q);
        close(got.K0,std::log(m0),3e-13,"enumerated local K");
        close(got.K1,m1/m0,3e-13,"enumerated local K1");
        close(got.K2,m2/m0-(m1/m0)*(m1/m0),3e-13,"enumerated local K2");
    }
}
void fullSpa() {
    // H=1 is a Bernoulli conditional genotype law; the symmetric score has
    // K_center(t)=16*log(cosh(t/2)), with an analytic saddlepoint for score 7.
    const int n=16;
    Eigen::VectorXd h=Eigen::VectorXd::Ones(n),r=h,e(n),g(n);
    for(int i=0;i<n;++i) {e[i]=i<8?-1:1;g[i]=i<8?0:1;}
    g[15]=0;
    // q is estimated at 7/16, rather than forced to .5. Build an independent
    // scalar long-double CGF and bisection oracle for the actual q.
    const auto got=eval(r,e,g,h,false,.01,true);
    const long double q=7.L/16, target=7.L;
    auto one=[&](long double t,int which) {
        long double total=0;
        for(int i=0;i<n;++i) {
            const long double wi=e[i],v=q*std::exp(t*wi),a=1-q+v;
            if(which==0) total+=std::log(a);
            if(which==1) total+=wi*v/a;
            if(which==2) total+=wi*wi*v*(1-q)/(a*a);
        }return total;
    };
    long double lo=0,hi=10;
    for(int k=0;k<150;++k) {const long double mid=(lo+hi)/2;if(one(mid,1)<target)lo=mid;else hi=mid;}
    const long double t=(lo+hi)/2,k0=one(t,0),k2=one(t,2);
    const long double w=std::sqrt(2*(t*target-k0)),v=t*std::sqrt(k2);
    const long double rstar=w+std::log(v/w)/w;
    const double want=-std::log10(std::erfc(static_cast<double>(rstar/std::sqrt(2.L))));
    require(got.values[11]<=2,"full local SPA converges");
    close(got.values[12],want,3e-10,"independent long-double full-CGF SPA");
    // Enumeration quantifies the lattice tail, not an equality claim: SPA is
    // an approximation and should lie between strict and inclusive tails.
    long double strict=0,inclusive=0;
    for(int code=0;code<(1<<n);++code) {
        int count=0;double score=0;
        for(int i=0;i<n;++i) if(code&(1<<i)) {++count;score+=e[i];}
        const long double p=std::pow(q,count)*std::pow(1-q,n-count);
        if(std::abs(score)>7) strict+=p;
        if(std::abs(score)>=7) inclusive+=p;
    }
    const double spaP=std::pow(10.,-got.values[12]);
    require(spaP>strict&&spaP<inclusive,"SPA within strict/inclusive enumerated lattice tails");
}
// Independent scalar CGF and bisection. This deliberately neither calls the
// production support helper nor its shared root solver/BN tail functions.
double scalarFiniteSpa(const Eigen::VectorXd &w,const Eigen::VectorXd &h,
                       long double q,long double observed) {
    long double mean=0,lower=0,upper=0,logLo=0,logHi=0;
    for(int i=0;i<w.size();++i) {
        mean+=q*h[i]*w[i];
        if(w[i]>0) {upper+=h[i]*w[i];logLo+=h[i]*std::log1p(-q);logHi+=h[i]*std::log(q);}
        if(w[i]<0) {lower+=h[i]*w[i];logLo+=h[i]*std::log(q);logHi+=h[i]*std::log1p(-q);}
    }
    auto cumulants=[&](long double t,int order) {
        long double v=0;
        for(int i=0;i<w.size();++i) {
            const long double x=t*w[i],a=1-q+q*std::exp(x),p=q*std::exp(x)/a;
            if(order==0)v+=h[i]*std::log(a);
            if(order==1)v+=h[i]*w[i]*p;
            if(order==2)v+=h[i]*w[i]*w[i]*p*(1-p);
        }return v;
    };
    long double p=0;
    for(bool low:{false,true}) {
        const long double s=mean+(low?-1:1)*std::abs(observed-mean);
        if((low&&s<lower-1e-12L)||(!low&&s>upper+1e-12L))continue;
        if(std::abs(s-(low?lower:upper))<1e-12L) {p+=std::exp(low?logLo:logHi);continue;}
        long double l=-1,u=1;
        while(cumulants(l,1)>s)l*=2;
        while(cumulants(u,1)<s)u*=2;
        for(int j=0;j<180;++j){const long double m=(l+u)/2;if(cumulants(m,1)<s)l=m;else u=m;}
        const long double t=(l+u)/2;
        const long double z=std::copysign(std::sqrt(2*(t*s-cumulants(t,0))),t);
        const long double v=t*std::sqrt(cumulants(t,2));
        const long double rs=z+std::log(v/z)/z;
        p+=.5L*std::erfc((low?-1:1)*rs/std::sqrt(2.L));
    }
    return -std::log10(p);
}
void finiteSupportTails() {
    // Highly asymmetric support: the reflected tail is impossible, while
    // the observed tail is strictly interior and has a well-defined SPA root.
    // Branch A: choose R orthogonal to the requested weight so lambda=0.
    const int n=20;
    Eigen::VectorXd h=Eigen::VectorXd::Ones(n),w(n),r=Eigen::VectorXd::Ones(n),e(n),g=Eigen::VectorXd::Zero(n);
    w[0]=10;
    for(int i=1;i<n;++i)w[i]=i<10?.1:-.1;
    r[0]=.01;e=w.array()/r.array();g[0]=1;
    for(bool reverse:{false,true}) {
        const Eigen::VectorXd env=reverse?-e:e;
        const Eigen::VectorXd weights=reverse?-w:w;
        const auto got=eval(r,env,g,h,false,.01,true);
        const double want=scalarFiniteSpa(weights,h,1.L/n,g.dot(weights));
        require(got.values[11]==0,"support-exterior tail does not cause normal fallback A");
        close(got.values[12],want,3e-9,"finite-support A matches independent scalar CGF");
        require(std::abs(got.values[12]+std::log10(std::erfc(std::abs(got.values[8])/std::sqrt(2.))))>.05,
                "strong-tail SPA differs materially from normal substitution");
    }
    // Branch B: R has zero sum separately among carriers/noncarriers, so the
    // OLS residual is exactly R; E*R has the same asymmetric-tail failure.
    g.setZero();g[0]=g[1]=1;
    for(int i=0;i<n;++i){r[i]=(i==0||(i>=2&&i<11))?1:-1;w[i]=i==0?10:(i==1?1:(i==19?-.1:.1));}
    e=w.array()/r.array();
    const auto b=eval(r,e,g,h,true,.01,true);
    require(b.values[11]==0,"support-exterior tail does not cause normal fallback B");
    close(b.values[12],scalarFiniteSpa(w,h,.1L,g.dot(w)),3e-9,"finite-support B independent scalar CGF");

    // Inclusive endpoint atoms are not zero. Enumerate all 2^6 states of a
    // small Bernoulli score to verify both endpoint masses exactly.
    h=Eigen::VectorXd::Ones(6);w.resize(6);w<<10,.1,.1,-.1,-.1,-.1;
    r=Eigen::VectorXd::Ones(6);r[0]=.01;e=w.array()/r.array();
    g.resize(6);g<<1,1,1,0,0,0;
    long double mass=0;
    for(int code=0;code<64;++code) {
        long double s=0;for(int i=0;i<6;++i)if(code&(1<<i))s+=w[i];
        if(s<=-.3L+1e-14L||s>=10.2L-1e-14L)mass+=1.L/64;
    }
    const auto endpoint=eval(r,e,g,h,false,.01,true);
    require(endpoint.values[11]==0,"endpoint atoms are resolved without a divergent root");
    close(endpoint.values[12],-std::log10(mass),1e-12,"inclusive endpoint matches exhaustive distribution");
}
void waldRank() {
    const int n=100;
    Eigen::VectorXd h(n),g(n),r(n),e(n),y(n);
    for(int i=0;i<n;++i) {
        h[i]=i%3;g[i]=h[i]>0&&i%5==0?1:0;e[i]=std::sin(i*.37);
        r[i]=std::cos(i*.23)+.3*g[i];
        y[i]=.2*h[i]+.4*e[i]+.3*g[i]*e[i]+std::sin(i*.31);
    }
    auto wd=std::make_shared<interaction_wald::WaldData>();
    wd->trait=interaction_wald::TraitType::Linear;wd->covar=e;wd->y=y;
    const auto got=eval(r,e,g,h,true,100,true,wd);
    auto independent=*wd;
    independent.covar.resize(n,2);independent.covar.col(0)=e;independent.covar.col(1)=h;
    const double want=interaction_wald::waldInteractionLog10P(independent,g,e);
    require(std::isfinite(got.values[6]),"local Wald rank-safe K-1 H adjustment");
    close(got.values[6],want,1e-12,"local Wald includes one independent H column");
    spagxe_local::LocalKernel kernel(r,{"E"},{e},wd,1,100);
    Eigen::MatrixXd repeated(n,4);repeated.col(0)=h;repeated.col(1)=2-h.array();repeated.col(2)=h;repeated.col(3).setConstant(2);
    const auto duplicate=kernel.evaluate(g,h,repeated);
    close(duplicate.values[6],want,1e-12,"repeated/constant H columns safely skipped");
}
void degeneracy() {
    const Eigen::VectorXd r=Eigen::VectorXd::LinSpaced(8,-1,1),e=r.array().square();
    const Eigen::VectorXd h=Eigen::VectorXd::Zero(8),g=h;
    const auto noH=eval(r,e,g,h,false);
    require(std::isnan(noH.values[5])&&noH.values[11]==8,"H0 no-test is explicit");
    Eigen::VectorXd bad=h;bad[1]=.5;
    bool thrown=false;
    try{eval(r,e,g,bad,false);}catch(const std::invalid_argument&){thrown=true;}
    require(thrown,"fractional H rejected");
    bad=g;bad[2]=1;thrown=false;
    try{eval(r,e,bad,h,false);}catch(const std::invalid_argument&){thrown=true;}
    require(thrown,"G above H rejected");
}
}
int main() {
    try{normalOracle();cgfEnumeration();fullSpa();finiteSupportTails();waldRank();degeneracy();
        std::cout<<checks<<" local kernel checks passed\n";return 0;}
    catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}
}
