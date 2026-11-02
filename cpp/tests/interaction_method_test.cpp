// SPAGxE v1 statistical regression tests.
// Copyright (c) 2026 Yuzhuo Ma, yuzhuoma@amss.ac.cn.
// Independent scalar oracles avoid production SPA / p-value helper calls.
#include "interaction/interaction_test.hpp"
#include "util/eigen_exp.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
int checks = 0;
void require(bool ok, const std::string &name) {
    ++checks;
    if (!ok) throw std::runtime_error(name);
}
void close(double got, double want, double tol, const std::string &name) {
    require(std::isfinite(got) && std::isfinite(want) &&
            std::abs(got-want) <= tol * (1.0+std::abs(want)), name);
}
double normalMagnitude(double z) {
    return -std::log10(std::erfc(std::abs(z)/std::sqrt(2.0)));
}
std::vector<double> evaluate(const Eigen::VectorXd &r, const Eigen::VectorXd &e,
                            Eigen::VectorXd g, double q, double epsilon,
                            double cutoff=100.0,
                            std::shared_ptr<const AFData> af=nullptr) {
    InteractionScan method(r, {"E"}, {e}, nullptr, nullptr, epsilon, cutoff, af);
    require(method.resultSize()==14, "nine-column environment schema");
    std::vector<double> out;
    method.getResultVec(g,q,0,out);
    require(out.size()==14, "result shape");
    return out;
}
Eigen::VectorXd originalRProjection(const Eigen::VectorXd &r,
                                    const Eigen::VectorXd &g) {
    // Explicit two-by-two inverse is independent of the centered-G code.
    const double n=g.size(), sg=g.sum(), sgg=g.squaredNorm();
    const double sr=r.sum(), sgr=g.dot(r), det=n*sgg-sg*sg;
    const double a=(sgg*sr-sg*sgr)/det;
    const double b=(n*sgr-sg*sr)/det;
    return (r.array()-a-b*g.array()).matrix();
}
void checkNormal(const Eigen::VectorXd &r, const Eigen::VectorXd &e,
                 const Eigen::VectorXd &g, double q, bool branchB) {
    const auto out=evaluate(r,e,g,q,branchB?1.0:0.0);
    const double f=2.0*q*(1.0-q);
    const Eigen::VectorXd centeredG=g.array()-2.0*q;
    const double zg=centeredG.dot(r)/std::sqrt(f*r.squaredNorm());
    close(out[1],zg,1e-12,"centered marginal score");
    close(out[0],normalMagnitude(zg),1e-11,"marginal p");
    Eigen::VectorXd w;
    if (branchB) {
        const Eigen::VectorXd r0=originalRProjection(r,g);
        w=e.array()*r0.array();
    } else {
        const double lambda=(e.array()*r.array().square()).sum()/r.squaredNorm();
        w=(e.array()-lambda)*r.array();
    }
    const double variance=f*w.squaredNorm();
    const double score=centeredG.dot(w);
    const double z=score/std::sqrt(variance);
    close(out[8],z,2e-11,"normal GxE score matches original R weight");
    close(out[12],normalMagnitude(z),2e-11,"separate score-test p");
    close(out[5],out[12],1e-13,"SPA-only final p");
    close(out[9],score/variance,2e-11,"centered score beta");
    close(out[9]/out[10],out[8],1e-13,"beta/SE equals score z");
    close(out[13],branchB?2.0:1.0,0.0,"branch code");
    close(out[11],1.0,0.0,"designed normal status");
}
void safeExpCompatibilityTest() {
    Eigen::ArrayXd values(13);
    values << -500.0,-100.0,-10.0,-2.0,-1.0,-0.1,0.0,0.1,1.0,2.0,10.0,100.0,500.0;
    const Eigen::ArrayXd result=spagxe_numeric::eigenSafeExp(-values);
    const Eigen::ArrayXd expected=values;
    values.setZero(); // returned result must own its storage
    for(Eigen::Index i=0;i<result.size();++i)
        close(result[i],std::exp(-expected[i]),2e-13,"safe exp, expression input and owning result");
}
void normalAndBranchTests() {
    Eigen::VectorXd g(12),r(12),e(12);
    g << 0,1,0,2,1,0,1,2,0,1,2,1;
    r << 1.4,-0.1,0.9,0.2,1.0,0.4,-0.4,1.8,0.8,-0.6,1.1,0.3;
    e << 0,1,0,1,0,1,0,1,1,0,1,0;
    const double q=g.mean()/2;
    require(std::abs(r.sum())>1.0,"nonzero-mean residual fixture");
    checkNormal(r,e,g,q,false);
    checkNormal(r,e,g,q,true);
    checkNormal(r,(e.array()+2.0).matrix(),g,q,true);
    // The restored paper/R weight differs from the upstream extra projection.
    const Eigen::VectorXd r0=originalRProjection(r,g);
    const Eigen::VectorXd original=e.array()*r0.array();
    const double a=original.dot(r0)/r0.squaredNorm();
    const Eigen::VectorXd upstream=original-a*r0;
    require(std::abs(original.squaredNorm()-upstream.squaredNorm())>1e-4,
            "fixture distinguishes original R and upstream Branch B");
    // Batched GEMM and marker-wise code must agree for two environments.
    InteractionScan method(r,{"E","E2"},{e,(e.array()+2.0).matrix()},
                       nullptr,nullptr,1.0,100.0);
    Eigen::MatrixXd batch(g.size(),2);
    batch.col(0)=g;
    batch.col(1)=2.0-g.array();
    std::vector<std::vector<double>> batched;
    method.getResultBatch(batch,{q,1.0-q},{0,1},batched);
    for (int k=0;k<2;++k) {
        Eigen::VectorXd one=batch.col(k);
        std::vector<double> single;
        method.getResultVec(one,k==0?q:1.0-q,k,single);
        for (size_t j=0;j<single.size();++j)
            if (std::isfinite(single[j]))
                close(batched[k][j],single[j],1e-11,"batch/single consistency");
    }
}
void exactSpaAndEnumerationTest() {
    const Eigen::VectorXd r=Eigen::VectorXd::Ones(8);
    Eigen::VectorXd e(8),g(8);
    e << -1,-1,-1,-1,1,1,1,1;
    g << 0,0,0,0,2,2,2,1;
    const auto out=evaluate(r,e,g,0.5,0.0,0.0);
    // Here K(t)=16 log cosh(t/2); solve K'(t)=7 analytically.
    const double t=2.0*std::atanh(7.0/8.0);
    const double k=16.0*std::log(std::cosh(t/2.0));
    const double k2=4.0*(1.0-49.0/64.0);
    const double w=std::sqrt(2.0*(7.0*t-k));
    const double v=t*std::sqrt(k2);
    const double rstar=w+std::log(v/w)/w;
    const double analyticP=std::erfc(rstar/std::sqrt(2.0));
    close(out[12],-std::log10(analyticP),2e-10,"full-CGF analytic saddlepoint");
    require(out[11]<=2.0,"full-CGF SPA converged");
    // Enumerate all 3^8 genotype configurations under independent Bin(2,.5).
    double strict=0.0,inclusive=0.0;
    for(int code=0;code<6561;++code) {
        int rem=code; double p=1.0,score=0.0;
        for(int i=0;i<8;++i) {
            const int gi=rem%3; rem/=3;
            p*=gi==1?0.5:0.25;
            score+=gi*e[i];
        }
        if(std::abs(score)>7.0)strict+=p;
        if(std::abs(score)>=7.0)inclusive+=p;
    }
    close(strict,2.0/65536.0,1e-14,"strict binomial enumeration");
    close(inclusive,34.0/65536.0,1e-14,"inclusive binomial enumeration");
    require(analyticP>strict && analyticP<inclusive,
            "SPA lies between strict and inclusive discrete tail probabilities");
    // An approximation is not an exact discrete p-value: this assertion only
    // validates the CGF/solver and the expected continuity bracket.
}
std::shared_ptr<AFData> makeAF(int n) {
    auto af=std::make_shared<AFData>();
    af->nPC=1;
    af->onePlusPCs.resize(n,2);
    af->onePlusPCs.col(0).setOnes();
    for(int i=0;i<n;++i)af->onePlusPCs(i,1)=(i-(n-1)/2.0)/n;
    const Eigen::MatrixXd inv=(af->onePlusPCs.transpose()*af->onePlusPCs).inverse();
    af->XtX_inv_Xt=inv*af->onePlusPCs.transpose();
    af->sqrt_XtX_inv_diag=inv.diagonal().cwiseSqrt();
    return af;
}
void highAltAndRankTests() {
    constexpr int n=100;
    Eigen::VectorXd r(n),e(n),g=Eigen::VectorXd::Zero(n);
    for(int i=0;i<n;++i) {
        r[i]=0.3+std::sin(0.37*i)+0.2*std::cos(0.11*i);
        e[i]=(i%3)==0?1.0:0.0;
    }
    for(int i=0;i<5;++i)g[i]=1.0;
    const double q=g.mean()/2.0;
    const auto af=makeAF(n);
    for(bool branchB:{false,true}) {
        const double eps=branchB?1.0:0.0;
        const auto minor=evaluate(r,e,g,q,eps,100.0,af);
        const auto major=evaluate(r,e,(2.0-g.array()).matrix(),1.0-q,eps,100.0,af);
        close(minor[0],major[0],2e-11,"high-ALT marginal p symmetry");
        close(minor[12],major[12],2e-11,"high-ALT GxE p symmetry");
        close(minor[8],-major[8],2e-11,"ALT-oriented score sign");
        close(minor[9],-major[9],2e-11,"ALT-oriented beta sign");
        close(minor[10],major[10],2e-11,"ALT-oriented SE");
        require(std::isfinite(major[12]),"high-ALT marker remains testable");
    }
    auto singular=makeAF(n);
    singular->onePlusPCs.col(1).setOnes();
    bool rejected=false;
    try { (void)evaluate(r,e,g,q,0.0,100.0,singular); }
    catch(const std::runtime_error &) {rejected=true;}
    require(rejected,"rank-deficient PC design rejected");
}
}
int main() {
    try {
        safeExpCompatibilityTest();
        normalAndBranchTests();
        exactSpaAndEnumerationTest();
        highAltAndRankTests();
        std::cout<<"SPAGxE method: "<<checks<<" checks passed.\n";
        return 0;
    } catch(const std::exception &e) {
        std::cerr<<"SPAGxE method test failed: "<<e.what()<<"\n";
        return 1;
    }
}

