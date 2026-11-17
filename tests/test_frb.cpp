#include "frb_sim.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace frbsim;
namespace {
void require(bool cond, const std::string& msg){ if(!cond) throw std::runtime_error(msg); }
std::vector<Cx> as_complex(const std::vector<double>& vals) {
    std::vector<Cx> out; out.reserve(vals.size());
    for (double v: vals) out.emplace_back(v,0.0);
    return out;
}
std::map<std::string,std::vector<Cx>> load_reference(const std::string& filename){
    std::ifstream in(filename); require(bool(in), "Cannot open golden TSV: " + filename);
    std::map<std::string,std::vector<Cx>> expected;
    std::string line;
    while (std::getline(in,line)) {
        if (line.empty() || line.front()=='#') continue;
        std::istringstream s(line);
        std::string name; std::size_t i; double re,im;
        require(bool(s>>name>>i>>re>>im), "Malformed golden TSV line");
        auto& vals=expected[name];
        require(i == vals.size(), "Unexpected fixture index for " + name);
        vals.emplace_back(re,im);
    }
    return expected;
}
void compare(const std::map<std::string,std::vector<Cx>>& expected,
             const std::string& key, const std::vector<Cx>& got,
             double atol=5e-12, double rtol=3e-13) {
    auto it=expected.find(key); require(it!=expected.end(), "Missing golden fixture: " + key);
    const auto& want=it->second;
    require(want.size()==got.size(), "Size mismatch for " + key);
    for(std::size_t i=0;i<got.size();++i){
        const double err=std::abs(want[i]-got[i]);
        if (!(err <= atol+rtol*std::abs(want[i]))) {
            std::ostringstream msg;
            msg << key << " at " << i << " expected=" << want[i]
                << " got=" << got[i] << " error=" << err;
            throw std::runtime_error(msg.str());
        }
    }
    std::cout << "PASS " << key << " (" << got.size() << " values)\n";
}
}
int main(int argc,char** argv){
    try {
        require(argc==2, "Usage: frb_sim_test golden.tsv");
        auto golden=load_reference(argv[1]);
        BoxCarBurst burst{0.2,0.35,3.2,15.0,165e6};
        Source src{burst,92.0,-25.0};
        std::vector<double> times{0.3,0.53,0.68,1.1};
        std::vector<double> freq{138e6,155e6,172e6};
        std::vector<UVW> uvw{{0.,0.,0.},{120.,-40.,10.},{-50.,80.,-5.}};
        std::vector<double> delays;
        for(auto f:freq) delays.push_back(burst.delay(f));
        compare(golden,"delay",as_complex(delays));
        std::vector<double> eval;
        for(auto t:times) for(auto f:freq) eval.push_back(burst.evaluate(t,f));
        compare(golden,"evaluate",as_complex(eval));
        compare(golden,"sample",as_complex(burst.sample(times,.2,freq,1.28e6,5,7)));
        LMN lmn=src.lmn(90.,-26.);
        compare(golden,"lmn",as_complex(std::vector<double>{lmn[0],lmn[1],lmn[2]}));
        std::vector<Cx> phase;
        for(auto row:uvw) for(auto f:freq) phase.push_back(phase_term(row,lmn,f));
        compare(golden,"phase",phase);
        compare(golden,"ds_1",src.dynamic_spectrum(times,.2,freq,1.28e6,Correlations::I,5,7));
        compare(golden,"ds_2",src.dynamic_spectrum(times,.2,freq,1.28e6,Correlations::XX_YY,5,7));
        compare(golden,"ds_4",src.dynamic_spectrum(times,.2,freq,1.28e6,Correlations::XX_XY_YX_YY,5,7));
        std::vector<double> t3(times.begin(), times.begin()+3);
        compare(golden,"tabrow",predict_tabrow(src,90.,-26.,t3,.2,
                           uvw,freq,1.28e6,Correlations::XX_XY_YX_YY,5,7));
        std::vector<UVW> fixed{{120.,-40.,10.},{-50.,80.,-5.}};
        compare(golden,"reg_fixed_legacy",predict_regular(src,90.,-26.,t3,.2,fixed,freq,1.28e6,
                           Correlations::XX_XY_YX_YY,5,7,true,true));
        std::vector<UVW> uvws{{120.,-40.,10.},{-50.,80.,-5.},
                              {130.,-45.,12.},{-55.,85.,-7.},
                              {125.,-43.,11.},{-53.,83.,-6.}};
        compare(golden,"reg_variable",predict_regular(src,90.,-26.,t3,.2,uvws,freq,1.28e6,
                           Correlations::XX_YY,5,7,false));
        // Physical invariants independent of Python implementation:
        require(std::abs(phase_term({0.,0.,0.},lmn,150e6)-Cx(1.,0.))<1e-14, "zero UVW phase");
        const Source central{burst,90.,-26.};
        require(std::abs(phase_term({100.,200.,30.},central.lmn(90.,-26.),150e6)-Cx(1.,0.))<1e-12,
                "phase centre must give unit complex visibility");
        const BoxCarBurst nodm{0.25,0.5,3.2,0.,165e6};
        require(nodm.evaluate(0.25,150e6)==nodm.peak,
                "pulse leading boundary inclusive");
        require(nodm.evaluate(0.75,150e6)==nodm.peak,
                "pulse trailing boundary inclusive");
        auto actual = predict_regular(src,90.,-26.,t3,.2,fixed,freq,1.28e6,
                                     Correlations::XX_XY_YX_YY,5,7,true,false);
        auto legacy = predict_regular(src,90.,-26.,t3,.2,fixed,freq,1.28e6,
                                     Correlations::XX_XY_YX_YY,5,7,true,true);
        require(actual != legacy, "test must demonstrate upstream fixed UVW row-order bug");
        // Independent row-order reference: expand the corrected baseline/time geometry.
        std::vector<double> expanded_times;
        std::vector<UVW> expanded_uvw;
        for (double time : t3) for (const UVW& baseline : fixed) {
            expanded_times.push_back(time);
            expanded_uvw.push_back(baseline);
        }
        auto independently_correct = predict_tabrow(src,90.,-26.,expanded_times,.2,
                        expanded_uvw,freq,1.28e6,Correlations::XX_XY_YX_YY,5,7);
        require(actual.size() == independently_correct.size(), "corrected row count");
        for(std::size_t i=0;i<actual.size();++i)
            require(std::abs(actual[i]-independently_correct[i])<5e-12,
                    "corrected time-major UVW must match tabrow prediction");
        std::vector<std::complex<float>> buffer(actual.size(),{1.f,0.f});
        add_inplace(buffer,actual);
        for(std::size_t i=0;i<buffer.size();++i)
            require(std::abs(buffer[i]-std::complex<float>(1.f,0.f)-static_cast<std::complex<float>>(actual[i]))<2e-6f,
                    "in-place injection mismatch");
        std::cout << "PASS physical invariants and in-place injection\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
