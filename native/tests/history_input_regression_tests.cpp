#include <dynsyslab/core/HistoryInputRegression.hpp>
#include <dynsyslab/core/HistoryKoopman.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
int main() {
    using namespace dynsyslab::core;
    try {
        const std::complex<double> m=std::exp(std::complex<double>(-0.03,0.1));
        std::vector<HistoryInputSample> rows;
        for(int k=0;k<900;++k) {
            const double x=std::sin(k*0.123),y=std::cos(k*0.071),u=std::sin(k*0.273);
            const std::complex<double> z(x,y),b=std::complex<double>(0.1,0.2)
                +x*std::complex<double>(0.03,-0.04)+y*std::complex<double>(-0.07,0.02);
            rows.push_back({{x,y},z,m*z+b*u,u});
        }
        const auto fitted=fitHistoryInputMap(rows,m,0.01,1,{1,1},0);
        require(fitted.rank==3,"regression rank changed");
        require(scoreHistoryInputMap(fitted,rows).relativeResidual<1e-12,"exact affine recovery failed");
        const auto constant=fitHistoryInputMap(rows,m,0.01,0,{1,1},0);
        require(scoreHistoryInputMap(constant,rows).relativeResidual>0.1,"constant ablation unexpectedly exact");
        const auto standard=std::get<BilinearControlMap>(fitted.canonicalRepresentation());
        const auto factor=(m-1.0)/(std::log(m)/0.01);
        const std::complex<double> z(0.4,-0.7);
        const auto reconstructed=factor*(standard.offsets[0]+standard.operators[0]*z
            +standard.conjugateOperators[0]*std::conj(z));
        require(std::abs(reconstructed-fitted.evaluate(z))<1e-13,"canonical map conversion failed");
        HankelLeftBasis basis; basis.rows=3; basis.rank=1;
        basis.leftSingularVectors={0.2,0.3,0.4}; basis.singularValues={2};
        CausalHistoryProjector online(basis);
        online.append(1); online.append(2); require(!online.projectLatest(),"premature projection");
        online.append(3); online.append(4);
        require(std::abs(online.projectLatest()->at(0)-1.45)<1e-14,"decimated online projection mismatch");
        bool rejected=false;
        try { (void)fitHistoryInputMap({},m,0.01,1,{1,1}); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"empty fit was accepted");
        std::cout<<"history input regression and canonical map tests passed\n";
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
