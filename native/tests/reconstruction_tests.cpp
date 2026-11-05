#include <dynsyslab/core/Reconstruction.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    try {
        using C=std::complex<double>;
        const std::vector<C> psi{{1,0},{0,1},{2,0},{0,2},{-1,0},{0,-1}};
        const std::vector<C> expectedC{{2,0},{0,-0.5},{-1,0},{0,1.5}};
        std::vector<C> values(6);
        for (std::size_t sample=0;sample<3;++sample)
            for (std::size_t output=0;output<2;++output)
                for (std::size_t mode=0;mode<2;++mode)
                    values[sample*2+output]+=psi[sample*2+mode]*expectedC[output*2+mode];
        const auto fit=dynsyslab::core::fitObservableReconstruction(psi,3,2,values,2);
        if (fit.matrix.size()!=4 || fit.relativeResidual>1e-12)
            throw std::runtime_error("projection residual");
        const auto reconstructed=dynsyslab::core::reconstructObservable(fit,psi,3);
        for (std::size_t i=0;i<values.size();++i)
            if (std::abs(reconstructed[i]-values[i])>1e-12)
                throw std::runtime_error("observable reconstruction");
        std::cout << "dynsyslab_reconstruction_tests: PASS\n"; return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_reconstruction_tests: FAIL: " << error.what() << '\n'; return 1;
    }
}
