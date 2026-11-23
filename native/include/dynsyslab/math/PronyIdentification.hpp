#pragma once
#include <dynsyslab/math/DiscreteResolvent.hpp>
#include <optional>
#include <string>
#include <vector>
namespace dynsyslab::math {
// Poles and residues of y_k=sum_j c_j mu_j^k from uniform samples. Exponential fitting goes back
// to Prony, J. Ecole Polytechnique 1(22) (1795) 24-76. This uses the SVD matrix pencil of Hua &
// Sarkar, IEEE Trans. Acoust. Speech Signal Process. 38(5) (1990) 814-824, with the truncation
// of Sarkar & Pereira, IEEE Antennas Propag. Mag. 37(1) (1995) 48-55. Simple poles only:
// repeated poles appear as a close pair and raise the cluster warning.
struct PronyOptions {
    double sampleInterval{1};
    std::optional<std::size_t> order;  // Empty: count singular values above rankTolerance*s_1.
    std::size_t pencilParameter{0};    // L; 0 selects min(N/2,128).
    double rankTolerance{1e-9};
    double clusterFactor{1};           // Cluster: separation below factor*(2 pi/N) per sample.
    double gapThreshold{10};           // s_M/s_{M+1} below this means no clear model order.
};
struct PronyPole {std::complex<double> multiplier,residue,eigenvalue;};  // eigenvalue is principal.
struct PronyOrderScan {std::size_t order{};double relativeResidual{};};
struct PronyResult {
    std::vector<PronyPole> poles;           // Sorted by descending |residue|.
    DiscretePoleModel model;
    std::vector<double> singularValues;
    std::size_t order{},pencilParameter{};
    std::string orderSource;                // "declared" or "rank".
    double relativeResidual{},singularGap{},minimumSeparation{},rayleighLimit{},vandermondeCondition{};
    std::vector<PronyOrderScan> orderScan;
    bool clustered{false},nearContinuous{false};
    std::vector<std::string> warnings;
};
PronyResult pronyIdentify(const ComplexVector& samples,const PronyOptions& ={});
}
