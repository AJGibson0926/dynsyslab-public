#pragma once
// Robust state feedback u = K z for a bilinear K-ROM with a residual sector:
//
//   z_plus = (A + sum_i u_i B_i) z + B0 u + r,   ||r|| <= c_x ||z|| + c_u ||u||.
//
// Variables X > 0, L, mu > 0, with K = L X^-1. For every vertex theta of the input
// box, with Phi = A X + sum_i theta_i B_i X + B0 L:
//
//   [ lambda X   Phi'      g_x X   g_u L' ]
//   [ Phi        X - mu I  0       0      ]  >= 0,    [ ubar_i^2  L_i ]
//   [ g_x X      0         mu I    0      ]           [ L_i'      X   ]  >= 0,
//   [ g_u L      0         0       mu I   ]           X <= R^2 I.
//
// Consequence, for the MODEL WITH ANY r IN THE SECTOR: on E = {z' X^-1 z <= 1},
// |K_i z| <= ubar_i, V(z) = z' X^-1 z obeys V(z_plus) <= lambda V(z), and E is
// forward invariant. E lies in the ball ||z|| <= R, where the sector is assumed.
// Nothing is claimed for the physical plant unless the sector holds for it there.
//
// Methods: quadratic stabilization via the change of variables L = K X (Bernussou,
// Peres & Geromel, Syst. Control Lett. 13 (1989) 65-72; Boyd, El Ghaoui, Feron &
// Balakrishnan, Linear Matrix Inequalities in System and Control Theory, SIAM
// 1994), S-procedure and invariant ellipsoids (same book, Ch. 2, 5 and 7), and the
// bilinear SafEDMD setting (Strasser, Worthmann et al., arXiv:2402.03145).
#include <dynsyslab/core/BilinearKrom.hpp>
#include <dynsyslab/core/DenseSdp.hpp>
#include <dynsyslab/core/ResidualSector.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct CertifiedDesignSpec {
    std::vector<double> inputBounds; // ubar_i > 0
    double stateRadius{1.0};         // R: region where the model and sector are trusted
    double contraction{0.99};        // lambda in (0, 1): V(z_plus) <= lambda V(z)
    double solverMargin{1e-4};       // LMIs are solved as F >= solverMargin I (normalized units)
    double certificateMargin{1e-6};  // the verifier requires lambda_min(F) >= this
    SdpOptions solver{};
};

enum class CertifiedDesignStatus { certified, infeasible, inconclusive, verificationFailed };
[[nodiscard]] const char* certifiedDesignStatusName(CertifiedDesignStatus status);

struct CertifiedDesignResult {
    CertifiedDesignStatus status{CertifiedDesignStatus::inconclusive};
    bool certified{false};
    // Populated only when certified; never a best-effort controller.
    std::vector<double> gain;   // m x n, u = K z (note the sign differs from LQR's u = -K x)
    std::vector<double> X;      // n x n
    double mu{};
    double inscribedRadius{};   // largest ball around 0 inside E
    LmiCertificate certificate;
    SdpStatus solverStatus{SdpStatus::iterationLimit};
    double primalResidual{}, dualResidual{}, gap{};
    std::size_t iterations{};
    std::string message;
};

[[nodiscard]] CertifiedDesignResult designCertifiedBilinearController(const BilinearKrom& model,
    const ResidualSectorBound& sector, const CertifiedDesignSpec& spec);

// Rebuilds every LMI from the deployed (K, X, mu), independent of the solver.
[[nodiscard]] LmiCertificate verifyCertifiedDesign(const BilinearKrom& model,
    const ResidualSectorBound& sector, const CertifiedDesignSpec& spec,
    const std::vector<double>& gain, const std::vector<double>& X, double mu);

}  // namespace dynsyslab::core
