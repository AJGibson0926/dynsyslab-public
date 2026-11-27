#pragma once
// Bridge from a KoopmanModel with a stored continuous bilinear control map
//
//   w' = G w + sum_i (b0_i + N_i w) u_i,    w = split(psi)
//
// to the discrete, target-centred BilinearKrom used by the certified design:
//
//   z = w - w*,  v = u - u*,
//   z_plus = A_eff z + B0_eff v + sum_i v_i B_i z + c,
//   A_eff  = A + sum_i u*_i B_i,   B0_eff = B0 + [B_1 w* ... B_m w*],
//   c      = A_eff w* + B0 u* - w*.
//
// A = exp(G dt) is the exact transition of the frozen generator. B0 is the exact
// zero-order-hold integral of b0. B_i = int_0^dt exp(G (dt-s)) N_i exp(G s) ds is the
// exact first-order-in-u term of exp((G + u N) dt), computed by Van Loan's block
// exponential (C. F. Van Loan, "Computing integrals involving the matrix exponential",
// IEEE Trans. Automat. Control 23 (1978) 395-404). The O(u^2) remainder is not
// representable in this model family; conversionDefect reports it and the data-based
// residual sector must cover it.
//
// A nonzero constant c means the target is not an equilibrium of the model at u*. It is
// reported and never dropped silently: a sector cannot hold a nonzero residual at z = 0.
//
// Bilinear Koopman form: Goswami & Paley, IEEE TAC 67 (2022) 2715-2728; Bruder, Fu &
// Vasudevan, IEEE RA-L 6 (2021) 4369-4376. Target-centred lifting: Strasser, Worthmann
// et al., SafEDMD, arXiv:2402.03145.
#include <dynsyslab/core/BilinearKrom.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

// Exact real-split form of a stored bilinear map (conjugate operators included), as the one real
// bilinear model: A = generator, B0 = offsets b0, B[i] = N_i; sampleTime 0 (continuous).
[[nodiscard]] BilinearKrom continuousBilinearSplit(const KoopmanModel& model,
    const std::string& representation);

struct KoopmanBilinearBridgeOptions {
    std::string representation{"bilinear"};
    double sampleTime{};                // dt > 0
    std::vector<double> targetLifted;   // w*, empty = 0
    std::vector<double> targetInput;    // u*, empty = 0
    double equilibriumTolerance{1e-9};  // ||c|| <= tol * max(1, ||w*||)
};

struct KoopmanBilinearBridge {
    BilinearKrom model;                  // centred (z, v); model.A == A_eff
    BilinearKrom continuous;
    std::vector<double> frozenTransition;              // A = exp(G dt), n x n
    std::vector<double> discreteOffsets;               // uncentred B0, n x m
    std::vector<std::vector<double>> discreteOperators;// uncentred B_i
    std::vector<double> targetLifted, targetInput;
    std::vector<double> constantTerm;    // c
    double constantNorm{}, constantTolerance{};
    bool targetIsEquilibrium{false};
    bool heldTargetInput{false};         // u* != 0, so model.A differs from A
};
[[nodiscard]] KoopmanBilinearBridge bridgeKoopmanBilinear(const KoopmanModel& model,
    const KoopmanBilinearBridgeOptions& options);

// Regresses B0 and B_i on centred held-input pairs with A frozen at bridge.model.A.
[[nodiscard]] BilinearFit refitBridgeOperators(const KoopmanBilinearBridge& bridge,
    const BilinearDataset& centred, const BilinearFitOptions& options = {});

struct BilinearConversionDefect {
    double transition{};  // max ||exp((G + theta N) dt) - (A + theta B)||_F over the box vertices
    double input{};       // max ||Gamma(theta) theta - B0 theta||_2 over the box vertices
};
// Vertices are u* +/- inputBounds. Zero when N = 0 for the transition part.
[[nodiscard]] BilinearConversionDefect conversionDefect(const KoopmanBilinearBridge& bridge,
    const std::vector<double>& inputBounds);

}  // namespace dynsyslab::core
