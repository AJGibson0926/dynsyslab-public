#pragma once
// Real elliptic integrals in the parameter convention m = k^2 (DLMF 19.2), through Carlson's
// symmetric integrals R_F and R_D computed by duplication (B. C. Carlson, "Numerical computation of
// real or complex elliptic integrals", Numer. Algorithms 10 (1995) 13-26; DLMF 19.36.1).
// The incomplete integrals are continued to every real amplitude by F(phi + n pi|m) =
// F(phi|m) + 2 n K(m) and E(phi + n pi|m) = E(phi|m) + 2 n E(m) (DLMF 19.2.10), so a phase built from
// them is continuous across turns. Valid for m < 1 (and m sin^2 phi < 1).

namespace dynsyslab::core {

[[nodiscard]] double carlsonRF(double x, double y, double z);
[[nodiscard]] double carlsonRD(double x, double y, double z);
[[nodiscard]] double ellipticK(double m);                     // K(m)
[[nodiscard]] double ellipticE(double m);                     // E(m)
[[nodiscard]] double ellipticF(double phi, double m);         // F(phi|m)
[[nodiscard]] double ellipticEIncomplete(double phi, double m);  // E(phi|m)

}  // namespace dynsyslab::core
