#pragma once
// Real elliptic integrals in the parameter convention m = k^2 (DLMF 19.2), through Carlson's
// symmetric integrals R_F, R_D, R_C and R_J computed by duplication (B. C. Carlson, "Numerical computation of
// real or complex elliptic integrals", Numer. Algorithms 10 (1995) 13-26; DLMF 19.36.1).
// The incomplete integrals are continued to every real amplitude by F(phi + n pi|m) =
// F(phi|m) + 2 n K(m) and E(phi + n pi|m) = E(phi|m) + 2 n E(m) (DLMF 19.2.10), so a phase built from
// them is continuous across turns; the third kind likewise, Pi(phi + n pi) = Pi(phi) + 2 n Pi(n|m).
// Valid for m < 1 (and m sin^2 phi < 1); the third kind for n < 1 (no Cauchy principal value).
// R_C uses duplication too, not its arccos/arccosh closed form, which loses precision as its
// arguments meet (as they do inside R_J's duplication).

namespace dynsyslab::core {

[[nodiscard]] double carlsonRF(double x, double y, double z);
[[nodiscard]] double carlsonRD(double x, double y, double z);
[[nodiscard]] double ellipticK(double m);                     // K(m)
[[nodiscard]] double ellipticE(double m);                     // E(m)
[[nodiscard]] double ellipticF(double phi, double m);         // F(phi|m)
[[nodiscard]] double ellipticEIncomplete(double phi, double m);  // E(phi|m)
[[nodiscard]] double carlsonRC(double x, double y);
[[nodiscard]] double carlsonRJ(double x, double y, double z, double p);
// Pi(n; phi|m) = int_0^phi dtheta / ((1 - n sin^2 theta) sqrt(1 - m sin^2 theta)) (DLMF 19.2.7, 19.25.14).
[[nodiscard]] double ellipticPi(double phi, double n, double m);

}  // namespace dynsyslab::core
