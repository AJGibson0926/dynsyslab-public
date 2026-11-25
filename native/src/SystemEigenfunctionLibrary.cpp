// The system eigenfunction library: hand-derived Koopman eigenfunctions of the catalog systems, as
// text in each system's state and parameter names (DynamicalKnowledge::eigenfunctions, loaded by
// libraryEigenfunctions). Every entry is checked against the configured equations by the
// generator identity grad(psi).f = lambda psi in native/tests/system_eigenfunction_library_tests.cpp.
//
// Kinds: "exact" (one eigenvalue on the stated region), "principal" (linearization at an equilibrium,
// accurate near it), "shell" (exact, with the eigenvalue i n Omega(H) of the orbit's energy shell, so
// psi(x(t)) = exp(i n Omega(H) t) psi(x(0)) along every orbit of the region; Mezic 2020, section 11).
//
// Integrable phases: for a conserved H and a phase phi with dphi/dt = Omega(H), exp(i phi) is such a
// shell eigenfunction, and H (any function of it) has eigenvalue 0. The phases below invert the
// closed-form solutions (Jacobi elliptic functions, DLMF 22.13 and 22.19; Brizard (arXiv:1108.4970) for the
// pendulum; Jaffe and Brumer 1984 for the pendulum eigenfunctions; DeMarcus 1978, Kuru and Negro 2008,
// Delisle-Doray et al. 2019 and Henrici and Kappeler 2008 for the Morse, isotonic, Poschl-Teller and
// cosh oscillators; the Kepler equation for the mean anomaly). Elliptic integrals use the parameter
// m = k^2.
#include "DynamicalKnowledgeInternal.hpp"

#include <cmath>
#include <string>

namespace dynsyslab::core::detail {
namespace {
using S = std::string;
const S mezic = "Mezic (2020), J. Nonlinear Sci. 30:2091-2145, section 11 (arXiv:1702.07597)";
const S brizard = "Brizard, Jacobi zeta function and action-angle coordinates for the pendulum, arXiv:1108.4970";
const S jaffe = "Jaffe and Brumer (1984), J. Phys. Chem. 88:4829-4839, section II.F";
const S dlmf22 = "https://dlmf.nist.gov/22.19";
const S kuru = "Kuru and Negro (2008), classical factorization, arXiv:0709.4649";
const S demarcus = "DeMarcus (1978), Am. J. Phys., https://doi.org/10.1119/1.11110";
const S delisle = "Delisle-Doray et al. (2019), arXiv:1812.11582";
const S henrici = "Henrici and Kappeler (2008), arXiv:0802.4252";
const S murray = "Murray and Dermott (1999), Solar System Dynamics, chapter 2";
const S mezicLinear = "Mezic (2005), Nonlinear Dyn. 41:309-325";
const S arnold = "https://doi.org/10.1007/978-1-4757-2063-1";

double parameter(const SystemConfiguration& c, const char* key, double fallback = 0) {
    const auto found = c.parameters.find(key);
    return found == c.parameters.end() ? fallback : found->second;
}
S paren(const S& s) { return "(" + s + ")"; }

struct Library {
    std::vector<EigenfunctionEntry>& out;
    void add(S id, S title, S kind, S real, S imag, S rate, S rateImag, S conditions, S derivation,
             std::vector<S> refs, bool applicable) {
        out.push_back({std::move(id), std::move(title), std::move(kind), std::move(real), std::move(imag), std::move(rate),
                       std::move(rateImag), std::move(conditions), std::move(derivation), std::move(refs), applicable});
    }
    void invariant(S id, S title, S expression, S conditions, std::vector<S> refs, bool applicable) {
        add(std::move(id), std::move(title), "exact", std::move(expression), "0", "0", "0", std::move(conditions),
            "conserved by the configured field: eigenvalue 0", std::move(refs), applicable);
    }
    // psi = (lambda + b) x + v for x' = v, v' = -a x - b v, lambda = -b/2 + i sqrt(a - b^2/4) (and its
    // conjugate): the left eigenvectors of the linear field.
    void linearModes(const S& prefix, const S& x, const S& v, const S& a, const S& b, const S& kind, const S& conditions,
                     const S& derivation, bool applicable) {
        const S w = "sqrt(" + a + " - " + b + "^2/4)";
        add(prefix, "Oscillatory mode (lambda = -b/2 + i w)", kind, b + "/2*" + x + " + " + v, w + "*" + x, "-" + b + "/2", w,
            conditions, derivation, {mezicLinear}, applicable);
        add(prefix + "_conjugate", "Conjugate oscillatory mode", kind, b + "/2*" + x + " + " + v, "-" + w + "*" + x,
            "-" + b + "/2", "-" + w, conditions, derivation, {mezicLinear}, applicable);
    }
};

void pendulum(Library& L, const SystemConfiguration& c, bool free) {
    const bool conservative = free && parameter(c, "damping_n_m_s") == 0;
    const S a = "(gravity_m_s2/length_m)", H = "(velocity^2/2 + " + a + "*(1 - cos(angle)))";
    L.invariant("energy", "Energy H = v^2/2 + a(1 - cos theta)", H, "Zero damping and driver", {brizard}, conservative);
    // Libration 0 < H < 2a: alpha from sin(theta/2) = k sin(alpha), v = 2k sqrt(a) cos(alpha),
    // phi = pi F(alpha|k^2) / (2K), Omega = pi sqrt(a) / (2K) (Brizard, arXiv:1108.4970).
    const S m = paren(H + "/(2*" + a + ")"), k = "sqrt" + m;
    const S alpha = "atan2(sin(angle/2)/" + k + ", velocity/(2*" + k + "*sqrt" + a + "))";
    const S phi = paren("pi*ellipf(" + alpha + ", " + m + ")/(2*ellipk(" + m + "))");
    L.add("libration_phase", "Libration phase exp(i phi)", "shell", "cos" + phi, "sin" + phi, "0",
          "pi*sqrt" + a + "/(2*ellipk(" + m + "))", "Zero damping and driver; libration 0 < H < 2a, |theta| < pi",
          "theta = 2 asin(k sn(u|k^2)), v = 2 k sqrt(a) cn(u|k^2), phi = pi u/(2K): dphi/dt = Omega(H)",
          {brizard, jaffe, mezic, dlmf22}, conservative);
    // Rotation H > 2a, m = 2a/H: phi_r = pi F(theta/2|m)/K, dphi_r/dt = sign(v) pi sqrt(H/2)/K.
    const S mr = paren("2*" + a + "/" + H), phir = paren("pi*ellipf(angle/2, " + mr + ")/ellipk(" + mr + ")");
    L.add("rotation_phase", "Rotation phase exp(i phi_r)", "shell", "cos" + phir, "sin" + phir, "0",
          "velocity/sqrt(velocity^2)*pi*sqrt(" + H + "/2)/ellipk(" + mr + ")", "Zero damping and driver; rotation H > 2a",
          "|v| = sqrt(2H) sqrt(1 - m sin^2(theta/2)), so dF(theta/2|m)/dt = sign(v) sqrt(H/2)",
          {brizard, dlmf22, mezic}, conservative);
    L.linearModes("principal", "angle", "velocity", a, "(damping_n_m_s/(mass_kg*length_m^2))", "principal",
                  "Zero driver; underdamped; near the hanging equilibrium",
                  "left eigenvectors of Df(0) = [[0, 1], [-a, -b]], b = damping/(m l^2)",
                  free && parameter(c, "gravity_m_s2") / parameter(c, "length_m", 1) >
                              0.25 * std::pow(parameter(c, "damping_n_m_s") / (parameter(c, "mass_kg", 1) *
                                                      std::pow(parameter(c, "length_m", 1), 2)), 2));
}

void duffing(Library& L, const SystemConfiguration& c, bool free) {
    const bool ok = free && parameter(c, "delta") == 0 && parameter(c, "alpha") > 0 && parameter(c, "beta") > 0;
    const S H = "(velocity^2/2 - alpha*position^2/2 + beta*position^4/4)", D = "sqrt(alpha^2 + 4*beta*" + H + ")";
    L.invariant("energy", "Energy H = v^2/2 - alpha x^2/2 + beta x^4/4", H, "Zero damping and driver", {dlmf22}, ok);
    // Outer orbits H > 0: x = A cn(nu t|m), v = -A nu sn dn, A^2 = (alpha + D)/beta, nu^2 = D,
    // m = (alpha + D)/(2D); u = F(atan2(sn, cn)|m), phi = pi u/(2K).
    const S A = "sqrt((alpha + " + D + ")/beta)", m = "((alpha + " + D + ")/(2*" + D + "))";
    const S d = "sqrt(1 - " + m + " + " + m + "*(position/" + A + ")^2)";
    const S am = "atan2(-velocity/(" + A + "*sqrt(" + D + ")*" + d + "), position/" + A + ")";
    const S phi = paren("pi*ellipf(" + am + ", " + m + ")/(2*ellipk(" + m + "))");
    L.add("outer_phase", "Outer-orbit phase exp(i phi)", "shell", "cos" + phi, "sin" + phi, "0",
          "pi*sqrt(" + D + ")/(2*ellipk(" + m + "))", "Zero damping and driver; outer orbits H > 0",
          "x = A cn(sqrt(D) t|m) solves x'' = alpha x - beta x^3; phase inverted through atan2(sn, cn)",
          {dlmf22, mezic}, ok);
    // Inner orbits -alpha^2/(4 beta) < H < 0 (either well, s = sign(x)): x = s A dn(nu t|m),
    // v = -s A nu m sn cn, m = 2D/(alpha + D), nu = A sqrt(beta/2); am = atan2(sn^2, sn cn) in [0, pi],
    // phi = pi F(am|m)/K.
    const S mi = "(2*" + D + "/(alpha + " + D + "))", nu = paren(A + "*sqrt(beta/2)");
    const S ami = "atan2(1 - position^2/" + paren(A + "^2") + ", -position/sqrt(position^2)*velocity/(" + A + "*" + nu + "))";
    const S phii = paren("pi*ellipf(" + ami + ", " + mi + ")/ellipk(" + mi + ")");
    L.add("well_phase", "Single-well phase exp(i phi)", "shell", "cos" + phii, "sin" + phii, "0",
          "pi*" + nu + "/ellipk(" + mi + ")", "Zero damping and driver; inner orbits -alpha^2/(4 beta) < H < 0, either well",
          "x = s A dn(nu t|m) solves x'' = alpha x - beta x^3; sn^2 = (1 - x^2/A^2)/m, sn cn = -s v/(A nu m)",
          {dlmf22, mezic}, ok);
}

void integrable(Library& L, const SystemDescriptor& d, const SystemConfiguration& c, bool free) {
    const S& id = d.identifier;
    const bool ok = free && parameter(c, "damping") == 0;
    if (id == "isotonic_oscillator") {
        const S H = "(momentum^2/2 + omega^2*position^2/2 + inverse_square_strength/(2*position^2))";
        L.invariant("energy", "Energy", H, "Zero damping and driver", {kuru}, ok);
        L.add("isochronous_mode", "Isochronous mode Z = q^2 - H/omega^2 - i q p/omega", "exact",
              "position^2 - " + H + "/omega^2", "-position*momentum/omega", "0", "2*omega", "Zero damping and driver; position > 0",
              "Z' = 2 i omega Z for every orbit (the frequency does not depend on the energy)", {kuru}, ok);
    } else if (id == "morse_oscillator") {
        const S H = "(momentum^2/(2*mass) + depth*(1 - exp(-inverse_length*position))^2)", k = "sqrt(" + H + "/depth)";
        const S s2 = "(1 - " + H + "/depth)", W = "(inverse_length*sqrt(2*(depth - " + H + ")/mass))";
        L.invariant("energy", "Energy", H, "Zero damping and driver", {demarcus}, ok);
        L.add("phase", "Bound-orbit phase exp(i phi)", "shell", "(" + s2 + "*exp(inverse_length*position) - 1)/" + k,
              "-inverse_length*momentum*" + s2 + "*exp(inverse_length*position)/(mass*" + k + "*" + W + ")", "0", W,
              "Zero damping and driver; bound orbits 0 < H < depth",
              "elementary bound solution of the Morse oscillator; Omega = a sqrt(2(D - H)/m)", {demarcus, kuru, mezic}, ok);
    } else if (id == "poschl_teller") {
        const S ch = "((exp(inverse_length*position) + exp(-inverse_length*position))/2)";
        const S sh = "((exp(inverse_length*position) - exp(-inverse_length*position))/2)";
        const S H = "(momentum^2/2 - depth/" + ch + "^2)", W = "(inverse_length*sqrt(-2*" + H + "))";
        const S A = "sqrt((" + H + " + depth)/(-" + H + "))";
        L.invariant("energy", "Energy", H, "Zero damping and driver", {delisle}, ok);
        L.add("phase", "Bound-orbit phase exp(i phi)", "shell", "inverse_length*momentum*" + ch + "/(" + A + "*" + W + ")",
              sh + "/" + A, "0", W, "Zero damping and driver; bound orbits -depth < H < 0",
              "sinh(a position) = A sin(Omega t) with Omega = a sqrt(-2H)", {delisle, mezic}, ok);
    } else if (id == "cosh_oscillator") {
        const S H = "(momentum^2/2 + strength*(exp(position) + exp(-position)))", m = "((" + H + " - 2*strength)/(" + H + " + 2*strength))";
        const S A = "sqrt((" + H + " - 2*strength)/(4*strength))", nu = "sqrt((" + H + " + 2*strength)/2)";
        const S u = "((exp(position/2) - exp(-position/2))/2)", ch = "((exp(position/2) + exp(-position/2))/2)";
        const S dn = "sqrt(1 - " + m + " + " + m + "*(" + u + "/" + A + ")^2)";
        const S am = "atan2(-" + ch + "*momentum/(2*" + A + "*" + nu + "*" + dn + "), " + u + "/" + A + ")";
        const S phi = paren("pi*ellipf(" + am + ", " + m + ")/(2*ellipk(" + m + "))");
        L.invariant("energy", "Energy", H, "Zero damping and driver", {henrici}, ok);
        L.add("phase", "Orbit phase exp(i phi)", "shell", "cos" + phi, "sin" + phi, "0",
              "pi*" + nu + "/(2*ellipk(" + m + "))", "Zero damping and driver; H > 2 strength",
              "sinh(position/2) = A cn(nu t|m); phase inverted through atan2(sn, cn)", {henrici, dlmf22, mezic}, ok);
    } else if (id == "radial_twist") {
        L.invariant("action", "Action I = (position^2 + momentum^2)/2", "(position^2 + momentum^2)/2", "Zero damping and driver", {arnold}, ok);
        L.add("rotation", "Rotation q - i p", "shell", "position", "-momentum", "0", "omega + twist*(position^2 + momentum^2)/2",
              "Zero damping and driver", "(q - i p)' = i Omega(I) (q - i p), Omega = omega + twist I (q = position, p = momentum)", {arnold, mezic}, ok);
    } else if (id == "kepler_planar") {
        const S H = "((p_1^2 + p_2^2)/2 - mu/sqrt(q_1^2 + q_2^2))", a = "(-mu/(2*" + H + "))";
        const S rv = "((q_1*p_1 + q_2*p_2)/sqrt(mu*" + a + "))";
        const S M = paren("atan2(" + rv + ", 1 - sqrt(q_1^2 + q_2^2)/" + a + ") - " + rv);
        L.invariant("energy", "Energy", H, "Zero damping and driver", {murray}, ok);
        L.invariant("angular_momentum", "Angular momentum", "q_1*p_2 - q_2*p_1", "Zero damping and driver", {murray}, ok);
        L.add("mean_anomaly", "Mean-anomaly phase exp(i M)", "shell", "cos" + M, "sin" + M, "0",
              "sqrt(mu/" + a + "^3)", "Zero damping and driver; bound noncollision orbits H < 0",
              "e cos E = 1 - r/a, e sin E = (q.p)/sqrt(mu a), M = E - e sin E (Kepler's equation)", {murray, mezic}, ok);
    } else if (id == "euler_top") {
        L.invariant("casimir", "Casimir |M|^2", "M_1^2 + M_2^2 + M_3^2", "Zero damping and driver", {arnold}, ok);
        L.invariant("energy", "Energy", "M_1^2/(2*inertia_1) + M_2^2/(2*inertia_2) + M_3^2/(2*inertia_3)",
                    "Zero damping and driver", {arnold}, ok);
    } else if (id == "coupled_duffing_hamiltonian") {
        const S H1 = "p_1^2/2 + omega_1^2*q_1^2/2 + beta_1*q_1^4/4", H2 = "p_2^2/2 + omega_2^2*q_2^2/2 + beta_2*q_2^4/4";
        L.invariant("energy", "Energy", H1 + " + " + H2 + " + coupling*q_1^2*q_2^2", "Zero damping and driver", {arnold}, ok);
        const bool separable = ok && parameter(c, "coupling") == 0;
        L.invariant("energy_1", "First oscillator energy", H1, "Zero coupling, damping and driver", {arnold}, separable);
        L.invariant("energy_2", "Second oscillator energy", H2, "Zero coupling, damping and driver", {arnold}, separable);
    }
}

}  // namespace

std::vector<EigenfunctionEntry> builtInEigenfunctions(const SystemDescriptor& d, const SystemConfiguration& c) {
    std::vector<EigenfunctionEntry> out;
    Library L{out};
    const bool free = c.driver.amplitude == 0;
    const S& id = d.identifier;
    if (id == "nonlinear_pendulum") pendulum(L, c, free);
    else if (id == "duffing") duffing(L, c, free);
    else if (id == "harmonic_oscillator")
        L.linearModes("mode", "position", "velocity", "(stiffness_n_m/mass_kg)", "(damping_kg_s/mass_kg)", "exact",
                      "Zero driver; underdamped (damping^2 < 4 mass stiffness)",
                      "left eigenvectors of [[0, 1], [-k/m, -c/m]]: exact for the linear field",
                      free && 4 * parameter(c, "mass_kg", 1) * parameter(c, "stiffness_n_m", 1) >
                                  parameter(c, "damping_kg_s") * parameter(c, "damping_kg_s"));
    else if (id == "triple_well_duffing")
        L.invariant("energy", "Energy", "velocity^2/2 + position^6/6 - (1 + inner_saddle^2)*position^4/4 + inner_saddle^2*position^2/2",
                    "Zero damping and driver", {dlmf22}, free && parameter(c, "damping") == 0);
    else if (id == "three_leaf_hamiltonian_clover")
        L.invariant("energy", "Hamiltonian", "(x_1^2 + x_2^2)^2 - x_1^3 + 3*x_1*x_2^2", "Zero driver", {arnold}, free);
    else integrable(L, d, c, free);
    return out;
}

}  // namespace dynsyslab::core::detail
