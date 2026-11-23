#include "DynamicalKnowledgeInternal.hpp"

namespace dynsyslab::core::detail {
namespace {
ResearchPaper paper(std::string id, std::string title, std::vector<std::string> authors, int year, std::string venue,
                    std::string doi, std::string findings, std::string applicability,
                    std::vector<std::string> topics) {
    ResearchPaper p;
    p.id = std::move(id);
    p.title = std::move(title);
    p.authors = std::move(authors);
    p.year = year;
    p.venue = std::move(venue);
    p.doi = std::move(doi);
    if (!p.doi.empty()) p.urls = {"https://doi.org/" + p.doi};
    p.findings = std::move(findings);
    p.applicability = std::move(applicability);
    p.verification = p.doi.empty() ? "Citation details from the published record; no DOI is registered here."
                                   : "Bibliographic record checked against Crossref.";
    p.topics = std::move(topics);
    return p;
}

ResearchPaper arnoldMechanics(std::string applicability) {
    return paper("arnold-1989", "Mathematical Methods of Classical Mechanics", {"V. I. Arnold"}, 1989,
                 "Graduate Texts in Mathematics 60, Springer, second edition", "10.1007/978-1-4757-2063-1",
                 "Reference for Hamiltonian flows, integrability and action-angle coordinates.",
                 std::move(applicability), {"Hamiltonian mechanics", "integrability", "action-angle coordinates"});
}
} // namespace

std::vector<ResearchPaper> oscillatorResearchPapers(const std::string &id) {
    if (id == "cosh_oscillator")
        return {paper("toda-1967", "Vibration of a Chain with Nonlinear Interaction", {"Morikazu Toda"}, 1967,
                      "Journal of the Physical Society of Japan 22(2), 431-436", "10.1143/JPSJ.22.431",
                      "Introduces the exponential-interaction chain and its exact periodic solutions.",
                      "The native cosh potential has the two-exponential form of a reduced Toda-type pair; "
                      "its parameters are native.",
                      {"Toda interaction", "exact periodic solutions"})};
    if (id == "coupled_duffing_hamiltonian")
        return {paper("kovacic-brennan-2011", "The Duffing Equation: Nonlinear Oscillators and their Behaviour",
                      {"Ivana Kovacic", "Michael J. Brennan"}, 2011, "Wiley", "10.1002/9780470977859",
                      "Reference for Duffing oscillators, their energy structure and frequency response.",
                      "Applies to each uncoupled oscillator; the coupled pair is integrable only at zero coupling.",
                      {"Duffing oscillator", "Hamiltonian coupling"})};
    if (id == "euler_top")
        return {arnoldMechanics("Covers free rigid-body motion; the native parameters are the principal moments.")};
    if (id == "kepler_planar")
        return {arnoldMechanics("Covers the Kepler problem; the native analytic scope is bound noncollision orbits.")};
    if (id == "radial_twist")
        return {paper("arnold-kozlov-neishtadt-2006", "Mathematical Aspects of Classical and Celestial Mechanics",
                      {"Vladimir I. Arnold", "Valery V. Kozlov", "Anatoly I. Neishtadt"}, 2006,
                      "Encyclopaedia of Mathematical Sciences 3, Springer, third edition", "10.1007/978-3-540-48926-9",
                      "Reference for twist, action-dependent frequency, averaging and perturbation theory.",
                      "The native system is an elementary integrable twist model.",
                      {"twist", "averaging", "perturbation theory"})};
    if (id == "isotonic_oscillator")
        return {paper("calogero-1969", "Solution of a Three-Body Problem in One Dimension", {"F. Calogero"}, 1969,
                      "Journal of Mathematical Physics 10(12), 2191-2196", "10.1063/1.1664820",
                      "Exact treatment of a one-dimensional system with harmonic and inverse-square interactions.",
                      "The native model is the single-coordinate harmonic-plus-inverse-square oscillator on q>0.",
                      {"isotonic oscillator", "inverse-square potential"})};
    if (id == "morse_oscillator")
        return {paper("morse-1929", "Diatomic Molecules According to the Wave Mechanics. II. Vibrational Levels",
                      {"Philip M. Morse"}, 1929, "Physical Review 34(1), 57-64", "10.1103/PhysRev.34.57",
                      "Introduces the Morse potential.", "The native model uses the classical bounded orbits, 0<H<D.",
                      {"Morse potential", "anharmonic oscillator"})};
    if (id == "poschl_teller")
        return {paper("poschl-teller-1933", "Bemerkungen zur Quantenmechanik des anharmonischen Oszillators",
                      {"G. Poschl", "E. Teller"}, 1933, "Zeitschrift fur Physik 83, 143-151", "",
                      "Introduces the solvable hyperbolic potential family.",
                      "The native model uses the classical bounded orbits of the sech-squared well.",
                      {"Poschl-Teller potential", "solvable potentials"})};
    if (id == "triple_well_duffing")
        return {paper("chacon-diaz-bejarano-1995", "Homoclinic and heteroclinic chaos in a triple-well oscillator",
                      {"R. Chacon", "J. Diaz Bejarano"}, 1995, "Journal of Sound and Vibration 186(2), 269-278",
                      "10.1006/jsvi.1995.0448", "Studies homoclinic and heteroclinic chaos in a triple-well oscillator.",
                      "Provides context; the native polynomial normalization is its own API definition.",
                      {"triple well", "homoclinic chaos", "heteroclinic chaos"})};
    if (id == "three_leaf_hamiltonian_clover") {
        auto p = paper("kirichuka-2012", "Periodic solutions for Hamiltonian systems", {"A. Kirichuka"}, 2012,
                       "Section 4.6, equations (23)-(24)", "",
                       "Source of the three-leaf Hamiltonian equations.",
                       "The native implementation uses the published polynomial with an input gain.",
                       {"Hamiltonian systems", "three-leaf clover"});
        p.urls = {"https://www.lumii.lv/uploads/sadirb./2012/LUMII_2012_A_Kirichuka.pdf"};
        return {p};
    }
    return {};
}
} // namespace dynsyslab::core::detail
