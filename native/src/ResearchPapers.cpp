#include "DynamicalKnowledgeInternal.hpp"

namespace dynsyslab::core::detail {
std::vector<ResearchPaper> builtInResearchPapers(const std::string &identifier) {
    if(identifier=="toda_lattice") {
        ResearchPaper p;p.id="henrici-kappeler-toda-2008";
        p.title="Global Birkhoff coordinates for the periodic Toda lattice";
        p.authors={"Andreas Henrici","Thomas Kappeler"};p.year=2008;
        p.doi="10.1088/0951-7715/21/12/001";p.urls={"https://arxiv.org/abs/0802.4252"};
        p.findings="Global Birkhoff coordinates for periodic Toda on reduced symplectic leaves; the open chain has different scattering data.";
        return {p};
    }
    if (identifier == "quadratic_triangular") return triangularResearchPapers();
    if (auto oscillators = oscillatorResearchPapers(identifier); !oscillators.empty())
        return oscillators;
    if (auto extra = additionalResearchPapers(identifier); !extra.empty())
        return extra;
    if (identifier != "rpe") {
        ResearchPaper paper;
        if (identifier == "lorenz") {
            paper.id = "lorenz-1963";
            paper.title = "Deterministic Nonperiodic Flow";
            paper.authors = {"Edward N. Lorenz"};
            paper.year = 1963;
            paper.venue = "Journal of the Atmospheric Sciences 20(2), 130-141";
            paper.doi = "10.1175/1520-0469(1963)020<0130:dnf>2.0.co;2";
            paper.urls = {"https://journals.ametsoc.org/view/journals/atsc/20/2/"
                          "1520-0469_1963_020_0130_dnf_2_0_co_2.xml"};
            paper.findings = "Original dissipative-flow model and analysis of nonperiodic behavior "
                             "and sensitivity.";
        } else if (identifier == "rossler") {
            paper.id = "rossler-1976";
            paper.title = "An equation for continuous chaos";
            paper.authors = {"O. E. Rossler"};
            paper.year = 1976;
            paper.venue = "Physics Letters A 57(5), 397-398";
            paper.doi = "10.1016/0375-9601(76)90101-8";
            paper.urls = {"https://www.cfm.brown.edu/people/dobrush/am34/tutorial/"
                          "AN-EQUATION-FOR-CONTINUOUS-CHAOS_ROSSLER.pdf"};
            paper.findings = "Introduces a continuous chaotic flow with one quadratic nonlinearity "
                             "and a folded return-map geometry.";
        } else if (identifier == "chua") {
            paper.id = "matsumoto-1984-report";
            paper.title = "A Chaotic Attractor From Chua's Circuit";
            paper.authors = {"T. Matsumoto"};
            paper.year = 1984;
            paper.venue = "University of California Berkeley technical report UCB/ERL M84/36";
            paper.urls = {"https://www2.eecs.berkeley.edu/Pubs/TechRpts/1984/307.html",
                          "https://digicoll.lib.berkeley.edu/record/138262/files/ERL-m-84-36.pdf"};
            paper.findings = "Original report on the circuit's chaotic attractor; background for "
                             "the piecewise-linear model.";
        } else if (identifier == "hindmarsh_rose") {
            paper.id = "hindmarsh-rose-1984";
            paper.title = "A model of neuronal bursting using three coupled first order "
                          "differential equations";
            paper.authors = {"J. L. Hindmarsh", "R. M. Rose"};
            paper.year = 1984;
            paper.venue = "Proceedings of the Royal Society of London B 221(1222), 87-102";
            paper.doi = "10.1098/rspb.1984.0024";
            paper.urls = {"https://pubmed.ncbi.nlm.nih.gov/6144106/",
                          "https://doi.org/10.1098/rspb.1984.0024"};
            paper.findings = "Bursting model developed through equilibrium, saddle, separatrix and "
                             "limit-cycle analysis of its component dynamics.";
        } else if (identifier == "henon_heiles") {
            paper.id = "henon-heiles-1964";
            paper.title =
                "The applicability of the third integral of motion: Some numerical experiments";
            paper.authors = {"Michel Henon", "Carl Heiles"};
            paper.year = 1964;
            paper.venue = "Astronomical Journal 69, 73-79";
            paper.doi = "10.1086/109234";
            paper.urls = {"https://doi.org/10.1086/109234"};
            paper.findings = "Original numerical study of the Hamiltonian model and the "
                             "applicability of an additional integral.";
        } else if (identifier == "fput") {
            paper.id = "fermi-pasta-ulam-1955";
            paper.title = "Studies of nonlinear problems";
            paper.authors = {"Enrico Fermi", "John R. Pasta", "Stanislaw M. Ulam"};
            paper.year = 1955;
            paper.venue = "Los Alamos Scientific Laboratory report LA-1940";
            paper.doi = "10.2172/4376203";
            paper.urls = {"https://cds.cern.ch/record/425600",
                          "https://digital.library.unt.edu/ark:/67531/metadc1018377/"};
            paper.findings = "Original nonlinear-chain numerical experiments; model background for "
                             "the FPUT catalogue.";
        } else
            return {};
        paper.topics = {"model definition", "phase-space dynamics"};
        paper.relatedFactIds = {identifier + "-equilibria"};
        paper.applicability =
            "Match the publication's equations, normalization and parameters before using a "
            "parameter-specific result. The registered equilibrium seeds are separately derived "
            "from DynSysLab's equations.";
        paper.verification = "Bibliographic record or primary-paper abstract inspected; no "
                             "complete full-text review claimed.";
        if(identifier=="henon_heiles"){
            ResearchPaper orbit;orbit.id="brack-mehta-tanaka-2001";
            orbit.title="Occurrence of periodic Lame functions at bifurcations in chaotic Hamiltonian systems";
            orbit.authors={"M. Brack","M. Mehta","K. Tanaka"};orbit.year=2001;
            orbit.venue="arXiv:nlin.CD/0105048v3";
            orbit.urls={"https://arxiv.org/abs/nlin/0105048","https://epub.uni-regensburg.de/1660/1/brack11.pdf"};
            orbit.findings="Equations 30-33 give the exact straight-line orbit along each of the three symmetry axes using Jacobi sn squared.";
            orbit.applicability="Native q=lambda*s gives q''=-q+q^2. Exact bounded symmetry-axis motion and its homoclinic limit only; bifurcated-orbit approximations are not exact native solutions.";
            orbit.verification="Author-hosted full text section 4, equations 27-33 inspected; native RHS, energy and numerical-refinement tests.";
            orbit.topics={"analytic solutions","symmetry axes","separatrix","periodic orbits"};return {paper,orbit};
        }
        if(identifier=="lorenz"){
            ResearchPaper method;method.id="krauskopf-etal-2005-manifold-survey";
            method.title="A survey of methods for computing (un)stable manifolds of vector fields";
            method.authors={"Bernd Krauskopf","Hinke M. Osinga","Eusebius J. Doedel","Michael E. Henderson","John Guckenheimer","Alexander Vladimirsky","Michael Dellnitz","Oliver Junge"};
            method.year=2005;method.venue="International Journal of Bifurcation and Chaos 15(3), 763-791";
            method.doi="10.1142/S0218127405012533";
            method.urls={"https://research-information.bris.ac.uk/ws/portalfiles/portal/2988488/2004r14.pdf"};
            method.findings="Section 2 describes boundary-value continuation of geodesic fronts and adaptive triangulation of two-dimensional invariant manifolds, including the Lorenz stable manifold.";
            method.applicability="Numerical finite-domain surface approximation; seed, rank, tangency and mesh refinement must be checked independently. No global completeness certificate.";
            method.verification="Primary author survey full text section 2 inspected; independent native implementation tested on exact curved graphs and Lorenz.";
            method.topics={"stable manifolds","unstable manifolds","surface continuation"};return {paper,method};
        }
        return {paper};
    }
    std::vector<ResearchPaper> papers;
    auto add = [&](const char *id, const char *title, int year, const char *venue, const char *doi,
                   std::vector<std::string> authors,
                   std::vector<std::string> urls) -> ResearchPaper & {
        ResearchPaper paper;
        paper.id = id;
        paper.title = title;
        paper.year = year;
        paper.venue = venue;
        paper.doi = doi;
        paper.authors = std::move(authors);
        paper.urls = std::move(urls);
        papers.push_back(std::move(paper));
        return papers.back();
    };
    auto &chang = add("chang-chen-1986", "Growth of a gas bubble in a viscous fluid", 1986,
                      "Physics of Fluids 29(11), 3530-3536", "10.1063/1.865830",
                      {"Hsueh-Chia Chang", "Liang-Heng Chen"},
                      {"https://doi.org/10.1063/1.865830",
                       "https://academicweb.nd.edu/~changlab/Publications.html",
                       "https://arxiv.org/pdf/physics/9911072"});
    chang.topics = {"viscosity", "Melnikov method", "cavitation threshold",
                    "homoclinic bifurcation"};
    chang.findings = "Relates dynamic cavitation pressure to viscosity using numerical and "
                     "Melnikov perturbation analysis, guided by a double-zero normal form.";
    chang.applicability =
        "Viscous spherical Rayleigh-Plesset model; perturbative result and specified "
        "pressure-change/initial-condition setting, not a universal threshold.";
    chang.coordinates = "radial bubble dynamics; normal-form reduction";
    chang.verification =
        "Bibliography verified against the author's publication list; method identified in "
        "Harkin-Nadim-Kaper introduction and reference 4. Original full text not yet inspected.";
    chang.relatedFactIds = {"rpe-equilibria", "rpe-energy"};

    auto &manifolds = add("hegedus-kullmann-2012",
                          "Basins of attraction in a harmonically excited spherical bubble model",
                          2012, "Periodica Polytechnica Mechanical Engineering 56(2), 125-132",
                          "10.3311/pp.me.2012-2.08", {"Ferenc Hegedus", "Laszlo Kullmann"},
                          {"https://pp.bme.hu/me/article/view/1243",
                           "https://pp.bme.hu/me/article/download/1243/647"});
    manifolds.topics = {"stable manifolds", "unstable manifolds", "basins of attraction",
                        "Melnikov method", "periodic orbits"};
    manifolds.findings =
        "Computes invariant manifolds of saddle periodic solutions to delimit attraction basins; "
        "examines transverse intersections and a Melnikov estimate.";
    manifolds.applicability =
        "Harmonic forcing, finite viscosity and isothermal gas; examples place ambient pressure "
        "below vapor pressure. Computed curves require matching parameters and forcing phase.";
    manifolds.coordinates = "stroboscopic radius-velocity Poincare section";
    manifolds.verification =
        "Primary publisher full text inspected, especially section 3 and manifold figures.";
    manifolds.relatedFactIds = {"rpe-periodic-manifolds", "rpe-equilibria"};

    auto &harkin =
        add("harkin-nadim-kaper-1999", "On acoustic cavitation of slightly subcritical bubbles",
            1999, "Physics of Fluids 11, 274-287", "10.1063/1.869878",
            {"Anthony Harkin", "Ali Nadim", "Tasso J. Kaper"},
            {"https://arxiv.org/abs/physics/9911072"});
    harkin.topics = {"Blake threshold", "normal form", "Melnikov method", "viscosity"};
    harkin.findings = "Derives a damped forced distinguished-limit equation near the Blake radius "
                      "and compares cavitation predictions with the full RPE.";
    harkin.applicability =
        "Slightly subcritical isothermal bubbles; reduced coordinates and distinguished scaling. "
        "Melnikov intersection and escape thresholds must be distinguished.";
    harkin.coordinates = "scaled near-Blake radius displacement and velocity";
    harkin.verification = "Primary arXiv full text inspected.";
    harkin.relatedFactIds = {"rpe-equilibria", "rpe-periodic-manifolds"};

    auto &prosperetti = add("prosperetti-1975",
                            "Nonlinear oscillations of gas bubbles in liquids: transient solutions "
                            "and the connection between subharmonic signal and cavitation",
                            1975, "Journal of the Acoustical Society of America 57(4), 810-821",
                            "10.1121/1.380523", {"Andrea Prosperetti"},
                            {"https://pages.jh.edu/aprospe1/publications/PapersPublished/"
                             "GasBubbles/JASA_NonlinOscilTrans.pdf"});
    prosperetti.topics = {"separatrices", "subharmonic resonance", "averaging", "transients"};
    prosperetti.findings = "Analyzes critical points and separatrices in the slow subharmonic "
                           "amplitude plane, including damped transients and attraction domains.";
    prosperetti.applicability = "Weakly nonlinear near-resonant reduction; not a direct global "
                                "radius-velocity manifold formula.";
    prosperetti.coordinates = "u=C cos(phi), v=C sin(phi); averaged amplitude plane";
    prosperetti.verification =
        "Author-hosted full text inspected, equations 9 and 12-13 and Figures 2, 4 and 6.";
    prosperetti.relatedFactIds = {"rpe-prosperetti-1975"};

    auto &kudryashov =
        add("kudryashov-sinelshchikov-2014",
            "Analytical solutions of the Rayleigh equation for empty and gas-filled bubble", 2014,
            "Journal of Physics A: Mathematical and Theoretical 47, 405202",
            "10.1088/1751-8113/47/40/405202", {"Nikolay A. Kudryashov", "Dmitry I. Sinelshchikov"},
            {"https://arxiv.org/abs/1409.6699"});
    kudryashov.topics = {"analytic solutions", "Weierstrass function", "Sundman transformation"};
    kudryashov.findings = "Constructs parametric special-function solutions for inviscid empty and "
                          "gas-filled bubble cases.";
    kudryashov.applicability = "Specific inviscid parameter families; match gas law and "
                               "capillarity assumptions before evaluation.";
    kudryashov.coordinates = "transformed radius with parametric physical-time reconstruction";
    kudryashov.verification =
        "Primary arXiv text inspected; native reductions independently checked against equations.";
    kudryashov.relatedFactIds = {"rpe-energy"};

    auto &mancas =
        add("mancas-rosu-2016",
            "Evolution of spherical cavitation bubbles: parametric and closed-form solutions", 2016,
            "Physics of Fluids 28, 022009", "10.1063/1.4942237", {"S. C. Mancas", "H. C. Rosu"},
            {"https://arxiv.org/abs/1508.01157"});
    mancas.topics = {"analytic solutions", "capillarity", "Weierstrass function", "Abel equation"};
    mancas.findings = "Studies parametric solutions of vacuous-bubble RPE with surface tension; "
                      "discusses numerical treatment when viscosity is nonzero.";
    mancas.applicability = "Vacuous-bubble assumptions; the viscous discussion does not supply a "
                           "general closed-form viscous trajectory.";
    mancas.coordinates = "radius with parametric time";
    mancas.verification =
        "Primary arXiv text inspected; first preprint 2015, journal publication 2016.";
    mancas.relatedFactIds = {"rpe-energy"};
    return papers;
}
} // namespace dynsyslab::core::detail
