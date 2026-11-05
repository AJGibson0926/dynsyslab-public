#include "DynamicalKnowledgeInternal.hpp"
namespace dynsyslab::core::detail {
std::vector<ResearchPaper> additionalResearchPapers(const std::string &id) {
    ResearchPaper p;
    p.verification = "Primary bibliographic record or abstract inspected; scope is stated below.";
    if (id == "quadratic_saddle_node" || id == "cubic_bifurcation") {
        p.id = "kuznetsov-2004";
        p.title = "Elements of Applied Bifurcation Theory";
        p.authors = {"Yuri A. Kuznetsov"};
        p.year = 2004;
        p.venue = "Applied Mathematical Sciences 112, Springer, third edition";
        p.doi = "10.1007/b98848";
        p.urls = {"https://link.springer.com/book/10.1007/b98848"};
        p.findings = "Reference for fold, pitchfork and cusp normal forms, local stability and "
                     "continuation of bifurcation branches.";
        p.applicability = "The native scalar systems are explicit normal-form fixtures. Their "
                          "computed branches and stability still use the configured native "
                          "coefficients rather than book-specific parameter values.";
        p.verification = "Springer bibliographic record and chapter scope inspected; native "
                         "equilibria and branch residuals are independently tested.";
        p.topics = {"fold bifurcation", "pitchfork bifurcation", "cusp", "normal forms"};
        return {p};
    }
    if (id == "radial_limit_cycle") {
        p.id = "guckenheimer-holmes-1983";
        p.title = "Nonlinear Oscillations, Dynamical Systems, and Bifurcations of Vector Fields";
        p.authors = {"John Guckenheimer", "Philip Holmes"};
        p.year = 1983;
        p.venue = "Applied Mathematical Sciences 42, Springer New York";
        p.doi = "10.1007/978-1-4612-1140-2";
        p.urls = {"https://link.springer.com/book/10.1007/978-1-4612-1140-2"};
        p.findings = "Reference for planar flows, local bifurcation normal forms, periodic orbits "
                     "and their stability.";
        p.applicability = "The native radial system is an exactly solvable planar normal-form "
                          "fixture. Its cycle radius and Floquet data are derived and checked "
                          "against the configured native coefficients.";
        p.verification = "Springer bibliographic record and chapter scope inspected; native "
                         "orbit, return-map and stability identities are independently tested.";
        p.topics = {"limit cycles", "normal forms", "periodic orbits", "Floquet stability"};
        return {p};
    }
    if (id == "harmonic_oscillator" || id == "duffing" || id == "van_der_pol" ||
        id == "nonlinear_pendulum" || id == "spring_coupled_pendulums") {
        p.id = "teschl-ode";
        p.title = "Ordinary Differential Equations and Dynamical Systems";
        p.authors = {"Gerald Teschl"};
        p.urls = {"https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode.pdf"};
        p.venue = "American Mathematical Society, Graduate Studies in Mathematics 140";
        p.year = 2012;
        p.findings = "Reference for linear systems, conservative planar dynamics, equilibria and "
                     "periodic-orbit methods.";
        p.applicability = "Use the stated hypotheses; native model-specific equilibria and energy "
                          "formulas are independently derived and checked.";
        p.topics = {"equilibria", "phase-space methods", "analytic solutions"};
        ResearchPaper special;
        special.id = "nist-dlmf-elliptic";
        special.title = "Digital Library of Mathematical Functions: Jacobian Elliptic Functions";
        special.authors = {"NIST DLMF"};
        special.urls = {"https://dlmf.nist.gov/22.19", "https://dlmf.nist.gov/22.13",
                        "https://dlmf.nist.gov/23.10"};
        special.findings = "Elliptic-function identities and pendulum applications used by "
                           "analytic trajectory evaluators.";
        special.applicability =
            "Conservative integrable limits only; distinguish elliptic parameter m from modulus k.";
        special.verification =
            "Primary DLMF formulas inspected and native residual/identity tests passed.";
        special.topics = {"Jacobi functions", "Weierstrass functions", "pendulum"};
        return {p, special};
    }
    if (id == "marmottant") {
        p.id = "marmottant-2005";
        p.title = "A model for large amplitude oscillations of coated bubbles accounting for "
                  "buckling and rupture";
        p.authors = {"Philippe Marmottant", "Sander van der Meer", "Marcia Emmer",
                     "Michel Versluis",     "Nico de Jong",        "Sascha Hilgenfeldt",
                     "Detlef Lohse"};
        p.year = 2005;
        p.venue = "Journal of the Acoustical Society of America 118(6), 3499-3505";
        p.doi = "10.1121/1.2109427";
        p.urls = {"https://liphy-annuaire.univ-grenoble-alpes.fr/pages_personnelles/"
                  "philippe_marmottant/Publications_files/Marmottant2005_JASACoatedBubble.pdf"};
        p.findings = "Defines buckled, elastic and ruptured shell regimes; equilibrium pressure "
                     "balance and a small-amplitude linearization are derived.";
        p.applicability = "Linearized oscillations are approximations, not exact full-model "
                          "solutions. Native rupture history is an additional state: a phase-space "
                          "sidecar must stay in one memory-independent smooth shell regime.";
        p.verification = "Author-hosted full text inspected, equations 2-6 and sections II.B-II.C; "
                         "native hidden rupture-state behavior inspected.";
        p.topics = {"shell regimes", "hysteresis", "equilibria", "linearization"};
    } else if (id == "reddy_szeri") {
        p.id = "reddy-szeri-2002-translation";
        p.title =
            "Coupled dynamics of translation and collapse of acoustically driven microbubbles";
        p.authors = {"Anil J. Reddy", "Andrew J. Szeri"};
        p.year = 2002;
        p.venue = "Journal of the Acoustical Society of America 112(4), 1346-1352";
        p.doi = "10.1121/1.1502899";
        p.urls = {"https://pubmed.ncbi.nlm.nih.gov/12398441/"};
        p.findings = "Derives coupled radial and translational microbubble dynamics and compares "
                     "with experiment and decoupled models.";
        p.applicability = "Full translating/collapsing dynamics require numerical solution. A separately verified stationary-radius shape family follows from the native closure; radial-only RPE formulas are not transferred.";
        p.topics = {"model definition", "radial-translation coupling"};
        ResearchPaper shape;
        shape.id="reddy-szeri-2002-shape";
        shape.title="Shape stability of unsteadily translating bubbles";
        shape.authors={"Anil J. Reddy","Andrew J. Szeri"};
        shape.year=2002;
        shape.venue="Physics of Fluids 14(7), 2216-2224";
        shape.doi="10.1063/1.1483840";
        shape.urls={"https://doi.org/10.1063/1.1483840","https://mech.ubc.ca/andrew-szeri/"};
        shape.findings="Shape-mode evolution for translating bubbles; provides the model provenance for shape-stability studies.";
        shape.applicability="The stationary-radius oscillator formula is independently derived from native equations. No claim that the paper supplies an exact nonlinear hydrodynamic trajectory.";
        shape.verification="Author publication list verifies title, authors, year, venue, pages and DOI; full paper was not accessible in this review. Native reduction is separately tested.";
        shape.topics={"shape modes","stability","translation"};
        shape.relatedFactIds={"reddy-szeri-shape-limit"};
        return {p,shape};
    } else if (id == "hodgkin_huxley") {
        p.id = "hodgkin-huxley-1952";
        p.title = "A quantitative description of membrane current and its application to "
                  "conduction and excitation in nerve";
        p.authors = {"A. L. Hodgkin", "A. F. Huxley"};
        p.year = 1952;
        p.venue = "Journal of Physiology 117(4), 500-544";
        p.doi = "10.1113/jphysiol.1952.sp004764";
        p.urls = {"https://physoc.onlinelibrary.wiley.com/doi/10.1113/jphysiol.1952.sp004764"};
        p.findings = "Original membrane-current and gating model.";
        p.applicability =
            "Fixed-voltage exponentials normally describe a clamp experiment. The native analytic family additionally verifies a constant-voltage invariant set of the free-voltage equations.";
        p.topics = {"model definition", "gating", "voltage clamp"};
    } else if (id == "double_pendulum") {
        return doublePendulumResearchPapers();
    } else if (id == "kuramoto") {
        p.id = "strogatz-2000";
        p.title = "From Kuramoto to Crawford: exploring the onset of synchronization in "
                  "populations of coupled oscillators";
        p.authors = {"Steven H. Strogatz"};
        p.year = 2000;
        p.venue = "Physica D 143, 1-20";
        p.doi = "10.1016/s0167-2789(00)00094-4";
        p.urls = {"https://www.sciencedirect.com/science/article/pii/S0167278900000944",
                  "https://web.mit.edu/~jadbabai/www/ESE680/Strogatz_Kuramoto.pdf"};
        p.findings = "Reviews the onset of collective synchronization and stability in oscillator "
                     "populations.";
        p.applicability = "Distinguish thermodynamic/population results from the finite native "
                          "oscillator array and its frequency distribution.";
        p.topics = {"synchronization", "stability", "Kuramoto model"};
        ResearchPaper pair;
        pair.id="ha-ko-park-zhang-2016";
        pair.title="Collective synchronization of classical and quantum oscillators";
        pair.authors={"Seung-Yeal Ha","Dongnam Ko","Jinyeong Park","Xiongtao Zhang"};
        pair.year=2016;pair.venue="EMS Surveys in Mathematical Sciences 3, 209-267";
        pair.doi="10.4171/EMSS/17";
        pair.urls={"https://ems.press/content/serial-article-files/36988"};
        pair.findings="Example 4.2 gives the exact two-oscillator Adler phase-difference equation and locking condition.";
        pair.applicability="The native symmetric frequency pair has detuning 2 times frequency_spread. Phase-shift neutrality is retained.";
        pair.verification="Publisher full text inspected, Example 4.2; independent projective-exponential and native residual tests.";
        pair.topics={"Adler equation","phase locking","exact two-oscillator reduction"};
        pair.relatedFactIds={"kuramoto-relative-locking"};
        ResearchPaper peer;
        peer.id="gupta-campa-ruffo-2014";
        peer.title="Kuramoto model of synchronization: Equilibrium and nonequilibrium aspects";
        peer.authors={"Shamik Gupta","Alessandro Campa","Stefano Ruffo"};
        peer.year=2014;peer.venue="Journal of Statistical Mechanics: Theory and Experiment, R08001";
        peer.doi="10.1088/1742-5468/14/08/R08001";
        peer.urls={"https://arxiv.org/abs/1403.2083"};
        peer.findings="Synchronization review spanning finite models, mean-field descriptions, noise and inertia.";
        peer.applicability="Provides context; continuum thresholds are not substituted for native finite-array stability.";
        peer.verification="DynSys Class Hub source record consulted read-only; primary arXiv abstract and journal DOI verified.";
        peer.topics={"synchronization","finite-size cautions","mean-field limits"};
        return {p,pair,peer};
    } else
        return {};
    return {p};
}
} // namespace dynsyslab::core::detail
