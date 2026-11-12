#include "LatticeResearch.hpp"
namespace dynsyslab::core::detail {
DynamicalKnowledge latticeResearch(const std::string &id) {
    DynamicalKnowledge k;
    k.researchStatus =
        "reviewed; thermodynamic exact results are distinct from stochastic sweep trajectories";
    DynamicalFact scope;
    scope.id = id + "-statistical-scope";
    scope.title = "Equilibrium statistics versus sweep dynamics";
    scope.kind = "applicability_limit";
    scope.statement =
        "The archived order parameter and energy are observations of a finite stochastic spin "
        "configuration. Exact equilibrium partition functions or thermodynamic critical points are "
        "not deterministic time trajectories for those two observations.";
    scope.conditions = "Check lattice dimension, coupling normalization, external field, boundary "
                       "conditions and finite-size effects. No ODE equilibrium/manifold analysis "
                       "is inferred from the two recorded observations.";
    scope.coordinates = "Monte Carlo sweeps and lattice order/energy observations; k_B=1";
    ResearchPaper p;
    if (id == "ising") {
        p.id = "onsager-1944";
        p.title =
            "Crystal Statistics. I. A Two-Dimensional Model with an Order-Disorder Transition";
        p.authors = {"Lars Onsager"};
        p.year = 1944;
        p.venue = "Physical Review 65, 117";
        p.doi = "10.1103/PhysRev.65.117";
        p.urls = {"https://journals.aps.org/pr/abstract/10.1103/PhysRev.65.117"};
        p.findings =
            "Exact zero-field square-lattice partition function and thermodynamic transition.";
        p.applicability = "Two-dimensional ferromagnetic equilibrium model, zero field; no "
                          "equality with a single finite Monte Carlo history.";
        p.verification = "Primary publisher abstract and bibliographic record inspected.";
        k.papers.push_back(p);
        p.id = "yang-1952";
        p.title = "The Spontaneous Magnetization of a Two-Dimensional Ising Model";
        p.authors = {"C. N. Yang"};
        p.year = 1952;
        p.venue = "Physical Review 85, 808";
        p.doi = "10.1103/PhysRev.85.808";
        p.urls = {"https://journals.aps.org/pr/abstract/10.1103/PhysRev.85.808"};
        p.findings = "Exact spontaneous magnetization and long-range order in the two-dimensional "
                     "equilibrium Ising model.";
        k.papers.push_back(p);
    } else if (id == "potts") {
        p.id = "beffara-duminil-copin-2012";
        p.title = "The self-dual point of the two-dimensional random-cluster model is critical";
        p.authors = {"V. Beffara", "H. Duminil-Copin"};
        p.year = 2012;
        p.urls = {"https://www.unige.ch/~duminil/publi/FK_critical_point.pdf"};
        p.findings = "Theorem 3 gives beta_c=log(1+sqrt(q)) for the zero-field square-lattice "
                     "Potts Hamiltonian -sum delta(s_i,s_j).";
        p.applicability = "Native Hamiltonian -J*sum delta(s_i,s_j), ferromagnetic J>0 and zero "
                          "external field: T_c=J/log(1+sqrt(q)) in the thermodynamic limit; finite "
                          "simulations have finite-size effects.";
        p.verification = "Primary author PDF, Hamiltonian convention and Theorem 3 inspected.";
        k.papers.push_back(p);
    } else {
        p.id = "viswanathan-etal-2022";
        p.title = "What does it take to solve the 3D Ising model? Minimal necessary conditions for "
                  "a valid solution";
        p.authors = {"G. M. Viswanathan", "M. A. G. Portillo", "E. P. Raposo", "M. G. E. da Luz"};
        p.year = 2022;
        p.urls = {"https://arxiv.org/abs/2205.12357"};
        p.findings = "Necessary consistency tests for proposed exact three-dimensional Ising "
                     "partition functions, with discussion of failures of claimed solutions.";
        p.applicability =
            "The native cubic-lattice model does not inherit Onsager's two-dimensional formula. No "
            "verified general analytic sweep trajectory is registered.";
        p.verification = "Primary arXiv abstract and bibliographic record inspected; no claimed "
                         "new exact solution is treated as established.";
        k.papers.push_back(p);
    }
    for (auto &paper : k.papers) {
        paper.topics = {"equilibrium statistical mechanics", "analytic applicability"};
        scope.references.push_back(paper.id);
    }
    k.facts.push_back(scope);
    return k;
}
} // namespace dynsyslab::core::detail
