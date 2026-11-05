#include <dynsyslab/core/Delay.hpp>
namespace dynsyslab::core {
DynamicalKnowledge mackeyGlassKnowledge() {
    DynamicalKnowledge result;
    result.researchStatus =
        "reviewed; constant-history analytic regimes and delay-state limitations recorded";
    ResearchPaper model;
    model.id = "mackey-glass-1977";
    model.title = "Oscillation and chaos in physiological control systems";
    model.authors = {"Michael C. Mackey", "Leon Glass"};
    model.year = 1977;
    model.venue = "Science 197(4300), 287-289";
    model.doi = "10.1126/science.267326";
    model.urls = {"https://pubmed.ncbi.nlm.nih.gov/267326/"};
    model.topics = {"delay equations", "oscillations", "physiological feedback"};
    model.findings = "Nonlinear physiological feedback with delay can produce oscillatory and "
                     "aperiodic dynamics.";
    model.applicability = "Native feedback convention beta*x(t-tau)/(1+x(t-tau)^n)-gamma*x(t); "
                          "consult parameter and history assumptions.";
    model.verification = "Primary publication record and abstract inspected; no claim that the "
                         "general nonlinear equation has a closed-form trajectory.";
    ResearchPaper delayed;
    delayed.id = "khusainov-pokojovy-azizbayov-2015";
    delayed.title = "On the Cauchy problem for a linear harmonic oscillator with pure delay";
    delayed.authors = {"Denys Ya. Khusainov", "Michael Pokojovy", "Elvin I. Azizbayov"};
    delayed.year = 2015;
    delayed.venue = "Advances in Difference Equations 2015, 197";
    delayed.doi = "10.1186/s13662-015-0538-z";
    delayed.urls = {"https://d-nb.info/1110771886/34"};
    delayed.topics = {"delayed exponential", "method of steps", "explicit solutions"};
    delayed.findings = "Section 3.2 defines the piecewise-polynomial delayed exponential and "
                       "verifies its delay identity.";
    delayed.applicability = "The native scalar exponent=0 limit uses the first-order delayed "
                            "fundamental kernel and its constant-history convolution, derived "
                            "independently; it is not the second-order oscillator in the paper.";
    delayed.verification =
        "Primary full text, delayed exponential definition and derivative relation inspected; "
        "native scalar residual and interval-boundary tests verify the reformulation.";
    result.papers = {model, delayed};
    DynamicalFact equilibria;
    equilibria.id = "mackey-glass-constant-histories";
    equilibria.title = "Constant-history equilibria";
    equilibria.kind = "equilibrium_relation";
    equilibria.statement =
        "Candidate histories satisfy beta*h/(1+h^n)=gamma*h. For positive n and beta>gamma>0 the "
        "positive branch is h=(beta/gamma-1)^(1/n). Delay stability is infinite-dimensional.";
    equilibria.formula =
        R"(h=0;\quad h=(\beta/\gamma-1)^{1/n};\quad\lambda+\gamma-f'(h)e^{-\lambda\tau}=0)";
    equilibria.conditions = "Check the real feedback domain, zero-state definition and signs. A "
                            "scalar instantaneous Jacobian does not determine delay stability.";
    equilibria.coordinates = "Native state and time; full constant history on [-tau,0]";
    equilibria.references = {model.id};
    result.facts = {equilibria};
    return result;
}
} // namespace dynsyslab::core
