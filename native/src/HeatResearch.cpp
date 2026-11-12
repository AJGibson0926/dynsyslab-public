#include <dynsyslab/core/PDE.hpp>
namespace dynsyslab::core {
DynamicalKnowledge heatEquationKnowledge() {
    DynamicalKnowledge result;
    result.researchStatus =
        "reviewed; continuum Fourier solutions and discrete-stencil identities distinguished";
    ResearchPaper paper;
    paper.id = "mit-dyatlov-heat-equation";
    paper.title = "18.03 Differential Equations: The heat equation (Section 10.1)";
    paper.authors = {"Semyon Dyatlov"};
    paper.venue = "Massachusetts Institute of Technology, course notes";
    paper.urls = {"https://math.mit.edu/~dyatlov/18.03/lectures/Section10-1.pdf",
                  "https://math.mit.edu/~dyatlov/18.03/rs/RS12s.pdf"};
    paper.topics = {"heat equation", "Fourier series", "Dirichlet and Neumann eigenfunctions"};
    paper.findings =
        "Boundary eigenfunctions evolve independently with exponential diffusion rates.";
    paper.applicability =
        "Constant coefficients on an interval; match boundary conditions and initial function. "
        "Native finite-difference eigenvalues differ from continuum eigenvalues.";
    paper.verification =
        "Primary course material and explicit sine-series worked solutions inspected. Native "
        "mode/rate and mesh-convergence tests verify the implementation.";
    result.papers = {paper};
    DynamicalFact modes;
    modes.id = "heat-fourier-semigroup";
    modes.title = "Diffusion modes and equilibrium space";
    modes.kind = "analytic_family";
    modes.statement = "Periodic modes have wave number 2*pi*k/L; homogeneous Dirichlet/Neumann "
                      "modes have pi*k/L. With zero decay the constant mode is an equilibrium for "
                      "periodic/Neumann conditions; zero-value boundaries exclude it. Zero "
                      "diffusivity and zero decay make every field stationary.";
    modes.formula = R"(\lambda_k=-\kappa q_k^2-c)";
    modes.conditions =
        "kappa>=0 and c>=0; finite domain; specified homogeneous boundary conditions. Stationary "
        "configurations require checking the actual mesh operator.";
    modes.coordinates = "Physical x in [0,L], time in seconds, native field units";
    modes.references = {paper.id};
    result.facts = {modes};
    return result;
}
} // namespace dynsyslab::core
