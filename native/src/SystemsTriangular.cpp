#include "SystemSupport.hpp"
#include "AnalyticSystems.hpp"
#include "DynamicalKnowledgeInternal.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {
constexpr auto brunton = "https://doi.org/10.1371/journal.pone.0150171";
constexpr auto kaiser = "https://arxiv.org/abs/1707.01146";
}
std::unique_ptr<ODESystem> makeTriangularSystem(SystemConfiguration config) {
    double a = -1, b = -3, c = 1;
    if (config.preset == "brunton_2016") a = -.05, b = -1;
    else if (config.preset == "kaiser_2021") a = -.1, b = -1;
    else if (config.preset == "infinite_closure_saddle") a = 1, b = -3;
    else if (!config.preset.empty() && config.preset != "gradient_closed_reference")
        throw std::invalid_argument("Unknown quadratic_triangular preset");
    if (config.preset.empty()) config.preset = "gradient_closed_reference";
    a = resolvedParameter(config, "a", a);
    b = resolvedParameter(config, "b", b);
    c = resolvedParameter(config, "c", c);
    const double g1 = resolvedParameter(config, "input_gain_1", 1);
    const double g2 = resolvedParameter(config, "input_gain_2", 0);
    for (double value : {a,b,c,g1,g2})
        if (!std::isfinite(value))
            throw std::invalid_argument("Triangular system parameters must be finite");
    auto initial = selectedInitial(config, {1,0});
    for (double value : initial)
        if (!std::isfinite(value))
            throw std::invalid_argument("Triangular initial state must be finite");
    const auto x1 = AnalyticExpression::state(2,0);
    const auto x2 = AnalyticExpression::state(2,1);
    const auto one = AnalyticExpression::constant(1,2);
    const auto f1 = a*x1, f2 = b*x2+c*x1.power(2);
    std::vector<Observable> observables;
    const auto add = [&](const char* name, const char* label, const AnalyticExpression& e) {
        observables.push_back(Observable::fromExpressions(
            {name,label,"dimensionless","dimensionless",{},false},{e}));
    };
    add("one","1",one);
    add("x_1_squared","x_1^2",x1.power(2));
    add("x_1_dot","dx_1/dt (zero input)",f1);
    add("x_2_dot","dx_2/dt (zero input)",f2);
    std::map<std::string,std::vector<std::string>> presets{
        {"state",{"x_1","x_2"}}};
    if (b != 2*a || c == 0) {
        const double k = c == 0 ? 0 : c/(b-2*a);
        if (std::isfinite(k)) {
            const auto eta = x2+k*x1.power(2);
            add("eta","eta = x_2 + c/(b-2a) x_1^2",eta);
            add("eta_dot","d eta/dt (zero input)",b*eta);
            presets["eigenfunctions"] = {"one","x_1","x_1_squared","eta"};
            if (a==1 && b==-3 && c==1) {
                const auto invariant = x1.power(3)*eta;
                add("mixed_invariant","I = x_1^3 eta",invariant);
                add("invariant_exponential","exp(I)",invariant.unary("exp"));
            }
        }
    }
    auto descriptor = descriptorFor("quadratic_triangular");
    if (!presets.count("eigenfunctions")) descriptor.observableTimeDerivatives.erase("eta");
    return std::make_unique<FunctionalSystem>(
        std::move(descriptor),config,std::move(initial),
        [=](double,const double* x,double u,double* f) {
            f[0]=a*x[0]+g1*u; f[1]=b*x[1]+c*x[0]*x[0]+g2*u;
        },
        [=](double,const double* x) { return std::vector<double>{a,0,2*c*x[0],b}; },
        [=](const double*) { return std::vector<double>{g1,g2}; },
        [](double,const double* x) { return !std::isfinite(x[0]) || !std::isfinite(x[1]); },
        std::vector<double>{0,0},std::move(observables),std::map<std::string,std::string>{},
        std::move(presets),std::vector<AnalyticExpression>{f1,f2},
        std::vector<AnalyticExpression>{AnalyticExpression::constant(g1,2),
                                        AnalyticExpression::constant(g2,2)});
}

std::vector<AnalyticSolutionFamily> triangularAnalyticSolutions(const SystemConfiguration& c) {
    const double a=c.parameters.at("a"), b=c.parameters.at("b"), q=c.parameters.at("c");
    return {{{"quadratic-triangular-exact","Exact triangular flow",
        R"(x_1=e^{a s}x_{1,0},\quad x_2=e^{b s}x_{2,0}+c x_{1,0}^2(e^{2as}-e^{bs})/(2a-b),\ s=t-t_0)",
        "No forcing; at 2a=b replace the quotient by s exp(b s)",
        "closed form, including quadratic resonance",{brunton,kaiser}},
        [=](const AnalyticState& x,double t0)->AnalyticTrajectory {
            if(c.driver.amplitude!=0)
                return {{},"Exact autonomous flow requires zero driver amplitude"};
            return {[=](double t) {
                const double s=t-t0, e1=std::exp(a*s), e2=std::exp(b*s), delta=2*a-b;
                // expm1 preserves the resonant limit; the separated exponentials
                // avoid 0*infinity in stable long-time nonresonant evaluations.
                const double h=delta==0 ? s*e2 : (std::abs(delta*s)<.5
                    ? e2*std::expm1(delta*s)/delta : (std::exp(2*a*s)-e2)/delta);
                return AnalyticState{e1*x[0],e2*x[1]+q*x[0]*x[0]*h};
            },{}};
        }}};
}

DynamicalKnowledge triangularKnowledge(const SystemConfiguration& config) {
    DynamicalKnowledge r;
    r.researchStatus="reviewed_exact_reduction";
    const auto& p=config.parameters;
    const double a=p.at("a"),b=p.at("b"),c=p.at("c");
    const bool free=config.driver.amplitude==0, nonresonant=b!=2*a || c==0;
    const auto add=[&](const char* id,const char* title,const char* kind,
                      const char* statement,const char* formula,bool applicable,
                      const char* conditions) {
        r.facts.push_back({id,title,kind,statement,formula,conditions,
            "dimensionless state x_bar=(x_1,x_2)",{brunton,kaiser},applicable,
            applicable?"Configured equations satisfy the stated conditions":
                       "Configured parameters or driver do not satisfy the stated conditions",
            "Derived for the native equations; cited papers supply benchmark context. "
            "The analytic benchmark report gives the separate proofs of the derived facts."});
    };
    add("quadratic_triangular-equilibria","Origin and linear rates","equilibria",
        "The origin is an equilibrium of the unforced field; it is unique when a*b!=0. "
        "It is globally asymptotically stable for a<0 and b<0.",
        R"(x_1=x_2=0,\quad \operatorname{spec} Df(0)=\{a,b\})",free,"Zero driver");
    r.equilibriumGuides=[free](const auto& lo,const auto& hi) {
        std::vector<EquilibriumGuide> out;
        if(free && lo.size()==2 && hi.size()==2 &&
           lo[0]<=0 && hi[0]>=0 && lo[1]<=0 && hi[1]>=0)
            out.push_back({"quadratic_triangular-equilibria",{0,0}});
        return out;
    };
    add("quadratic_triangular-eigenfunctions","Exact gradient-closed eigenfunctions","eigenfunctions",
        "The constant, x_1, x_1^2 and eta have exact rates 0,a,2a,b. "
        "Their partial derivatives close in the same span; derivatives above order two vanish.",
        R"(z=(1,x_1,x_1^2,\eta)^T,\quad \eta=x_2+c x_1^2/(b-2a))",
        free&&nonresonant,"Zero driver; b!=2a, or c=0 with eta=x_2");
    // Library eigenfunctions (loaded by libraryEigenfunctions). Brunton et al. (2016), eqs. 25-27;
    // Kaiser et al. (2021), eqs. 50-52. Df(0)=diag(a,b), so the linearized principal
    // eigenfunctions at the origin are x_1 and x_2; eta is the exact one that x_2 linearizes.
    const std::string eta=c==0?"x_2":"x_2 + c/(b - 2*a)*x_1^2";
    r.eigenfunctions={
        {"x_1","Exact eigenfunction x_1","exact","x_1","0","a","0","Zero driver",
         "x_1' = a x_1",{brunton,kaiser},free},
        {"eta","Exact eigenfunction eta","exact",eta,"0","b","0","Zero driver; b!=2a, or c=0",
         "eta' = x_2' + 2c x_1 x_1'/(b-2a) = b eta",{brunton,kaiser},free&&nonresonant},
        {"principal_a","Linearized principal eigenfunction at the origin (rate a)","principal","x_1","0","a","0",
         "Zero driver; near the origin","left eigenvector (1,0) of Df(0)",{brunton},free},
        {"principal_b","Linearized principal eigenfunction at the origin (rate b)","principal","x_2","0","b","0",
         "Zero driver; near the origin","left eigenvector (0,1) of Df(0); linearization of eta",{brunton},free}};
    add("quadratic_triangular-slow-manifold","Invariant parabola and fast deviation","invariant_manifold",
        "The zero set eta=0 is invariant; eta evolves at rate b. "
        "It is the distinguished slow graph when b<a<0. At the reference parameters "
        "it differs from the x_2 nullcline.",
        R"(x_2=-c x_1^2/(b-2a),\quad \eta(t)=e^{bt}\eta(0))",
        free&&nonresonant,"Zero driver; b!=2a, or c=0 with graph x_2=0; slow interpretation additionally b<a<0");
    add("quadratic_triangular-cubic-resonance","Nonunique tangent invariant graphs","invariant_manifold",
        "At the reference parameters, every graph x_2=x_1^2+C*x_1^3 is invariant. "
        "Tangency to the slow eigendirection does not uniquely select the parabola.",
        R"(x_2=x_1^2+C x_1^3,\quad C\in\mathbb R)",
        free&&a==-1&&b==-3&&c==1,"Reference a=-1,b=-3,c=1; zero driver");
    add("quadratic_triangular-quadratic-resonance","Jordan coupling at quadratic resonance","resonance",
        "For b=2a and c!=0 the natural state-containing polynomial lift is invariant "
        "but not diagonalizable. The exact flow remains available.",
        R"(x_2(t)=e^{bt}(x_{2,0}+ct x_{1,0}^2))",
        free&&b==2*a&&c!=0,"Zero driver; b=2a and c!=0");
    add("quadratic_triangular-smooth-eigenfunctions","Global smooth eigenfunction classification","eigenfunctions",
        "Every nonzero global C-infinity eigenfunction at the stable reference is a finite "
        "sum of x_1^p*eta^q with p+3q=n and eigenvalue -n. Any finite smooth extension "
        "admits a finite derivative-closed enlargement. Smoothness at the origin is essential.",
        R"(\psi=\sum_{p+3q=n}c_{pq}x_1^p\eta^q,\quad \lambda=-n)",
        free&&a==-1&&b==-3&&c==1,"Reference a=-1,b=-3,c=1; global C-infinity; zero driver");
    add("quadratic_triangular-infinite-extension","Global entire seed forces infinite extension","eigenfunctions",
        "The polynomial eigenfunction bank remains derivative-closed, but adding exp(I), "
        "I=x_1^3*(x_2-x_1^2/5), forces the independent family x_1^(3n)*exp(I). "
        "Each member is an eigenfunction; no finite derivative closure contains the added seed.",
        R"(I=x_1^3(x_2-x_1^2/5),\quad \partial_2^n e^I=x_1^{3n}e^I)",
        free&&a==1&&b==-3&&c==1,"Saddle a=1,b=-3,c=1; global entire functions; zero driver");
    return r;
}

std::vector<ResearchPaper> triangularResearchPapers() {
    ResearchPaper b;
    b.id="brunton-etal-2016-finite-koopman";
    b.title="Koopman Invariant Subspaces and Finite Linear Representations of Nonlinear Dynamical Systems for Control";
    b.authors={"Steven L. Brunton","Bingni W. Brunton","Joshua L. Proctor","J. Nathan Kutz"};
    b.year=2016; b.venue="PLOS ONE 11(2), e0150171"; b.doi="10.1371/journal.pone.0150171";
    b.urls={brunton}; b.topics={"finite Koopman embedding","polynomial eigenfunctions","control"};
    b.verification="Primary full text, equations 21-27 inspected";
    b.applicability="Native a=mu, b=lambda, c=-lambda gives the paper's quadratic form; "
                    "the default case is equivalent via X_1=x_1, X_2=3*x_2, mu=-1, lambda=-3";
    b.relatedFactIds={"quadratic_triangular-eigenfunctions","quadratic_triangular-slow-manifold"};
    ResearchPaper k=b;
    k.id="kaiser-kutz-brunton-2021-kronic";
    k.title="Data-driven discovery of Koopman eigenfunctions for control";
    k.authors={"Eurika Kaiser","J. Nathan Kutz","Steven L. Brunton"};
    k.year=2021; k.venue="arXiv:1707.01146v4"; k.doi="10.48550/arXiv.1707.01146";
    k.urls={kaiser}; k.verification="Primary full text, Section 6, equations 50-52 inspected";
    return {b,k};
}
} // namespace dynsyslab::core::detail
