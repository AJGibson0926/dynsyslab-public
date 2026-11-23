#include <dynsyslab/core/GlobalJetClosure.hpp>
#include <dynsyslab/core/System.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

using namespace dynsyslab::core;
namespace {
void require(bool ok,const char* message) {
    if (!ok) throw std::runtime_error(message);
}
Observable map(std::vector<AnalyticExpression> terms) {
    return Observable::fromExpressions({"test","Test map","","native",{},false},std::move(terms));
}
Observable rhs(const json::Value& result) {
    return Observable::fromExpressionDefinition(result.at("orders").array().back().at("rhs_observable"));
}
}

int main() {
    try {
        using E=AnalyticExpression;
        const auto x=E::state(2,0), y=E::state(2,1);
        const auto one=E::constant(1,2);
        const double q[]{2,3};
        auto result=globalJetClosure({y,-x-y/2},map({x}));
        require(result.at("minimality_proved").boolean(),"oscillator minimal order proof");
        require(rhs(result).evaluateReal(q)==std::vector<double>({3,-3.5}),"companion dynamics");
        result=globalJetClosure({y,-x},map({x.power(2)+y.power(2)}));
        require(result.at("certified_order").number(-1)==0 && rhs(result).evaluateReal(q)[0]==0,
            "closed conserved quantity");
        require(!result.at("state_recovery_obstruction").isNull(),"all-order sign symmetry");
        result=globalJetClosure({x,y},map({x+y,x-y}));
        require(result.at("orders").array()[0].at("state_recovery").string()=="certified_inverse",
            "invertible mixed coordinate chart");
        require(rhs(result).evaluateReal(q)==std::vector<double>({2,3}),"mixed vector dynamics");

        result=globalJetClosure({y,-x.unary("sin")},map({x}));
        require(result.at("certified_order").number(-1)==1,"analytic position jet");
        require(std::abs(rhs(result).evaluateReal(q)[1]+std::sin(2.0))<1e-14,"analytic rate renaming");
        result=globalJetClosure({one,y},map({x.unary("sin")}));
        require(result.at("certified_order").number(-1)==1 && rhs(result).evaluateReal(q)==std::vector<double>({3,-2}),
            "sine measurement oscillator equation by polynomial lift");
        require(!result.at("minimality_proved").boolean(),"rank on a lifted relaxation used as source proof");
        require(result.at("state_recovery_obstruction").at("method").string()=="closed_subsystem_dependencies",
            "unmeasured state cannot affect the sine history");
        result=globalJetClosure({one,y},map({x.unary("exp")*x.unary("sin")}));
        require(result.at("certified_order").number(-1)==1 && rhs(result).evaluateReal(q)==std::vector<double>({3,2}),
            "mixed exponential/trigonometric sensor equation");
        result=globalJetClosure({one,y},map({x.unary("sech"),x.unary("tanh")}),GlobalJetOptions{0,2});
        require(result.at("certified_order").number(-1)==0,"sech/tanh closed coordinates");
        const double hyperbolic[]{1/std::cosh(.7),std::tanh(.7)};
        const auto hyperbolicRate=rhs(result).evaluateReal(hyperbolic);
        require(std::abs(hyperbolicRate[0]+hyperbolic[0]*hyperbolic[1])<1e-13
            &&std::abs(hyperbolicRate[1]-(1-hyperbolic[1]*hyperbolic[1]))<1e-13,
            "sech polynomial lift dynamics");
        require(!result.at("polynomial_lift").at("real_algebraic_graph_exact").boolean(),
            "transcendental graph must remain a relaxation");
        result=globalJetClosure({one,y},map({x.power(0.5)}));
        require(result.at("certified_order").isNull(),"unknown branch must not invent an order");
        GlobalJetOptions positive;
        positive.positiveAxes={0};
        result=globalJetClosure({x,y},map({x.power(0.5)}),positive);
        require(result.at("polynomial_lift").at("real_algebraic_graph_exact").boolean() &&
            result.at("finite_order_existence").string()=="guaranteed_on_regular_domain","positive square-root lift");
        require(result.at("certified_order").number(-1)==0 && rhs(result).evaluateReal(q)[0]==1.,
            "root sensor dynamics after reduction on its graph");

        const auto t=E::state(1,0);
        result=globalJetClosure({E::constant(1,1)/(t.power(2)+1)},map({t.power(2)}),GlobalJetOptions{0,1});
        require(result.at("guaranteed_order_dynamics").at("order").string()=="72",
            "general rational degree bound and implicit dynamics");
        result=globalJetClosure({E::constant(1,1)},map({t.power(3)}));
        require(result.at("certified_order").number(-1)==0,"odd-root closure");
        require(std::abs(evaluateJetRhs(result.at("orders").array()[0],{-8})[0]-12.)<1e-12,"real cube-root evaluation");
        for (double value:{1e-200,1e200}) {
            const double root=std::cbrt(value);
            const double actual=evaluateJetRhs(result.at("orders").array()[0],{value})[0];
            require(std::abs(actual/(3*root*root)-1)<1e-13,"representable root of extreme radicand");
        }
        GlobalJetOptions all{2,2}; all.searchAllOrders=true;
        result=globalJetClosure({E::constant(1,1)},map({t.power(3)}),all);
        require(result.at("certified_order").number(-1)==0 && result.at("first_smooth_recovering_order").number(-1)==2,
            "singular closed law followed by smooth reconstructing jet");
        require(evaluateJetRhs(result.at("orders").array()[2],{-8,12,-12})==std::vector<double>({12,-12,6}),
            "smooth third derivative from the cubic measurement");
        result=globalJetClosure({t.power(2)},map({t.power(3)}),GlobalJetOptions{0,4});
        require(std::abs(evaluateJetRhs(result.at("orders").array()[0],{-8})[0]-48.)<1e-12,"second cube-root law");

        GlobalJetOptions options;
        options.highestOrder=0; options.polynomialDegree=1;
        options.pairSeeds={{{1,1},{-1,1}}}; options.emitSolverProblems=true;
        result=globalJetClosure({one,y},map({x.power(2)}),options);
        const auto& obstruction=result.at("orders").array()[0];
        require(obstruction.at("certificate").at("method").string()=="exact_pair","exact unequal successors");
        require(!obstruction.at("closure_obstruction_smt2").string().empty(),"solver problem missing");
        options.positiveAxes={0};
        result=globalJetClosure({one,y},map({x.power(2)}),options);
        require(result.at("orders").array()[0].at("certificate").isNull(),"seed outside positive domain");

        SystemConfiguration config; config.identifier="rpe"; config.parameters["kappa"]=4./3.;
        auto system=createSystem(config);
        options.highestOrder=3; options.emitSolverProblems=false;
        options.pairSeeds={{{0.9457041294916282,14.559306001295163},{0.7810263280328975,-14.828898931901062}},
            {{0.9087071031,-5.8121728919},{1.1002768956,-6.4529294201}},
            {{0.7500373016335908,5.627534424684314},{1.1887083801676752,7.575355314915713}}};
        result=globalJetClosure(system->driftExpressions(),system->observableCatalog().at("echo"),options);
        require(result.at("order_lower_bound").number()==3 && result.at("certified_order").isNull(),
            "RPE lower and upper evidence conflated");
        require(result.at("orders").array()[2].at("certificate").at("method").string()=="rational_interval_contraction",
            "RPE order-two interval validation");
        require(result.at("finite_order_bound").at("sufficient_order").string()=="2672672","RPE contact bound");
        require(result.at("guaranteed_order_dynamics").at("order").string()=="2672672","implicit companion law");
        require(result.at("component_state_recovery").at("finite_injective_order_exists").boolean(),"RPE recovery theorem");
        config.parameters["mu_pa_s"]=0;
        system=createSystem(config);
        result=globalJetClosure(system->driftExpressions(),system->observableCatalog().at("echo"),
            GlobalJetOptions{0,1});
        require(result.at("component_state_recovery").isNull(),"inviscid model misclassified as dissipative");
        if (std::getenv("DYNSYSLAB_Z3_LIBRARY")) {
            GlobalJetOptions decide{0,1};
            decide.solverTimeoutMs=5000;
            const auto translation=E::constant(1,1);
            result=globalJetClosure({translation},map({t.power(2)}),decide);
            require(result.at("orders").array()[0].at("closure_solver").at("status").string()=="sat",
                "solver failed to detect the two branches");
            result=globalJetClosure({translation},map({t.power(3)+t}),decide);
            require(result.at("status").string()=="solver_closed_implicit_graph" &&
                result.at("orders").array()[0].at("state_recovery").string()=="solver_injective",
                "solver failed to certify the strictly monotone cubic");
            result=globalJetClosure({translation/(t.power(2)+1)},map({t.power(2)}),decide);
            require(result.at("order_lower_bound").number()==1,"rational lifted obstruction");
            result=globalJetClosure({translation/(t.power(2)+1)},map({t.power(3)+t}),decide);
            require(result.at("status").string()=="solver_closed_implicit_graph","rational lifted global graph");
            result=globalJetClosure({translation},map({t.unary("exp")}),decide);
            require(result.at("orders").array()[0].at("state_recovery").string()=="undetermined",
                "a collision in a relaxed exponential graph cannot refute injectivity");
        }
        std::cout<<"Global jet proofs and executable equations passed\n";
        return 0;
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n'; return 1;}
}
