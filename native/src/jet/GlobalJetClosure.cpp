#include <dynsyslab/core/GlobalJetClosure.hpp>
#include "JetEquationBasis.hpp"
#include "JetFiberCertificate.hpp"
#include "JetGlobalBounds.hpp"
#include "JetRationalBound.hpp"
#include "JetPolynomialLift.hpp"
#include "JetRank.hpp"
#include "JetDecision.hpp"
#include "JetAnalyticChart.hpp"

namespace dynsyslab::core {
namespace {
using namespace jet;
using List = json::Value::Array;
using Record = json::Value::Object;
using Polys = std::vector<Polynomial>;

void recordOrder(Record& report,json::Value order,std::size_t lower) {
    const auto& facts=std::as_const(order);
    const auto k=facts.at("order");
    const auto status=facts.at("status").string();
    const bool explicitLaw=status=="certified_closed";
    const auto recovery=facts.at("state_recovery").string();
    const bool injective=recovery=="certified_inverse" || recovery=="solver_injective";
    const bool smooth=facts.at("ode_regularity").string()=="smooth near each realizable jet";
    if (injective && !report.count("first_tested_injective_order")) report["first_tested_injective_order"]=k;
    if (smooth && !report.count("first_smooth_rhs_order")) report["first_smooth_rhs_order"]=k;
    if (smooth && injective && !report.count("first_smooth_recovering_order")) report["first_smooth_recovering_order"]=k;
    if ((explicitLaw || status=="solver_closed_implicit_graph") && report.at("certified_order").isNull()) {
        report["certified_order"]=k;
        report["minimality_proved"]=double(lower)==k.number();
        report["finite_order_existence"]="guaranteed_on_regular_domain";
        report["status"]=explicitLaw?"certified_equations":"solver_closed_implicit_graph";
    }
    report["orders"].array().push_back(std::move(order));
}

bool rationalTree(const json::Value& node) {
    const auto op=node.at("operation").string();
    if (op=="power" && std::floor(node.at("scalar").number())!=node.at("scalar").number()) return false;
    if (op!="state" && op!="constant" && op!="negative" && op!="add" && op!="subtract"
        && op!="multiply" && op!="divide" && op!="power") return false;
    for (const auto& arg:node.at("arguments").array()) if (!rationalTree(arg)) return false;
    return true;
}

List encodePolys(const Polys& rows) {
    List result;
    for (const auto& row:rows) result.emplace_back(encode(row));
    return result;
}

json::Value runtimeEquations(const Polys& equations) {
    const auto n=equations.size();
    List expressions;
    for (const auto& p:equations) {
        auto sum=AnalyticExpression::constant(0,n);
        for (const auto& [powers,c]:p) {
            const auto rounded=c.convert_to<double>();
            if (!std::isfinite(rounded) || (rounded==0 && c!=0)) return {};
            auto term=AnalyticExpression::constant(rounded,n);
            for (std::size_t i=0;i<n;++i) if (powers[i])
                term=term*AnalyticExpression::state(n,i).power(powers[i]);
            sum=sum+term;
        }
        expressions.emplace_back(sum.toJson());
    }
    return Record{{"schema","dynsyslab.observable-expression/v1"},{"name","jet_rhs"},
        {"expressions",expressions}};
}

std::optional<Polys> expressAll(const EquationBasis& basis,const Polys& targets) {
    Polys results;
    for (const auto& target:targets) {
        auto p=basis.express(target);
        if (!p) return {};
        results.push_back(std::move(*p));
    }
    return results;
}

// In particular r' = v turns the radius jet into a literal state chart.
// This path retains negative powers, which a polynomial dictionary cannot span.
std::optional<Polys> inStateChart(const Polys& jet,const Polys& states,const Polys& targets) {
    std::vector<std::size_t> slots;
    for (const auto& x:states) {
        const auto found=std::find(jet.begin(),jet.end(),x);
        if (found==jet.end()) return {};
        slots.push_back(std::size_t(found-jet.begin()));
    }
    Polys output;
    for (const auto& p:targets) {
        Polynomial result;
        for (const auto& [powers,c]:p) {
            Powers q(jet.size());
            for(std::size_t j=0;j<slots.size();++j) q[slots[j]]=powers[j];
            addTerm(result,q,c);
        }
        output.push_back(std::move(result));
    }
    return output;
}
bool equalAt(const Polys& polys,const Box& a,const Box& b) {
    for (const auto& p:polys) if(evaluate(p,a).low!=evaluate(p,b).low) return false;
    return true;
}

struct Analysis {
    Polys f, jet, level, states;
    std::set<std::size_t> poles;
    const GlobalJetOptions& options;
    std::size_t width;
    const PolynomialLift* lift{};

    json::Value checkOrder(std::size_t k) {
        Polys next;
        const auto* graph=lift?&lift->graph:nullptr;
        for (const auto& p:level) next.push_back(onGraph(lie(p,f),graph));
        jet.insert(jet.end(),level.begin(),level.end());
        Record out{{"order",double(k)},{"coordinate_count",double(jet.size())}};
        std::string closureProblem, inverseProblem;
        if (options.emitSolverProblems || options.solverTimeoutMs) {
            const Polys empty;
            const auto& constraints=lift?lift->graph:empty;
            const auto& positive=lift?lift->positiveAxes:options.positiveAxes;
            closureProblem=fiberProblem(jet,next,f.size(),poles,positive,constraints);
            inverseProblem=fiberProblem(jet,states,f.size(),poles,positive,constraints);
            if (lift && !lift->exactGraph) out["solver_domain"]="algebraic_superset_of_source_graph";
            if (options.emitSolverProblems) {
                out["closure_obstruction_smt2"]=closureProblem;
                out["noninjectivity_smt2"]=inverseProblem;
            }
        }
        std::optional<std::size_t> rank;
        bool rankGrowth=false;
        if (!lift) try {
            rank=genericRank(jet,f.size());
            auto extended=jet;
            extended.insert(extended.end(),next.begin(),next.end());
            const auto nextRank=*rank==f.size()?*rank:genericRank(extended,f.size());
            rankGrowth=nextRank>*rank;
            out["generic_rank"]=double(*rank);
            out["generic_successor_rank"]=double(nextRank);
        } catch (const SearchLimit& e) {out["rank_limit"]=e.what();}
        if (!lift) for (const auto& [a,b]:options.pairSeeds) {
            auto certificate=certifyFiber(jet,next,a,b,poles,options.positiveAxes);
            if (!certificate.isNull()) {
                out["status"]="certified_fiber_obstruction";
                out["certificate"]=std::move(certificate);
                level=std::move(next);
                return out;
            }
        }
        if (rankGrowth) {
            out["status"]="certified_fiber_obstruction";
            out["certificate"]=Record{{"method","generic_rank_increase_on_regular_fibers"}};
            level=std::move(next);
            return out;
        }
        auto law=lift?std::optional<Polys>{}:inStateChart(jet,states,next);
        auto inverse=lift?std::optional<Polys>{}:inStateChart(jet,states,states);
        std::vector<int> rootDegree(width,1);
        if (!law) {
            EquationBasis basis(jet,graph);
            for (std::size_t d=0;d<=options.polynomialDegree;++d) {
                try {basis.addDegree(d,f.size());}
                catch (const SearchLimit& e) {out["library_limit"]=e.what(); break;}
                law=expressAll(basis,next);
                if (law) {
                    inverse=expressAll(basis,states);
                    // A low-degree law may need a higher-degree inverse.
                    for (++d; !inverse && d<=options.polynomialDegree; ++d) {
                        try { basis.addDegree(d,f.size()); inverse=expressAll(basis,states); }
                        catch (const SearchLimit&) {break;}
                    }
                    break;
                }
            }
            // g^3 = P(q) has a unique real solution even at P=0.
            if (!law) {
                Polys radicands;
                for (std::size_t i=0;i<width;++i) {
                    auto value=basis.express(next[i]);
                    if (!value && next[i].size()<=32) {
                        try {value=basis.express(power(next[i],3,f.size()));}
                        catch (const SearchLimit&) {break;}
                        rootDegree[i]=3;
                    }
                    if (!value) break;
                    radicands.push_back(std::move(*value));
                }
                if (radicands.size()==width) {
                    law=std::move(radicands);
                    inverse=expressAll(basis,states);
                }
            }
        }
        if (law) {
            Polys rhs;
            for (std::size_t i=width;i<jet.size();++i) rhs.push_back({{unit(jet.size(),i),1}});
            rhs.insert(rhs.end(),law->begin(),law->end());
            out["status"]="certified_closed";
            out["method"]="exact_laurent_identity";
            out["rhs_in_jet"]=encodePolys(rhs);
            out["rhs_observable"]=runtimeEquations(rhs);
            out["ode_regularity"]="smooth near each realizable jet";
            const bool radical=std::find(rootDegree.begin(),rootDegree.end(),3)!=rootDegree.end();
            if (radical) {
                List degrees(jet.size()-width,json::Value(1.));
                for (int degree:rootDegree) degrees.emplace_back(double(degree));
                out["rhs_root_degrees"]=std::move(degrees);
                out["rhs_radicands_in_jet"]=std::move(out["rhs_in_jet"]);
                out.erase("rhs_in_jet");
                out.erase("rhs_observable");
                out["method"]="exact_real_odd_root_identity";
                out["ode_regularity"]="smooth away from zeros of nontrivial radicands; uniqueness at zeros undetermined";
            }
            out[lift?"jet_in_lifted_state":"jet_in_state"]=encodePolys(jet);
            out["state_recovery"]=inverse?"certified_inverse":"undetermined";
            if (inverse) out["state_in_jet"]=encodePolys(*inverse);
            else if (!radical && rank && *rank<f.size())
                out["state_recovery"]="impossible_on_open_state_space_by_rank_and_closed_law";
            if (!lift && !radical && !inverse) for (const auto& [l,r]:options.pairSeeds) {
                if (l==r) continue;
                const Box a=pointBox(l), b=pointBox(r);
                bool regular=true;
                for (auto i:poles) regular=regular && excludesZero(a[i]) && excludesZero(b[i]);
                for (auto i:options.positiveAxes) regular=regular && l[i]>0 && r[i]>0;
                if (regular && equalAt(jet,a,b)) {
                    out["state_recovery"]="indistinguishable_at_every_order";
                    out["indistinguishable_pair"]=List{encodeBox(a),encodeBox(b)};
                    break;
                }
            }
        } else out["status"]=out.count("library_limit")?"explicit_library_limit":"no_identity_in_tested_library";
        if (!law && options.solverTimeoutMs) {
            auto decision=decideFiber(closureProblem,options.solverTimeoutMs);
            const auto status=decision.at("status").string();
            if (status=="sat") out["status"]=lift && !lift->exactGraph
                ?"collision_on_algebraic_superset":"solver_fiber_obstruction";
            if (status=="unsat") {
                out["status"]="solver_closed_implicit_graph";
                out[lift?"jet_in_lifted_state":"jet_in_state"]=encodePolys(jet);
                out[lift?"successor_in_lifted_state":"successor_in_state"]=encodePolys(next);
                out["ode_regularity"]="single_valued_on_image; regularity and ODE uniqueness undetermined";
            }
            out["closure_solver"]=std::move(decision);
        }
        if (options.solverTimeoutMs && !inverse &&
            (law || out.at("status").string()=="solver_closed_implicit_graph") &&
            (!out.count("state_recovery") || out.at("state_recovery").string()=="undetermined")) {
            auto decision=decideFiber(inverseProblem,options.solverTimeoutMs);
            const auto status=decision.at("status").string();
            out["state_recovery"]=status=="unsat"?"solver_injective":
                status=="sat" && (!lift || lift->exactGraph)?"solver_noninjective_at_this_order":"undetermined";
            out["injectivity_solver"]=std::move(decision);
        }
        level=std::move(next);
        return out;
    }
};
} // namespace

json::Value globalJetClosure(const std::vector<AnalyticExpression>& drift,
    const Observable& coordinates,const GlobalJetOptions& options) {
    using namespace jet;
    const auto n=coordinates.stateDimension();
    if (options.solverTimeoutMs>300000) throw std::invalid_argument("solver_timeout_ms exceeds 300000");
    List positive;
    for (auto axis:options.positiveAxes) {
        if (axis>=n) throw std::invalid_argument("positive axis out of range");
        positive.emplace_back(double(axis));
    }
    if (options.highestOrder>12 || options.polynomialDegree>6 || options.pairSeeds.size()>32)
        throw std::invalid_argument("global search supports order 0..12, degree 0..6 and at most 32 seed pairs");
    for (const auto& [a,b]:options.pairSeeds) {
        if (a.size()!=n || b.size()!=n) throw std::invalid_argument("global seed dimension mismatch");
        for (const auto* state:{&a,&b}) for (double x:*state)
            if (!std::isfinite(x)) throw std::invalid_argument("nonfinite global seed");
    }
    Record report{{"schema","dynsyslab.global-jet-closure/v1"},
        {"domain","regular real domain of the autonomous expressions intersected with positive_axes"},
        {"positive_axes",std::move(positive)},
        {"coefficient_convention","exact rational value of each stored binary64 constant"},
        {"finite_order_existence","undetermined"},{"status","undetermined"},
        {"certified_order",json::Value{}},{"orders",List{}},
        {"search_all_orders",options.searchAllOrders},
        {"max_order",double(options.highestOrder)},{"max_polynomial_degree",double(options.polynomialDegree)}};
    if (drift.empty() || coordinates.expressions().empty()) {
        report["reason"]="analytic expressions unavailable for drift or coordinate map";
        return report;
    }
    if (drift.size()!=n) throw std::invalid_argument("global drift dimension mismatch");
    std::vector<json::Value> field,map;
    for (const auto& p:drift) {
        if (p.stateDimension()!=n) throw std::invalid_argument("global drift axis dimension mismatch");
        field.push_back(p.toJson());
    }
    for (const auto& p:coordinates.expressions()) map.push_back(p.toJson());
    auto hidden=hiddenSubsystem(field,map);
    if (!hidden.isNull()) report["state_recovery_obstruction"]=std::move(hidden);
    const bool rational=std::all_of(field.begin(),field.end(),rationalTree) && std::all_of(map.begin(),map.end(),rationalTree);
    if (rational) {
        report["finite_order_existence"]="guaranteed_on_regular_domain";
        report["existence_theorem"]="noetherian_pair_difference_ideals";
    }
    Analysis analysis{{},{},{},{},{},options,map.size()};
    std::optional<PolynomialLift> lift;
    try {
        if (n>32 || map.size()>32) throw SearchLimit("symbolic analysis supports up to 32 source and observable components");
        try {
            for (const auto& node:field) analysis.f.push_back(parseLaurent(node,analysis.poles));
            for (const auto& node:map) analysis.level.push_back(parseLaurent(node,analysis.poles));
        } catch (const Unsupported&) {
            auto chart=derivativeCoordinateChart(field,map,options.highestOrder);
            if (!chart.isNull() && !options.searchAllOrders) {
                recordOrder(report,std::move(chart),0);
                return report;
            }
            if (!chart.isNull()) {
                // This direct source identity may also cover expressions which
                // the lifted algebra cannot represent. Keep it as a fallback.
                report["direct_axis_chart"]=std::move(chart);
            }
            lift.emplace(field,map,options.positiveAxes);
            analysis.lift=&*lift;
            analysis.f=std::move(lift->field);
            analysis.level=std::move(lift->observable);
            analysis.poles.clear();
            report["polynomial_lift"]=Record{{"auxiliaries",std::move(lift->definitions)},
                {"real_algebraic_graph_exact",lift->exactGraph},{"graph_zero_polynomials",encodePolys(lift->graph)}};
            List positiveLift;
            for (auto axis:lift->positiveAxes) positiveLift.emplace_back(double(axis));
            report["polynomial_lift"].at("positive_axes")=std::move(positiveLift);
            report["finite_order_existence"]="guaranteed_on_regular_domain";
            report["existence_theorem"]="noetherian_pair_difference_ideals_on_polynomial_lift";
        }
        List poles;
        for (auto axis:analysis.poles) poles.emplace_back(double(axis));
        report["nonzero_axes"]=std::move(poles);
        report["finite_order_bound"]=finiteOrderBound(analysis.f,analysis.level,analysis.poles);
        if (lift && rational) {
            auto bound=rationalOrderBound(field,map);
            const auto old=report.at("finite_order_bound").at("sufficient_order").string();
            const auto candidate=bound.at("sufficient_order").string();
            if (Degree(candidate)<Degree(old)) report["finite_order_bound"]=std::move(bound);
        }
        if (!lift) {
            auto symmetry=signSymmetry(analysis.f,analysis.level,options.positiveAxes);
            if (!symmetry.isNull() && !report.count("state_recovery_obstruction"))
                report["state_recovery_obstruction"]=std::move(symmetry);
            auto recovery=volumeRecovery(analysis.f,analysis.level,analysis.poles);
            if (!recovery.isNull()) report["component_state_recovery"]=std::move(recovery);
        }
        for (std::size_t i=0;i<n;++i) analysis.states.push_back({{unit(analysis.f.size(),i),1}});
        std::size_t lowerBound=0;
        for (std::size_t k=0;k<=options.highestOrder;++k) {
            auto order=analysis.checkOrder(k);
            const auto status=order.at("status").string();
            if (k==lowerBound && (status=="certified_fiber_obstruction" || status=="solver_fiber_obstruction"))
                ++lowerBound;
            report["order_lower_bound"]=double(lowerBound);
            if (report.count("direct_axis_chart") && report.at("direct_axis_chart").at("order").number()==double(k)) {
                order=std::move(report["direct_axis_chart"]); report.erase("direct_axis_chart");
            }
            recordOrder(report,std::move(order),lowerBound);
            if (!options.searchAllOrders && !report.at("certified_order").isNull()) return report;
        }
        if (report.at("certified_order").isNull())
            report["reason"]="The tested orders and candidate equations did not certify a closure.";
    } catch (const SearchLimit& e) { report["reason"]=e.what(); }
      catch (const Unsupported& e) {report["reason"]=e.what();}
    if (report.count("direct_axis_chart")) {
        recordOrder(report,std::move(report["direct_axis_chart"]),0); report.erase("direct_axis_chart");
    }
    if (report["status"].string()=="undetermined" &&
        report["finite_order_existence"].string()=="guaranteed_on_regular_domain") {
        report["status"]="finite_order_proved_equations_unresolved";
        if (!report.count("finite_order_bound") && rational && n<=32)
            report["finite_order_bound"]=rationalOrderBound(field,map);
    }
    if (report["status"].string()=="finite_order_proved_equations_unresolved" && report.count("finite_order_bound")) {
        report["guaranteed_order_dynamics"]=Record{
            {"representation","implicit_lie_recurrence_on_realizable_jets"},
            {"order",report.at("finite_order_bound").at("sufficient_order")},
            {"source_drift",List(field.begin(),field.end())},{"source_observable",List(map.begin(),map.end())},
            {"lie_operator","D p(x) = grad(p)(x) times f(x), componentwise"},
            {"fiber","choose any regular source state x in the query domain with D^j h(x) = q[j] for j = 0,...,K"},
            {"shift_rhs","dq[j]/dt = q[j+1] for j < K"},
            {"last_rhs","dq[K]/dt = D^(K+1)h(x), independent of the chosen state"}};
    }
    return report;
}
std::vector<double> evaluateJetRhs(const json::Value& order,const std::vector<double>& q) {
    if (order.at("coordinate_count").number(-1)!=double(q.size()) || q.empty())
        throw std::invalid_argument("jet RHS evaluation dimension mismatch");
    for (double x:q) if (!std::isfinite(x)) throw std::invalid_argument("jet must be finite");
    const auto& ordinary=order.at("rhs_observable");
    if (!ordinary.isNull()) {
        const auto observable=Observable::fromExpressionDefinition(ordinary);
        if (observable.stateDimension()!=q.size()) throw std::invalid_argument("RHS expression dimension mismatch");
        return observable.evaluateReal(q.data());
    }
    const auto& degrees=order.at("rhs_root_degrees").array();
    const auto& rows=order.at("rhs_radicands_in_jet").array();
    if (rows.size()!=q.size() || degrees.size()!=q.size())
        throw std::invalid_argument("order has no numerically evaluable explicit equations");
    std::vector<double> result;
    for (std::size_t i=0;i<rows.size();++i) {
        Rational sum=0;
        for (const auto& term:rows[i].array()) {
            Rational value(term.at("coefficient").string());
            const auto& powers=term.at("powers").array();
            if (powers.size()!=q.size()) throw std::invalid_argument("invalid radicand dimension");
            for (std::size_t j=0;j<q.size();++j) {
                const double p=powers[j].number(-1);
                if (p<0 || p>6 || p!=std::floor(p)) throw std::invalid_argument("invalid radicand exponent");
                for (int k=0;k<int(p);++k) value*=Rational(q[j]);
            }
            sum+=value;
        }
        const double degree=degrees[i].number();
        if (degree!=1 && degree!=3) throw std::invalid_argument("unsupported real root");
        double value;
        if (degree==3 && sum!=0) {
            // Scale before rounding. The real root can be representable even
            // when the rational radicand overflows or underflows binary64.
            const int exponent=(int(boost::multiprecision::msb(abs(numerator(sum))))-
                int(boost::multiprecision::msb(denominator(sum))))/3;
            const Degree scale=Degree(1)<<(3*std::abs(exponent));
            if (exponent>=0) sum/=scale; else sum*=scale;
            value=std::ldexp(std::cbrt(sum.convert_to<double>()),exponent);
        } else value=sum.convert_to<double>();
        result.push_back(value);
        if (!std::isfinite(result.back())) throw std::overflow_error("jet RHS overflow");
    }
    return result;
}
} // namespace dynsyslab::core
