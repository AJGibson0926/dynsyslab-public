#include <dynsyslab/workflow/CoordinateClosureQuery.hpp>
#include <dynsyslab/core/CoordinateClosure.hpp>
#include <dynsyslab/core/GlobalJetClosure.hpp>
#include <dynsyslab/core/System.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {

json::Value numbers(const std::vector<double>& values) {
    json::Value::Array result;
    for(const auto value:values)result.emplace_back(value);
    return result;
}

std::vector<double> parseNumbers(const json::Value& value) {
    std::vector<double> result;
    for(const auto& item:value.array()) {
        const double number=item.number(std::numeric_limits<double>::quiet_NaN());
        if(!std::isfinite(number))throw std::invalid_argument("closure query arrays must contain finite numbers");
        result.push_back(number);
    }
    return result;
}

std::size_t count(const json::Value& value,const char* name,std::size_t fallback) {
    if(value.isNull())return fallback;
    const double number=value.number(-1);
    if(!std::isfinite(number)||number<0||number!=std::floor(number)
        ||number>=double(std::numeric_limits<std::size_t>::max()))
        throw std::invalid_argument(std::string(name)+" must be a nonnegative integer");
    return static_cast<std::size_t>(number);
}

json::Value witnessJson(const dynsyslab::core::ObservableCollisionWitness& witness) {
    if(!witness.found)return json::Value{};
    json::Value::Array mismatch;
    for(std::size_t index=0;index<witness.successorDerivativeLeft.size();++index)
        mismatch.emplace_back(witness.successorDerivativeLeft[index]-witness.successorDerivativeRight[index]);
    return json::Value::Object{
        {"left_state",numbers(witness.leftState)},
        {"right_state",numbers(witness.rightState)},
        {"left_coordinates",numbers(witness.leftCoordinates)},
        {"right_coordinates",numbers(witness.rightCoordinates)},
        {"successor_derivative_left",numbers(witness.successorDerivativeLeft)},
        {"successor_derivative_right",numbers(witness.successorDerivativeRight)},
        {"successor_derivative_difference",std::move(mismatch)},
        {"normalized_residual",witness.normalizedResidual},
        {"normalized_successor_difference",witness.normalizedSuccessorDifference},
        {"normalized_state_separation",witness.normalizedStateSeparation},
        {"refinement_iterations",double(witness.refinementIterations)}};
}

}  // namespace

json::Value evaluateCoordinateClosureQuery(const json::Value& request) {
    const auto schema=request.at("schema").string();
    const bool global=schema=="dynsyslab.observable-closure-query/v2";
    if(!global && schema!="dynsyslab.observable-closure-query/v1")
        throw std::invalid_argument("unsupported observable closure query schema");
    const auto& specification=request.at("system");
    if(global) {
        const auto sources=unsigned(!request.at("expression").isNull())+unsigned(!request.at("coordinate_axes").isNull())
            +unsigned(!request.at("coordinates").isNull())+unsigned(!request.at("observable").isNull());
        if(sources!=1) throw std::invalid_argument("v2 requires exactly one measured coordinate specification");
        for(const auto* key:{"positive_axes","coordinate_axes","coordinates","pair_seeds","evaluate_at_jet"}) {
            const auto& value=request.at(key);
            if(!value.isNull() && !std::holds_alternative<json::Value::Array>(value.storage()))
                throw std::invalid_argument(std::string(key)+" must be an array");
        }
    }
    if(global && (!request.at("domain").isNull() || !specification.at("driver").isNull()
        || !request.at("input").isNull() || request.at("include_input_terms").boolean()))
        throw std::invalid_argument("global queries use autonomous zero-input equations on their regular domain");
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier=specification.at("name").string(specification.at("identifier").string());
    configuration.preset=specification.at("preset").string();
    for(const auto& item:specification.at("parameters").object())
        configuration.parameters[item.first]=item.second.number();
    const auto& inlineField=specification.at("equations");
    std::unique_ptr<dynsyslab::core::ODESystem> system;
    std::vector<dynsyslab::core::AnalyticExpression> field;
    if(inlineField.isNull()) {
        system=dynsyslab::core::createSystem(configuration);
        if(global) field=system->driftExpressions();
    } else {
        if(!global || !configuration.identifier.empty() || !configuration.parameters.empty()
            || !configuration.preset.empty() || inlineField.array().empty())
            throw std::invalid_argument("supply either a registered system name or a nonempty v2 equations array");
        for(const auto& item:inlineField.array())
            field.push_back(dynsyslab::core::AnalyticExpression::fromJson(item));
    }
    const auto name=request.at("observable").string();
    const auto dimension=system?system->descriptor().stateDimension:field.size();
    auto observable=[&]() -> dynsyslab::core::Observable {
        if(!request.at("expression").isNull())
            return dynsyslab::core::Observable::fromExpressionDefinition(request.at("expression"));
        if(global && !request.at("coordinate_axes").isNull()) {
            std::vector<dynsyslab::core::Observable> components;
            for(const auto& axis:request.at("coordinate_axes").array()) {
                auto i=count(axis,"coordinate axis",dimension);
                components.push_back(dynsyslab::core::Observable::stateComponent(dimension,i,
                    "x"+std::to_string(i),"state component"));
            }
            return dynsyslab::core::Observable::stack("coordinates","Coordinates",components);
        }
        if(!system) throw std::invalid_argument("custom equations require expression or coordinate_axes");
        const auto catalog=system->observableCatalog();
        if(global && !request.at("coordinates").isNull()) {
            std::vector<std::string> names;
            for(const auto& item:request.at("coordinates").array()) names.push_back(item.string());
            return catalog.coordinates(names);
        }
        return catalog.coordinatePresets().count(name)?catalog.coordinatePreset(name):catalog.at(name);
    }();
    if(global) {
        dynsyslab::core::GlobalJetOptions exact;
        exact.highestOrder=count(request.at("highest_derivative_order"),"highest_derivative_order",4);
        exact.polynomialDegree=count(request.at("polynomial_degree"),"polynomial_degree",2);
        exact.emitSolverProblems=request.at("emit_smtlib").boolean();
        exact.searchAllOrders=request.at("search_all_orders").boolean();
        exact.solverTimeoutMs=count(request.at("solver_timeout_ms"),"solver_timeout_ms",0);
        for(const auto& axis:request.at("positive_axes").array())
            exact.positiveAxes.push_back(count(axis,"positive axis",dimension));
        for(const auto& pair:request.at("pair_seeds").array())
            exact.pairSeeds.emplace_back(parseNumbers(pair.at("left")),parseNumbers(pair.at("right")));
        auto result=dynsyslab::core::globalJetClosure(field,observable,exact);
        result.at("system")=system?system->descriptor().identifier:"supplied_equations";
        result.at("observable")=observable.descriptor().name;
        if(!request.at("evaluate_at_jet").isNull()) {
            const auto q=parseNumbers(request.at("evaluate_at_jet"));
            for(const auto& order:result.at("orders").array()) {
                if(order.at("status").string()=="certified_closed" && order.at("coordinate_count").number()==double(q.size())) {
                    result.at("evaluated_rhs")=numbers(dynsyslab::core::evaluateJetRhs(order,q));
                    result.at("evaluated_order")=order.at("order");
                    break;
                }
            }
            if(result.at("evaluated_rhs").isNull())
                throw std::invalid_argument("evaluate_at_jet requires explicit certified equations matching the jet dimension");
        }
        return result;
    }
    dynsyslab::core::ObservableClosureOptions options;
    const auto& domain=request.at("domain");
    options.stateMinimum=parseNumbers(domain.at("minimum").isNull()?domain.at("lower"):domain.at("minimum"));
    options.stateMaximum=parseNumbers(domain.at("maximum").isNull()?domain.at("upper"):domain.at("maximum"));
    options.highestDerivativeOrder=count(request.at("highest_derivative_order"),
        "highest_derivative_order",4);
    options.sampleCount=count(request.at("sample_count"),"sample_count",512);
    options.pairCandidateCount=count(request.at("pair_candidate_count"),"pair_candidate_count",64);
    options.maximumRefinementIterations=count(request.at("maximum_refinement_iterations"),
        "maximum_refinement_iterations",80);
    options.seed=static_cast<std::uint64_t>(count(request.at("seed"),"seed",1));
    options.minimumNormalizedStateSeparation=request.at("minimum_normalized_state_separation").number(0.05);
    options.collisionTolerance=request.at("collision_tolerance").number(1.0e-9);
    options.successorDifferenceTolerance=
        request.at("successor_difference_tolerance").number(1.0e-8);
    options.rankTolerance=request.at("rank_tolerance").number(1.0e-10);
    for(const auto& pair:request.at("pair_seeds").array()) {
        const auto left=parseNumbers(pair.at("left"));
        const auto right=parseNumbers(pair.at("right"));
        if(left.size()!=dimension||right.size()!=dimension)
            throw std::invalid_argument("closure pair seed state dimension mismatch");
        options.pairSeeds.insert(options.pairSeeds.end(),left.begin(),left.end());
        options.pairSeeds.insert(options.pairSeeds.end(),right.begin(),right.end());
    }
    const auto report=dynsyslab::core::analyzeObservableClosure(*system,observable,options);
    json::Value::Array orders;
    std::optional<std::size_t> firstCandidate;
    for(const auto& order:report.orders) {
        const bool fullRank=order.minimumSampledRank==report.stateDimension;
        if(!order.coordinateCollision.found&&fullRank&&!firstCandidate)
            firstCandidate=order.highestDerivativeOrder;
        const char* status=order.closureObstruction.found?"closure_obstruction_found":
            (order.coordinateCollision.found?"coordinate_collision_found_no_obstruction_in_search":
                (fullRank?"no_coordinate_collision_found_in_search":
                    "rank_deficient_no_collision_found_in_search"));
        orders.emplace_back(json::Value::Object{
            {"highest_derivative_order",double(order.highestDerivativeOrder)},
            {"coordinate_count",double(order.coordinateCount)},
            {"minimum_sampled_rank",double(order.minimumSampledRank)},
            {"minimum_sampled_singular_value",order.minimumSampledSingularValue},
            {"minimum_rank_state",numbers(order.minimumRankState)},
            {"locally_full_rank_on_samples",fullRank},
            {"status",status},
            {"coordinate_collision_witness",witnessJson(order.coordinateCollision)},
            {"closure_obstruction_witness",witnessJson(order.closureObstruction)}});
    }
    json::Value::Object result{
        {"schema","dynsyslab.observable-closure-result/v1"},
        {"system",system->descriptor().identifier},
        {"observable",observable.descriptor().name},
        {"state_dimension",double(report.stateDimension)},
        {"observable_dimension",double(report.observableDimension)},
        {"requested_samples",double(report.requestedSamples)},
        {"accepted_samples",double(report.acceptedSamples)},
        {"domain",json::Value::Object{{"minimum",numbers(options.stateMinimum)},
            {"maximum",numbers(options.stateMaximum)},
            {"coordinate_system",domain.at("coordinate_system").string("native")}}},
        {"orders",std::move(orders)},
        {"analysis_scope","bounded_domain_numerical"},
        {"search_is_proof",false},
        {"global_guarantee",false},
        {"proof_certificate",json::Value{}},
        {"interpretation","A numerical closure obstruction is a pair of distinct states with matching observable jets within tolerance but a separated successor derivative. It is strong counterexample evidence, not an exact certificate. A coordinate collision with the same successor only challenges injectivity. Failure to find either in a finite bounded search is inconclusive globally."}};
    result["first_full_rank_no_collision_found_order"]=
        firstCandidate?json::Value(double(*firstCandidate)):json::Value{};
    return result;
}

}  // namespace dynsyslab::workflow
