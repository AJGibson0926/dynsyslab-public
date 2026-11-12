#include <dynsyslab/workflow/CoordinateClosureQuery.hpp>

#include <iostream>
#include <stdexcept>

namespace {
void check(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}
}

int main() {
    try {
        const auto query=dynsyslab::core::json::parse(R"({
          "schema":"dynsyslab.observable-closure-query/v1",
          "system":{"name":"rpe","parameters":{"kappa":1.3333333333333333}},
          "observable":"echo",
          "domain":{"minimum":[0.5,-10],"maximum":[1.5,10]},
          "highest_derivative_order":2,
          "sample_count":16,
          "pair_candidate_count":2,
          "maximum_refinement_iterations":20,
          "collision_tolerance":1e-10,
          "successor_difference_tolerance":1e-8,
          "pair_seeds":[{"left":[0.7500373016335908,5.627534424684314],
            "right":[1.1887083801676752,7.575355314915713]}]
        })");
        const auto result=dynsyslab::workflow::evaluateCoordinateClosureQuery(query);
        check(result.at("schema").string()=="dynsyslab.observable-closure-result/v1",
            "closure result schema");
        check(result.at("analysis_scope").string()=="bounded_domain_numerical",
            "closure analysis scope");
        check(!result.at("search_is_proof").boolean(true)
            &&!result.at("global_guarantee").boolean(true),
            "numerical closure search was mislabeled as proof");
        const auto& order=result.at("orders").array()[2];
        check(order.at("status").string()=="closure_obstruction_found",
            "known order-two closure obstruction status");
        check(order.at("closure_obstruction_witness")
                .at("normalized_successor_difference").number()>1.0e-8,
            "order-two successor derivative was not separated");
        const auto global=dynsyslab::core::json::parse(R"({
          "schema":"dynsyslab.observable-closure-query/v2",
          "system":{"equations":[{"encoding":"analytic-state-expression/v1","operation":"state",
            "state_dimension":1,"state_index":0}]},"coordinate_axes":[0]})");
        const auto exact=dynsyslab::workflow::evaluateCoordinateClosureQuery(global);
        check(exact.at("status").string()=="certified_equations" &&
            exact.at("certified_order").number(-1)==0,"global custom-equation query");
        check(exact.at("orders").array()[0].at("rhs_observable").at("expressions").array().size()==1,
            "reusable RHS expression missing");
        std::cout<<"observable closure query contract passed\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
