#include <chrono>
#include <cmath>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/KromAuditQuery.hpp>
#include <dynsyslab/workflow_c_api.h>
#include <filesystem>
#include <iostream>
#include <stdexcept>
namespace c = dynsyslab::core;
namespace j = c::json;
void check(bool v, const char *message) {
  if (!v)
    throw std::runtime_error(message);
}
j::Value abi(const j::Value &q) {
  dsl_workflow_result_handle *result = nullptr;
  check(dsl_workflow_krom_audit_json(j::dump(q).c_str(), &result) ==
            DSL_STATUS_OK,
        "C ABI audit failed");
  std::size_t n = 0;
  dsl_workflow_result_json(result, nullptr, 0, &n);
  std::string s(n, '\0');
  const auto status = dsl_workflow_result_json(result, s.data(), n, &n);
  dsl_workflow_result_destroy(result);
  check(status == DSL_STATUS_OK, "Result read failed");
  s.resize(n - 1);
  return j::parse(s);
}
int main() {
  try {
    auto root =
        std::filesystem::temp_directory_path() /
        ("dynsyslab-krom-audit-" +
         std::to_string(
             std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root);
    c::SystemConfiguration config;
    config.identifier = "quadratic_triangular";
    config.parameters = {{"a", -1}, {"b", -2}, {"c", 0}};
    auto model = c::KoopmanModel::fromPolynomialCoefficients(
        config, 2, {1, 0, 0, 1}, {{1, 0}, {0, 0}, {0, 0}, {1, 0}},
        {{-1, 0}, {-2, 0}});
    model.setReconstruction("state", {{1., 0., 0., 1.}, 2, 2});
    c::saveKoopmanModel(root / "model.npz", model);
    const auto system = c::createSystem(config);
    j::Value::Object parameters;
    for (auto [key, value] : system->configuration().parameters)
      parameters[key] = value;
    c::NpzRecord record;
    std::vector<double> time, x;
    for (int k = 0; k <= 10; ++k) {
      double t = k * .01;
      time.push_back(t);
      x.insert(x.end(), {std::exp(-t), 2 * std::exp(-2 * t)});
    }
    record.arrays["time"] = c::NpyArray::fromDoubles({11}, time);
    record.arrays["values"] = c::NpyArray::fromDoubles({11, 2}, x);
    record.arrays["metadata_json"] = c::NpyArray::fromString(
        j::dump(j::Value::Object{{"system", "quadratic_triangular"},
                                 {"system_parameters", parameters},
                                 {"coordinate_system", "dimensionless"}}));
    c::writeNpzAtomic(root / "source.npz", record);
    j::Value q =
        j::Value::Object{{"schema", "dynsyslab.krom-audit-query/v1"},
                         {"action", "trajectory"},
                         {"krom", (root / "model.npz").generic_u8string()},
                         {"source", (root / "source.npz").generic_u8string()},
                         {"input_map", "none"}};
    const auto direct = dynsyslab::workflow::evaluateKromAuditQuery(q),
               via = abi(q);
    check(j::dump(direct) == j::dump(via), "C++/C parity");
    check(
        direct.at("trajectory").at("step_defect_relative_to_change").number(1) <
            1e-11,
        "Exact linear replay defect");
    const auto &replay = direct.at("trajectory").at("replay_state").array();
    for (std::size_t i = 0; i < x.size(); ++i)
      check(std::abs(replay[i].number() - x[i]) < 1e-12,
            "Free replay mismatch");
    q.at("action") = "derivative_closure";
    q.at("extend") = false;
    check(abi(q).at("closure").at("status").string() == "closed",
          "Linear bank derivative closure");
    dsl_workflow_result_handle *out =
        reinterpret_cast<dsl_workflow_result_handle *>(1);
    check(dsl_workflow_krom_audit_json(nullptr, &out) ==
                  DSL_STATUS_INVALID_ARGUMENT &&
              !out,
          "Null request");
    check(dsl_workflow_krom_audit_json("{}", nullptr) ==
              DSL_STATUS_INVALID_ARGUMENT,
          "Null output");
    check(dsl_workflow_krom_audit_json("{}", &out) ==
                  DSL_STATUS_INVALID_ARGUMENT &&
              !out,
          "Bad schema");
    std::filesystem::remove_all(root);
    std::cout << "krom_audit_query_tests: PASS\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
