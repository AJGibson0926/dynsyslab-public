#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/LieAlgebraQuery.hpp>
#include <dynsyslab/workflow_c_api.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
void check(bool ok, const char *why) {
  if (!ok)
    throw std::runtime_error(why);
}
std::string invoke(const std::string &request) {
  dsl_workflow_result_handle *result = nullptr;
  check(dsl_workflow_lie_algebra_json(request.c_str(), &result) ==
            DSL_STATUS_OK,
        "C ABI query failed");
  size_t required = 0;
  auto status = dsl_workflow_result_json(result, nullptr, 0, &required);
  if (status != DSL_STATUS_BUFFER_TOO_SMALL) {
    dsl_workflow_result_destroy(result);
    throw std::runtime_error("C ABI buffer sizing failed");
  }
  std::string text(required, '\0');
  status =
      dsl_workflow_result_json(result, text.data(), text.size(), &required);
  dsl_workflow_result_destroy(result);
  check(status == DSL_STATUS_OK, "C ABI result read failed");
  text.resize(required - 1);
  return text;
}
} // namespace
int main() {
  try {
    namespace json = dynsyslab::core::json;
    size_t count = 0;
    for (const auto &entry :
         std::filesystem::directory_iterator(DYNSYSLAB_LIE_EXAMPLES)) {
      if (entry.path().extension() != ".json")
        continue;
      std::ifstream stream(entry.path());
      std::ostringstream text;
      text << stream.rdbuf();
      auto request = json::parse(text.str());
      const auto direct = dynsyslab::workflow::evaluateLieAlgebraQuery(request);
      const auto viaAbi = json::parse(invoke(text.str()));
      check(json::dump(direct) == json::dump(viaAbi), "C++/C ABI mismatch");
      check(direct.at("schema").string() == "dynsyslab.lie-algebra-result/v1",
            "Schema mismatch");
      ++count;
    }
    check(count >= 18, "Case studies missing");
    dsl_workflow_result_handle *output =
        reinterpret_cast<dsl_workflow_result_handle *>(1);
    check(dsl_workflow_lie_algebra_json(nullptr, &output) ==
                  DSL_STATUS_INVALID_ARGUMENT &&
              !output,
          "Null input contract");
    check(dsl_workflow_lie_algebra_json("{}", nullptr) ==
              DSL_STATUS_INVALID_ARGUMENT,
          "Null output contract");
    check(dsl_workflow_lie_algebra_json("{}", &output) ==
                  DSL_STATUS_INVALID_ARGUMENT &&
              !output,
          "Bad schema contract");
    std::cout << count << " Lie queries matched through C++ and C ABI\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
