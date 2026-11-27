#include <dynsyslab/workflow/JobTemplates.hpp>

#include <stdexcept>

namespace dynsyslab::workflow {
namespace {

constexpr const char* solveTemplate = R"({
  "schema": "dynsyslab.compute-job/v1", "kind": "solve", "label": "interactive-solve",
  "system": {"name": "duffing", "preset": "double_well_reference", "parameters": {}},
  "time": {"start": 0.0, "stop": 10.0, "output_step": 0.01, "coordinate_system": "automatic"},
  "initial_conditions": {"mode": "default", "count": 1, "clusters": 1, "seed": 0,
                         "values": [], "ranges": []},
  "controller": {"kind": "none", "settings": {}, "target": {"kind": "constant", "state": [0.0, 0.0]}},
  "campaign": {"enabled": false, "mode": "cartesian", "axes": []},
  "solver": {"method": "RK45", "rtol": 1e-8, "atol": 1e-10},
  "execution": {"max_threads": 0, "max_memory_mib": 0},
  "queue": {"retry_limit": 0, "depends_on": [], "hold": false, "reuse": "none"},
  "archive": {"enabled": true}
})";

constexpr const char* postprocessTemplate = R"({
  "schema": "dynsyslab.compute-job/v1", "kind": "postprocess", "label": "interactive-postprocess",
  "sources": [],
  "methods": [
    {"name": "hdmd", "enabled": true,
     "settings": {"observable": 0, "matrix_size": 200, "sampling_window": 1000}},
    {"name": "sindy", "enabled": true, "settings": {"archive_storage": "coefficient-first"}}
  ],
  "execution": {"max_threads": 0, "max_memory_mib": 0},
  "archive": {"enabled": true}
})";

constexpr const char* estimateTemplate = R"({
  "schema": "dynsyslab.compute-job/v1", "kind": "postprocess", "label": "interactive-estimation",
  "sources": [],
  "methods": [
    {"name": "observability", "enabled": true,
     "settings": {"krom": "", "measurement": "echo", "state_reconstruction": "state",
                  "observability_method": "continuous-gramian", "observability_horizon": 10.0}}
  ],
  "archive": {"enabled": true}
})";

}  // namespace

std::vector<std::string> jobTemplateNames() { return {"solve", "postprocess", "estimate"}; }

json::Value jobTemplate(const std::string& name) {
    if (name == "solve") return json::parse(solveTemplate);
    if (name == "postprocess") return json::parse(postprocessTemplate);
    if (name == "estimate") return json::parse(estimateTemplate);
    throw std::invalid_argument("template must be solve, postprocess or estimate");
}

}  // namespace dynsyslab::workflow
