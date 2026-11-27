#pragma once
// Dictionary terms written in a job, shared by the steps that take a dictionary (eigenfunction
// expansions, SINDy). A term is one of
//   {"monomial" | "chebyshev" | "laguerre" | "hermite" | "fourier_hermite": [orders...]},
//   {"exponential" | "sine" | "cosine": [w...], "phase": p},
//   {"gaussian" | "laplacian" | "lorentzian" | "matern": {"node": [...], "width", "smoothness", "amplitude"}},
//   {"expression": "sin(r)*v"} (in `variables`), or a term in the archived form ({"kind": ...});
// each optionally with "center" and "scale" (core/ExpansionDictionary.hpp).
#include <dynsyslab/core/ExpansionDictionary.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <memory>
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

namespace dictionary_terms_detail {
[[nodiscard]] std::vector<double> numbers(const json::Value& value);
}

[[nodiscard]] dynsyslab::core::DictionaryTerm dictionaryTerm(const json::Value& spec, std::size_t dimension,
                                                             const std::vector<std::string>& variables);
// A whole dictionary: a list of terms, an archived dictionary ({"schema": ...}), or radial kernels
// on centers ({"kind": "gaussian", "centers": [[...]], "width": w, "include_constant": false}).
// Null gives no dictionary.
[[nodiscard]] std::shared_ptr<const dynsyslab::core::ExpansionDictionary> jobDictionary(
    const json::Value& spec, std::size_t dimension, const std::vector<std::string>& variables = {});

}  // namespace dynsyslab::workflow::postprocess_detail
