#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/SystemLabels.hpp>
namespace dynsyslab::core {
std::vector<std::string> systemStateLabels(const SystemDescriptor &descriptor,
                                           std::size_t dimension) {
    std::vector<std::string> result;
    if (descriptor.identifier == "kuramoto") {
        for (std::size_t j = 0; j < dimension; ++j)
            result.push_back(R"(\theta_{)" + std::to_string(j + 1) + "}");
        return result;
    }
    result = descriptor.stateLabels;
    while (result.size() < dimension)
        result.push_back("state " + std::to_string(result.size()));
    result.resize(dimension);
    return result;
}
} // namespace dynsyslab::core
