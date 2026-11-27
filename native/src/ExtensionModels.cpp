#include "ExtensionModels.hpp"

#include <dynsyslab/core/Json.hpp>

#include <array>
#include <stdexcept>
#include <string>

namespace dynsyslab::core::detail {
namespace {

void check(int status, const std::array<char, 1024>& error, const char* operation) {
    if (status == 0) return;
    const std::string message(error.data());
    throw std::runtime_error(std::string("extension ") + operation + " failed"
                             + (message.empty() ? std::string() : ": " + message));
}

template <class Table>
ExtensionModel described(const Table& table, const char* family) {
    if (table.struct_size < sizeof(Table) || !table.descriptor_json)
        throw std::runtime_error(std::string("extension ") + family + " table is invalid");
    const auto descriptor = json::parse(table.descriptor_json);
    ExtensionModel model;
    model.family = family;
    model.identifier = descriptor.at("identifier").string();
    model.label = descriptor.at("label").string(model.identifier);
    model.descriptorJson = table.descriptor_json;
    if (model.identifier.empty()) throw std::runtime_error("extension model has no identifier");
    return model;
}

}  // namespace

std::vector<ExtensionModel> extensionModels(const dsl_extension_v1& extension) {
    std::vector<ExtensionModel> result;
    if (extension.struct_size < sizeof(dsl_extension_v1)) return result;
    for (std::size_t index = 0; index < extension.delay_count; ++index) {
        const auto* table = &extension.delays[index];
        auto model = described(*table, "dde");
        if (!table->rhs) throw std::runtime_error("extension delay model has no rhs");
        model.delayRhs = [table](const double* window, std::size_t windowCount, double* derivative,
                                 std::size_t stateCount) {
            std::array<char, 1024> error{};
            check(table->rhs(table->user_data, window, windowCount, derivative, stateCount,
                             error.data(), error.size()), error, "delay rhs");
        };
        result.push_back(std::move(model));
    }
    for (std::size_t index = 0; index < extension.field_count; ++index) {
        const auto* table = &extension.fields[index];
        auto model = described(*table, "pde");
        if (!table->rhs || !table->initial_field) throw std::runtime_error("extension field model is incomplete");
        model.initialField = [table](const std::string& configuration, std::size_t count) {
            std::vector<double> field(count);
            std::array<char, 1024> error{};
            check(table->initial_field(table->user_data, configuration.c_str(), field.data(), count,
                                       error.data(), error.size()), error, "initial field");
            return field;
        };
        model.fieldRhs = [table](double time, const double* field, std::size_t count, double* derivative) {
            std::array<char, 1024> error{};
            check(table->rhs(table->user_data, time, field, count, derivative, error.data(), error.size()),
                  error, "field rhs");
        };
        result.push_back(std::move(model));
    }
    for (std::size_t index = 0; index < extension.lattice_count; ++index) {
        const auto* table = &extension.lattices[index];
        auto model = described(*table, "lattice");
        if (!table->initial || !table->sweep) throw std::runtime_error("extension lattice model is incomplete");
        model.latticeInitial = [table](std::uint64_t seed, std::int32_t* sites, std::size_t count) {
            std::array<char, 1024> error{};
            check(table->initial(table->user_data, seed, sites, count, error.data(), error.size()),
                  error, "lattice initial state");
        };
        model.latticeSweep = [table](std::uint64_t seed, std::uint64_t sweep, std::int32_t* sites,
                                     std::size_t count, double* observables, std::size_t observableCount) {
            std::array<char, 1024> error{};
            check(table->sweep(table->user_data, seed, sweep, sites, count, observables, observableCount,
                               error.data(), error.size()), error, "lattice sweep");
        };
        result.push_back(std::move(model));
    }
    return result;
}

}  // namespace dynsyslab::core::detail
