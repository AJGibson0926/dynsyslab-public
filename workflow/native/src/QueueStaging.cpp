#include "QueueStaging.hpp"

#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/Text.hpp>

#include <fstream>
#include <stdexcept>
#include <utility>

namespace dynsyslab::workflow::detail {

json::Value stageInputs(const json::Value& job, const std::filesystem::path& root, const std::string& id) {
    const auto& plan = job.at("queue").at("staging");
    const auto& inputs = plan.at("inputs").array();
    if (inputs.empty()) return {};
    const auto directory = std::filesystem::absolute(root / L"job_manager" / L"staging" / core::wide(id));
    std::filesystem::create_directories(core::extendedLengthPath(directory));
    json::Value::Array copied;
    for (const auto& item : inputs) {
        const auto source = std::filesystem::absolute(std::filesystem::u8path(item.string())).lexically_normal();
        if (!std::filesystem::is_regular_file(core::extendedLengthPath(source)))
            throw std::invalid_argument("staging input is missing: " + source.generic_u8string());
        const auto staged = directory / source.filename();
        std::filesystem::copy_file(core::extendedLengthPath(source), core::extendedLengthPath(staged),
                                   std::filesystem::copy_options::overwrite_existing);
        const auto sourceDigest = core::sha256File(core::extendedLengthPath(source));
        const auto stagedDigest = core::sha256File(core::extendedLengthPath(staged));
        if (sourceDigest != stagedDigest)
            throw std::runtime_error("checksum mismatch while staging " + source.generic_u8string());
        copied.emplace_back(json::Value::Object{{"source", source.generic_u8string()},
            {"staged", staged.generic_u8string()}, {"sha256", sourceDigest},
            {"bytes", static_cast<double>(std::filesystem::file_size(core::extendedLengthPath(staged)))}});
    }
    json::Value manifest = json::Value::Object{{"schema", "dynsyslab.staging-manifest/v1"},
        {"job_id", id}, {"staging_directory", directory.generic_u8string()},
        {"inputs", json::Value(std::move(copied))},
        {"output_location", plan.at("output_location").string("primary")},
        {"stage_out", plan.at("stage_out").boolean(true)},
        {"cleanup_scratch", plan.at("cleanup_scratch").boolean(true)}};
    std::ofstream out(core::extendedLengthPath(directory / L"staging-manifest.json"), std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("unable to write the staging manifest");
    out << json::dump(manifest);
    return manifest;
}

void cleanupStaging(const json::Value& manifest) {
    if (manifest.isNull() || !manifest.at("cleanup_scratch").boolean(true)) return;
    const auto directory = std::filesystem::u8path(manifest.at("staging_directory").string());
    std::error_code ignored;
    std::filesystem::remove_all(core::extendedLengthPath(directory), ignored);
}

}  // namespace dynsyslab::workflow::detail
