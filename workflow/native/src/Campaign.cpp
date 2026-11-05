#include <dynsyslab/workflow/Campaign.hpp>

#include <dynsyslab/core/Hierarchy.hpp>
#include <dynsyslab/core/Text.hpp>

#include "SolveDispatch.hpp"
#include "WorkflowSupport.hpp"

#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::workflow {
namespace {

struct Axis {
    std::string path;
    std::vector<json::Value> values;
};

bool supportsAxis(const std::string& path) {
    const auto begins = [&path](const char* prefix) {
        return path.rfind(prefix, 0) == 0;
    };
    return path == "system.preset" || begins("system.parameters.")
        || begins("controller.settings.") || begins("controller.target.")
        || begins("time.") || begins("solver.")
        || begins("initial_conditions.");
}

std::vector<Axis> axes(const json::Value& campaign) {
    if (campaign.at("mode").string("cartesian") != "cartesian")
        throw std::invalid_argument("native campaigns currently require cartesian mode");
    std::vector<Axis> result;
    for (const auto& source : campaign.at("axes").array()) {
        Axis axis;
        axis.path = source.at("path").string();
        if (!supportsAxis(axis.path))
            throw std::invalid_argument(
                "campaign axis must address system parameters, controller settings/target, "
                "time, solver, or initial-condition fields");
        axis.values = source.at("values").array();
        if (axis.values.empty()) throw std::invalid_argument("campaign axis has no values");
        result.push_back(std::move(axis));
    }
    if (result.empty()) throw std::invalid_argument("enabled campaign requires at least one axis");
    return result;
}

std::size_t pointCount(const std::vector<Axis>& axes) {
    std::size_t result = 1;
    for (const auto& axis : axes) {
        if (result > 100000 / axis.values.size())
            throw std::invalid_argument("campaign contains more than 100000 Cartesian points");
        result *= axis.values.size();
    }
    return result;
}

void atomicText(const std::filesystem::path& path, const std::string& contents) {
    std::filesystem::create_directories(path.parent_path());
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    { std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
      if (!stream) throw std::runtime_error("unable to create campaign manifest");
      stream << contents << '\n'; }
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(),
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish campaign manifest");
#else
    std::filesystem::rename(temporary, path);
#endif
}

dynsyslab::core::ExecutionContext pointExecution(
    const dynsyslab::core::ExecutionContext& parent,
    const std::filesystem::path& progressPath, std::size_t point,
    std::size_t count) {
    dynsyslab::core::ExecutionContext child;
    child.limits = parent.limits;
    child.cancellation = [&parent] { return parent.cancelled(); };
    child.progress = [&parent, progressPath, point, count](
                         dynsyslab::core::ExecutionProgress event) {
        const double local = std::clamp(event.fraction, 0.0, 1.0);
        const double global = (static_cast<double>(point) + local) / count;
        const std::string message = "Campaign point " + std::to_string(point + 1)
            + " of " + std::to_string(count) + ": " + event.message;
        detail::atomicProgress(progressPath, global,
            point + static_cast<std::size_t>(local >= 1.0), count, message);
        parent.report({"campaign/" + event.stage, message, global,
            point + static_cast<std::size_t>(local >= 1.0), count});
    };
    return child;
}

}  // namespace

bool hasEnabledCampaign(const json::Value& job) {
    return job.at("campaign").at("enabled").boolean(false);
}

std::size_t campaignPointCount(const json::Value& job) {
    return hasEnabledCampaign(job) ? pointCount(axes(job.at("campaign"))) : 1;
}

std::wstring CampaignResult::summary() const {
    std::wostringstream out;
    out << L"Native campaign completed\r\n\r\nCampaign: " << dynsyslab::core::wide(label)
        << L"\r\nPoints: " << points.size() << L"\r\nPython required: no\r\n\r\nManifest:\r\n"
        << manifest.wstring();
    return out.str();
}

CampaignResult executeCampaign(const json::Value& job,
                                           const std::filesystem::path& archiveRoot,
                                           const std::filesystem::path& progressPath,
                                           const dynsyslab::core::ExecutionContext* execution) {
    const auto selectedAxes = axes(job.at("campaign"));
    const std::size_t count = campaignPointCount(job);
    CampaignResult result;
    result.label = job.at("label").string("campaign");
    const std::string stamp = detail::timestamp();
    result.id = job.at("campaign").at("id").string(
        dynsyslab::core::stableIdentifier(result.label + "|" + stamp + "|campaign"));
    result.manifest = archiveRoot / L"campaigns"
        / dynsyslab::core::wide(
            detail::safeName(result.label, "campaign") + "__" + stamp + ".campaign.json");
    result.points.reserve(count);
    detail::atomicProgress(progressPath, 0.0, 0, count, "Preparing campaign");
    for (std::size_t point = 0; point < count; ++point) {
        if (execution) execution->report({"campaign", "Running campaign point",
            static_cast<double>(point) / count, point, count});
        json::Value child = job;
        std::size_t divisor = 1;
        CampaignPoint record;
        for (const auto& axis : selectedAxes) {
            const std::size_t selected = (point / divisor) % axis.values.size();
            json::setPath(child, axis.path, axis.values[selected]);
            record.coordinates.push_back(json::dump(axis.values[selected]));
            divisor *= axis.values.size();
        }
        child.at("label") = result.label + "-point-" + std::to_string(point);
        child.at("campaign").at("enabled") = json::Value(false);
        child.at("campaign").at("id") = result.id;
        child.at("campaign").at("kind") = job.at("campaign").at("kind").string("control_sweep");
        child.at("campaign").at("point_index") = json::Value(static_cast<double>(point));
        dynsyslab::core::ExecutionContext childExecution;
        const dynsyslab::core::ExecutionContext* active = nullptr;
        if (execution) {
            childExecution = pointExecution(*execution, progressPath, point, count);
            active = &childExecution;
        }
        record.solve = detail::executeSingleSolve(child, archiveRoot,
            progressPath.empty() ? std::filesystem::path{} :
                progressPath.parent_path() / (progressPath.stem().wstring()
                    + L"-point-" + std::to_wstring(point) + L".json"), active);
        result.points.push_back(std::move(record));
        detail::atomicProgress(progressPath,
            static_cast<double>(point + 1) / count, point + 1, count,
            "Completed campaign point " + std::to_string(point + 1)
                + " of " + std::to_string(count));
    }
    json::Value::Object manifest;
    manifest["format"] = json::Value("dynsyslab.campaign/v3");
    manifest["id"] = json::Value(result.id); manifest["label"] = json::Value(result.label);
    manifest["created_utc"] = json::Value(stamp);
    manifest["kind"] = json::Value(job.at("campaign").at("kind").string("control_sweep"));
    manifest["axes"] = job.at("campaign").at("axes");
    json::Value::Array points;
    for (std::size_t index = 0; index < result.points.size(); ++index) {
        json::Value::Object point;
        point["index"] = json::Value(static_cast<double>(index));
        json::Value::Array coordinates;
        for (const auto& coordinate : result.points[index].coordinates)
            coordinates.push_back(json::parse(coordinate));
        point["coordinates"] = json::Value(std::move(coordinates));
        point["hierarchy_manifest"] = json::Value(
            result.points[index].solve.hierarchyManifest.generic_u8string());
        json::Value::Array archives;
        for (const auto& path : result.points[index].solve.archives)
            archives.emplace_back(path.generic_u8string());
        point["solution_archives"] = json::Value(std::move(archives));
        points.emplace_back(std::move(point));
    }
    manifest["points"] = json::Value(std::move(points));
    atomicText(result.manifest, json::dump(json::Value(std::move(manifest))));
    if (execution) execution->report({"campaign", "Campaign archived", 1.0, count, count});
    return result;
}

}  // namespace dynsyslab::workflow
