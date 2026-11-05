#include <dynsyslab/core/Hierarchy.hpp>

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::core {
namespace {

std::string quote(const std::string& value) {
    std::ostringstream out;
    out << '"';
    for (const unsigned char item : value) {
        switch (item) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (item < 0x20) out << "\\u" << std::hex << std::setw(4)
                                 << std::setfill('0') << static_cast<int>(item) << std::dec;
            else out << static_cast<char>(item);
        }
    }
    return out.str() + '"';
}

template<class Items, class Writer>
void array(std::ostringstream& out, const Items& items, Writer writer) {
    out << '[';
    bool first = true;
    for (const auto& item : items) {
        if (!first) out << ',';
        first = false;
        writer(out, item);
    }
    out << ']';
}

void stringArray(std::ostringstream& out, const std::vector<std::string>& items) {
    array(out, items, [](auto& stream, const auto& item) { stream << quote(item); });
}

void numberArray(std::ostringstream& out, const std::vector<double>& items) {
    array(out, items, [](auto& stream, double item) {
        stream << std::setprecision(17) << item;
    });
}

std::string genericPath(const std::filesystem::path& value) {
    return value.generic_u8string();
}

std::filesystem::path portableReference(const std::filesystem::path& target,
                                        const std::filesystem::path& manifestDirectory) {
    std::error_code targetError, baseError, relativeError;
    const auto absoluteTarget = target.is_absolute()
        ? target : std::filesystem::absolute(target, targetError);
    const auto absoluteBase = manifestDirectory.is_absolute()
        ? manifestDirectory : std::filesystem::absolute(manifestDirectory, baseError);
    if (!targetError && !baseError) {
        const auto relative = std::filesystem::relative(
            absoluteTarget, absoluteBase, relativeError);
        if (!relativeError && !relative.empty()) return relative;
    }
    return target;
}

}  // namespace

std::string stableIdentifier(const std::string& value) {
    std::uint64_t hash = 1469598103934665603ull;
    for (const unsigned char item : value) {
        hash ^= item;
        hash *= 1099511628211ull;
    }
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}

std::string parameterSetKey(const SystemConfiguration& configuration) {
    return equationsFingerprint(configuration);
}

std::string hierarchyJson(const BatchRecord& batch) {
    std::ostringstream out;
    out << "{\"format\":\"dynsyslab.hierarchy/v2\",\"tier\":\"batch\",\"id\":"
        << quote(batch.id) << ",\"label\":" << quote(batch.label)
        << ",\"created_utc\":" << quote(batch.createdUtc)
        << ",\"data_role\":" << quote(batch.dataRole)
        << ",\"campaign_id\":" << quote(batch.campaignId)
        << ",\"campaign_kind\":" << quote(batch.campaignKind) << ",\"tags\":";
    stringArray(out, batch.tags);
    out << ",\"properties\":{";
    bool firstProperty = true;
    for (const auto& entry : batch.properties) {
        if (!firstProperty) out << ',';
        firstProperty = false;
        out << quote(entry.first) << ':' << quote(entry.second);
    }
    out << "},\"systems\":";
    array(out, batch.systems, [](auto& stream, const SystemRunRecord& system) {
        stream << "{\"id\":" << quote(system.id) << ",\"label\":" << quote(system.label)
               << ",\"system\":" << quote(system.systemIdentifier)
               << ",\"system_label\":" << quote(system.systemLabel)
               << ",\"system_family\":" << quote(system.systemFamily)
               << ",\"preset\":" << quote(system.preset)
               << ",\"preset_label\":" << quote(system.presetLabel)
               << ",\"parameter_set_key\":" << quote(system.parameterSetKey)
               << ",\"coordinate_system\":" << quote(system.coordinateSystem)
               << ",\"equations_version\":" << quote(system.equationsVersion)
               << ",\"equations_fingerprint\":" << quote(system.equationsFingerprint)
               << ",\"parameter_overrides\":{";
        bool firstParameter = true;
        for (const auto& entry : system.parameterOverrides) {
            if (!firstParameter) stream << ',';
            firstParameter = false;
            stream << quote(entry.first) << ':' << std::setprecision(17) << entry.second;
        }
        stream << "},\"resolved_parameters\":{";
        firstParameter = true;
        for (const auto& entry : system.resolvedParameters) {
            if (!firstParameter) stream << ',';
            firstParameter = false;
            stream << quote(entry.first) << ':' << std::setprecision(17) << entry.second;
        }
        stream << "},\"parameter_units\":{";
        firstParameter = true;
        for (const auto& entry : system.parameterUnits) {
            if (!firstParameter) stream << ',';
            firstParameter = false;
            stream << quote(entry.first) << ':' << quote(entry.second);
        }
        stream << "},\"dimensionless_groups\":{";
        firstParameter = true;
        for (const auto& entry : system.dimensionlessGroups) {
            if (!firstParameter) stream << ',';
            firstParameter = false;
            stream << quote(entry.first) << ':' << std::setprecision(17) << entry.second;
        }
        stream << "},\"dimensionless_group_units\":{";
        firstParameter = true;
        for (const auto& entry : system.dimensionlessGroupUnits) {
            if (!firstParameter) stream << ',';
            firstParameter = false;
            stream << quote(entry.first) << ':' << quote(entry.second);
        }
        stream << "},\"state_names\":";
        stringArray(stream, system.stateNames);
        stream << ",\"clusters\":";
        array(stream, system.clusters, [](auto& clusterStream, const ClusterRecord& cluster) {
            clusterStream << "{\"id\":" << quote(cluster.id) << ",\"label\":"
                          << quote(cluster.label) << ",\"trajectories\":";
            array(clusterStream, cluster.trajectories,
                [](auto& trajectoryStream, const TrajectoryRecord& trajectory) {
                    trajectoryStream << "{\"id\":" << quote(trajectory.id)
                                     << ",\"label\":" << quote(trajectory.label)
                                     << ",\"initial_state\":";
                    numberArray(trajectoryStream, trajectory.initialState);
                    trajectoryStream << ",\"solution_archive\":"
                                     << quote(genericPath(trajectory.solutionArchive))
                                     << ",\"tags\":";
                    stringArray(trajectoryStream, trajectory.tags);
                    trajectoryStream << ",\"time_series\":";
                    array(trajectoryStream, trajectory.timeSeries,
                        [](auto& seriesStream, const TimeSeriesReference& series) {
                            seriesStream << "{\"role\":" << quote(series.role)
                                         << ",\"time_array\":" << quote(series.timeArray)
                                         << ",\"value_array\":" << quote(series.valueArray)
                                         << '}';
                        });
                    trajectoryStream << '}';
                });
            clusterStream << '}';
        });
        stream << '}';
    });
    out << '}';
    return out.str();
}

void writeHierarchyAtomic(const std::filesystem::path& path, const BatchRecord& batch) {
    if (path.empty()) throw std::invalid_argument("hierarchy path is empty");
    std::filesystem::create_directories(path.parent_path());
    const auto temporary = path.parent_path() / (path.filename().wstring() + L".partial");
    BatchRecord portable = batch;
    for (auto& system : portable.systems)
        for (auto& cluster : system.clusters)
            for (auto& trajectory : cluster.trajectories)
                trajectory.solutionArchive = portableReference(
                    trajectory.solutionArchive, path.parent_path());
    {
        std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
        if (!stream) throw std::runtime_error("unable to create hierarchy manifest");
        stream << hierarchyJson(portable) << '\n';
        if (!stream) throw std::runtime_error("unable to write hierarchy manifest");
    }
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(),
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("unable to publish hierarchy manifest");
#else
    std::filesystem::rename(temporary, path);
#endif
}

}  // namespace dynsyslab::core
