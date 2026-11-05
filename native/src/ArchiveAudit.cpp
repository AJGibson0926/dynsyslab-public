#include <dynsyslab/core/ArchiveAudit.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <map>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
namespace json = dynsyslab::core::json;

struct Record {
    std::filesystem::path path;
    std::uintmax_t bytes{};
    std::string exactKey;
    std::string planKey;
};

std::uint64_t hashBytes(std::uint64_t hash, const void* data, std::size_t size) {
    const auto* bytes = static_cast<const std::uint8_t*>(data);
    for (std::size_t index = 0; index < size; ++index) {
        hash ^= bytes[index];
        hash *= 1099511628211ull;
    }
    return hash;
}

std::string arrayKey(const NpyArray& array) {
    std::uint64_t hash = 1469598103934665603ull;
    hash = hashBytes(hash, array.descriptor.data(), array.descriptor.size());
    hash = hashBytes(hash, array.shape.data(), array.shape.size() * sizeof(std::size_t));
    hash = hashBytes(hash, array.data.data(), array.data.size());
    std::ostringstream out;
    out << array.descriptor << ':';
    for (const auto extent : array.shape) out << extent << 'x';
    out << ':' << std::hex << hash;
    return out.str();
}

json::Value planIdentity(const json::Value& metadata, const NpyArray& values,
                         const NpyArray* time) {
    json::Value::Object plan;
    for (const auto* key : {"system", "system_preset", "system_parameters",
                            "resolved_system_parameters", "initial_state_override",
                            "initial_field", "integration_form", "coordinate_system",
                            "time_span", "options", "solver", "forcing",
                            "controller_configuration"}) {
        const auto& value = metadata.at(key);
        if (!value.isNull()) plan[key] = value;
    }
    json::Value::Array valueShape;
    for (const auto extent : values.shape) valueShape.emplace_back(static_cast<double>(extent));
    plan["values_shape"] = json::Value(std::move(valueShape));
    if (time) {
        json::Value::Array timeShape;
        for (const auto extent : time->shape) timeShape.emplace_back(static_cast<double>(extent));
        plan["time_shape"] = json::Value(std::move(timeShape));
    }
    return json::Value(std::move(plan));
}

Record inspect(const std::filesystem::path& path) {
    const auto archive = readNpz(path);
    if (!archive.contains("values"))
        throw std::runtime_error("not a solution-like archive with a values array");
    const auto& values = archive.at("values");
    const auto* time = archive.contains("time") ? &archive.at("time") : nullptr;
    const auto metadataText = archive.metadataJson();
    const auto metadata = metadataText.empty() ? json::Value(json::Value::Object{})
                                                : json::parse(metadataText);
    std::string exact = arrayKey(values);
    if (time) exact += '|' + arrayKey(*time);
    return {path, std::filesystem::file_size(path), std::move(exact),
            json::dump(planIdentity(metadata, values, time))};
}

bool containsPart(const std::filesystem::path& path, const std::string& part) {
    return std::any_of(path.begin(), path.end(), [&](const auto& item) {
        return item.u8string() == part;
    });
}

bool selected(const std::filesystem::path& path, const ArchiveAuditScope& scope,
              const std::chrono::system_clock::time_point now) {
    if (path.filename().u8string().find(".partial.npz") != std::string::npos) return false;
    if (!scope.includeDerived && containsPart(path, "derived")) return false;
    if (!scope.systems.empty() && std::none_of(scope.systems.begin(), scope.systems.end(),
        [&](const auto& name) { return containsPart(path, name); })) return false;
    if (scope.modifiedWithinDays >= 0.0) {
        const auto modified = std::filesystem::last_write_time(path);
        const auto translated = now + (modified - std::filesystem::file_time_type::clock::now());
        if (now - translated > std::chrono::duration<double>(scope.modifiedWithinDays * 86400.0))
            return false;
    }
    return true;
}

std::vector<ArchiveDuplicateGroup> groups(const std::map<std::string, std::vector<Record>>& map,
                                          const std::string& kind) {
    std::vector<ArchiveDuplicateGroup> result;
    for (const auto& [key, records] : map) {
        (void)key;
        if (records.size() < 2) continue;
        ArchiveDuplicateGroup group; group.kind = kind;
        for (const auto& record : records) group.files.push_back(record.path);
        std::sort(group.files.begin(), group.files.end());
        if (kind == "exact_solution") {
            std::vector<std::uintmax_t> sizes;
            for (const auto& record : records) sizes.push_back(record.bytes);
            std::sort(sizes.begin(), sizes.end(), std::greater<>());
            for (std::size_t index = 1; index < sizes.size(); ++index)
                group.redundantBytes += sizes[index];
        }
        result.push_back(std::move(group));
    }
    std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) {
        return left.files.front() < right.files.front();
    });
    return result;
}

bool sameArray(const NpzRecord& left, const NpzRecord& right, const std::string& name) {
    if (left.contains(name) != right.contains(name)) return false;
    if (!left.contains(name)) return true;
    const auto& a = left.at(name); const auto& b = right.at(name);
    return a.descriptor == b.descriptor && a.shape == b.shape && a.data == b.data;
}

bool sameSolution(const std::filesystem::path& left,
                  const std::filesystem::path& right) {
    const auto a = readNpz(left), b = readNpz(right);
    return sameArray(a, b, "time") && sameArray(a, b, "values");
}

std::vector<ArchiveDuplicateGroup> exactGroups(
        const std::map<std::string, std::vector<Record>>& buckets) {
    std::vector<ArchiveDuplicateGroup> result;
    for (const auto& [key, candidates] : buckets) {
        (void)key;
        std::vector<std::vector<const Record*>> partitions;
        for (const auto& candidate : candidates) {
            const auto match = std::find_if(partitions.begin(), partitions.end(),
                [&](const auto& group) { return sameSolution(group.front()->path, candidate.path); });
            if (match == partitions.end()) partitions.push_back({&candidate});
            else match->push_back(&candidate);
        }
        for (const auto& partition : partitions) {
            if (partition.size() < 2) continue;
            ArchiveDuplicateGroup group; group.kind = "exact_solution";
            std::vector<std::uintmax_t> sizes;
            for (const auto* record : partition) {
                group.files.push_back(record->path); sizes.push_back(record->bytes);
            }
            std::sort(group.files.begin(), group.files.end());
            std::sort(sizes.begin(), sizes.end(), std::greater<>());
            for (std::size_t index = 1; index < sizes.size(); ++index)
                group.redundantBytes += sizes[index];
            result.push_back(std::move(group));
        }
    }
    std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) {
        return left.files.front() < right.files.front();
    });
    return result;
}

json::Value stringArray(const std::vector<std::string>& values) {
    json::Value::Array result;
    for (const auto& value : values) result.emplace_back(value);
    return json::Value(std::move(result));
}

json::Value duplicateJson(const ArchiveDuplicateGroup& group) {
    json::Value::Array files;
    for (const auto& path : group.files) files.emplace_back(path.generic_u8string());
    json::Value::Object result;
    result["kind"] = group.kind;
    result["count"] = static_cast<double>(group.files.size());
    result["files"] = json::Value(std::move(files));
    result["redundant_bytes"] = static_cast<double>(group.redundantBytes);
    return json::Value(std::move(result));
}

std::filesystem::path safeSection(const std::string& name) {
    const auto path = std::filesystem::u8path(name);
    if (name.empty() || path.is_absolute() || path.has_parent_path()
            || name == "." || name == "..")
        throw std::invalid_argument("archive section must be one relative directory name");
    return path;
}

}  // namespace

ArchiveAuditReport auditArchives(const std::filesystem::path& root,
                                 const ArchiveAuditScope& scope) {
    ArchiveAuditReport report; report.root = std::filesystem::absolute(root);
    std::map<std::string, std::vector<Record>> exact, plans;
    const auto now = std::chrono::system_clock::now();
    const auto sections = scope.sections.empty() ? std::vector<std::string>{"solutions"}
                                                  : scope.sections;
    for (const auto& section : sections) {
        const auto base = report.root / safeSection(section);
        if (!std::filesystem::exists(base)) continue;
        for (const auto& entry : std::filesystem::recursive_directory_iterator(base)) {
            if (!entry.is_regular_file() || entry.path().extension() != ".npz") continue;
            try {
                if (!selected(entry.path(), scope, now)) continue;
                auto record = inspect(entry.path());
                report.totalBytes += record.bytes; ++report.scannedFiles;
                exact[record.exactKey].push_back(record);
                plans[record.planKey].push_back(std::move(record));
            } catch (const std::exception& error) {
                report.unreadableFiles.push_back({entry.path(), error.what()});
            }
        }
    }
    report.exactDuplicateGroups = exactGroups(exact);
    report.samePlanGroups = groups(plans, "same_plan");
    return report;
}

std::string archiveAuditJson(const ArchiveAuditReport& report,
                             const ArchiveAuditScope& scope) {
    json::Value::Object scopeJson;
    scopeJson["modified_within_days"] = scope.modifiedWithinDays < 0.0
        ? json::Value{} : json::Value(scope.modifiedWithinDays);
    scopeJson["sections"] = stringArray(scope.sections);
    scopeJson["systems"] = stringArray(scope.systems);
    scopeJson["include_derived"] = scope.includeDerived;
    json::Value::Array exact, plans, unreadable;
    for (const auto& group : report.exactDuplicateGroups) exact.push_back(duplicateJson(group));
    for (const auto& group : report.samePlanGroups) plans.push_back(duplicateJson(group));
    for (const auto& issue : report.unreadableFiles) {
        json::Value::Object value;
        value["path"] = issue.path.generic_u8string(); value["error"] = issue.message;
        unreadable.emplace_back(std::move(value));
    }
    std::uintmax_t redundantBytes{}; std::size_t redundantFiles{};
    for (const auto& group : report.exactDuplicateGroups) {
        redundantBytes += group.redundantBytes;
        redundantFiles += group.files.size() - 1;
    }
    json::Value::Object result;
    result["schema"] = "dynsyslab.archive-audit/v1";
    result["root"] = report.root.generic_u8string();
    result["scope"] = json::Value(std::move(scopeJson));
    result["scanned_files"] = static_cast<double>(report.scannedFiles);
    result["total_bytes"] = static_cast<double>(report.totalBytes);
    result["redundant_file_count"] = static_cast<double>(redundantFiles);
    result["redundant_bytes"] = static_cast<double>(redundantBytes);
    result["exact_duplicate_groups"] = json::Value(std::move(exact));
    result["same_plan_groups"] = json::Value(std::move(plans));
    result["unreadable_files"] = json::Value(std::move(unreadable));
    result["read_only"] = true;
    return json::dump(json::Value(std::move(result)));
}

}  // namespace dynsyslab::core
