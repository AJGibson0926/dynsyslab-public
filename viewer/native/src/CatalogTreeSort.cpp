#include "dynsyslab_viewer/CatalogTreeSort.hpp"

#include <algorithm>
#include <charconv>
#include <cctype>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>

namespace dsl {
namespace {

struct AggregateBuild {
    CatalogTreeAggregate value;
    std::map<std::string, double> sources;
};

std::string sourceKey(std::string source) {
    std::transform(source.begin(), source.end(), source.begin(), [](unsigned char c) {
        if (c == '\\') return '/';
        return static_cast<char>(std::tolower(c));
    });
    return source;
}

bool trajectoryBearing(const CatalogItem& item) {
    return item.kind == "trajectory" || item.kind == "control_attempt";
}

std::size_t trajectoryCount(const CatalogItem& item) {
    const auto found = item.properties.find("trajectory count");
    if (found != item.properties.end()) {
        std::size_t value = 0;
        const char* begin = found->second.data();
        const char* end = begin + found->second.size();
        const auto parsed = std::from_chars(begin, end, value);
        if (parsed.ec == std::errc{} && parsed.ptr == end) return value;
    }
    return trajectoryBearing(item) ? 1 : 0;
}

bool naturalLess(const std::string& left, const std::string& right) {
    std::size_t a = 0, b = 0;
    while (a < left.size() && b < right.size()) {
        const unsigned char ca = left[a], cb = right[b];
        if (std::isdigit(ca) && std::isdigit(cb)) {
            std::size_t ae = a, be = b;
            while (ae < left.size() && std::isdigit(static_cast<unsigned char>(left[ae]))) ++ae;
            while (be < right.size() && std::isdigit(static_cast<unsigned char>(right[be]))) ++be;
            std::size_t az = a, bz = b;
            while (az + 1 < ae && left[az] == '0') ++az;
            while (bz + 1 < be && right[bz] == '0') ++bz;
            if (ae - az != be - bz) return ae - az < be - bz;
            const int numeric = left.compare(az, ae - az, right, bz, be - bz);
            if (numeric != 0) return numeric < 0;
            a = ae; b = be; continue;
        }
        const auto al = static_cast<unsigned char>(std::tolower(ca));
        const auto bl = static_cast<unsigned char>(std::tolower(cb));
        if (al != bl) return al < bl;
        ++a; ++b;
    }
    return left.size() < right.size();
}

bool trainingDataItem(const Catalog& catalog, int index,
                      const std::unordered_map<std::string, int>& byId,
                      bool presented) {
    if (index < 0 || static_cast<std::size_t>(index) >= catalog.items.size())
        return false;
    std::vector<unsigned char> visited(catalog.items.size());
    while (index >= 0 && static_cast<std::size_t>(index) < catalog.items.size() &&
           !visited[static_cast<std::size_t>(index)]) {
        visited[static_cast<std::size_t>(index)] = 1;
        const auto& item = catalog.items[static_cast<std::size_t>(index)];
        if (item.kind == "presentation_group" && item.label == "Training data")
            return true;
        if (!presented) {
            if (item.category == "Training data") return true;
            const auto archiveCategory = item.properties.find("archive category");
            if (archiveCategory != item.properties.end() &&
                archiveCategory->second == "Training data") return true;
        }
        const auto parent = byId.find(item.parentId);
        if (parent == byId.end()) break;
        index = parent->second;
    }
    return false;
}

}  // namespace

OrderedCatalogTree orderedCatalogTree(const Catalog& catalog,
                                      CatalogTreeSort sort) {
    OrderedCatalogTree result;
    result.children.resize(catalog.items.size());
    result.aggregates.resize(catalog.items.size());
    std::unordered_map<std::string, int> byId;
    for (std::size_t i = 0; i < catalog.items.size(); ++i)
        byId.emplace(catalog.items[i].id, static_cast<int>(i));
    for (std::size_t i = 0; i < catalog.items.size(); ++i) {
        const auto parent = byId.find(catalog.items[i].parentId);
        if (!catalog.items[i].parentId.empty() && parent != byId.end() &&
            parent->second != static_cast<int>(i))
            result.children[static_cast<std::size_t>(parent->second)].push_back(
                static_cast<int>(i));
        else
            result.roots.push_back(static_cast<int>(i));
    }

    std::vector<unsigned char> state(catalog.items.size());
    std::vector<AggregateBuild> built(catalog.items.size());
    std::function<const AggregateBuild&(int)> aggregate = [&](int index)
        -> const AggregateBuild& {
        const auto position = static_cast<std::size_t>(index);
        if (state[position] == 2) return built[position];
        if (state[position] == 1) return built[position];
        state[position] = 1;
        const auto& item = catalog.items[position];
        auto& current = built[position];
        current.value.newest = item.modified;
        current.value.trajectoryCount = trajectoryCount(item);
        if (!item.source.empty()) current.sources[sourceKey(item.source)] = item.bytes;
        for (const int child : result.children[position]) {
            const auto& nested = aggregate(child);
            current.value.newest = std::max(current.value.newest, nested.value.newest);
            current.value.trajectoryCount += nested.value.trajectoryCount;
            for (const auto& source : nested.sources) {
                auto& bytes = current.sources[source.first];
                bytes = std::max(bytes, source.second);
            }
        }
        for (const auto& source : current.sources)
            current.value.totalBytes += source.second;
        result.aggregates[position] = current.value;
        state[position] = 2;
        return current;
    };
    for (std::size_t i = 0; i < catalog.items.size(); ++i)
        (void)aggregate(static_cast<int>(i));

    const auto before = [&](int left, int right) {
        const auto& a = result.aggregates[static_cast<std::size_t>(left)];
        const auto& b = result.aggregates[static_cast<std::size_t>(right)];
        if (sort == CatalogTreeSort::date && a.newest != b.newest)
            return a.newest > b.newest;
        if (sort == CatalogTreeSort::trajectoryCount &&
            a.trajectoryCount != b.trajectoryCount)
            return a.trajectoryCount > b.trajectoryCount;
        if (sort == CatalogTreeSort::totalBytes && a.totalBytes != b.totalBytes)
            return a.totalBytes > b.totalBytes;
        return naturalLess(catalog.items[static_cast<std::size_t>(left)].label,
                           catalog.items[static_cast<std::size_t>(right)].label);
    };
    const auto sortLeaves = [&](std::vector<int>& siblings) {
        std::vector<std::size_t> slots;
        std::vector<int> leaves;
        for (std::size_t i = 0; i < siblings.size(); ++i) {
            const int index = siblings[i];
            if (result.children[static_cast<std::size_t>(index)].empty()) {
                slots.push_back(i);
                leaves.push_back(index);
            }
        }
        std::stable_sort(leaves.begin(), leaves.end(), before);
        for (std::size_t i = 0; i < slots.size(); ++i)
            siblings[slots[i]] = leaves[i];
    };
    sortLeaves(result.roots);
    for (auto& siblings : result.children) sortLeaves(siblings);
    return result;
}

bool isTrainingDataItem(const Catalog& catalog, int index) {
    std::unordered_map<std::string, int> byId;
    for (std::size_t i = 0; i < catalog.items.size(); ++i)
        byId.emplace(catalog.items[i].id, static_cast<int>(i));
    const bool presented = std::any_of(catalog.items.begin(), catalog.items.end(),
        [](const auto& item) { return item.kind == "presentation_group"; });
    return trainingDataItem(catalog, index, byId, presented);
}

int preferredInitialCatalogItem(const Catalog& catalog,
                                CatalogTreeSort sort) {
    const auto tree = orderedCatalogTree(catalog, sort);
    std::unordered_map<std::string, int> byId;
    for (std::size_t i = 0; i < catalog.items.size(); ++i)
        byId.emplace(catalog.items[i].id, static_cast<int>(i));
    const bool presented = std::any_of(catalog.items.begin(), catalog.items.end(),
        [](const auto& item) { return item.kind == "presentation_group"; });
    int selected = -1;
    const auto before = [&](int left, int right) {
        const auto& leftItem = catalog.items[static_cast<std::size_t>(left)];
        const auto& rightItem = catalog.items[static_cast<std::size_t>(right)];
        const bool leftSet = leftItem.label.rfind("Training set #", 0) == 0;
        const bool rightSet = rightItem.label.rfind("Training set #", 0) == 0;
        if (leftSet != rightSet) return leftSet;
        const auto contextual = [](const CatalogItem& item) {
            const auto found = item.properties.find("physics");
            return found != item.properties.end() && !found->second.empty() &&
                found->second != "General parameters";
        };
        const bool leftContext = contextual(leftItem);
        const bool rightContext = contextual(rightItem);
        if (leftContext != rightContext) return leftContext;
        const auto& a = tree.aggregates[static_cast<std::size_t>(left)];
        const auto& b = tree.aggregates[static_cast<std::size_t>(right)];
        if (sort == CatalogTreeSort::date && a.newest != b.newest)
            return a.newest > b.newest;
        if (sort == CatalogTreeSort::trajectoryCount &&
            a.trajectoryCount != b.trajectoryCount)
            return a.trajectoryCount > b.trajectoryCount;
        if (sort == CatalogTreeSort::totalBytes && a.totalBytes != b.totalBytes)
            return a.totalBytes > b.totalBytes;
        return naturalLess(catalog.items[static_cast<std::size_t>(left)].label,
                           catalog.items[static_cast<std::size_t>(right)].label);
    };
    for (std::size_t i = 0; i < catalog.items.size(); ++i) {
        const auto& item = catalog.items[i];
        if (item.kind == "presentation_group" || item.views.empty() ||
            !trainingDataItem(catalog, static_cast<int>(i), byId, presented)) continue;
        if (selected < 0 || before(static_cast<int>(i), selected))
            selected = static_cast<int>(i);
    }
    return selected;
}

}  // namespace dsl
