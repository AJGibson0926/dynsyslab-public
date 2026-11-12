#include "CatalogPresentation.hpp"
#include "CatalogPresentationText.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <sstream>

namespace dsl {
namespace {

std::string evidence(const CatalogItem& item, const std::string& originalLabel) {
    std::string result = originalLabel + " " + item.kind + " " + item.category + " " + item.role;
    for (const auto& tag : item.tags) result += " " + tag;
    return lower(result);
}

std::string activity(const CatalogItem& item, const std::string& originalLabel,
                     bool trainsKrom) {
    const auto text = evidence(item, originalLabel);
    if (trainsKrom && item.kind != "krom") return "Training data";
    if (item.kind == "control_attempt" || has(property(item, "controller"), "lqr"))
        return "KLQR controls";
    if (item.kind == "control_map" || item.kind == "control_error" ||
        item.kind == "control_data" || has(text, "bpsi")) return "Bψ identification";
    if (item.kind == "krom") return "KROM models";
    if (item.kind == "koopman_validation") return "KROM validation";
    if (item.kind == "koopman_spectrum") return "Koopman spectra";
    if (item.kind == "krom_control_campaign") return "Control searches";
    if (item.kind == "observability" || item.kind == "estimation") return "Estimation";
    if (item.kind == "decomposition") return "Post-processing";
    if (item.kind == "batch") {
        if (has(text, "control")) return "Control runs";
        if (has(text, "periodic")) return "Periodic forcing";
        if (has(text, "campaign") || has(text, "sweep")) return "Parameter sweeps";
        if (has(text, "forced")) return "Driven responses";
        return "Free responses";
    }
    if (has(text, "persistent excitation") || has(text, "persistent-excitation") ||
        has(text, "symmetric pulses") || has(text, "symmetric-pulses"))
        return "Excitation data";
    if (has(item.category, "training")) return "Training data";
    return "Other data";
}

std::string majorLane(const std::string& study) {
    if (study == "Training data") return "Training data";
    if (study == "KROM models" || study == "KROM validation" ||
        study == "Koopman spectra" || study == "Bψ identification" ||
        study == "Post-processing") return "Models & assets";
    if (study == "KLQR controls" || study == "Control searches" ||
        study == "Estimation") return "Control & estimation";
    return "Other data";
}

std::string trainingMode(const CatalogItem& item,
                         const std::string& originalLabel) {
    const auto forcing = lower(property(item, "forcing type"));
    if (forcing == "free" || forcing == "unforced" || forcing == "autonomous")
        return "Free";
    if (forcing == "forced" || forcing == "driven" || forcing == "controlled" ||
        forcing == "periodic") return "Forced";
    const auto role = lower(property(item, "data role") + " " + item.role);
    if (has(role, "unforced") || has(role, "free")) return "Free";
    if (has(role, "forced") || has(role, "control") || has(role, "driven"))
        return "Forced";
    const auto text = evidence(item, originalLabel);
    if (has(text, "unforced") || has(text, "free")) return "Free";
    if (has(text, "forced") || has(text, "periodic") || has(text, "driven") ||
        has(text, "control")) return "Forced";
    return "Unspecified";
}

std::string targetLabel(const std::string& raw) {
    const auto text = lower(raw);
    if (text.empty()) return {};
    if (text.find("equilibrium") != std::string::npos) return "Target · EC";
    const auto offset = text.find("offset");
    if (offset != std::string::npos) {
        std::string number;
        for (std::size_t i = offset + 6; i < text.size(); ++i) {
            if (std::isdigit(static_cast<unsigned char>(text[i]))) number += text[i];
            else if (!number.empty()) break;
        }
        return number.empty() ? "Target · offset" : "Target · offset #" + number;
    }
    if (text.find("echo") != std::string::npos) return "Target · echo";
    if (text.find("radius") != std::string::npos) return "Target · radius";
    auto result = titleWords(raw);
    if (result.size() > 36) result.resize(36);
    return "Target · " + result;
}

std::vector<std::string> strategyLabels(const CatalogItem& item,
                                        const std::string& originalLabel,
                                        const std::string& study) {
    if (study != "KLQR controls" && study != "Bψ identification" &&
        study != "Excitation data") return {};
    const auto representation = lower(property(item, "control representation"));
    const auto text = lower(originalLabel + " " + representation);
    std::vector<std::string> result;
    if (has(text, "control-map-errors") || has(text, "error_fields") ||
        has(text, "error-fields")) result.push_back("Error fields");
    else if (has(text, "control-representations")) result.push_back("Representation comparison");
    else if (has(text, "exact") || has(text, "analytic")) result.push_back("Analytic reference");
    else if (has(text, "persistent")) result.push_back("Persistent-excitation data");
    else if (has(text, "symmetric")) result.push_back("Symmetric-pulse data");
    else if (has(text, "pulse")) result.push_back("Pulse data");
    if (has(text, "conjugate-complete")) result.push_back("Conjugate-complete bilinear fit");
    else if (has(text, "complex-bilinear")) result.push_back("Complex bilinear fit");
    else if (has(text, "bilinear")) result.push_back("Bilinear fit");
    else if (has(text, "constant")) result.push_back("Constant fit");
    return result;
}

std::string shortDescriptor(const CatalogItem& item, const std::string& originalLabel) {
    const auto text = lower(originalLabel);
    if (item.kind == "batch") {
        if (has(text, "phase-portrait") || has(text, "periodic")) return "Periodic forcing";
        if (has(text, "derivative")) return "Echo derivative";
        if (has(text, "echo")) return "Echo control";
        if (has(text, "baseline")) return "Baseline";
    }
    if (item.kind == "krom") {
        if (has(text, "hdmd2000")) return "H-DMD 2000";
        if (has(text, "hdmd500")) return "H-DMD 500";
        if (has(text, "hdmd")) return "H-DMD";
        if (has(text, "sindy")) return "SINDy";
        if (has(text, "extract")) {
            const auto modes = property(item, "KROM modes");
            return modes.empty() ? "Extract" : modes + "-mode extract";
        }
    }
    if (has(text, "persistent-excitation")) return "Persistent excitation";
    if (has(text, "symmetric-pulses")) return "Symmetric pulses";
    if (has(text, "phase-portrait")) return "Phase portrait";
    if (has(text, "baseline")) return "Baseline";
    return {};
}

std::string noun(const CatalogItem& item, const std::string& originalParentKind,
                 const std::string& originalLabel, bool trainsKrom) {
    if (originalParentKind == "cluster") return "IC";
    if (trainsKrom && originalParentKind.empty()) return "Training set";
    if (item.kind == "batch") return "Run";
    if (item.kind == "system_run") return "Parameter set";
    if (item.kind == "cluster") return "IC set";
    if (item.kind == "control_attempt") return "Attempt";
    if (item.kind == "krom") return "KROM";
    if (item.kind == "control_map") return "Fit";
    if (item.kind == "control_error") return "Error";
    if (item.kind == "control_data" && has(originalLabel, "error")) return "Error";
    if (item.kind == "control_data" && has(originalLabel, "representation"))
        return "Comparison";
    if (item.kind == "control_data") return "Dataset";
    if (item.kind == "koopman_validation") return "Validation";
    if (item.kind == "koopman_spectrum") return "Spectrum";
    if (item.kind == "decomposition") return "Result";
    return "Data";
}

std::string stableId(const std::string& key) {
    std::uint64_t hash = 1469598103934665603ull;
    for (const unsigned char value : key) { hash ^= value; hash *= 1099511628211ull; }
    std::ostringstream out; out << "presentation-" << std::hex << hash;
    return out.str();
}

std::string orderKey(const CatalogItem& item, const std::string& originalLabel) {
    const auto created = property(item, "created UTC");
    if (!created.empty()) return created + item.id;
    return item.source + originalLabel + item.id;
}

bool orderedBefore(const CatalogItem& left, const std::string& leftLabel,
                   const CatalogItem& right, const std::string& rightLabel) {
    const auto leftCreated = property(left, "created UTC");
    const auto rightCreated = property(right, "created UTC");
    if (!leftCreated.empty() && !rightCreated.empty() && leftCreated != rightCreated)
        return leftCreated < rightCreated;
    if (left.modified > 0.0 && right.modified > 0.0 && left.modified != right.modified)
        return left.modified < right.modified;
    return orderKey(left, leftLabel) < orderKey(right, rightLabel);
}

std::string archiveAttempt(const CatalogItem& item) {
    const auto text = lower(item.source);
    const auto marker = text.rfind("attempt-");
    if (marker == std::string::npos) return {};
    std::string number;
    for (std::size_t i = marker + 8; i < text.size(); ++i) {
        if (std::isdigit(static_cast<unsigned char>(text[i]))) number += text[i];
        else break;
    }
    return number.empty() ? std::string{} : "Search attempt #" + number;
}

}  // namespace

void applyCatalogPresentation(Catalog& catalog) {
    if (catalog.items.empty() || std::any_of(catalog.items.begin(), catalog.items.end(),
        [](const auto& item) { return item.kind == "presentation_group"; })) return;
    const bool sourceBacked = std::any_of(catalog.items.begin(), catalog.items.end(),
        [](const auto& item) { return !item.system.empty() &&
            (item.category == "Archived solves" || item.properties.count("native catalog")); });
    if (!sourceBacked) return;

    const std::size_t originalSize = catalog.items.size();
    std::map<std::string, std::size_t> byId;
    std::map<std::string, std::vector<std::size_t>> children;
    std::map<std::string, std::string> originalLabels, originalParents, originalKinds;
    for (std::size_t i = 0; i < originalSize; ++i) {
        auto& item = catalog.items[i];
        byId[item.id] = i; children[item.parentId].push_back(i);
        originalLabels[item.id] = item.label;
        originalParents[item.id] = item.parentId;
        originalKinds[item.id] = item.kind;
        item.properties.emplace("archive label", item.label);
        item.properties.emplace("archive category", item.category);
        if (!item.parentId.empty()) item.properties.emplace("archive parent", item.parentId);
    }

    std::set<std::string> trainingItems;
    static const std::set<std::string> directTrainingRelations{
        "derived from", "fit from trajectory", "extracted from"};
    static const std::set<std::string> reverseTrainingRelations{
        "produced derived result", "used by spectral fit", "has extracted KROM"};
    for (std::size_t i = 0; i < originalSize; ++i) {
        const auto& item = catalog.items[i];
        for (const auto& relation : item.relations) {
            const auto target = byId.find(relation.target);
            if (target == byId.end()) continue;
            const auto& targetItem = catalog.items[target->second];
            if (item.kind == "krom" && targetItem.kind != "krom" &&
                directTrainingRelations.count(relation.kind))
                trainingItems.insert(targetItem.id);
            if (item.kind != "krom" && targetItem.kind == "krom" &&
                reverseTrainingRelations.count(relation.kind))
                trainingItems.insert(item.id);
        }
    }
    const auto directTrainingItems = trainingItems;
    for (const auto& id : directTrainingItems) {
        auto parent = originalParents[id];
        std::set<std::string> visited;
        while (!parent.empty() && visited.insert(parent).second) {
            trainingItems.insert(parent);
            parent = originalParents[parent];
        }
    }

    std::vector<std::string> systems(originalSize), physics(originalSize);
    for (std::size_t i = 0; i < originalSize; ++i) {
        systems[i] = catalog.items[i].system;
        physics[i] = compactPhysics(catalog.items[i]);
        if (physics[i] == "General parameters" &&
            property(catalog.items[i], "physics preset").empty() &&
            property(catalog.items[i], "system preset").empty() &&
            property(catalog.items[i], "preset").empty()) physics[i].clear();
    }
    for (std::size_t pass = 0; pass < originalSize + 1; ++pass) {
        bool changed = false;
        for (std::size_t i = 0; i < originalSize; ++i) {
            const auto parent = byId.find(originalParents[catalog.items[i].id]);
            if (parent != byId.end()) {
                if (systems[i].empty() && !systems[parent->second].empty()) {
                    systems[i] = systems[parent->second]; changed = true;
                }
                if (physics[i].empty() && !physics[parent->second].empty()) {
                    physics[i] = physics[parent->second]; changed = true;
                }
            }
            std::set<std::string> childSystems, childPhysics;
            for (const auto child : children[catalog.items[i].id]) {
                if (!systems[child].empty()) childSystems.insert(systems[child]);
                if (!physics[child].empty()) childPhysics.insert(physics[child]);
            }
            if (systems[i].empty() && childSystems.size() == 1) {
                systems[i] = *childSystems.begin(); changed = true;
            }
            if (physics[i].empty() && childPhysics.size() == 1) {
                physics[i] = *childPhysics.begin(); changed = true;
            }
        }
        if (!changed) break;
    }
    for (std::size_t pass = 0; pass < originalSize + 1; ++pass) {
        bool changed = false;
        for (std::size_t i = 0; i < originalSize; ++i) {
            if (!physics[i].empty()) continue;
            std::set<std::string> candidates;
            for (const auto& relation : catalog.items[i].relations) {
                static const std::set<std::string> inheritsFromTarget{
                    "controlled by KROM", "derived from", "evaluates KROM",
                    "expands training domain of", "extracted from",
                    "fit from trajectory", "initialized from spectrum", "uses KROM"};
                if (inheritsFromTarget.count(relation.kind) == 0) continue;
                const auto target = byId.find(relation.target);
                if (target != byId.end() && !physics[target->second].empty())
                    candidates.insert(physics[target->second]);
            }
            if (candidates.size() == 1) {
                physics[i] = *candidates.begin();
                changed = true;
            }
        }
        if (!changed) break;
    }

    std::map<std::string, int> kromNumber;
    std::map<std::string, std::vector<std::size_t>> kromFamilies;
    for (std::size_t i = 0; i < originalSize; ++i)
        if (catalog.items[i].kind == "krom")
            kromFamilies[compactSystem(systems[i]) + "|" + (physics[i].empty()
                ? "General parameters" : physics[i])].push_back(i);
    for (auto& family : kromFamilies) {
        std::sort(family.second.begin(), family.second.end(), [&](auto a, auto b) {
            return orderedBefore(catalog.items[a], originalLabels[catalog.items[a].id],
                                 catalog.items[b], originalLabels[catalog.items[b].id]);
        });
        for (std::size_t i = 0; i < family.second.size(); ++i)
            kromNumber[catalog.items[family.second[i]].id] = static_cast<int>(i + 1);
    }

    std::map<std::string, std::string> groups;
    auto ensureGroup = [&](const std::string& key, const std::string& label,
                           const std::string& parent, const std::string& system) {
        const auto found = groups.find(key);
        if (found != groups.end()) return found->second;
        CatalogItem group;
        group.id = stableId(key); group.kind = "presentation_group";
        group.category = "Data"; group.label = label; group.system = system;
        group.subtitle = "Derived navigation group; archive records remain immutable";
        group.parentId = parent; group.primary = false;
        group.properties["presentation grouping"] = label;
        groups.emplace(key, group.id); catalog.items.push_back(std::move(group));
        return groups.at(key);
    };

    std::vector<std::string> studies(originalSize), lanes(originalSize);
    std::map<std::string, std::set<std::string>> lanesByPhysics;
    for (std::size_t i = 0; i < originalSize; ++i) {
        const auto itemId = catalog.items[i].id;
        if (!originalParents[itemId].empty()) continue;
        const auto itemSystem = compactSystem(systems[i]);
        const auto physical = physics[i].empty() ? "General parameters" : physics[i];
        studies[i] = activity(catalog.items[i], originalLabels[itemId],
            trainingItems.count(itemId) != 0);
        lanes[i] = majorLane(studies[i]);
        lanesByPhysics[itemSystem + "|" + physical].insert(lanes[i]);
    }

    const auto kromBranchLabel = [&](const std::string& id) {
        const auto number = kromNumber.find(id);
        std::string label = number == kromNumber.end() ? "KROM source not recorded"
            : "KROM #" + std::to_string(number->second);
        const auto source = byId.find(id);
        const auto original = originalLabels.find(id);
        if (source != byId.end() && original != originalLabels.end()) {
            const auto context = shortDescriptor(catalog.items[source->second], original->second);
            if (!context.empty()) label += " · " + context;
        }
        return label;
    };

    for (std::size_t i = 0; i < originalSize; ++i) {
        catalog.items[i].properties["dynamical system"] = catalog.items[i].system;
        catalog.items[i].system = compactSystem(systems[i]);
        catalog.items[i].category = "Data";
        const auto itemId = catalog.items[i].id;
        const auto physical = physics[i].empty() ? "General parameters" : physics[i];
        catalog.items[i].properties["physics"] = physical;
        if (!originalParents[itemId].empty()) continue;
        const auto itemSystem = catalog.items[i].system;
        const auto base = itemSystem + "|" + physical;
        const auto& study = studies[i];
        const auto& lane = lanes[i];
        std::string sourceKrom;
        if (study == "KLQR controls" || study == "Bψ identification") {
            for (const auto& relation : catalog.items[i].relations) {
                const auto target = byId.find(relation.target);
                if (target != byId.end() && catalog.items[target->second].kind == "krom" &&
                    (relation.kind == "controlled by KROM" || relation.kind == "uses KROM" ||
                     relation.kind == "expands training domain of" || relation.kind == "derived from")) {
                    sourceKrom = relation.target; break;
                }
            }
        }
        const auto controlTarget = property(catalog.items[i], "control target");
        const auto strategies = strategyLabels(catalog.items[i], originalLabels[itemId], study);

        const auto physicsParent = ensureGroup(base + "|physics", physical, {}, itemSystem);
        static const std::vector<std::string> laneOrder{
            "Training data", "Models & assets", "Control & estimation", "Other data"};
        for (const auto& candidate : laneOrder)
            if (lanesByPhysics[base].count(candidate))
                (void)ensureGroup(base + "|lane|" + candidate,
                    candidate, physicsParent, itemSystem);
        auto parent = groups.at(base + "|lane|" + lane);
        if (lane == "Training data") {
            const auto mode = trainingMode(catalog.items[i], originalLabels[itemId]);
            parent = ensureGroup(base + "|lane|" + lane + "|mode|" + mode,
                mode, parent, itemSystem);
        } else {
            parent = ensureGroup(base + "|study|" + study, study, parent, itemSystem);
        }

        if (catalog.items[i].kind == "krom") {
            const auto context = shortDescriptor(catalog.items[i], originalLabels[itemId]);
            if (!context.empty()) parent = ensureGroup(base + "|krom-method|" + context,
                context, parent, itemSystem);
        }

        if (study == "KLQR controls" || study == "Bψ identification") {
            parent = ensureGroup(base + "|study|" + study + "|krom|" + sourceKrom,
                kromBranchLabel(sourceKrom), parent, itemSystem);
        }
        if (study == "Bψ identification") {
            const auto attempt = archiveAttempt(catalog.items[i]);
            if (!attempt.empty()) parent = ensureGroup(parent + "|attempt|" + attempt,
                attempt, parent, itemSystem);
        }
        if (study == "KLQR controls") {
            const auto target = targetLabel(controlTarget);
            if (!target.empty()) parent = ensureGroup(parent + "|target|" +
                controlTarget, target, parent, itemSystem);
        }
        for (const auto& strategy : strategies)
            parent = ensureGroup(base + "|strategy|" + parent + "|" + strategy,
                strategy, parent, itemSystem);
        catalog.items[i].parentId = parent;
    }

    std::map<std::string, std::vector<std::size_t>> siblings;
    for (std::size_t i = 0; i < originalSize; ++i)
        siblings[catalog.items[i].parentId].push_back(i);
    for (auto& family : siblings) {
        std::sort(family.second.begin(), family.second.end(), [&](auto a, auto b) {
            return orderedBefore(catalog.items[a], originalLabels[catalog.items[a].id],
                                 catalog.items[b], originalLabels[catalog.items[b].id]);
        });
        for (std::size_t order = 0; order < family.second.size(); ++order) {
            auto& item = catalog.items[family.second[order]];
            const auto parentKind = originalKinds[originalParents[item.id]];
            const auto labelNoun = noun(item, parentKind, originalLabels[item.id],
                trainingItems.count(item.id) != 0 || lanes[family.second[order]] == "Training data");
            const auto knownKrom = kromNumber.find(item.id);
            const auto sequence = item.kind == "krom" && knownKrom != kromNumber.end()
                ? knownKrom->second : static_cast<int>(order + 1);
            item.label = labelNoun + " #" + std::to_string(sequence);
            item.properties["presentation order"] = std::to_string(sequence);
        }
    }

    std::map<std::string, CatalogItem*> finalById;
    for (auto& item : catalog.items) finalById[item.id] = &item;
    for (auto& item : catalog.items)
        for (auto& relation : item.relations) {
            const auto target = finalById.find(relation.target);
            if (target != finalById.end()) relation.label = target->second->label;
        }
    for (auto& item : catalog.items) if (item.kind == "presentation_group") {
        std::size_t count = 0;
        for (std::size_t i = 0; i < originalSize; ++i) {
            std::string parent = catalog.items[i].parentId;
            std::set<std::string> visited;
            while (!parent.empty() && visited.insert(parent).second) {
                if (parent == item.id) { ++count; break; }
                const auto found = finalById.find(parent);
                parent = found == finalById.end() ? std::string{} : found->second->parentId;
            }
        }
        item.properties["item count"] = std::to_string(count);
    }
}

}  // namespace dsl
