#include "dynsyslab_viewer/Workbench.hpp"
namespace dsl {
void Workbench::requestPhaseContext() {
    if (phaseCancellation_) {
        phaseCancellation_->store(true);
        return;
    }
    if (taskActive_ || observablesActive_ || plots_.isKrom() || mainSources_.empty())
        return;
    auto state = plots_.captureState();
    if (!state.payload)
        return;
    const auto domain = plots_.phaseSpaceDomain();
    if (!domain)
        return;
    const auto sources = mainSources_;
    const auto selected = selected_;
    phaseCancellation_ = std::make_shared<std::atomic_bool>(false);
    const auto cancelled = phaseCancellation_;
    taskActive_ = true;
    updateExportButtons();
    task_ = std::async(std::launch::async, [this, state, domain, sources, selected, cancelled]() {
        AsyncResult result;
        result.itemIndex = selected;
        result.phaseContext = true;
        result.sources = sources;
        try {
            result.payload = std::make_shared<Payload>(backend_.applyPhaseSpace(
                *state.payload, sources, "compute", domain, cancelled.get()));
        } catch (const std::exception &error) {
            result.payload = std::make_shared<Payload>(*state.payload);
            if (!cancelled->load())
                result.payload->warnings.push_back(std::string("Phase context: ") + error.what());
        }
        return result;
    });
}
} // namespace dsl
