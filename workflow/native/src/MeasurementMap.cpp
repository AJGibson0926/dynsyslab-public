#include "MeasurementMap.hpp"
#include <dynsyslab/core/HistoryObservationArchive.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;

c::FrozenHistoryTransform historyTransform(const c::HankelBasis& basis, double sampleInterval,
                                           c::HankelCoordinateScaling scaling) {
    const auto& layout = basis.layout();
    if (layout.channels.size() != 1 || !layout.alignEnd || !basis.mean().empty() || layout.window.delayStride != 1
        || layout.channels.front().lag != 0)
        throw std::invalid_argument("a causal measurement map needs one channel, end-aligned unit-stride windows "
                                    "and no centering");
    c::PooledHankelBasis pooled;
    pooled.delays = layout.window.delays;
    pooled.rank = basis.rank();
    pooled.leftSingularVectors = basis.modes();  // delays x rank, oldest delay first
    const double scale = layout.channels.front().scale;  // the window was divided by it before projection
    for (std::size_t k = 0; k < basis.rank(); ++k)
        pooled.singularValues.push_back(scale * (scaling == c::HankelCoordinateScaling::whitened ? basis.singularValues()[k] : 1.0));
    return c::FrozenHistoryTransform(std::move(pooled), sampleInterval, {},
                                     layout.channelObservables.empty() ? "" : layout.channelObservables.front());
}

c::FrozenHistoryTransform readMeasurementMap(const std::filesystem::path& path) {
    const auto archive = c::readNpz(path);
    const auto meta = archive.metadataJson().empty() ? json::Value{} : json::parse(archive.metadataJson());
    const auto& block = meta.at("hankel_basis");
    if (block.isNull()) return c::readFrozenHistoryTransform(archive);
    const auto basis = c::hankelBasisFromArchive(block, archive.at("basis_modes").doubles(),
                                                 archive.at("singular_values").doubles(), archive.at("mean").doubles());
    return historyTransform(basis, block.at("sample_interval").number(),
                            c::hankelCoordinateScaling(block.at("coordinate_scaling").string()));
}

MeasurementMap loadMeasurementMap(const c::KoopmanModel& model, const std::filesystem::path& kromPath,
                                  const json::Value& settings) {
    const auto& domain = model.coordinates();
    if (!domain.computableFromMeasurements())
        throw std::invalid_argument("this K-ROM's domain has no measurement map");
    const auto path = settings.at("history_transform").string(
        domain.providerArtifact.empty() ? kromPath.u8string() : domain.providerArtifact);
    MeasurementMap out{path, readMeasurementMap(std::filesystem::u8path(path))};
    if (domain.providerIdentity.rfind("sha256:", 0) == 0) {
        if (domain.providerIdentity != out.transform.identity())
            throw std::invalid_argument("measurement map differs from the K-ROM's training coordinates");
    } else if (settings.at("unverified_transform_reason").string().empty()) {
        throw std::invalid_argument("legacy K-ROM has no measurement-map identity; provide unverified_transform_reason");
    }
    return out;
}

}  // namespace dynsyslab::workflow::postprocess_detail
