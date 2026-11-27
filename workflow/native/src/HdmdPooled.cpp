// Pooled (ensemble) H-DMD: one Hankel DMD over many trajectories of one system. Each trajectory
// contributes its own window pairs (X, Y); the pooled matrices are their columns side by side, so no
// window and no shifted pair ever crosses from one trajectory into the next (Arbabi & Mezic, SIAM J.
// Appl. Dyn. Syst. 16:2096, 2017; ensemble DMD over independent runs: Tu et al., J. Comput. Dyn.
// 1:391, 2014, Sec. 3). Every trajectory uses the same settings and channels.
//
// The krom's samples are the windows of all trajectories in order: training_sample_indices are rows
// of the input records concatenated in origin order, training_time each sample's own time, and the
// interpolation is on the physical states at the samples (each record's state source).
#include "HankelAnalysis.hpp"
#include "HdmdObjects.hpp"
#include "SampledDomain.hpp"
#include "TimeSeriesSource.hpp"
#include <dynsyslab/core/FunctionInterpolation.hpp>
#include <dynsyslab/core/HankelRepresentations.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Reconstruction.hpp>
#include <dynsyslab/core/Spectral.hpp>

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

struct Member {
    SourceData view, windowed;
    std::size_t first{};      // first kept row of the view
    HankelAnalysis hankel;
};

Member prepare(const SourceData& source, const J& settings) {
    Member m;
    m.view = joinedTimeSeries(source, settings);
    m.windowed = hankelTimeWindow(m.view, settings, m.first);
    m.hankel = prepareHankel(m.windowed, settings, true, true);
    return m;
}

}  // namespace

std::filesystem::path runPooledHdmd(const SourceData& full, const J& settings, const std::filesystem::path& root,
                                    const std::string& label) {
    std::vector<Member> members;
    members.push_back(prepare(full, settings));
    for (const auto& path : settings.at("additional_sources").array())
        members.push_back(prepare(loadSource(std::filesystem::u8path(path.string())), settings));
    const auto& lead = members.front().hankel;
    // Pool the pairs: columns of each trajectory, starts/ends as rows of the concatenated records.
    c::HankelMatrix pooled;
    pooled.window = lead.matrix.window;
    pooled.channels = lead.matrix.channels;
    pooled.rows = lead.matrix.rows;
    std::vector<std::size_t> recordOffset;  // first concatenated row of each member's view
    std::size_t offset = 0;
    for (const auto& m : members) {
        if (m.hankel.matrix.rows != pooled.rows || m.hankel.names != lead.names || m.hankel.alignEnd != lead.alignEnd)
            throw std::invalid_argument("pooled hdmd: every trajectory must give the same Hankel layout and channels");
        recordOffset.push_back(offset);
        pooled.columns += m.hankel.matrix.columns;
        for (std::size_t j = 0; j < m.hankel.matrix.columns; ++j) {
            pooled.starts.push_back(offset + m.first + m.hankel.matrix.starts[j]);
            pooled.ends.push_back(offset + m.first + m.hankel.matrix.ends[j]);
        }
        offset += m.view.samples();
    }
    pooled.sourceSamples = offset;
    pooled.x.resize(pooled.rows * pooled.columns);
    pooled.y.resize(pooled.rows * pooled.columns);
    for (std::size_t row = 0, at = 0; row < pooled.rows; ++row) {
        for (const auto& m : members) {
            const auto& h = m.hankel.matrix;
            std::copy(h.x.begin() + std::ptrdiff_t(row * h.columns), h.x.begin() + std::ptrdiff_t((row + 1) * h.columns),
                      pooled.x.begin() + std::ptrdiff_t(at));
            std::copy(h.y.begin() + std::ptrdiff_t(row * h.columns), h.y.begin() + std::ptrdiff_t((row + 1) * h.columns),
                      pooled.y.begin() + std::ptrdiff_t(at));
            at += h.columns;
        }
    }
    const double dt = uniformDeltaTime(members.front().windowed);
    for (const auto& m : members)
        if (std::abs(uniformDeltaTime(m.windowed) - dt) > 1e-9 * dt)
            throw std::invalid_argument("pooled hdmd: every trajectory must share one sample interval");
    c::HankelDMDOptions options;
    options.matrixSize = pooled.columns;
    options.sampleStride = pooled.window.sampleStride;
    options.snapshotStride = pooled.window.snapshotStride;
    options.alignToWindowEnd = lead.alignEnd;
    options.channelScales = lead.scales;
    options.rank = std::size_t(settings.at("rank").number(0));
    options.singularValueTolerance = settings.at("svd_tolerance").number(1.0e-10);
    options.optimalHardThreshold = settings.at("rank_rule").string("tolerance") == "optimal_hard_threshold";
    // The terminal amplitude belongs to one trajectory's last window; pooled modes keep DMD order unless asked.
    options.terminalAmplitudeEnergySort = settings.at("mode_sort").string("none") != "none";
    const auto result = c::fitHankelDMD(std::move(pooled), lead.names, dt, options);

    // Samples: each pooled column's label row, mapped back to its trajectory.
    c::EigenfunctionSamples samples;
    samples.modes = result.rank;
    samples.values = result.coordinates;
    std::vector<double> states;
    std::size_t stateDimension = 0;
    std::vector<std::vector<double>> times(members.size());
    std::vector<std::size_t> owner;
    for (const auto row : result.coordinateSampleOffsets) {
        std::size_t member = std::upper_bound(recordOffset.begin(), recordOffset.end(), row) - recordOffset.begin() - 1;
        const auto local = row - recordOffset[member];
        samples.sampleIndices.push_back(row);
        samples.time.push_back(members[member].view.time().at(local));
        times[member].push_back(samples.time.back());
        owner.push_back(member);
    }
    std::vector<std::vector<double>> memberStates(members.size());
    for (std::size_t i = 0; i < members.size(); ++i) {
        if (times[i].empty()) continue;
        const auto found = statesAtSamples(members[i].view.record, times[i]);
        if (!found) throw std::invalid_argument("pooled hdmd interpolates on the physical state; a trajectory has no producing system");
        stateDimension = found->dimension;
        memberStates[i] = found->values;
    }
    std::vector<std::size_t> taken(members.size(), 0);
    for (const auto member : owner) {
        const auto k = taken[member]++;
        states.insert(states.end(), memberStates[member].begin() + std::ptrdiff_t(k * stateDimension),
                      memberStates[member].begin() + std::ptrdiff_t((k + 1) * stateDimension));
    }
    std::size_t neighbors = std::size_t(settings.at("interpolation_neighbors").number(8));
    const auto method = settings.at("interpolation_method").isNull() ? c::InterpolationMethod{}
        : c::interpolationMethodFromJson(settings.at("interpolation_method"), &neighbors);
    std::shared_ptr<const c::FunctionRepresentation> interpolation = c::interpolateSamples(states, samples.samples(),
        c::FunctionDomain::physicalState(stateDimension), samples.values, method, neighbors, "pooled hdmd eigenfunctions at their samples");
    const std::size_t labelOffset = lead.alignEnd ? lead.matrix.ends.at(0) - lead.matrix.starts.at(0) : 0;
    std::shared_ptr<const c::FunctionRepresentation> window =
        c::windowMapFromModes(recordWindowLayout(lead, members.front().view.record, labelOffset), result.modes, result.modeGram, result.rank);
    c::DataOrigin origin{"hdmd", settings, {}};
    for (const auto& m : members)
        origin.inputs.push_back(dataReference(m.view.path, "time_series"));
    auto spectrum = result.spectrum;
    spectrum.withOrigin(origin);
    c::KoopmanModel krom(std::move(spectrum), samples, origin, {interpolation, window});
    const auto channels = result.observableNames.size();
    for (std::size_t channel = 0; channel < channels; ++channel)  // C of each channel the matrix read
        krom.setReconstruction(result.observableNames[channel],
            {std::vector<std::complex<double>>(result.reconstructionMatrix.begin() + channel * result.rank,
                                               result.reconstructionMatrix.begin() + (channel + 1) * result.rank), 1, result.rank});
    if (!settings.at("reconstruction_observables").array().empty()) {
        auto system = c::createSystem(systemConfiguration(members.front().view));
        const auto catalog = system->observableCatalog();
        for (const auto& item : settings.at("reconstruction_observables").array()) {
            const auto& observable = catalog.at(item.string());
            const auto real = observable.evaluateRealBatch(states, samples.samples());
            krom.setReconstruction(item.string(), c::fitObservableReconstruction(result.coordinates, samples.samples(), result.rank,
                std::vector<std::complex<double>>(real.begin(), real.end()), observable.outputDimension(),
                settings.at("reconstruction_ridge").number()));
        }
    }
    const auto errors = krom.errors();
    // Each eigenpair carries H-DMD's one-step residual of the evaluated representation along the training
    // trajectory, ||psi(x_k+1) - mu psi(x_k)|| / ||psi(x_k)|| (core/Eigenpair.hpp, "one-step"; Arbabi &
    // Mezic, SIAM J. Appl. Dyn. Syst. 16, 2017; the residual of Colbrook & Townsend 2024).
    for(std::size_t mode=0;mode<errors.koopmanResidual.size();++mode)
        krom.setResidual(mode,c::EigenpairResidual{errors.koopmanResidual[mode],c::residual_definition::oneStep,
            c::residual_scope::empiricalTraining});
    std::map<std::string, c::NpyArray> arrays;
    arrays["sample_interval"] = c::NpyArray::fromDoubles({1}, {result.deltaTime});
    arrays["eigenfunction_value_error"] = c::NpyArray::fromDoubles({result.rank}, errors.valueRelative);
    arrays["reconstruction_matrix"] = c::NpyArray::fromComplex({channels, result.rank}, result.reconstructionMatrix);
    arrays["hankel_singular_values"] = c::NpyArray::fromDoubles({result.singularValues.size()}, result.singularValues);
    std::vector<double> counts;
    for (const auto& m : members) counts.push_back(double(m.hankel.matrix.columns));
    arrays["trajectory_columns"] = c::NpyArray::fromDoubles({counts.size()}, counts);
    auto meta = json::parse(metadata(members.front().view, "hankel-dmd", label, ""));
    meta.object()["pooled"] = J::Object{{"trajectories", double(members.size())}, {"columns", double(result.matrixColumns)},
        {"boundaries_preserved", true}, {"sample_indices", "rows of the input records concatenated in origin order"}};
    meta.object()["hankel_window"] = lead.provenance;
    meta.object()["rank"] = double(result.rank);
    meta.object()["delta_t"] = result.deltaTime;
    meta.object()["eigenvalue_domain"] = "discrete-time";
    meta.object()["observables"] = json::parse(jsonStringArray(result.observableNames));
    meta.object()["eigenfunction_errors"] = J::Object{{"representation", errors.representation}, {"domain", errors.domain},
        {"samples", double(errors.samples)}};
    return publishKrom(root, members.front().view, krom, "koopman-hankel-dmd-pooled", json::dump(meta), std::move(arrays));
}

}  // namespace dynsyslab::workflow::postprocess_detail
