// Post-processing methods that project orbit records onto a system's known spectrum (OrbitAverages.hpp).
//
// band_projection (pooled over its sources, one orbit record each): each record's Fourier coefficients
// in the library band's exact phase, and from them the family's spectral measure of each observable,
// smoothed by the rational circle kernel of Rigged DMD (Colbrook, Drysdale & Horning, arXiv:2405.00782),
//   rho_g(theta) = sum_j w_j sum_n |a_{j,n}|^2 K_eps(theta - n Omega_j dt)   (math::smoothOrbitSpectra),
// and the wave-packet amplitude of each record at each packet centre,
//   A_g(Omega_j) = [sum_n |a_{j,n}|^2 |K_eps(theta - n Omega_j dt)|^2]^(1/2)  (math::spectralPacket).
//   {"name": "band_projection", "settings": {
//      "band": "phase", "observables": ["echo"], "harmonics": 255,
//      "weights": "uniform" | [..] | {"canonical_area": {"invariant": "energy", "band": "phase"}},
//      "invariant": "energy",                        recorded per orbit (optional)
//      "sample_interval": 0.1, "angles": [...] | {"count": 1024}, "bandwidth": 0.05, "order": 6,
//      "packet_centers": [0.5],
//      "krom": "...", "measures": [{"observable": "echo", "as": "coordinates", "output": 150}]}}
// A report: orbit_frequency, orbit_weight, orbit_invariant, fourier_harmonics, fourier_coefficients
// (orbits x observables x harmonics), angles, reference_density__<obs>, and, with orbits ordered by
// frequency, reference_packet_rms__<obs>__<k>. With "krom", its child also gains the reference as a band
// (method "orbit-fourier-reference") of measures named as that krom names the observables, and the
// report's metadata "density_errors" gives the relative L2 error, on the common angle grid, of each
// density that krom already has of those outputs (e.g. Rigged DMD's) against the reference.
//
// packet_validation (pooled; held-out orbit records): a fitted krom's Rigged DMD wave packets evaluated
// along each record's delay windows (the fit's own Hankel settings, checked against its coordinate
// contract), against the band's reference packets.
//   {"name": "packet_validation", "settings": {
//      "krom": "...", "observables": [0], "window": {...},   the fit's Hankel sampling
//      "band": "phase", "harmonics": 255,
//      "packets": [{"observable": "echo", "measure": {"band": "packets", "observable": "coordinates", "output": 150}}]}}
// A report: validation_frequency (sorted), reference_packet_rms__<obs>__<k>, packet_rms__<obs>__<k>.
#include "OrbitAverages.hpp"
#include "HankelAnalysis.hpp"
#include "KoopmanSpectralSettings.hpp"
#include "ObjectPublishing.hpp"

#include <dynsyslab/core/FunctionExpansion.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/math/SpectralFilters.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

std::filesystem::path writeReport(const std::filesystem::path& root, c::NpzRecord record, const std::string& stem,
                                  const std::string& method, const std::vector<std::filesystem::path>& inputs) {
    const auto directory = kindFolder(root, c::object_kind::report);
    std::filesystem::create_directories(directory);
    for (std::size_t suffix = 0; suffix < 10000; ++suffix) {
        const auto path = directory / c::wide(safeName(stem) + (suffix ? "__" + std::to_string(suffix) : "") + ".npz");
        if (std::filesystem::exists(path)) continue;
        stampDataObject(record, c::object_kind::report, method, inputs);
        c::writeNpzAtomicNew(path, record);
        return path;
    }
    throw std::runtime_error("unable to allocate a " + method + " archive path");
}

std::vector<std::string> strings(const J& list) {
    std::vector<std::string> out;
    for (const auto& v : list.array()) out.push_back(v.string());
    return out;
}

double packetPower(const math::FourierSeries& s) {
    double total = 0;
    for (const auto& [n, a] : s.coefficients) total += std::norm(a);
    return total;
}

// Orbit indices ordered by frequency.
std::vector<std::size_t> byFrequency(const std::vector<OrbitSeries>& orbits) {
    std::vector<std::size_t> order(orbits.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return orbits[a].frequency < orbits[b].frequency; });
    return order;
}

}  // namespace

std::filesystem::path runBandProjection(const std::vector<std::filesystem::path>& paths, const J& settings,
                                        const std::filesystem::path& root, const std::string& label) {
    if (paths.empty()) throw std::invalid_argument("band_projection needs orbit records");
    std::vector<SourceData> records;
    for (const auto& path : paths) records.push_back(loadSource(path));
    const auto band = recordBand(records.front(), settings.at("band").string());
    const auto observables = strings(settings.at("observables"));
    if (observables.empty()) throw std::invalid_argument("band_projection needs observables");
    const int harmonics = int(settings.at("harmonics").number(64));
    std::vector<OrbitSeries> orbits;
    for (const auto& r : records) orbits.push_back(orbitSeries(r, band, observables, harmonics, settings.at("invariant").string()));
    const auto weights = recordWeights(settings.at("weights"), records);
    const auto angles = spectral_settings::angleGrid(settings.at("angles"));
    const double dt = settings.at("sample_interval").number(), epsilon = settings.at("bandwidth").number();
    const auto order = spectral_settings::count(settings.at("order"), 1);
    if (!(dt > 0) || !(epsilon > 0) || !order || order > 8) throw std::invalid_argument("band_projection needs a sample interval, bandwidth and order 1..8");
    const math::RationalCircleKernel kernel(epsilon, order);
    std::vector<double> centers;
    if (!settings.at("packet_centers").isNull()) centers = spectral_settings::numbers(settings.at("packet_centers"));

    const std::size_t n = orbits.size(), m = observables.size(), H = 2 * std::size_t(harmonics) + 1;
    c::NpzRecord report;
    std::vector<double> frequency, invariant, harmonicsList;
    std::vector<std::complex<double>> coefficients;
    for (const auto& o : orbits) {
        frequency.push_back(o.frequency);
        invariant.push_back(o.invariant);
        for (const auto& s : o.channels) for (const auto& [k, a] : s.coefficients) coefficients.push_back(a);
    }
    for (int k = -harmonics; k <= harmonics; ++k) harmonicsList.push_back(k);
    report.arrays["orbit_frequency"] = c::NpyArray::fromDoubles({n}, frequency);
    report.arrays["orbit_weight"] = c::NpyArray::fromDoubles({n}, weights);
    report.arrays["orbit_invariant"] = c::NpyArray::fromDoubles({n}, invariant);
    report.arrays["fourier_harmonics"] = c::NpyArray::fromDoubles({H}, harmonicsList);
    report.arrays["fourier_coefficients"] = c::NpyArray::fromComplex({n, m, H}, coefficients);
    report.arrays["angles"] = c::NpyArray::fromDoubles({angles.size()}, angles);
    const auto ordered = byFrequency(orbits);
    std::vector<double> sortedFrequency;
    for (const auto j : ordered) sortedFrequency.push_back(orbits[j].frequency);
    report.arrays["validation_frequency"] = c::NpyArray::fromDoubles({n}, sortedFrequency);
    J::Array masses;
    std::vector<c::SpectralMeasure> measures;
    for (std::size_t i = 0; i < m; ++i) {
        std::vector<math::WeightedOrbitSpectrum> spectra;
        for (std::size_t j = 0; j < n; ++j) spectra.push_back({orbits[j].channels[i], weights[j]});
        const auto smoothed = math::smoothOrbitSpectra(spectra, dt, angles, kernel);
        const auto key = safeArrayKey(observables[i]);
        report.arrays["reference_density__" + key] = c::NpyArray::fromDoubles({angles.size()}, smoothed.density);
        masses.emplace_back(J::Object{{"observable", observables[i]}, {"total_mass", smoothed.totalMass},
                                      {"zero_harmonic_mass", smoothed.zeroHarmonicMass}});
        for (std::size_t k = 0; k < centers.size(); ++k) {
            std::vector<double> rms;
            for (const auto j : ordered) rms.push_back(std::sqrt(packetPower(math::spectralPacket(orbits[j].channels[i], dt, centers[k], kernel))));
            report.arrays["reference_packet_rms__" + key + "__" + std::to_string(k)] = c::NpyArray::fromDoubles({n}, rms);
        }
        c::SpectralMeasure measure;
        measure.observable = observables[i];
        measure.density = smoothed.density;
        measure.totalMass = smoothed.totalMass;
        measure.assumptions = "Zero-harmonic (invariant) mass " + std::to_string(smoothed.zeroHarmonicMass) + ".";
        measures.push_back(std::move(measure));
    }
    std::vector<std::filesystem::path> inputs(paths.begin(), paths.end());
    J::Array densityErrors;  // each fitted band density of a compared observable against the reference
    // The reference as a band of a krom's spectrum, beside the krom's own measures of the same observables.
    if (const auto kromPath = std::filesystem::u8path(settings.at("krom").string()); !kromPath.empty()) {
        auto model = c::loadKoopmanModel(kromPath);
        c::SpectralBand reference;
        reference.name = settings.at("reference_band").string("fourier-reference");
        reference.method = "orbit-fourier-reference";
        reference.grid = angles;
        reference.bandwidth = epsilon;
        reference.order = order;
        reference.assumptions = "The family's spectral measure from weighted Birkhoff averages along its orbit records in the "
                                "catalogue band " + band.name + "'s exact phase; not an estimate from this krom's operator.";
        auto spectrum = model.spectrum();
        spectrum.addBand(reference);
        for (const auto& item : settings.at("measures").array()) {
            const auto name = item.at("observable").string();
            const auto at = std::find(observables.begin(), observables.end(), name);
            if (at == observables.end()) throw std::invalid_argument("band_projection: measure of an observable not projected: " + name);
            auto measure = measures[std::size_t(at - observables.begin())];
            measure.observable = item.at("as").string(name);
            measure.output = std::size_t(item.at("output").number(0));
            measure.band = reference.name;
            // Relative L2 error on the common grid of every measure the krom already has of this output.
            for (const auto& fitted : spectrum.measures()) {
                const auto* other = spectrum.band(fitted.band);
                if (fitted.observable != measure.observable || fitted.output != measure.output || !other
                    || other->grid.size() != angles.size() || fitted.density.size() != angles.size()) continue;
                bool same = true;
                for (std::size_t k = 0; k < angles.size() && same; ++k) same = std::abs(other->grid[k] - angles[k]) <= 1e-12;
                if (!same) continue;
                double error = 0, norm = 0;
                for (std::size_t k = 0; k < angles.size(); ++k) {
                    error += std::pow(fitted.density[k] - measure.density[k], 2); norm += std::pow(measure.density[k], 2);
                }
                densityErrors.emplace_back(J::Object{{"observable", name}, {"band", fitted.band},
                    {"density_relative_l2", norm > 0 ? J(std::sqrt(error / norm)) : J()},
                    {"mass", fitted.totalMass}, {"reference_mass", measure.totalMass}});
            }
            spectrum.addMeasure(std::move(measure));
        }
        model.setSpectralEvidence(spectrum);
        const J meta = J::Object{{"schema", "dynsyslab.derived/v1"}, {"method", "band_projection"}, {"label", label},
                                 {"source_krom", kromPath.generic_u8string()}};
        const auto directory = kindFolder(root, c::object_kind::krom);
        std::filesystem::create_directories(directory);
        const auto path = directory / c::wide(safeName(kromPath.stem().u8string() + "__band-projection__" + timestamp()) + ".npz");
        auto withInputs = inputs;
        withInputs.push_back(kromPath);
        c::saveKoopmanModel(path, model, withDataObject(json::dump(meta), c::object_kind::krom, withInputs));
        registerAdditionalObject(path);
    }
    const J meta = J::Object{{"schema", "dynsyslab.derived/v1"}, {"method", "band_projection"}, {"label", label},
        {"band", band.name}, {"observables", settings.at("observables")}, {"harmonics", double(harmonics)},
        {"sample_interval", dt}, {"bandwidth", epsilon}, {"order", double(order)}, {"masses", J(std::move(masses))},
        {"density_errors", J(std::move(densityErrors))},
        {"coefficient_source", "weighted Birkhoff averages along the records in the band's exact phase (Das & Yorke 2018)"}};
    report.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(meta));
    return writeReport(root, std::move(report), paths.front().stem().u8string() + "__band-projection__" + timestamp(), "band-projection", inputs);
}

std::filesystem::path runPacketValidation(const std::vector<std::filesystem::path>& paths, const J& settings,
                                          const std::filesystem::path& root, const std::string& label) {
    if (paths.empty()) throw std::invalid_argument("packet_validation needs held-out orbit records");
    const auto kromPath = std::filesystem::u8path(settings.at("krom").string());
    const auto model = c::loadKoopmanModel(kromPath);
    if (!model.compression()) throw std::invalid_argument("packet_validation needs a spectral fit's krom");
    const auto contract = json::parse(model.compression()->coordinateContract);
    std::vector<SourceData> records;
    for (const auto& path : paths) records.push_back(loadSource(path));
    const auto band = recordBand(records.front(), settings.at("band").string());
    const int harmonics = int(settings.at("harmonics").number(64));
    std::vector<std::string> names;
    for (const auto& item : settings.at("packets").array()) names.push_back(item.at("observable").string());
    std::vector<OrbitSeries> orbits;
    for (const auto& r : records) orbits.push_back(orbitSeries(r, band, names, harmonics));
    const auto ordered = byFrequency(orbits);
    // The fit's delay windows along each record, from the fit's own Hankel settings.
    auto sampling = settings;
    sampling.object().erase("krom");
    std::vector<HankelAnalysis> windows;
    for (const auto& r : records) {
        windows.push_back(prepareHankel(r, sampling, false));
        const auto& p = windows.back().provenance;
        if (p.at("delays").number() != contract.at("delays").number() || p.at("rows").number() != contract.at("rows").number()
            || json::dump(p.at("channels")) != json::dump(contract.at("channels")))
            throw std::invalid_argument("packet_validation: these Hankel settings do not rebuild the fit's delay coordinates");
    }
    c::NpzRecord report;
    std::vector<double> frequency;
    for (const auto j : ordered) frequency.push_back(orbits[j].frequency);
    report.arrays["validation_frequency"] = c::NpyArray::fromDoubles({frequency.size()}, frequency);
    J::Array errors;
    for (std::size_t i = 0; i < names.size(); ++i) {
        const auto& spec = settings.at("packets").array()[i].at("measure");
        const auto* measure = model.spectrum().measure(spec.at("band").string(), spec.at("observable").string(),
                                                       std::size_t(spec.at("output").number(0)));
        if (!measure || !measure->packets) throw std::invalid_argument("packet_validation: the krom has no Rigged packets for " + names[i]);
        const auto* grid = model.spectrum().band(measure->band);
        const math::RationalCircleKernel kernel(grid->bandwidth, grid->order);
        const double dt = model.compression()->sampleInterval;
        const auto key = safeArrayKey(names[i]);
        for (std::size_t k = 0; k < grid->grid.size(); ++k) {
            std::vector<double> reference, fitted;
            for (const auto j : ordered) {
                reference.push_back(std::sqrt(packetPower(math::spectralPacket(orbits[j].channels[i], dt, grid->grid[k], kernel))));
                const auto& h = windows[j].matrix;
                double total = 0;
                std::vector<double> column(h.rows);
                for (std::size_t col = 0; col < h.columns; ++col) {
                    for (std::size_t row = 0; row < h.rows; ++row) column[row] = h.x[row * h.columns + col];
                    total += std::norm(measure->packets->evaluate(column.data())[k]);
                }
                fitted.push_back(std::sqrt(total / double(h.columns)));
            }
            double error = 0, norm = 0;
            for (std::size_t j = 0; j < reference.size(); ++j) { error += std::pow(fitted[j] - reference[j], 2); norm += reference[j] * reference[j]; }
            errors.emplace_back(J::Object{{"observable", names[i]}, {"center", grid->grid[k]},
                                          {"amplitude_relative_l2", norm > 0 ? J(std::sqrt(error / norm)) : J()}});
            report.arrays["reference_packet_rms__" + key + "__" + std::to_string(k)] = c::NpyArray::fromDoubles({reference.size()}, reference);
            report.arrays["packet_rms__" + key + "__" + std::to_string(k)] = c::NpyArray::fromDoubles({fitted.size()}, fitted);
        }
    }
    const J meta = J::Object{{"schema", "dynsyslab.derived/v1"}, {"method", "packet_validation"}, {"label", label},
        {"source_krom", kromPath.generic_u8string()}, {"band", band.name}, {"errors", J(std::move(errors))}};
    report.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(meta));
    std::vector<std::filesystem::path> inputs(paths.begin(), paths.end());
    inputs.push_back(kromPath);
    return writeReport(root, std::move(report), paths.front().stem().u8string() + "__packet-validation__" + timestamp(), "packet-validation", inputs);
}

}  // namespace dynsyslab::workflow::postprocess_detail
