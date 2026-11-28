// A K-ROM's eigenfunctions as data (KoopmanModel): H-DMD window map, the shared interpolation (the one
// every observable uses, each with its own interpolant), errors, extraction, krom archive round trip,
// child kroms that add a representation by reference to their parent, and time-series joining.
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/HankelRepresentations.hpp>
#include <dynsyslab/core/Spectral.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>
#include <dynsyslab/core/DataObject.hpp>
#include <cstdlib>
#include <filesystem>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::size_t samples30(const c::KoopmanModel& set) { return set.samples().sampleIndices[30]; }

double worst(const std::vector<double>& values) {
    double out = 0;
    for (double value : values) out = std::max(out, value);
    return out;
}

c::TimeSeriesRecord echoRecord(std::size_t samples, double dt, double shift = 0.0) {
    std::vector<double> time(samples), values(samples);
    for (std::size_t k = 0; k < samples; ++k) {
        time[k] = shift + dt * double(k);
        values[k] = std::exp(-0.02 * double(k)) * std::cos(0.3 * double(k))
                  + 0.5 * std::exp(-0.01 * double(k)) * std::cos(1.1 * double(k) + 0.4);
    }
    c::TimeSeriesRecord record{c::TimeSeries(std::move(time), std::move(values), 1), {{"echo", "echo", ""}},
                               {"coordinate-transform", {}, {{"echo.npz", "abc", "time_series"}}}, {}, std::nullopt};
    return record;
}

c::KoopmanModel fit(const c::TimeSeriesRecord& record, c::HankelDMDResult& result, c::HankelWindowLayout& layout) {
    c::HankelDMDOptions options;
    options.matrixSize = 60;
    options.samplingWindow = 12;
    options.alignToWindowEnd = true;
    std::vector<double> values(record.series.valuesData(), record.series.valuesData() + record.series.sampleCount());
    result = c::hankelDynamicModeDecomposition(values, record.series.sampleCount(), std::vector<std::string>{"echo"},
                                               record.series.deltaTime(), options);
    require(result.rank == 4, "two damped tones give four modes");
    layout.window.delays = result.delayCount;
    layout.channels = {{0, 0, 1.0}};
    layout.channelObservables = {"echo"};
    layout.alignEnd = true;
    layout.labelOffset = result.delayCount - 1;
    auto samples = std::make_shared<c::EigenfunctionSamples>();
    samples->modes = result.rank;
    samples->values = result.coordinates;
    // Measured echo, no producing system: the sampled values are the base points.
    auto input = c::FunctionDomain::supplied(1);
    input.names = {"echo"};
    std::vector<double> points;
    for (const auto row : result.coordinateSampleOffsets) {
        samples->sampleIndices.push_back(row);
        samples->time.push_back(record.series.time(row));
        points.push_back(record.series.value(row, 0));
    }
    c::DataOrigin origin{"hdmd", c::json::Value::Object{{"matrix_size", 60.0}}, record.origin.inputs};
    std::shared_ptr<const c::FunctionRepresentation> interpolation =
        c::interpolateSamples(points, points.size(), input, samples->values, {}, 4);
    std::shared_ptr<const c::FunctionRepresentation> window =
        c::windowMapFromModes(layout, result.modes, result.modeGram, result.rank);
    return c::KoopmanModel(c::Spectrum(result.spectrum).withOrigin(origin),
                           *samples, origin, {interpolation, window});
}

}  // namespace

int main() {
    try {
        const auto record = echoRecord(200, 0.1);
        c::HankelDMDResult result;
        c::HankelWindowLayout layout;
        const auto set = fit(record, result, layout);

        // The spectrum pairs conjugates.
        std::size_t paired = 0;
        for (const auto& value : set.spectrum().eigenvalues()) paired += value.conjugate ? 1 : 0;
        require(paired == 4, "spectrum pairs every oscillatory eigenvalue");

        // The window map reproduces H-DMD's own eigenfunction values on every training window.
        const auto window = set.errors("time-delay", &record);
        require(worst(window.valueRelative) < 1e-7, "window map reproduces H-DMD values");
        require(worst(window.koopmanResidual) < 1e-6, "window-map eigenfunctions advance by their eigenvalues");

        // Errors default to the first representation, on its domain.
        const auto defaults = set.errors();
        require(defaults.representation == "interpolation" && defaults.domain == "echo", "errors default to the first");
        require(worst(defaults.valueRelative) < 1e-12, "interpolation is exact at its own samples");

        // Domain recognition by the observables of the data handed in.
        require(set.recognizeDomain({"echo"}) == "echo", "echo data selects the echo domain");
        require(set.recognizeDomain({"echo"}, "time-delay") == "window", "echo windows are the window map's domain");
        bool rejected = false;
        try { (void)set.recognizeDomain({"state", "state"}); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "a domain the representation lacks is refused");

        // Extraction keeps the window map rows of the kept modes.
        const auto pair = set.extract({0, 1});
        const auto h = layout.window_at(record.series, set.samples().sampleIndices[30]);
        const auto full = set.evaluateRepresentation(h.data(), "time-delay");
        const auto kept = pair.evaluateRepresentation(h.data(), "time-delay");
        require(std::abs(kept[0] - full[0]) < 1e-9 * (1 + std::abs(full[0])) && std::abs(kept[1] - full[1]) < 1e-9 * (1 + std::abs(full[1])),
                "extracted window map matches the full one");

        {   // An eigenfunction is an observable: evaluated along a time series it gives a time series. On
            // the window map's domain each sample's window is cut from the series; samples without a
            // complete window are skipped, so row k of the output is window label k + first label.
            // One eigenfunction carries every representation; `as` picks the default.
            const auto all = set.eigenfunction(1);
            require(all.representationNames() == std::vector<std::string>{"interpolation", "time-delay"}
                    && all.representationName() == "interpolation", "an eigenfunction has all its representations");
            const auto psi = all.as("time-delay");
            require(psi.representationName() == "time-delay" && psi.representationNames().size() == 2
                    && psi.eigenvalue().continuous == all.eigenvalue().continuous, "as() keeps the eigenvalue and the others");
            const auto along = psi.evaluate(record.series);
            require(along.componentCount() == 2, "a complex observable gives real and imaginary columns");
            const auto first = layout.labelOffset;
            require(std::abs(along.time(0) - record.series.time(first)) < 1e-12, "the first complete window labels the first row");
            for (std::size_t k = 0; k < set.samples().samples(); k += 7) {
                const auto label = set.samples().sampleIndices[k];
                const auto stored = set.samples().values[k * set.liftedDimension() + 1];
                const std::complex<double> value(along.value(label - first, 0), along.value(label - first, 1));
                require(std::abs(value - stored) < 1e-7 * (1 + std::abs(stored)), "the window eigenfunction along the series reproduces H-DMD");
            }
        }

        // Interpolant methods: a list, the current one, switching it in place, errors per method.
        auto sweep = set;
        require(sweep.interpolation().data().method.name == "inverse_distance", "default interpolant");
        require(c::interpolationMethods().size() >= 7, "interpolant list");
        const double between = 0.5 * (record.series.value(samples30(set), 0) + record.series.value(samples30(set) + 1, 0));
        const auto idw = sweep.evaluateRepresentation(&between, "interpolation");
        for (const auto& name : c::interpolationMethods()) {
            c::InterpolationMethod method;
            method.name = name;
            sweep.setInterpolationMethod(method, 6);
            require(sweep.interpolation().data().method.name == name, "switched interpolant");
            const auto errors = sweep.errors();
            require(errors.valueRelative.size() == set.liftedDimension(), name + ": errors per mode");
            const auto value = sweep.evaluateRepresentation(&between, "interpolation");
            require(std::isfinite(value[0].real()), name + ": finite value off the samples");
        }
        c::InterpolationMethod thin;
        thin.name = "rbf_thin_plate";
        sweep.setInterpolationMethod(thin, 6);
        require(std::abs(sweep.evaluateRepresentation(&between, "interpolation")[0] - idw[0]) > 0, "methods differ off the samples");

        // Krom round trip keeps origin, representations, order and values; a method's arrays stay beside them.
        std::map<std::string, c::NpyArray> properties;
        properties["sample_interval"] = c::NpyArray::fromDoubles({1}, {result.deltaTime});
        properties["modes"] = c::NpyArray::fromComplex({result.featureDimension, result.rank}, result.modes);
        properties["mode_gram"] = c::NpyArray::fromComplex({result.rank, result.rank}, result.modeGram);
        const auto archive = c::makeKoopmanArchive(set, {}, properties, record.channels);
        require(c::isKrom(archive) && c::dataObjectKind(archive.metadataJson()) == c::object_kind::krom, "a krom declares its kind");
        const auto loaded = c::readKoopmanModel(archive);
        require(loaded.representations().front().name == "interpolation", "representation order survives");
        require(loaded.origin().method == "hdmd" && loaded.origin().inputs.front().path == "echo.npz", "origin survives");
        require(loaded.samples().values == set.samples().values, "sample values survive");
        {
            const auto again = c::readKoopmanModel(c::makeKoopmanArchive(sweep));
            require(again.interpolation().data().method.name == "rbf_thin_plate" && again.interpolation().data().neighbors == 6,
                    "interpolant method survives the archive");
        }
        require(worst(loaded.errors("time-delay", &record).valueRelative) < 1e-7, "loaded window map still exact");

        {   // Names: unique per set; a kind finds a representation only when one has it.
            require(set.representationNames() == std::vector<std::string>{"interpolation", "time-delay"}, "derived names");
            auto combined = set;
            const auto dictionary = std::make_shared<const c::ExpansionDictionary>(c::ExpansionDictionary::fromJson(
                c::json::parse(R"({"schema":"dynsyslab.expansion-dictionary/v1","dimension":1,"terms":[{"kind":"monomial","orders":[0]},{"kind":"monomial","orders":[1]}]})")));
            std::vector<std::complex<double>> coefficients;
            for (std::size_t j = 0; j < set.liftedDimension(); ++j) coefficients.insert(coefficients.end(), {{0.1 * j, 0}, {1, 0.2}});
            auto echo = c::FunctionDomain::supplied(1);
            echo.names = {"echo"};
            combined.addRepresentation(std::make_shared<const c::FunctionExpansion>(dictionary, coefficients, echo), true);
            require(combined.representationNames().front() == "polynomial", "a new default comes first");
            require(&combined.representation("time-delay") == &combined.representations()[2], "found by name");
            bool ambiguous = false, duplicate = false;
            try { (void)combined.representation("expansion"); } catch (const std::invalid_argument&) { ambiguous = true; }
            try { combined.addRepresentation(combined.representations()[1].function, false); } catch (const std::invalid_argument&) { duplicate = true; }
            require(ambiguous && duplicate, "an ambiguous kind and a repeated name are refused");

            // A child krom keeps its parent's parts by reference: none of their arrays are copied.
            const auto folder = std::filesystem::temp_directory_path() / ("dsl-extension-" + std::to_string(std::rand()));
            std::filesystem::create_directories(folder);
            const auto basePath = folder / "base.npz";
            c::writeNpzAtomic(basePath, archive);
            auto base = c::loadKoopmanModel(basePath);
            require(base.parent() && base.parent()->sha256 == c::sha256File(basePath), "a loaded krom is its own parent");
            base.addRepresentation(combined.representation("polynomial").function, true);
            const std::vector<std::size_t> pick{1, 2};
            const auto subset = base.extract(pick);
            const auto extension = c::makeKoopmanArchive(subset, {}, {}, record.channels);
            require(!extension.contains("window_map") && !extension.contains("interpolation_points")
                    && !extension.contains("eigenfunction_values") && !extension.contains("continuous_eigenvalues")
                    && extension.contains("expansion_coefficients"), "only the new representation is stored");
            const auto block = c::json::parse(extension.metadataJson()).at("krom");
            require(block.at("default_representation").string() == "polynomial", "its default is the new one");
            require(block.at("parent").at("mode_ids").array().size() == 2 && block.at("parent").at("mode_ids").array()[1].number() == 2.0,
                    "the parent's modes, in order");
            const auto read = c::readKoopmanModel(extension);
            require(read.representationNames() == std::vector<std::string>{"polynomial", "interpolation", "time-delay"},
                    "the parent's representations come back, in order");
            require(read.samples().samples() == set.samples().samples() && read.liftedDimension() == 2, "the parent's samples, its modes");
            const auto onWindow = read.evaluateRepresentation(h.data(), "time-delay");
            const auto direct = set.extract(pick).evaluateRepresentation(h.data(), "time-delay");
            require(std::abs(onWindow[0] - direct[0]) < 1e-9 * (1 + std::abs(direct[0])), "the referenced window map evaluates");
            std::filesystem::remove(basePath);
            bool missing = false;
            try { (void)c::readKoopmanModel(extension); } catch (const std::invalid_argument&) { missing = true; }
            require(missing, "a missing parent is refused");
            std::filesystem::remove_all(folder);
        }

        {   // The same framework for a general observable on a mapped domain: q = x_1 of a 2-D state.
            // Two observables of the same samples choose different interpolants.
            std::vector<double> states;
            std::vector<std::complex<double>> g;
            for (int k = 0; k < 40; ++k) {
                const double x0 = 0.1 * k, x1 = std::sin(0.3 * k);
                states.insert(states.end(), {x0, x1});
                g.emplace_back(x1 * x1);
            }
            const auto domain = c::FunctionDomain::fromStateMap(c::Observable::stateComponent(2, 1, "x1", "x1"));
            const auto nearest = c::interpolateSamples(states, 40, domain, g, {}, 1);
            c::InterpolationMethod linear;
            linear.name = "local_linear";
            const auto smooth = c::interpolateSamples(states, 40, domain, g, linear, 6);
            require(nearest->points().size() == 40 && nearest->points()[3] == states[7], "base points are the mapped states");
            c::ObservableDescriptor descriptor;
            descriptor.name = descriptor.label = "g";
            const auto a = c::Observable::fromRepresentation(descriptor, nearest);
            const auto b = c::Observable::fromRepresentation(descriptor, smooth);
            const double x[2] = {5.0, 0.37};
            require(a.representationKind() == "interpolation" && smooth->data().method.name == "local_linear",
                    "each observable keeps its own interpolant");
            require(std::abs(b.evaluateReal(x)[0] - 0.37 * 0.37) < 0.05, "evaluates through the state map");
            require(std::isfinite(a.evaluateReal(x)[0]), "nearest-sample observable evaluates");
        }

        // Time series join only on a shared clock.
        const auto joined = c::joinTimeSeriesRecords({record, echoRecord(200, 0.1)});
        require(joined.series.componentCount() == 2 && joined.channels.size() == 2, "join concatenates channels");
        bool clockRejected = false;
        try { (void)c::joinTimeSeriesRecords({record, echoRecord(200, 0.1, 0.05)}); }
        catch (const std::invalid_argument&) { clockRejected = true; }
        require(clockRejected, "join refuses different clocks");
        std::cout << "krom eigenfunction tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
