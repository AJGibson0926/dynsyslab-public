// Spectrum / eigenfunction objects: H-DMD window map, the shared interpolation (the one every
// observable uses, each with its own interpolant), errors, extraction, archive round trip, and
// time-series joining.
#include <dynsyslab/core/EigenfunctionArchive.hpp>
#include <dynsyslab/core/Spectral.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::size_t samples30(const c::EigenfunctionSet& set) { return set.samples().sampleIndices[30]; }

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

c::EigenfunctionSet fit(const c::TimeSeriesRecord& record, c::HankelDMDResult& result, c::HankelWindowLayout& layout) {
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
    auto interpolation = std::make_shared<c::SharedRepresentation>(
        c::interpolateSamples(points, points.size(), input, samples->values, {}, 4));
    auto window = c::WindowMapRepresentation::fromModes(layout,
        std::make_shared<const std::vector<std::complex<double>>>(result.modes),
        std::make_shared<const std::vector<std::complex<double>>>(result.modeGram), result.rank);
    return c::EigenfunctionSet(c::Spectrum::fromDiscrete(result.discreteEigenvalues, result.continuousEigenvalues,
                                                         result.deltaTime, origin),
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
        const auto window = set.errors("expansion", {}, &record);
        require(worst(window.valueRelative) < 1e-7, "window map reproduces H-DMD values");
        require(worst(window.koopmanResidual) < 1e-6, "window-map eigenfunctions advance by their eigenvalues");

        // Errors default to the first representation and its first domain.
        const auto defaults = set.errors();
        require(defaults.representation == "interpolation" && defaults.domain == "echo", "errors default to the first");
        require(worst(defaults.valueRelative) < 1e-12, "interpolation is exact at its own samples");

        // Domain recognition by the observables of the data handed in.
        require(set.recognizeDomain({"echo"}) == "echo", "echo data selects the echo domain");
        require(set.recognizeDomain({"echo"}, "expansion") == "window", "echo windows select the window domain");
        bool rejected = false;
        try { (void)set.recognizeDomain({"state", "state"}); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "a domain the representation lacks is refused");

        // Extraction materializes the window map rows for the kept modes.
        const auto pair = set.extract({0, 1});
        const auto h = layout.window_at(record.series, set.samples().sampleIndices[30]);
        const auto full = set.evaluate(h.data(), "window", "expansion");
        const auto kept = pair.evaluate(h.data(), "window", "expansion");
        require(std::abs(kept[0] - full[0]) < 1e-9 * (1 + std::abs(full[0])) && std::abs(kept[1] - full[1]) < 1e-9 * (1 + std::abs(full[1])),
                "extracted window map matches the full one");

        {   // An eigenfunction is an observable: evaluated along a time series it gives a time series. On
            // the window map's domain each sample's window is cut from the series; samples without a
            // complete window are skipped, so row k of the output is window label k + first label.
            const auto psi = set.eigenfunction(1, "window", "expansion");
            const auto along = psi.evaluate(record.series);
            require(along.componentCount() == 2, "a complex observable gives real and imaginary columns");
            const auto first = layout.labelOffset;
            require(std::abs(along.time(0) - record.series.time(first)) < 1e-12, "the first complete window labels the first row");
            for (std::size_t k = 0; k < set.samples().samples(); k += 7) {
                const auto label = set.samples().sampleIndices[k];
                const auto stored = set.samples().values[k * set.size() + 1];
                const std::complex<double> value(along.value(label - first, 0), along.value(label - first, 1));
                require(std::abs(value - stored) < 1e-7 * (1 + std::abs(stored)), "the window eigenfunction along the series reproduces H-DMD");
            }
        }

        // Interpolant methods: a list, the current one, switching it in place, errors per method.
        auto sweep = set;
        require(c::interpolationOf(sweep).data().method.name == "inverse_distance", "default interpolant");
        require(c::interpolationMethods().size() >= 7, "interpolant list");
        const double between = 0.5 * (record.series.value(samples30(set), 0) + record.series.value(samples30(set) + 1, 0));
        const auto idw = sweep.evaluate(&between, "echo");
        for (const auto& name : c::interpolationMethods()) {
            c::InterpolationMethod method;
            method.name = name;
            c::setInterpolationMethod(sweep, method, 6);
            require(c::interpolationOf(sweep).data().method.name == name, "switched interpolant");
            const auto errors = sweep.errors();
            require(errors.valueRelative.size() == set.size(), name + ": errors per mode");
            const auto value = sweep.evaluate(&between, "echo");
            require(std::isfinite(value[0].real()), name + ": finite value off the samples");
        }
        c::InterpolationMethod thin;
        thin.name = "rbf_thin_plate";
        c::setInterpolationMethod(sweep, thin, 6);
        require(std::abs(sweep.evaluate(&between, "echo")[0] - idw[0]) > 0, "methods differ off the samples");

        // Archive round trip keeps origin, representations, order and values.
        c::NpzRecord archive;
        auto objects = c::writeEigenfunctionSet(archive, set, record.channels);
        archive.arrays["sample_interval"] = c::NpyArray::fromDoubles({1}, {result.deltaTime});
        archive.arrays["modes"] = c::NpyArray::fromComplex({result.featureDimension, result.rank}, result.modes);
        archive.arrays["mode_gram"] = c::NpyArray::fromComplex({result.rank, result.rank}, result.modeGram);
        archive.arrays["metadata_json"] = c::NpyArray::fromString(
            c::json::dump(c::json::Value::Object{{"eigenfunction_objects", objects}}));
        const auto loaded = c::readEigenfunctionSet(archive);
        require(loaded.representations().front()->kind() == "interpolation", "representation order survives");
        require(loaded.origin().method == "hdmd" && loaded.origin().inputs.front().path == "echo.npz", "origin survives");
        {
            c::NpzRecord thinArchive;
            auto thinObjects = c::writeEigenfunctionSet(thinArchive, sweep, record.channels);
            thinArchive.arrays["metadata_json"] = c::NpyArray::fromString(
                c::json::dump(c::json::Value::Object{{"eigenfunction_objects", thinObjects}}));
            const auto again = c::readEigenfunctionSet(thinArchive, {"interpolation"});
            require(c::interpolationOf(again).data().method.name == "rbf_thin_plate" && c::interpolationOf(again).data().neighbors == 6,
                    "interpolant method survives the archive");
        }
        require(worst(loaded.errors("expansion", {}, &record).valueRelative) < 1e-7, "loaded window map still exact");

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
        std::cout << "eigenfunction objects tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
