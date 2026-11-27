#include <chrono>
#include <cmath>
#include <dynsyslab/core/EdmdControl.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/EigenfunctionValidation.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <iostream>
#include <stdexcept>
namespace {
namespace c = dynsyslab::core;
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool caught = false;
    try {
        f();
    } catch (const std::logic_error &) {
        caught = true;
    }
    require(caught, "Invalid contract accepted");
}
c::EigenfunctionFitResult polynomial() {
    c::FunctionDomain q;
    q.dimension = 2;
    q.names = {"xi_1", "xi_2"};
    q.providerArtifact = "frozen-projector.npz";
    q.providerIdentity = "causal-pod/v1";
    return {c::FunctionExpansion({2, {0, 0, 1, 0, 0, 1}, {1, 2}, {2, 3}}, {2., 1., {0, 1}, 0., 3., 0.}, q),
            c::Spectrum::fromDiscreteMultipliers({{.9, .2}, 0.}, .1), {}};
}
} // namespace
int main() {
    std::filesystem::path path;
    try {
        const auto fit = polynomial();
        const auto& bank = fit.functions;
        const double q[]{5, 8};
        require(bank.evaluate(0, q) == std::complex<double>(4, 2), "Normalized polynomial value");
        require(bank.gradient(0, q) == std::vector<std::complex<double>>({.5, {0, 1. / 3}}),
                "Normalized polynomial gradient");
        // An Eigenfunction object (Observable + Eigenvalue) shares the expansion's storage.
        const auto psi = c::KoopmanModel(polynomial(), 1).eigenfunction(0);
        require(psi.evaluate(q).front() == bank.evaluate(0, q), "Eigenfunction failed to retain storage");
        require(psi.representationKind() == "expansion", "Eigenfunction representation property");
        require(std::abs(psi.eigenvalue().discrete - std::complex<double>(.9, .2)) < 1e-15, "Eigenfunction eigenvalue");
        require(!fit.spectrum.continuousRate(1), "Zero mu acquired a generator");
        require(fit.spectrum.discreteMultipliers()[1] == std::complex<double>{},
                "Zero multiplier changed");
        rejects([&] { (void)fit.spectrum.continuousRates(); });
        rejects([&] { (void)fit.spectrum.transition(1, .05); });
        rejects([&] { (void)c::KoopmanModel(polynomial(), 1).extract({0, 0}); });
        rejects([&] { (void)bank.select({2}); });
        c::KoopmanModel model(fit, 1);
        require(!model.hasPhysicalSystem() && !model.physicalStateDimension(),
                "Invented physical system/dimension");
        require(model.eigenfunctionCount() == 2 && model.realLiftedDimension() == 3,
                "Wrong dimensions");
        require(!model.spectrum().hasResiduals(),
                "Manual coefficients invented accuracy");
        require(model.evaluateFromCoordinates(q) == bank.evaluate(q),
                "Runtime not using polynomial contract");
        rejects([&] { (void)model.evaluateFromState(q); });
        rejects([&] { (void)model.stateDimension(); });
        rejects([&] { (void)model.systemConfiguration(); });
        rejects([&] { (void)model.stateJacobianControlMap(q); });
        rejects([&] { (void)model.continuousGeneratorReal(); });
        model.setDiscreteInputMap("edmdc", {{1., 2.}, 2, 1, .1, "regression"});
        require(!model.hasControlMap("edmdc"), "Discrete input map masquerades as continuous Bpsi");
        rejects([&] { model.setDiscreteInputMap("bad", {{1., 2.}, 2, 1, .2, {}}); });
        model.setControlMap("regressed-a",
                            c::ConstantControlMap{{3., 4.}, 2, 1, "continuous regression"});
        model.addDatasetReference({"training.npz", "training", "forced", {}});
        model.addDatasetReference({"validation.npz", "validation", "unforced", {}});
        const auto reduced = model.extract({1});
        require(reduced.discreteInputMaps().at("edmdc").matrix[0] == 2.,
                "Discrete map row extraction");
        require(reduced.datasetReferences().size() == 2, "Extraction lost provenance");
        path =
            std::filesystem::temp_directory_path() /
            ("dsl-eigenfunction-contract-" +
             std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".npz");
        c::saveKoopmanModel(path, model);
        const auto loaded = c::loadKoopmanModel(path);
        require(!loaded.hasPhysicalSystem() && !loaded.physicalStateDimension(),
                "Archive invented physical state");
        require(loaded.evaluateCoordinates(q) == model.evaluateCoordinates(q),
                "Archive changed coefficients");
        require(loaded.spectrum().values() == model.spectrum().values(),
                "Archive changed primary spectrum");
        require(loaded.spectrum().continuousRate(0) == model.spectrum().continuousRate(0),
                "Archive changed rates");
        require(loaded.coordinates().providerArtifact == "frozen-projector.npz",
                "Archive lost frozen provider");
        require(loaded.discreteInputMaps().at("edmdc").matrix ==
                    model.discreteInputMaps().at("edmdc").matrix,
                "Archive lost discrete map");
        require(loaded.hasControlMap("regressed-a"), "Archive lost named continuous map");
        require(loaded.datasetReferences().size() == 2, "Archive lost dataset refs");
        const auto zero = bank.select({1});
        const auto validation =
            c::validateEigenfunctions(zero, fit.spectrum.extract({1}), {0, .1, .2, .3}, {5, 8, 1, 8, 1, 8, 1, 8}, {.1, .2});
        require(validation.modes[0].oneStep.absoluteRms == 0 &&
                    !validation.modes[0].generatorAvailable,
                "Discrete zero mode defect/generator availability");
        // A sampled linear model with a zero eigenvalue exports natively, without log(0).
        c::DmdControlResult dmd;
        dmd.dimension = 2;
        dmd.rank = 2;
        dmd.deltaTime = .1;
        dmd.hasEigenfunctions = true;
        dmd.spectrum = c::Spectrum::fromDiscreteMultipliers({.8, 0.}, dmd.deltaTime);
        dmd.eigenfunctionCoefficients = {1., 0., 0., 1.};
        require(dmd.eigenfunctionSet().spectrum.values() == dmd.spectrum.values(),
                "DMDc spectrum export");
        c::EdmdControlResult edmd;
        edmd.dictionary.dimension = 2;
        edmd.dictionary.degree = 1;
        edmd.dictionary.includeConstant = false;
        edmd.dictionary.exponents = c::polynomialExponentTable(2, 1, false);
        edmd.dictionary.center = {0, 0};
        edmd.dictionary.scale = {1, 1};
        edmd.lifted = dmd;
        require(edmd.eigenfunctionSet().spectrum.values() == dmd.spectrum.values(),
                "EDMDc spectrum export");
        c::FunctionDomain measured;
        measured.dimension = 1;
        std::vector<double> t, x, dx;
        for (std::size_t k = 0; k < 100; ++k) {
            t.push_back(.01 * k);
            x.push_back(std::exp(-.01 * k));
            dx.push_back(-x.back());
        }
        c::KoopmanSINDyOptions options;
        options.polynomialDegree = 1;
        options.normalizeStates = false;
        options.includeConstant = false;
        const auto sindyFit = c::fitSindyEigenfunctions(x, x.size(), t, {-1.}, measured, dx, options);
        require(sindyFit.diagnostics.algorithm == "sindy-generator" &&
                    sindyFit.diagnostics.sampleCount == 100,
                "SINDy fit evidence separated");
        require(sindyFit.spectrum.residual(0)->value < 1e-10, "SINDy wrapper changed mathematics");
        const auto high=c::Spectrum::fromStoredRates({{-.2,80}},.1,c::EigenvalueDomain::discreteTime);
        require(high.logarithmBranches()[0]==1,"Legacy high-frequency branch was not recovered");
        const auto highCopy=c::Spectrum::fromStoredDiscrete(high.values(),.1,
            high.logarithmBranches(),high.continuousRates());
        require(highCopy.continuousRates()==high.continuousRates(),"Archived high-frequency rate changed");
        rejects([&]{(void)c::Spectrum::fromStoredDiscrete(high.values(),.1,{0},high.continuousRates());});
        std::filesystem::remove(path);
        std::cout
            << "Eigenfunction contracts, compatibility, archives and discrete validation passed\n";
        return 0;
    } catch (const std::exception &e) {
        if (!path.empty())
            std::filesystem::remove(path);
        std::cerr << e.what() << '\n';
        return 1;
    }
}
