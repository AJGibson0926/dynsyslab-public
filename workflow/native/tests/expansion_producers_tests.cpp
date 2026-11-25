// DMD, EDMD and hand-written eigenfunctions produce K-ROMs that evaluate in their own bases.
#include <dynsyslab/core/EigenfunctionArchive.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace core = dynsyslab::core;
namespace json = core::json;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

json::Value job(const std::filesystem::path& source, const char* method, json::Value::Object settings) {
    return json::Value::Object{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"}, {"label", method},
        {"sources", json::Value::Array{source.generic_u8string()}},
        {"methods", json::Value::Array{json::Value::Object{{"name", method}, {"enabled", true}, {"settings", settings}}}},
        {"archive", json::Value::Object{{"enabled", true}}}};
}

// The K-ROM evaluated at the source states reproduces the method's stored eigenfunction values.
void reproducesStoredValues(const std::filesystem::path& archivePath, const std::vector<double>& states,
                            const std::string& what) {
    const auto archive = core::readNpz(archivePath);
    require(archive.contains("expansion_coefficients") && core::hasEigenfunctionObjects(archive),
            what + " stored no eigenfunctions object with its expansion");
    const auto model = core::loadKoopmanModel(archivePath);
    const auto stored = archive.at("eigenfunction_values").complexDoubles();
    const auto modes = archive.at("eigenfunction_values").shape.at(1);
    for (std::size_t sample : {std::size_t(3), std::size_t(17)}) {
        const auto psi = model.evaluate(states.data() + 2 * sample);
        for (std::size_t j = 0; j < modes; ++j)
            require(std::abs(psi[j] - stored[sample * modes + j]) < 1e-6 * (1 + std::abs(stored[sample * modes + j])),
                    what + " K-ROM does not reproduce its eigenfunction values");
    }
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root = std::filesystem::temp_directory_path()
             / ("dynsyslab-expansion-producers-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto sourcePath = root / "trajectory.npz";
        constexpr std::size_t samples = 60;
        std::vector<double> time(samples), states;
        for (std::size_t k = 0; k < samples; ++k) {
            const double t = 0.05 * double(k);
            time[k] = t;
            states.push_back(std::exp(-0.1 * t) * std::cos(1.7 * t));
            states.push_back(std::exp(-0.1 * t) * (-0.1 * std::cos(1.7 * t) - 1.7 * std::sin(1.7 * t)));
        }
        core::NpzRecord source;
        source.arrays["time"] = core::NpyArray::fromDoubles({samples}, time);
        source.arrays["values"] = core::NpyArray::fromDoubles({samples, 2}, states);
        source.arrays["metadata_json"] = core::NpyArray::fromString(
            R"({"format":"dynsyslab.trajectory/v1","label":"source","system":"duffing","system_preset":"double_well_reference","resolved_system_parameters":{},"integration_form":"dimensionless"})");
        core::writeNpzAtomic(sourcePath, source);

        const auto dmd = dynsyslab::workflow::executeJob(job(sourcePath, "dmd", {{"rank", 2.0}}), root / "dmd");
        reproducesStoredValues(dmd.archives.front(), states, "DMD");
        const auto edmd = dynsyslab::workflow::executeJob(job(sourcePath, "edmd", {
            {"dictionary", "polynomial"}, {"polynomial_degree", 1.0}, {"include_constant", false}}), root / "edmd");
        reproducesStoredValues(edmd.archives.front(), states, "EDMD");

        // Hand-written: psi = (r - i v / 1.7), an approximate eigenfunction of the damped rotation.
        const auto hand = dynsyslab::workflow::executeJob(job(sourcePath, "eigenfunction_expansion", {
            {"variables", json::Value::Array{"r", "v"}},
            {"eigenvalues", json::Value::Array{json::Value(json::Value::Array{-0.1, 1.7})}},
            {"eigenfunctions", json::Value::Array{json::Value::Object{{"real", "r"}, {"imag", "-v/1.7"}}}}}), root / "hand");
        const auto model = core::loadKoopmanModel(hand.archives.front());
        // r - i v/1.7 decomposes into monomials, so it shares terms with fitted polynomial banks.
        require(model.eigenfunctionRepresentation() == "analytic-polynomial", "hand-written text decomposes into monomials");
        const double x[]{0.8, -0.34};
        const auto psi = model.evaluate(x);
        require(std::abs(psi[0] - std::complex<double>(0.8, 0.2)) < 1e-12, "hand-written eigenfunction value");

        // The same delivery with explicit dictionary terms and coefficients: 2 r - i exp(-v).
        const auto terms = dynsyslab::workflow::executeJob(job(sourcePath, "eigenfunction_expansion", {
            {"variables", json::Value::Array{"r", "v"}},
            {"eigenvalues", json::Value::Array{json::Value(json::Value::Array{-0.1, 1.7})}},
            {"dictionary", json::Value::Array{json::Value::Object{{"monomial", json::Value::Array{1.0, 0.0}}},
                                              json::Value::Object{{"exponential", json::Value::Array{0.0, -1.0}}}}},
            {"coefficients", json::Value::Array{json::Value(json::Value::Array{2.0, json::Value(json::Value::Array{0.0, -1.0})})}}}),
            root / "terms");
        const auto termPsi = core::loadKoopmanModel(terms.archives.front()).evaluate(x);
        require(std::abs(termPsi[0] - std::complex<double>(1.6, -std::exp(0.34))) < 1e-12, "dictionary-term delivery");
        // The system's own library: eta = x_2 + c x_1^2/(b - 2a) of the triangular benchmark
        // (a=-1, b=-3, c=1 at gradient_closed_reference, so eta = x_2 - x_1^2), rate b.
        const auto triangularPath = root / "triangular.npz";
        core::NpzRecord triangular = source;
        triangular.arrays["metadata_json"] = core::NpyArray::fromString(
            R"({"format":"dynsyslab.trajectory/v1","label":"source","system":"quadratic_triangular","system_preset":"gradient_closed_reference","resolved_system_parameters":{},"integration_form":"dimensionless"})");
        core::writeNpzAtomic(triangularPath, triangular);
        const auto library = dynsyslab::workflow::executeJob(job(triangularPath, "eigenfunction_expansion", {
            {"library", json::Value::Array{"eta"}}}), root / "library");
        const auto libraryModel = core::loadKoopmanModel(library.archives.front());
        require(std::abs(libraryModel.evaluate(x)[0] - std::complex<double>(-0.34 - 0.64)) < 1e-12, "library eigenfunction value");
        require(std::abs(libraryModel.continuousEigenvalues().at(0) - std::complex<double>(-3.0)) < 1e-12, "library eigenvalue");
        std::filesystem::remove_all(root);
        std::cout << "expansion producer tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (!root.empty()) std::filesystem::remove_all(root);
        return 1;
    }
}
