#include <dynsyslab/core/SpectrumArchive.hpp>
#include <dynsyslab/core/SpectralKROM.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {

void writeSpectralEvidence(NpzRecord& archive, const Spectrum& spectrum) {
    const std::size_t m = spectrum.size();
    if (spectrum.hasResiduals()) {
        std::vector<double> values(m, std::numeric_limits<double>::quiet_NaN());
        json::Value::Array labels;
        for (std::size_t j = 0; j < m; ++j) {
            const auto& r = spectrum.residual(j);
            if (r) values[j] = r->value;
            labels.push_back(r ? json::Value(json::Value::Object{{"definition", r->definition}, {"scope", r->scope}}) : json::Value());
        }
        archive.arrays["eigenpair_residual"] = NpyArray::fromDoubles({m}, values);
        archive.arrays["eigenpair_residual_json"] = NpyArray::fromString(json::dump(json::Value(std::move(labels))));
    }
    if (const auto& p = spectrum.pseudospectrum(); p && !p->empty()) {
        archive.arrays["pseudospectrum_grid"] = NpyArray::fromComplex({p->grid.size()}, p->grid);
        archive.arrays["pseudospectrum_minimum"] = NpyArray::fromDoubles({p->minimum.size()}, p->minimum);
        if (!p->minimizingCoefficients.empty())
            archive.arrays["pseudospectrum_minimizing_coefficients"] =
                NpyArray::fromComplex({p->grid.size(), p->features}, p->minimizingCoefficients);
        archive.arrays["pseudospectrum_json"] = NpyArray::fromString(json::dump(json::Value::Object{
            {"definition", p->definition}, {"scope", p->scope}, {"features", double(p->features)}}));
    }
    if (!spectrum.continuous().empty()) {
        json::Value::Array components;
        for (const auto& component : spectrum.continuous()) components.push_back(spectralMeasureJson(component));
        archive.arrays["continuous_spectrum_json"] = NpyArray::fromString(json::dump(json::Value(std::move(components))));
    }
}

bool hasSpectralEvidence(const NpzRecord& archive) {
    return archive.contains("eigenpair_residual") || archive.contains("pseudospectrum_grid")
        || archive.contains("continuous_spectrum_json");
}

void readSpectralEvidence(const NpzRecord& archive, Spectrum& spectrum) {
    if (archive.contains("eigenpair_residual")) {
        const auto values = archive.at("eigenpair_residual").doubles();
        const auto labels = json::parse(archive.at("eigenpair_residual_json").string());
        if (values.size() != spectrum.size() || labels.array().size() != values.size())
            throw std::invalid_argument("eigenpair residuals do not match the spectrum");
        for (std::size_t j = 0; j < values.size(); ++j) {
            const auto& label = labels.array()[j];
            spectrum.setResidual(j, label.isNull() ? std::nullopt : std::optional<EigenpairResidual>(EigenpairResidual{
                values[j], label.at("definition").string(), label.at("scope").string()}));
        }
    }
    if (archive.contains("pseudospectrum_grid")) {
        const auto block = json::parse(archive.at("pseudospectrum_json").string());
        Pseudospectrum p{block.at("definition").string(), block.at("scope").string(),
                         archive.at("pseudospectrum_grid").complexDoubles(), archive.at("pseudospectrum_minimum").doubles(),
                         {}, std::size_t(block.at("features").number(0))};
        if (archive.contains("pseudospectrum_minimizing_coefficients"))
            p.minimizingCoefficients = archive.at("pseudospectrum_minimizing_coefficients").complexDoubles();
        spectrum.setPseudospectrum(std::move(p));
    }
    if (archive.contains("continuous_spectrum_json")) {
        std::vector<SpectralMeasure> components;
        const auto parsed = json::parse(archive.at("continuous_spectrum_json").string());
        for (const auto& component : parsed.array())
            components.push_back(spectralMeasureFromJson(component));
        spectrum.setContinuous(std::move(components));
    }
}

void writeSpectrumArrays(NpzRecord& archive, const Spectrum& spectrum) {
    const std::size_t m = spectrum.size();
    const bool sampled = spectrum.sampleInterval() > 0;
    if (spectrum.hasContinuousRates()) {
        archive.arrays["continuous_eigenvalues"] = NpyArray::fromComplex({m}, spectrum.continuousRates());
        // Derived summaries (Eigenvalue.hpp), for readers that want rates and frequencies by name.
        std::vector<double> rate(m), omega(m), frequency(m), quality(m);
        for (std::size_t j = 0; j < m; ++j) {
            const auto value = spectrum.eigenvalue(j);
            rate[j] = value.growthRate(); omega[j] = value.angularFrequency();
            frequency[j] = value.cyclicFrequency(); quality[j] = value.qualityFactor();
        }
        archive.arrays["growth_rate"] = NpyArray::fromDoubles({m}, rate);
        archive.arrays["angular_frequency"] = NpyArray::fromDoubles({m}, omega);
        archive.arrays["cyclic_frequency"] = NpyArray::fromDoubles({m}, frequency);
        archive.arrays["quality_factor"] = NpyArray::fromDoubles({m}, quality);
    }
    if (sampled) {
        archive.arrays["discrete_eigenvalues"] = NpyArray::fromComplex({m}, spectrum.discreteMultipliers());
        archive.arrays["sample_interval"] = NpyArray::fromDoubles({1}, {spectrum.sampleInterval()});
        std::vector<double> branches;
        for (const auto branch : spectrum.logarithmBranches()) branches.push_back(double(branch));
        if (!branches.empty()) archive.arrays["logarithm_branches"] = NpyArray::fromDoubles({branches.size()}, branches);
    }
    archive.arrays["eigenvalue_domain"] =
        NpyArray::fromString(spectrum.domain() == EigenvalueDomain::discreteTime ? "discrete" : "continuous");
    writeSpectralEvidence(archive, spectrum);
}

bool hasSpectrumArrays(const NpzRecord& archive) {
    return archive.contains("continuous_eigenvalues") || archive.contains("discrete_eigenvalues");
}

Spectrum readSpectrumArrays(const NpzRecord& archive, DataOrigin origin) {
    const bool discrete = archive.contains("eigenvalue_domain")
        ? archive.at("eigenvalue_domain").string() == "discrete"
        : !archive.contains("continuous_eigenvalues") && archive.contains("discrete_eigenvalues");
    const double interval = archive.contains("sample_interval") ? archive.at("sample_interval").doubles().at(0) : 0.0;
    std::vector<int> branches;
    if (archive.contains("logarithm_branches"))
        for (const double branch : archive.at("logarithm_branches").doubles()) branches.push_back(int(branch));
    const bool rates = archive.contains("continuous_eigenvalues");
    Spectrum spectrum;
    if (discrete) {
        const auto multipliers = archive.at("discrete_eigenvalues").complexDoubles();
        spectrum = rates ? Spectrum::fromStoredDiscrete(multipliers, interval, branches,
                                                        archive.at("continuous_eigenvalues").complexDoubles())
                         : Spectrum::fromDiscreteMultipliers(multipliers, interval, branches);
    } else {
        if (!rates) throw std::invalid_argument("a continuous spectrum needs continuous_eigenvalues");
        spectrum = Spectrum::fromStoredRates(archive.at("continuous_eigenvalues").complexDoubles(), interval,
                                             EigenvalueDomain::continuousTime, branches);
    }
    spectrum.withOrigin(std::move(origin));
    readSpectralEvidence(archive, spectrum);
    return spectrum;
}

}  // namespace dynsyslab::core
