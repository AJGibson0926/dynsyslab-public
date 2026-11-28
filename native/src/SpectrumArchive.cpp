#include <dynsyslab/core/SpectrumArchive.hpp>
#include <dynsyslab/core/KoopmanCompression.hpp>

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
        archive.arrays["pseudospectrum_json"] = NpyArray::fromString(json::dump(json::Value::Object{
            {"definition", p->definition}, {"scope", p->scope}}));
    }
    // The spectral bands and the observables' measures on them, as values; the measures' wave
    // packets and the pseudospectrum's minimizers are spectral functions (KromSpectralParts.cpp).
    if (!spectrum.bands().empty()) {
        json::Value::Array bands, measures;
        for (const auto& band : spectrum.bands()) bands.push_back(spectralBandJson(band));
        for (const auto& measure : spectrum.measures()) measures.push_back(spectralMeasureJson(measure));
        archive.arrays["spectral_bands_json"] = NpyArray::fromString(json::dump(json::Value(std::move(bands))));
        archive.arrays["spectral_measures_json"] = NpyArray::fromString(json::dump(json::Value(std::move(measures))));
    }
}

bool hasSpectralEvidence(const NpzRecord& archive) {
    return archive.contains("eigenpair_residual") || archive.contains("pseudospectrum_grid")
        || archive.contains("spectral_bands_json") || archive.contains("continuous_spectrum_json");
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
                         archive.at("pseudospectrum_grid").complexDoubles(), archive.at("pseudospectrum_minimum").doubles()};
        if (const auto& kept = spectrum.pseudospectrum(); kept && kept->grid == p.grid) p.minimizers = kept->minimizers;
        spectrum.setPseudospectrum(std::move(p));
    }
    // Bands (spectral_bands_json; continuous_spectrum_json before they were named bands).
    const char* bandsKey = archive.contains("spectral_bands_json") ? "spectral_bands_json" : "continuous_spectrum_json";
    if (archive.contains(bandsKey)) {
        std::vector<SpectralBand> bands;
        std::vector<SpectralMeasure> measures;
        const auto parsed = json::parse(archive.at(bandsKey).string());  // outlives the loop
        for (const auto& entry : parsed.array()) {
            if (entry.at("density").isNull()) { bands.push_back(spectralBandFromJson(entry)); continue; }
            // An archive written before bands and measures were separated: one measure per entry,
            // of the observable the entry was named for.
            SpectralBand c;
            c.name = entry.at("name").string(); c.support = entry.at("domain").string("unit-circle");
            c.method = entry.at("method").string(); c.assumptions = entry.at("assumptions").string();
            c.bandwidth = entry.at("bandwidth").number(0); c.order = std::size_t(entry.at("order").number(0));
            SpectralMeasure m;
            m.observable = m.band = c.name;
            for (const auto& value : entry.at("grid").array()) c.grid.push_back(value.number());
            for (const auto& value : entry.at("density").array()) m.density.push_back(value.number());
            m.totalMass = entry.at("total_mass").number(0);
            bands.push_back(std::move(c));
            measures.push_back(std::move(m));
        }
        if (archive.contains("spectral_measures_json")) {
            const auto listed = json::parse(archive.at("spectral_measures_json").string());
            for (const auto& entry : listed.array()) measures.push_back(spectralMeasureFromJson(entry));
        }
        // A measure keeps the wave packets the spectrum already has for it (a child's, from its parent).
        for (auto& m : measures)
            if (const auto* kept = spectrum.measure(m.band, m.observable, m.output)) m.packets = kept->packets;
        spectrum.setMeasures({});
        // A band keeps the generalized eigenfunctions the spectrum already has for it (a child's, from its parent).
        for (auto& b : bands)
            if (const auto* kept = spectrum.band(b.name); kept && b.eigenfunctions.empty()) b.eigenfunctions = kept->eigenfunctions;
        spectrum.setBands(std::move(bands));
        spectrum.setMeasures(std::move(measures));
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
    if (const auto ids = spectrum.modeIds(); [&] { for (std::size_t j = 0; j < ids.size(); ++j) if (ids[j] != j) return true; return false; }()) {
        std::vector<double> values(ids.begin(), ids.end());
        archive.arrays["mode_ids"] = NpyArray::fromDoubles({values.size()}, values);  // stable pair ids (Spectrum::modeId)
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
    if (archive.contains("mode_ids")) {
        std::vector<std::size_t> ids;
        for (const double id : archive.at("mode_ids").doubles()) ids.push_back(std::size_t(id));
        spectrum.setModeIds(std::move(ids));
    }
    readSpectralEvidence(archive, spectrum);
    return spectrum;
}

}  // namespace dynsyslab::core
