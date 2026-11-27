#include <algorithm>
#include <cmath>
#include <climits>
#include <dynsyslab/core/Spectrum.hpp>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
bool finite(std::complex<double> x) {
    return std::isfinite(x.real()) && std::isfinite(x.imag());
}
void check(const std::vector<std::complex<double>> &v, double dt, bool sampled) {
    if (v.empty() || !std::all_of(v.begin(), v.end(), finite))
        throw std::invalid_argument("Spectrum needs finite nonempty eigenvalues");
    if (!std::isfinite(dt) || dt < 0 || (sampled && dt == 0))
        throw std::invalid_argument("Sampled spectrum needs a positive finite interval");
}
// The branch k with rate*dt = log(exp(rate*dt)) + 2 pi i k.
std::vector<int> branchesOf(const std::vector<std::complex<double>> &rates, double dt) {
    std::vector<int> out;
    const double turn = 2 * std::acos(-1.);
    for (const auto rate : rates) {
        const double branch = std::round((rate.imag() * dt - std::arg(std::exp(rate * dt))) / turn);
        if (!std::isfinite(branch) || branch < INT_MIN || branch > INT_MAX)
            throw std::invalid_argument("Stored logarithm branch exceeds supported range");
        out.push_back(static_cast<int>(branch));
    }
    return out;
}
} // namespace

void Spectrum::pairConjugates() {
    conjugates_.assign(values_.size(), std::nullopt);
    for (std::size_t j = 0; j < values_.size(); ++j) {
        const auto v = values_[j];
        if (std::abs(v.imag()) <= 1e-12 * std::max(1.0, std::abs(v))) continue;
        for (std::size_t k = 0; k < values_.size(); ++k)
            if (k != j && std::abs(values_[k] - std::conj(v)) <= 1e-9 + 1e-7 * std::abs(v)) {
                conjugates_[j] = k;
                break;
            }
    }
}
Spectrum Spectrum::fromContinuousRates(std::vector<std::complex<double>> v, double dt) {
    return fromStoredRates(std::move(v), dt, EigenvalueDomain::continuousTime);
}
Spectrum Spectrum::fromStoredRates(std::vector<std::complex<double>> v, double dt, EigenvalueDomain origin,
                                   std::vector<int> branches) {
    check(v, dt, origin == EigenvalueDomain::discreteTime);
    if (!branches.empty() && branches.size() != v.size())
        throw std::invalid_argument("Spectrum branches must be mode aligned");
    Spectrum out;
    out.domain_ = origin;
    out.sampleInterval_ = dt;
    out.rates_ = v;
    out.rateAvailable_.assign(v.size(), true);
    out.branches_ = std::move(branches);
    if (origin == EigenvalueDomain::discreteTime) {
        if (out.branches_.empty()) out.branches_ = branchesOf(v, dt);
        for (auto &x : v) {
            x = std::exp(x * dt);
            if (!finite(x) || x == std::complex<double>{})
                throw std::invalid_argument(
                    "Stored rate conversion overflow/underflow; supply primary discrete multipliers instead");
        }
    }
    out.values_ = std::move(v);
    out.pairConjugates();
    return out;
}
Spectrum Spectrum::fromDiscreteMultipliers(std::vector<std::complex<double>> v, double dt,
                                           std::vector<int> branches) {
    check(v, dt, true);
    if (branches.empty())
        branches.assign(v.size(), 0);
    if (branches.size() != v.size())
        throw std::invalid_argument("Spectrum branches must be mode aligned");
    Spectrum out;
    out.domain_ = EigenvalueDomain::discreteTime;
    out.sampleInterval_ = dt;
    out.values_ = std::move(v);
    out.branches_ = std::move(branches);
    for (std::size_t i = 0; i < out.values_.size(); ++i) {
        const auto mu = out.values_[i];
        out.rateAvailable_.push_back(mu != std::complex<double>{});
        out.rates_.push_back(
            mu == std::complex<double>{}
                ? std::complex<double>{}
                : (std::log(mu) + std::complex<double>(0, 2 * std::acos(-1.) * out.branches_[i])) / dt);
        if (!finite(out.rates_.back()))
            throw std::invalid_argument("Continuous spectrum conversion overflow");
    }
    out.pairConjugates();
    return out;
}
Spectrum Spectrum::fromDiscrete(std::vector<std::complex<double>> multipliers,
                                const std::vector<std::complex<double>> &rates, double dt, DataOrigin origin) {
    if (multipliers.size() != rates.size())
        throw std::invalid_argument("Spectrum needs one rate per multiplier");
    check(rates, dt, true);
    auto out = fromStoredDiscrete(std::move(multipliers), dt, branchesOf(rates, dt), rates);
    out.origin_ = std::move(origin);
    return out;
}
bool Spectrum::hasContinuousRates() const noexcept {
    return std::all_of(rateAvailable_.begin(), rateAvailable_.end(), [](bool x) { return x; });
}
Spectrum Spectrum::fromStoredDiscrete(std::vector<std::complex<double>> mu, double dt, std::vector<int> branches,
                                      const std::vector<std::complex<double>> &rates) {
    auto out = fromDiscreteMultipliers(std::move(mu), dt, std::move(branches));
    if (rates.size() != out.size())
        throw std::invalid_argument("Archived rates must be mode aligned");
    for (std::size_t j = 0; j < out.size(); ++j) {
        if (!finite(rates[j]) || (!out.rateAvailable_[j] && rates[j] != std::complex<double>{}))
            throw std::invalid_argument("Invalid archived continuous rate");
        if (out.rateAvailable_[j] && (std::abs(std::exp(rates[j] * dt) - out.values_[j]) >
                                         1e-12 * std::max(1., std::abs(out.values_[j]))
            || std::abs((rates[j] - out.rates_[j]) * dt) > 1e-11 * std::max(1., std::abs(rates[j] * dt))))
            throw std::invalid_argument("Archived rate does not match its discrete multiplier");
        out.rates_[j] = rates[j];
    }
    return out;
}
std::optional<std::complex<double>> Spectrum::continuousRate(std::size_t mode) const {
    if (!rateAvailable_.at(mode))
        return {};
    return rates_.at(mode);
}
const std::vector<std::complex<double>> &Spectrum::continuousRates() const {
    if (!hasContinuousRates())
        throw std::invalid_argument("Zero discrete multiplier has no finite continuous generator; "
                                    "use the discrete spectrum");
    return rates_;
}
std::vector<std::complex<double>> Spectrum::discreteMultipliers(double dt) const {
    if (dt == 0)
        dt = sampleInterval_;
    if (!std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("Discrete spectrum requires a positive finite interval");
    if (domain_ == EigenvalueDomain::discreteTime && dt == sampleInterval_)
        return values_;
    std::vector<std::complex<double>> result;
    for (std::size_t j = 0; j < size(); ++j)
        result.push_back(transition(j, dt));
    return result;
}
std::complex<double> Spectrum::transition(std::size_t mode, double elapsed) const {
    if (!std::isfinite(elapsed))
        throw std::invalid_argument("Spectrum transition requires finite elapsed time");
    if (!rateAvailable_.at(mode)) {
        const double steps = elapsed / sampleInterval_;
        if (!std::isfinite(steps) || steps < 0 ||
            std::abs(steps - std::round(steps)) > 1e-9 * std::max(1., std::abs(steps)))
            throw std::invalid_argument("Zero sampled mode only supports nonnegative integer sample horizons");
        return std::round(steps) == 0 ? 1. : 0.;
    }
    const auto result = std::exp(rates_.at(mode) * elapsed);
    if (!finite(result))
        throw std::invalid_argument("Spectrum transition overflow");
    return result;
}
Spectrum Spectrum::extract(const std::vector<std::size_t> &modes) const {
    if (modes.empty() || std::set<std::size_t>(modes.begin(), modes.end()).size() != modes.size())
        throw std::invalid_argument("Spectrum selection must be nonempty and unique");
    Spectrum result;
    result.domain_ = domain_;
    result.sampleInterval_ = sampleInterval_;
    result.origin_ = origin_;
    for (auto j : modes) {
        result.values_.push_back(values_.at(j));
        result.rates_.push_back(rates_.at(j));
        result.rateAvailable_.push_back(rateAvailable_.at(j));
        if (!branches_.empty())
            result.branches_.push_back(branches_.at(j));
        if (!residuals_.empty())
            result.residuals_.push_back(residuals_.at(j));
    }
    result.continuous_ = continuous_;
    result.pseudospectrum_ = pseudospectrum_;
    result.pairConjugates();
    return result;
}
Eigenvalue Spectrum::eigenvalue(std::size_t mode) const {
    if (mode >= size()) throw std::out_of_range("eigenvalue index");
    Eigenvalue e;
    e.index = mode;
    e.sampleInterval = sampleInterval_;
    e.continuous = continuousRate(mode).value_or(std::complex<double>{});
    e.discrete = domain_ == EigenvalueDomain::discreteTime ? values_[mode]
               : sampleInterval_ > 0 ? std::exp(e.continuous * sampleInterval_) : std::complex<double>{};
    if (mode < conjugates_.size()) e.conjugate = conjugates_[mode];
    return e;
}
Eigenpair Spectrum::pair(std::size_t mode) const {
    return {eigenvalue(mode), std::nullopt, residual(mode)};
}
const std::optional<EigenpairResidual>& Spectrum::residual(std::size_t mode) const {
    static const std::optional<EigenpairResidual> none;
    if (mode >= size()) throw std::out_of_range("eigenpair index");
    return residuals_.empty() ? none : residuals_[mode];
}
bool Spectrum::hasResiduals() const noexcept {
    return std::any_of(residuals_.begin(), residuals_.end(), [](const auto& r) { return r.has_value(); });
}
Spectrum& Spectrum::setResidual(std::size_t mode, std::optional<EigenpairResidual> value) {
    if (mode >= size()) throw std::out_of_range("eigenpair index");
    if (value && (!std::isfinite(value->value) || value->value < 0 || value->definition.empty()))
        throw std::invalid_argument("an eigenpair residual is a finite nonnegative value with a definition");
    if (residuals_.empty()) residuals_.resize(size());
    residuals_[mode] = std::move(value);
    return *this;
}
Spectrum& Spectrum::setResiduals(const std::vector<double>& values, const std::string& definition, const std::string& scope) {
    if (!values.empty() && values.size() != size()) throw std::invalid_argument("one residual per eigenpair");
    residuals_.clear();
    for (std::size_t j = 0; j < values.size(); ++j) setResidual(j, EigenpairResidual{values[j], definition, scope});
    return *this;
}
Spectrum& Spectrum::setPseudospectrum(std::optional<Pseudospectrum> pseudospectrum) {
    if (pseudospectrum) pseudospectrum->validate();
    pseudospectrum_ = std::move(pseudospectrum);
    return *this;
}
Spectrum& Spectrum::adoptEvidence(const Spectrum& other) {
    if (other.size() != size()) throw std::invalid_argument("spectral evidence belongs to a spectrum of another size");
    residuals_ = other.residuals_;
    continuous_ = other.continuous_;
    pseudospectrum_ = other.pseudospectrum_;
    return *this;
}
std::vector<Eigenvalue> Spectrum::eigenvalues() const {
    std::vector<Eigenvalue> out;
    for (std::size_t j = 0; j < size(); ++j) out.push_back(eigenvalue(j));
    return out;
}
} // namespace dynsyslab::core
