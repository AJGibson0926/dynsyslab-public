#pragma once
// The eigenfunction representation specific to Hankel DMD. (Its interpolation is not specific: it is
// the shared interpolateSamples / SharedRepresentation every observable uses.)
//
//   window_map     the linear coordinate change from a Hankel window to psi:
//                  psi(h) = W^-1 U^T h = (Phi^H Phi)^-1 Phi^H h, with Phi = U W the Hankel modes
//                  (Arbabi & Mezic, SIAM J. Appl. Dyn. Syst. 16, 2017; Brunton et al., Nat. Commun. 8, 2017).
#include <dynsyslab/core/FunctionDomain.hpp>
#include <dynsyslab/core/EigenfunctionObjects.hpp>
#include <dynsyslab/core/Hankel.hpp>
#include <dynsyslab/core/FunctionInterpolation.hpp>

#include <mutex>

namespace dynsyslab::core {

class WindowMapRepresentation final : public EigenfunctionRepresentation {
public:
    // Full mode set: Hankel modes Phi (rows x rank) and their Gram matrix G = Phi^H Phi = W^H W.
    static std::shared_ptr<WindowMapRepresentation> fromModes(HankelWindowLayout layout,
        std::shared_ptr<const std::vector<std::complex<double>>> modes,
        std::shared_ptr<const std::vector<std::complex<double>>> gram, std::size_t rank);
    // Materialized map T (modes x rows).
    static std::shared_ptr<WindowMapRepresentation> fromMap(HankelWindowLayout layout,
                                                            std::vector<std::complex<double>> map, std::size_t modes);
    std::string kind() const override { return "expansion"; }  // linear basis on the window domain
    std::vector<std::string> domains() const override { return {"window"}; }
    std::size_t domainDimension(const std::string&) const override { return layout_.rows(); }
    std::vector<std::string> domainObservables(const std::string&) const override { return layout_.channelObservables; }
    // `point` is a raw window as produced by HankelWindowLayout::window_at; channel scales are applied here.
    std::vector<std::complex<double>> evaluate(const double* point, const std::string& domain) const override;
    std::vector<std::complex<double>> evaluateTraining(const EigenfunctionSamples&, const std::string& domain,
                                                       const TimeSeriesRecord* input) const override;
    std::shared_ptr<const EigenfunctionRepresentation> extract(const std::vector<std::size_t>& modes) const override;
    json::Value describe() const override;
    // The window map as a FunctionExpansion in degree-one monomials of the raw window.
    std::shared_ptr<const FunctionRepresentation> function(const std::string& domain) const override;
    const HankelWindowLayout& layout() const noexcept { return layout_; }
    std::size_t modes() const noexcept { return modeCount_; }
    // T rows for the given modes (modes x rows); materializes from Phi and G when needed.
    std::vector<std::complex<double>> mapRows(const std::vector<std::size_t>& modes) const;
    bool materialized() const noexcept { return !map_.empty(); }
    const std::vector<std::complex<double>>& map() const noexcept { return map_; }

private:
    WindowMapRepresentation() = default;
    HankelWindowLayout layout_;
    std::size_t modeCount_{};
    std::vector<std::complex<double>> map_;  // modes x rows, when materialized
    std::shared_ptr<const std::vector<std::complex<double>>> modes_, gram_;
    // G = L L^H, factored once; psi = G^-1 Phi^H h by two triangular solves.
    mutable std::vector<std::complex<double>> cholesky_;
    mutable std::once_flag choleskyOnce_;
    const std::vector<std::complex<double>>& factor() const;
    std::vector<std::complex<double>> solveGram(std::vector<std::complex<double>> b) const;
};

}  // namespace dynsyslab::core
