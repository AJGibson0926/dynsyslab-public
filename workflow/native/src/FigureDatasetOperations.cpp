// figure_dataset operations (FigureDatasetOperations.hpp).
#include "FigureDatasetOperations.hpp"

#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <cmath>
#include <complex>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace json = core::json;

FigureSeries psiSeries(const json::Value& spec, const core::NpzRecord& source) {
    const auto kromPath = spec.at("krom").string();
    if (kromPath.empty()) throw std::invalid_argument("figure_dataset psi needs a krom path");
    const auto model = core::loadKoopmanModel(std::filesystem::u8path(kromPath));
    const auto& states = source.at(spec.at("array").string("values"));
    const std::size_t columns = states.shape.size() == 2 ? states.shape[1] : 1;
    if (columns != model.stateDimension())
        throw std::invalid_argument("figure_dataset psi: the source rows are not states of the K-ROM");
    std::vector<std::size_t> modes;
    for (const auto& mode : spec.at("modes").array()) modes.push_back(static_cast<std::size_t>(mode.number()));
    if (modes.empty()) for (std::size_t k = 0; k < model.eigenfunctionCount(); ++k) modes.push_back(k);
    const auto part = spec.at("part").string("real");
    if (part != "real" && part != "imag" && part != "abs" && part != "phase")
        throw std::invalid_argument("figure_dataset psi part must be real, imag, abs or phase");
    const auto values = states.doubles();
    const std::size_t rows = values.size() / columns;
    FigureSeries series(modes.size(), std::vector<double>(rows));
    for (std::size_t r = 0; r < rows; ++r) {
        const auto psi = model.evaluate(values.data() + r * columns);
        for (std::size_t k = 0; k < modes.size(); ++k) {
            if (modes[k] >= psi.size()) throw std::invalid_argument("figure_dataset psi mode is outside the K-ROM");
            const auto value = psi[modes[k]];
            series[k][r] = part == "real" ? value.real() : part == "imag" ? value.imag()
                         : part == "abs" ? std::abs(value) : std::arg(value);
        }
    }
    return series;
}

void applyOperation(const std::string& operation, const json::Value& spec, FigureSeries& series,
                    const core::NpzRecord& source, const std::function<FigureSeries(const json::Value&)>& seriesOf) {
    if (series.empty() || series.front().empty()) throw std::invalid_argument("figure_dataset: empty series");
    const std::size_t rows = series.front().size();
    if (operation == "column" || operation == "psi") return;
    if (operation == "norm") {
        std::vector<double> norm(rows, 0.0);
        for (const auto& column : series)
            for (std::size_t r = 0; r < rows; ++r) norm[r] += column[r] * column[r];
        for (auto& value : norm) value = std::sqrt(value);
        series = {std::move(norm)};
        return;
    }
    if (operation == "difference") {
        if (spec.at("minus").isNull()) throw std::invalid_argument("figure_dataset difference needs \"minus\"");
        const auto other = seriesOf(spec.at("minus"));
        if (other.front().size() != rows || (other.size() != 1 && other.size() != series.size()))
            throw std::invalid_argument("figure_dataset difference: the two specs differ in shape");
        for (std::size_t k = 0; k < series.size(); ++k)
            for (std::size_t r = 0; r < rows; ++r) series[k][r] -= other[other.size() == 1 ? 0 : k][r];
        return;
    }
    if (operation == "diff") {
        if (rows < 2) throw std::invalid_argument("figure_dataset diff needs two rows");
        std::vector<double> step(rows - 1, 1.0);
        if (spec.at("per_time").boolean(false)) {
            const auto clock = source.at(spec.at("clock").string("time")).doubles();
            if (clock.size() != rows) throw std::invalid_argument("figure_dataset diff per_time needs a matching clock");
            for (std::size_t r = 0; r + 1 < rows; ++r) step[r] = clock[r + 1] - clock[r];
        }
        for (auto& column : series) {
            std::vector<double> differences(rows - 1);
            for (std::size_t r = 0; r + 1 < rows; ++r) differences[r] = (column[r + 1] - column[r]) / step[r];
            column = std::move(differences);
        }
        return;
    }
    throw std::invalid_argument("figure_dataset: unknown operation '" + operation + "'");
}

}  // namespace dynsyslab::workflow::postprocess_detail
