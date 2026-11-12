#pragma once
#include "JetPolynomial.hpp"

namespace dynsyslab::core::jet {
// Generic rank over the rational function field. Row elimination uses the
// identity pivot * row - entry * pivotRow, with an algebraically nonzero pivot.
inline std::size_t genericRank(const std::vector<Polynomial>& outputs,std::size_t n) {
    if (n>8 || outputs.size()>32) throw SearchLimit("symbolic rank matrix limit");
    std::vector<std::vector<Polynomial>> rows;
    for (const auto& p:outputs) {
        std::vector<Polynomial> row;
        for (std::size_t j=0;j<n;++j) row.push_back(partial(p,j));
        rows.push_back(std::move(row));
    }
    std::size_t rank=0;
    for (std::size_t j=0;j<n && rank<rows.size();++j) {
        std::size_t pivot=rank;
        while (pivot<rows.size() && rows[pivot][j].empty()) ++pivot;
        if (pivot==rows.size()) continue;
        std::swap(rows[rank],rows[pivot]);
        // Last column needs a pivot only, with no trailing elimination.
        for(std::size_t i=rank+1;i<rows.size();++i) if (!rows[i][j].empty())
            for(std::size_t k=j+1;k<n;++k) {
                auto value=multiply(rows[rank][j],rows[i][k]);
                addScaled(value,multiply(rows[i][j],rows[rank][k]),-1);
                rows[i][k]=std::move(value);
            }
        ++rank;
    }
    return rank;
}
} // namespace dynsyslab::core::jet
