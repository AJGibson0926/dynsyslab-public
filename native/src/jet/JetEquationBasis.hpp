#pragma once
#include "JetPolynomial.hpp"
#include <optional>

namespace dynsyslab::core::jet {
// Row reduction takes place in the source monomial basis. Each pivot stores
// its identical expression in q, so a zero remainder constructs a global law.
class EquationBasis {
    struct Row { Polynomial source, output; };
    std::map<Powers,Row> pivots_;
    const std::vector<Polynomial>& jet_;
    const std::vector<Polynomial>* graph_;
    std::size_t columns_{};

    void reduce(Polynomial& p, Polynomial& output) const {
        while (!p.empty()) {
            const auto pivot=pivots_.find(p.rbegin()->first);
            if (pivot==pivots_.end()) break;
            const Rational factor=p.rbegin()->second;
            addScaled(p,pivot->second.source,-factor);
            addScaled(output,pivot->second.output,factor);
        }
    }
    void insert(const Polynomial& p,const Powers& powers) {
        if (++columns_>256) throw SearchLimit("candidate monomial limit (256)");
        Polynomial source=onGraph(p,graph_), removed;
        reduce(source,removed);
        if (source.empty()) return;
        Polynomial expression{{powers,1}};
        addScaled(expression,removed,-1);
        const auto key=source.rbegin()->first;
        const Rational scale=source.rbegin()->second;
        for (auto& term:source) term.second/=scale;
        for (auto& term:expression) term.second/=scale;
        pivots_.emplace(key,Row{std::move(source),std::move(expression)});
    }
    void enumerate(const Polynomial& p,Powers& powers,std::size_t first,std::size_t remaining) {
        if (!remaining) {insert(p,powers); return;}
        for (std::size_t i=first;i<jet_.size();++i) {
            ++powers[i];
            enumerate(multiply(p,jet_[i]),powers,i,remaining-1);
            --powers[i];
        }
    }
public:
    explicit EquationBasis(const std::vector<Polynomial>& jet,const std::vector<Polynomial>* graph=nullptr)
        :jet_(jet),graph_(graph) {}
    void addDegree(std::size_t degree,std::size_t stateDimension) {
        Powers powers(jet_.size());
        enumerate({{Powers(stateDimension),1}},powers,0,degree);
    }
    std::optional<Polynomial> express(Polynomial p) const {
        p=onGraph(std::move(p),graph_);
        Polynomial output;
        reduce(p,output);
        if (!p.empty()) return {};
        return output;
    }
};
} // namespace dynsyslab::core::jet
