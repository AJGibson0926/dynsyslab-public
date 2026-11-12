#pragma once
#include "JetPolynomial.hpp"

namespace dynsyslab::core::jet {

// Closed intervals with exact rational endpoints. No floating-point rounding
// decisions enter the inclusion or the successor separation tests.
struct Interval {
    Rational low{0}, high{0};
};

inline Interval operator+(const Interval& a,const Interval& b) {
    return {a.low+b.low,a.high+b.high};
}
inline Interval operator-(const Interval& a,const Interval& b) {
    return {a.low-b.high,a.high-b.low};
}
inline Interval operator*(const Interval& a,const Interval& b) {
    const Rational products[]{a.low*b.low,a.low*b.high,a.high*b.low,a.high*b.high};
    return {*std::min_element(products,products+4),*std::max_element(products,products+4)};
}
inline Rational magnitude(const Interval& a) {return std::max(Rational(abs(a.low)),Rational(abs(a.high)));}
inline bool excludesZero(const Interval& a) {return a.low>0 || a.high<0;}

inline Interval intervalPower(Interval base,int exponent) {
    Interval result{1,1};
    if (exponent<0) {
        if (!excludesZero(base)) throw std::domain_error("interval contains a pole");
        base={1/base.high,1/base.low};
        exponent=-exponent;
    }
    while (exponent) {
        if (exponent%2) result=result*base;
        exponent/=2;
        if (exponent) base=base*base;
    }
    return result;
}
using Box = std::vector<Interval>;
inline Interval evaluate(const Polynomial& polynomial,const Box& box) {
    Interval result;
    for (const auto& [powers,c] : polynomial) {
        Interval term{c,c};
        for (std::size_t j=0;j<powers.size();++j)
            if (powers[j]) term=term*intervalPower(box[j],powers[j]);
        result=result+term;
    }
    return result;
}
inline Box pointBox(const std::vector<double>& seed) {
    Box point;
    for (double x : seed) point.push_back({Rational(x),Rational(x)});
    return point;
}
inline json::Value encodeBox(const Box& box) {
    json::Value::Array result;
    for (const auto& x : box) result.emplace_back(json::Value::Array{x.low.str(),x.high.str()});
    return result;
}
using Matrix = std::vector<std::vector<Rational>>;

// Choose a nonsingular square slice of the two-state system and form its
// inverse by reducing [J | I]. Unselected state components stay fixed.
inline bool inverseSlice(Matrix& augmented,std::size_t columns,std::vector<std::size_t>& pivots) {
    const auto m=augmented.size();
    for (std::size_t col=0;col<columns && pivots.size()<m;++col) {
        const auto row=pivots.size();
        auto pivot=row;
        while (pivot<m && augmented[pivot][col]==0) ++pivot;
        if (pivot==m) continue;
        std::swap(augmented[row],augmented[pivot]);
        const Rational divisor=augmented[row][col];
        for (auto& entry : augmented[row]) entry/=divisor;
        for (std::size_t i=0;i<m;++i) if (i!=row) {
            const Rational scale=augmented[i][col];
            for (std::size_t j=0;j<columns+m;++j) augmented[i][j]-=scale*augmented[row][j];
        }
        pivots.push_back(col);
    }
    return pivots.size()==m;
}
inline json::Value certifyFiber(const std::vector<Polynomial>& jet,
    const std::vector<Polynomial>& next,const std::vector<double>& left,
    const std::vector<double>& right,const std::set<std::size_t>& poles,
    const std::vector<std::size_t>& positiveAxes) {
    const auto n=left.size(), m=jet.size();
    const auto a=pointBox(left), b=pointBox(right);
    for (auto j : poles) if (!excludesZero(a[j]) || !excludesZero(b[j])) return {};
    std::vector<Rational> residual;
    bool exact=true;
    for (const auto& p : jet) {
        residual.push_back(evaluate(p,a).low-evaluate(p,b).low);
        exact=exact && residual.back()==0;
    }
    Box A=a,B=b;
    if (!exact) {
        if (m>2*n || m>12) return {};
        Matrix J(m,std::vector<Rational>(2*n+m));
        std::vector<std::vector<Polynomial>> derivatives(m);
        for (std::size_t i=0;i<m;++i) {
            for (std::size_t j=0;j<n;++j) {
                derivatives[i].push_back(partial(jet[i],j));
                J[i][j]=evaluate(derivatives[i][j],a).low;
                J[i][j+n]=-evaluate(derivatives[i][j],b).low;
            }
            J[i][2*n+i]=1;
        }
        std::vector<std::size_t> pivots;
        if (!inverseSlice(J,2*n,pivots)) return {};
        const Rational radius=Rational(1)/100000000;
        for (auto j : pivots) {
            auto& x=j<n?A[j]:B[j-n];
            x.low-=radius; x.high+=radius;
        }
        for (auto j : poles) if (!excludesZero(A[j]) || !excludesZero(B[j])) return {};
        std::vector<std::vector<Interval>> slopes(m,std::vector<Interval>(m));
        for (std::size_t i=0;i<m;++i) for (std::size_t k=0;k<m;++k) {
            const auto j=pivots[k];
            slopes[i][k]=j<n ? evaluate(derivatives[i][j],A) :
                Interval{}-evaluate(derivatives[i][j-n],B);
        }
        for (std::size_t i=0;i<m;++i) {
            Rational offset=0, norm=0;
            for (std::size_t l=0;l<m;++l) offset+=J[i][2*n+l]*residual[l];
            for (std::size_t k=0;k<m;++k) {
                Interval entry{i==k?1:0,i==k?1:0};
                for (std::size_t l=0;l<m;++l) {
                    const auto& y=J[i][2*n+l];
                    entry=entry-Interval{y,y}*slopes[l][k];
                }
                norm+=magnitude(entry);
            }
            // The preconditioned fixed-point map is a contraction into the box.
            if (abs(offset)+radius*norm>=radius) return {};
        }
    }
    for (auto j : positiveAxes) if (A[j].low<=0 || B[j].low<=0) return {};
    for (std::size_t i=0;i<next.size();++i) {
        const Interval gap=evaluate(next[i],A)-evaluate(next[i],B);
        if (excludesZero(gap)) return json::Value::Object{
            {"method",exact?"exact_pair":"rational_interval_contraction"},
            {"left_box",encodeBox(A)},{"right_box",encodeBox(B)},
            {"successor_component",double(i)},{"successor_gap",encodeBox({gap}).array().front()}};
    }
    return {};
}
} // namespace dynsyslab::core::jet
