#include "lookup.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
using receivers::Lookup;
void require(bool b){if(!b)throw std::runtime_error("lookup invariant failed");}
int main() {
    std::vector<double> x,y;
    for(int i=0;i<81;++i){x.insert(x.end(),{std::sin(i*.4),std::cos(i*.13)});y.insert(y.end(),{double(i),double(i*i)});}
    Lookup model(x,y,{.2,3.},2);double q[]{.17,-.31};
    const auto result=model.predict(q,8);
    std::vector<std::pair<double,int>> brute;
    for(int i=0;i<81;++i){
        const double a=std::asinh(x[2*i]/.2)-std::asinh(q[0]/.2), b=std::asinh(x[2*i+1]/3.)-std::asinh(q[1]/3.);
        brute.emplace_back(a*a+b*b,i);
    }
    std::sort(brute.begin(),brute.end());
    double sum=0,expected=0;
    for(int j=0;j<8;++j){require(result.neighbors[j].row==static_cast<std::size_t>(brute[j].second));
        sum+=1/brute[j].first;expected+=brute[j].second/brute[j].first;}
    require(std::abs(result.values[0]-expected/sum)<1e-12);
    Lookup duplicates({0,0,1},{-1,1,10},{1},1);double zero=0;
    const auto tied=duplicates.predict(&zero,2);require(tied.neighbors[0].row==0 && tied.values[0]==0.);
    bool caught=false;try{duplicates.predict(&zero,4);}catch(const std::invalid_argument&){caught=true;}require(caught);
}
