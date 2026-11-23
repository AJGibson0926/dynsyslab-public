#include <dynsyslab/math/DiscreteCycles.hpp>
#include <dynsyslab/math/Query.hpp>
#include <dynsyslab/math/HarmonicAnalysis.hpp>
#include <dynsyslab/math/CyclicCocycle.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace { void require(bool x,const char* message){if(!x)throw std::runtime_error(message);} }
int main(){try {
    using namespace dynsyslab;const double pi=std::acos(-1.);
    std::vector<double> x;for(int i=0;i<24;++i)x.push_back(double(i%3));
    auto r=math::diagnoseSampledCycles(x,1,{1,2,3,6},{1},{0},1e-12);
    require(r.detectedPeriod==3&&r.foldingRms==0,"arbitrary non-dyadic period");
    x[13]+=.1;r=math::diagnoseSampledCycles(x,1,{1,2,3,6},{1},{0},1e-4);
    require(r.detectedPeriod==0,"reject corrupted recurrence");
    x.clear();for(int i=0;i<16;++i)x.push_back(.3+2*pi*i);
    r=math::diagnoseSampledCycles(x,1,{1,2,4},{1},{2*pi},1e-12);
    require(r.detectedPeriod==1,"circular coordinate distance");
    bool rejected=false;try{math::diagnoseSampledCycles(x,1,{8},{1},{0},1e-6);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"insufficient repetition must be rejected");
    const auto s=math::discreteCycleSpectrum({1,0,-1,0},1);
    require(std::abs(s.at("total_mass").number()-.5)<1e-12,"Parseval mass");
    require(std::abs(s.at("masses").array()[1].number()-.25)<1e-12,"cyclic atom amplitude");
    require(s.at("reduced_denominators").array()[1].number()==4,"dyadic atom order");
    require(std::abs(s.at("conductor_masses").array().back().at("mass").number()-.5)<1e-12,"conductor mass groups both conjugate characters");
    auto q=core::json::parse(R"({"schema":"dynsyslab.math-query/v1","operation":"sampled_hankel_pairs","values":[0,1,2,3,4,5,6,7],"dimension":1,"delays":3,"snapshot_stride":2})");
    const auto h=math::evaluateQuery(q).at("result");
    require(h.at("samples").number()==3&&h.at("x").array()[3].number()==2,"Hankel block layout");
    require(h.at("y").array()[3].number()==3,"pair shift is not snapshot stride");
    x.clear();for(int j=0;j<6;++j)x.push_back(std::cos(2*pi*j/6)+2*std::cos(4*pi*j/6));
    const auto projection=math::projectDiscreteCycle(x,1,3);
    for(int j=0;j<6;++j)require(std::abs(projection.at("projected_values").array()[j].number()-2*std::cos(4*pi*j/6))<1e-12,"quotient character selection");
    require(std::abs(projection.at("rms_error").number()-std::sqrt(.5))<1e-12,"projection Parseval error");
    rejected=false;try{math::projectDiscreteCycle(x,1,4);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"reject non-divisor quotient");
    const auto alignment=math::alignDiscreteCycles({1,2,4},{4,1,2},1);
    require(alignment.at("shift").number()==1&&alignment.at("rms_error").number()==0,"recover forward cyclic phase");
    auto b=math::evaluateQuery(core::json::parse(R"({"schema":"dynsyslab.math-query/v1","operation":"fourier","action":"evaluate_many","series":{"angular_frequency":2,"coefficients":[[1,0.5,0],[-1,0.5,0]]},"phases":[0,1,2]})")).at("result");
    for(int j=0;j<3;++j)require(std::abs(b.at("values").array()[j].array()[0].number()-std::cos(j))<1e-12,"batched Fourier analytic reference");
    const auto evolution=math::evolveFourierSeries({2,{{1,.5},{-1,.5}}},{3,4,5},.2,3);
    for(int j=0;j<3;++j)require(std::abs(evolution[j].real()-std::cos(.2+2*j))<1e-12,"physical clock and phase anchoring");
    auto cocycle=math::cyclicCocycle({1,1,0,1, 1,0,2,1},2,2,{1,0},2);
    const auto& product=cocycle.at("monodromy").array();
    require(product[0].number()==1&&product[1].number()==1&&product[2].number()==2&&product[3].number()==3,"noncommuting cocycle product order");
    require(cocycle.at("deviations").array().back().number()==2,"cocycle rollout matches product");
    require(cocycle.at("character_synthesis_maximum_error").number()<1e-12,"matrix character synthesis convention");
    auto constant=math::cyclicCocycle({1,1,0,1, 1,0,2,1},2,1,{1,0},2);
    require(std::abs(constant.at("monodromy").array()[0].number()-1.5)<1e-12,"quotient product spans original cycle");
    auto shifted=math::cyclicCocycle({1,1,0,1, 1,0,2,1},2,2,{1,0},2,1);
    require(std::abs(shifted.at("monodromy").array()[0].number()-3)<1e-12,"nonzero phase preserves ordered product");
    auto partition=math::cyclicPartitionGeometry({0,1,.1,1.1},1,2);
    require(std::abs(partition.at("levels").array()[1].at("maximum_diameter").number()-.1)<1e-12,"nested residue cell geometry");
    require(!partition.at("certified_infinite_limit").boolean(),"finite geometry must never certify infinity");
    partition=math::cyclicPartitionGeometry({0,1,.1,1.1},1,2,{1,.1,1.1,.01});
    require(std::abs(partition.at("sample_transition").at("maximum_error").number()-.01)<1e-12,"independent cyclic image defect");
    rejected=false;try{math::cyclicPartitionGeometry({0,1,2},1,2);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"reject incompatible base-power hierarchy");
    x.clear();for(int j=0;j<12;++j)x.push_back(j%3+.01*j);
    auto tower=math::cyclicTowerPartitionGeometry(x,1,{1,3,6,12});
    require(tower.at("levels").array()[1].at("branching_factor").number()==3,"odd-base initial branch");
    require(std::abs(tower.at("levels").array()[1].at("maximum_diameter").number()-.09)<1e-12,"mixed tower residue geometry");
    require(tower.at("levels").array()[3].at("singleton_terminal_level").boolean(),"mixed tower singleton disclaimer");
    tower=math::cyclicTowerPartitionGeometry(x,1,{3,6});
    require(!tower.at("levels").array().back().at("singleton_terminal_level").boolean(),"partial hierarchy retains nontrivial cells");
    rejected=false;try{math::cyclicTowerPartitionGeometry(x,1,{2,3,12});}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"non-nested divisors cannot define a tower");
    auto mixed=math::evaluateQuery(core::json::parse(R"({"schema":"dynsyslab.math-query/v1","operation":"cycle_partition_geometry","dimension":1,"values":[0,1,2,0.1,1.1,2.1],"partition_periods":[1,3,6]})")).at("result");
    require(mixed.at("partition_periods").array().size()==3&&!mixed.at("certified_infinite_limit").boolean(),"mixed tower query contract");
    rejected=false;try{math::evaluateQuery(core::json::parse(R"({"schema":"dynsyslab.math-query/v1","operation":"cycle_partition_geometry","dimension":1,"values":[0,1],"base":2,"partition_periods":[1,2]})"));}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"ambiguous tower configuration rejected");
    std::cout<<"discrete cycles: PASS\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
