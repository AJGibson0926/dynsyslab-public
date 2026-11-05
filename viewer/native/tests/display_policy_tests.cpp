#include "dynsyslab_viewer/NativePayload.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}

void highDimensionalTrajectoryViews(){const auto path=std::filesystem::temp_directory_path()
    /"dynsyslab-native-high-dimensional-view-test.npz";dynsyslab::core::NpzRecord archive;
    archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles({8},
        {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5});std::vector<double> values(8*12);
    for(std::size_t index=0;index<values.size();++index)values[index]=static_cast<double>(index);
    archive.arrays["values"]=dynsyslab::core::NpyArray::fromDoubles({8,12},values);
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"Reddy-Szeri trajectory","system":"reddy_szeri","integration_form":"dimensionless"})");
    dynsyslab::core::writeNpzAtomic(path,archive);const auto time=dsl::buildNativePayload(path,"time");
    require(time.figures.size()==12,"all high-dimensional state components are visible");
    require(time.figures[0].title=="r = \\frac{R}{R_{0}} versus time","native state-component labels");
    const auto three=dsl::buildNativePayload(path,"three_d");require(three.figures.size()==1
        &&!three.figures[0].traces[0].z.empty()
        &&three.figures[0].xLabel.find("r = \\frac{R}{R_{0}}")!=std::string::npos
        &&three.figures[0].yLabel.find("\\frac{dr}{d\\tau}")!=std::string::npos,
        "native three-state trajectory projection and scientific labels");
    const auto fourier=dsl::buildNativePayload(path,"fourier");require(fourier.figures.size()==12
        &&fourier.figures[0].xLabel.find("tau")!=std::string::npos
        &&fourier.figures[0].yLabel=="amplitude",
        "native all-component Fourier view");std::error_code ignored;std::filesystem::remove(path,ignored);}

void automaticDetailInset(){dsl::Figure transient{"time","control","time","control",false};
    dsl::Trace trace;trace.label="u";for(int index=0;index<1000;++index){trace.x.push_back(index*.01);
        trace.y.push_back(index<5?1000.0:std::sin(index*.03));}transient.traces.push_back(trace);
    const auto inset=dsl::automaticDetailInset(transient);require(inset&&inset->automatic
        &&inset->yMax<5.0&&inset->yMin>-5.0,"native automatic transient inset");
    dsl::Figure ordinary=transient;for(std::size_t index=0;index<ordinary.traces[0].y.size();++index)
        ordinary.traces[0].y[index]=std::sin(index*.03);require(!dsl::automaticDetailInset(ordinary),
        "ordinary signal should not get an inset");ordinary.kind="phase";ordinary.traces.push_back(trace);
    require(!dsl::automaticDetailInset(ordinary),"phase figures should not get time insets");}

void allDecompositionModesRemainVisible(){const auto path=std::filesystem::temp_directory_path()
    /"dynsyslab-native-decomposition-display-test.npz";dynsyslab::core::NpzRecord archive;
    archive.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles({3},{0.0,1.0,2.0});
    archive.arrays["singular_values"]=dynsyslab::core::NpyArray::fromDoubles(
        {6},{6.0,5.0,4.0,3.0,2.0,1.0});std::vector<double> scores(18);
    for(std::size_t index=0;index<scores.size();++index)scores[index]=static_cast<double>(index);
    archive.arrays["scores"]=dynsyslab::core::NpyArray::fromDoubles({3,6},scores);
    archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        R"({"label":"six-mode PCA","method":"pca"})");
    dynsyslab::core::writeNpzAtomic(path,archive);
    const auto payload=dsl::buildNativePayload(path,"data_overview");
    require(payload.figures.size()>=2&&payload.figures.back().traces.size()==6,
        "decomposition display silently truncated retained modes");
    std::error_code ignored;std::filesystem::remove(path,ignored);}
}

int main(){try{highDimensionalTrajectoryViews();automaticDetailInset();
        allDecompositionModesRemainVisible();
        std::cout<<"dynsyslab_viewer_display_policy_tests: PASS\n";return 0;}
    catch(const std::exception& error){std::cerr<<"dynsyslab_viewer_display_policy_tests: FAIL: "
        <<error.what()<<'\n';return 1;}}
