#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace core=dynsyslab::core;
namespace json=core::json;

struct Solution {
    std::vector<double> time,values,control,target;
    std::size_t samples{},states{};
    json::Value metadata;
};

void close(double actual,double expected,double tolerance,const std::string& label) {
    if (!std::isfinite(actual) || std::abs(actual-expected)>tolerance)
        throw std::runtime_error(label+" differs: actual="+std::to_string(actual)
            +", expected="+std::to_string(expected));
}

Solution run(const std::filesystem::path& repository,
             const std::filesystem::path& output,const std::string& jobName) {
    const auto job=json::load((repository/"examples"/"native"/"case-studies"/jobName).wstring());
    const auto result=dynsyslab::workflow::executeJob(job,output,{});
    if (result.kind!="solve" || result.archives.size()!=1)
        throw std::runtime_error(jobName+" did not produce one native solve archive");
    const auto archive=core::readNpz(result.archives.front());
    const auto& values=archive.at("values");
    if (values.shape.size()!=2)
        throw std::runtime_error(jobName+" values array is not sample-by-state");
    Solution solution;
    solution.samples=values.shape[0]; solution.states=values.shape[1];
    solution.time=archive.at("time").doubles(); solution.values=values.doubles();
    solution.control=archive.at("control_values").doubles();
    solution.target=archive.at("control_target_values").doubles();
    solution.metadata=json::parse(archive.metadataJson());
    if (solution.time.size()!=solution.samples
        || solution.values.size()!=solution.samples*solution.states
        || solution.control.size()!=solution.samples
        || solution.target.size()!=solution.values.size())
        throw std::runtime_error(jobName+" archive arrays are not aligned");
    if (solution.metadata.at("execution_backend").at("python_required").boolean(true)
        || solution.metadata.at("control_evaluation").string()!="continuous_rhs"
        || solution.metadata.at("solver_configuration").at("method").string()!="BDF")
        throw std::runtime_error(jobName+" did not preserve the native continuous-BDF contract");
    return solution;
}

std::size_t nearest(const std::vector<double>& time,double target) {
    const auto found=std::lower_bound(time.begin(),time.end(),target);
    if (found==time.end()) return time.size()-1;
    if (found==time.begin()) return 0;
    const auto right=static_cast<std::size_t>(found-time.begin()),left=right-1;
    return std::abs(time[left]-target)<=std::abs(time[right]-target)?left:right;
}

double percentile99(std::vector<double> values) {
    for (double& value:values) value=std::abs(value);
    std::sort(values.begin(),values.end());
    const double position=.99*static_cast<double>(values.size()-1);
    const auto left=static_cast<std::size_t>(std::floor(position));
    const auto right=static_cast<std::size_t>(std::ceil(position));
    return values[left]+(position-left)*(values[right]-values[left]);
}

void verifyDuffing(const Solution& solution) {
    if (solution.samples!=200001 || solution.states!=2)
        throw std::runtime_error("Duffing parity sample geometry changed");
    const auto release=nearest(solution.time,10.0),last=solution.samples-1;
    close(solution.values[release*2],-.2501115763410269,2e-8,"Duffing release x");
    close(solution.values[release*2+1],-1.9367425822554373e-6,2e-8,
          "Duffing release velocity");
    close(solution.values[last*2],-.4837885510,2e-8,"Duffing final x");
    close(solution.values[last*2+1],-.0120155367,2e-8,"Duffing final velocity");
    std::vector<double> active;
    for (std::size_t index=0;index<solution.samples;++index)
        if (solution.time[index]>=1.0 && solution.time[index]<=10.0)
            active.push_back(solution.control[index]);
    close(percentile99(active),.3123375429,2e-5,"Duffing active-control p99");
    close(*std::max_element(solution.control.begin(),solution.control.end(),
        [](double left,double right){return std::abs(left)<std::abs(right);}),
        -431.8161068471,.01,"Duffing engagement impulse");
}

void verifyDuffingApproximation(const Solution& solution,double releaseX,
                                double releaseVelocity,double finalX,
                                double finalVelocity,double peak,
                                const std::string& label) {
    if (solution.samples!=200001 || solution.states!=2)
        throw std::runtime_error(label+" Duffing sample geometry changed");
    const auto release=nearest(solution.time,10.0),last=solution.samples-1;
    close(solution.values[release*2],releaseX,2e-6,label+" release x");
    close(solution.values[release*2+1],releaseVelocity,2e-6,label+" release velocity");
    close(solution.values[last*2],finalX,2e-6,label+" final x");
    close(solution.values[last*2+1],finalVelocity,2e-6,label+" final velocity");
    double actualPeak=0.0;
    for (const auto value:solution.control) actualPeak=std::max(actualPeak,std::abs(value));
    close(actualPeak,peak,.15,label+" control peak");
}

std::pair<double,double> secondHalfRmse(const Solution& solution) {
    const std::size_t first=solution.samples/2; double sum0=0.0,sum1=0.0;
    for (std::size_t index=first;index<solution.samples;++index) {
        const double e0=solution.values[index*2]-solution.target[index*2];
        const double e1=solution.values[index*2+1]-solution.target[index*2+1];
        sum0+=e0*e0; sum1+=e1*e1;
    }
    const double count=static_cast<double>(solution.samples-first);
    return {std::sqrt(sum0/count),std::sqrt(sum1/count)};
}

void verifyRpe(const Solution& solution,double radius,double velocity,
               const std::string& label) {
    if (solution.samples!=10001 || solution.states!=2)
        throw std::runtime_error(label+" RPE sample geometry changed");
    if (solution.metadata.at("control_units").string()!="p_control/p0"
        || solution.metadata.at("integration_form").string()!="dimensionless")
        throw std::runtime_error(label+" lost dimensionless dynamics/control metadata");
    const auto error=secondHalfRmse(solution);
    close(error.first,radius,5e-6,label+" radius RMSE");
    close(error.second,velocity,5e-6,label+" velocity RMSE");
}

}  // namespace

int main(int argc,char** argv) {
    std::filesystem::path output;
    try {
        if (argc!=2) throw std::invalid_argument("repository root argument is required");
        const std::filesystem::path repository=argv[1];
        output=std::filesystem::temp_directory_path()
            /("dynsyslab-native-klqr-parity-"+std::to_string(
                std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        verifyDuffing(run(repository,output,
            "duffing-minimal-exact-klqr-transfer-release.json"));
        verifyDuffingApproximation(run(repository,output,
            "duffing-minimal-constant-klqr-transfer-release.json"),
            -.2452138856,-.0026970207,-.4835196200,-.0119669181,
            1057.1714082,"Duffing constant");
        verifyDuffingApproximation(run(repository,output,
            "duffing-minimal-bilinear-dynamics-klqr-transfer-release.json"),
            .8938695394,0.0,-.5404897285,.0168130656,
            929.8443386,"Duffing bilinear dynamics");
        verifyDuffingApproximation(run(repository,output,
            "duffing-minimal-bilinear-map-klqr-transfer-release.json"),
            .8956283868,0.0,-.5397237012,.0212784038,
            984.9279435,"Duffing bilinear map");
        verifyRpe(run(repository,output,"rpe-minimal-exact-klqr-moderate.json"),
                  .00121273,.00412142,"exact");
        verifyRpe(run(repository,output,"rpe-minimal-constant-klqr-moderate.json"),
                  .00118033,.00397416,"constant");
        verifyRpe(run(repository,output,"rpe-minimal-bilinear-klqr-moderate.json"),
                  .00118979,.00400489,"bilinear");
        std::filesystem::remove_all(output);
        std::cout<<"dynsyslab_klqr_case_study_parity_tests: PASS\n"; return 0;
    } catch (const std::exception& error) {
        if (!output.empty()) { std::error_code ignored; std::filesystem::remove_all(output,ignored); }
        std::cerr<<"dynsyslab_klqr_case_study_parity_tests: FAIL: "<<error.what()<<'\n';
        return 1;
    }
}
