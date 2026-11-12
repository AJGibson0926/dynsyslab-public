#include <dynsyslab/workflow/DecompositionPostprocess.hpp>

#include <dynsyslab/core/Text.hpp>

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Decomposition.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {
std::string escape(const std::string& value) {
    std::string result;
    for (const char item:value) {
        if (item=='"'||item=='\\') result.push_back('\\');
        if (item=='\n') result+="\\n"; else if (item!='\r') result.push_back(item);
    }
    return result;
}
std::string timestamp() {
    const auto now=std::chrono::system_clock::now();
    const std::time_t raw=std::chrono::system_clock::to_time_t(now); std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc,&raw);
#else
    gmtime_r(&raw,&utc);
#endif
    std::ostringstream out; out<<std::put_time(&utc,"%Y%m%dT%H%M%S")<<'Z'; return out.str();
}
std::string safe(std::string value) {
    for (char& item:value) if (!std::isalnum(static_cast<unsigned char>(item))
        && item!='-'&&item!='_') item='-';
    return value;
}
std::string nativeSystem(std::string value) {
    std::transform(value.begin(),value.end(),value.begin(),[](unsigned char item) {
        return static_cast<char>(std::tolower(item));
    });
    if (value.find("rayleigh")!=std::string::npos||value=="rpe") return "rpe";
    if (value.find("duffing")!=std::string::npos) return "duffing";
    if (value.find("marmottant")!=std::string::npos) return "marmottant";
    if (value.find("reddy")!=std::string::npos) return "reddy_szeri";
    for (const auto& descriptor:dynsyslab::core::systemCatalog())
        if (value==descriptor.identifier) return value;
    return {};
}
dynsyslab::core::SystemConfiguration configuration(const std::string& text) {
    if (text.empty()) throw std::invalid_argument("named decomposition observable requires system metadata");
    const auto metadata=json::parse(text); dynsyslab::core::SystemConfiguration result;
    result.identifier=nativeSystem(metadata.at("system").string());
    result.preset=metadata.at("system_preset").string(metadata.at("preset").string());
    for (const auto& item:metadata.at("system_parameters").object())
        result.parameters.emplace(item.first,item.second.number());
    if (result.identifier.empty()) throw std::invalid_argument("source system is not native");
    return result;
}
std::filesystem::path destination(const std::filesystem::path& root,
    const std::filesystem::path& source,const std::string& method) {
    return root/L"derived"/dynsyslab::core::wide(method)/dynsyslab::core::wide(safe(source.stem().u8string()
        +"__"+method+"__"+timestamp())+".npz");
}
std::string metadata(const std::filesystem::path& source,const std::string& method,
    const std::string& label,const std::string& details) {
    return "{\"format\":\"dynsyslab.derived/v1\",\"label\":\""+escape(label)
        +"\",\"method\":\""+method+"\",\"source_archive\":\""
        +escape(source.u8string())+"\",\"execution_backend\":{\"name\":\"native-cpp\","
         "\"compiled\":true,\"python_required\":false},"+details+"}";
}
double deltaTime(const std::vector<double>& time) {
    if (time.size()<3) throw std::invalid_argument("decomposition requires three samples");
    const double dt=time[1]-time[0];
    for (std::size_t i=2;i<time.size();++i)
        if (!(dt>0.0)||std::abs((time[i]-time[i-1])-dt)>1.0e-10*std::max(1.0,std::abs(dt)))
            throw std::invalid_argument("decomposition requires uniform time steps");
    return dt;
}
}  // namespace

std::filesystem::path runDecompositionPostprocess(
    const std::string& method,const std::filesystem::path& sourcePath,
    const json::Value& settings,const std::filesystem::path& root,const std::string& label) {
    const auto source=dynsyslab::core::readNpz(sourcePath);
    const auto& stateArray=source.contains("values")?source.at("values"):source.at("states");
    const auto& timeArray=source.at("time");
    if (stateArray.shape.size()!=2||timeArray.shape.size()!=1
        ||stateArray.shape[0]!=timeArray.shape[0])
        throw std::invalid_argument("decomposition source must be samples by state");
    const auto states=stateArray.doubles(),time=timeArray.doubles();
    const std::size_t samples=stateArray.shape[0],dimension=stateArray.shape[1];
    dynsyslab::core::NpzRecord output; std::string details;
    if (method=="edmd") {
        dynsyslab::core::EDMDOptions options;
        options.regularization=settings.at("regularization").number(0.0);
        const std::string kind=settings.at("dictionary").string("polynomial");
        if (kind=="polynomial") {
            options.dictionary.kind=dynsyslab::core::DictionaryKind::polynomial;
            options.dictionary.polynomialDegree=static_cast<std::size_t>(
                settings.at("polynomial_degree").number(2));
            options.dictionary.includeConstant=settings.at("include_constant").boolean(true);
        } else if (kind=="rbf") {
            options.dictionary.kind=dynsyslab::core::DictionaryKind::radialBasis;
            options.dictionary.radialBasisWidth=settings.at("width").number(1.0);
            for (const auto& center:settings.at("centers").array()) {
                if (center.array().size()!=dimension) throw std::invalid_argument("EDMD RBF center dimension mismatch");
                for (const auto& value:center.array()) options.dictionary.radialBasisCenters.push_back(value.number());
                ++options.dictionary.centerCount;
            }
        } else throw std::invalid_argument("EDMD dictionary must be polynomial or rbf");
        const auto result=dynsyslab::core::extendedDynamicModeDecomposition(
            states,samples,dimension,deltaTime(time),options);
        output.arrays["eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
            {result.featureCount},result.discreteEigenvalues);
        output.arrays["discrete_eigenvalues"]=output.arrays["eigenvalues"];
        output.arrays["continuous_eigenvalues"]=dynsyslab::core::NpyArray::fromComplex(
            {result.featureCount},result.continuousEigenvalues);
        output.arrays["modes"]=dynsyslab::core::NpyArray::fromComplex(
            {dimension,result.featureCount},result.modes);
        output.arrays["eigenfunction_values"]=dynsyslab::core::NpyArray::fromComplex(
            {samples,result.featureCount},result.eigenfunctionValues);
        output.arrays["eigenfunction_coefficients"]=dynsyslab::core::NpyArray::fromComplex(
            {result.featureCount,result.featureCount},result.eigenfunctionCoefficients);
        output.arrays["operator"]=dynsyslab::core::NpyArray::fromDoubles(
            {result.featureCount,result.featureCount},result.operatorMatrix);
        details="\"dictionary\":\""+kind+"\",\"feature_count\":"
            +std::to_string(result.featureCount)+",\"regularization\":"
            +std::to_string(options.regularization)+",\"eigenvalue_domain\":\"discrete-time\"";
    } else if (method=="havok") {
        dynsyslab::core::HAVOKOptions options;
        options.delays=static_cast<std::size_t>(settings.at("delays").number(100));
        options.rank=static_cast<std::size_t>(settings.at("rank").number(15));
        options.observableIndex=static_cast<std::size_t>(settings.at("observable").number());
        const std::string observableName=settings.at("observable").string();
        dynsyslab::core::HAVOKResult result;
        if (observableName.empty())
            result=dynsyslab::core::havokDecomposition(states,time,samples,dimension,options);
        else {
            auto system=dynsyslab::core::createSystem(configuration(source.metadataJson()));
            const auto catalog=system->observableCatalog();
            result=dynsyslab::core::havokDecomposition(
                states,time,samples,dimension,catalog.at(observableName),options);
        }
        const std::size_t autonomous=result.rank-1;
        output.arrays["A"]=dynsyslab::core::NpyArray::fromDoubles(
            {autonomous,autonomous},result.generator);
        output.arrays["B"]=dynsyslab::core::NpyArray::fromDoubles({autonomous,1},result.input);
        output.arrays["forcing"]=dynsyslab::core::NpyArray::fromDoubles(
            {result.embeddedSamples},result.forcing);
        output.arrays["coordinates"]=dynsyslab::core::NpyArray::fromDoubles(
            {result.embeddedSamples,result.rank},result.coordinates);
        output.arrays["values"]=output.arrays["coordinates"];
        output.arrays["time"]=dynsyslab::core::NpyArray::fromDoubles(
            {result.embeddedSamples},result.coordinateTime);
        output.arrays["left_singular_vectors"]=dynsyslab::core::NpyArray::fromDoubles(
            {result.delays,result.rank},result.leftSingularVectors);
        output.arrays["singular_values"]=dynsyslab::core::NpyArray::fromDoubles(
            {result.rank},result.singularValues);
        output.arrays["right_singular_vectors"]=dynsyslab::core::NpyArray::fromDoubles(
            {result.rank,result.embeddedSamples},result.rightSingularVectors);
        output.arrays["observable"]=dynsyslab::core::NpyArray::fromDoubles({samples},result.observable);
        details="\"delays\":"+std::to_string(result.delays)+",\"rank\":"
            +std::to_string(result.rank)+",\"observable\":\""
            +escape(observableName.empty()?std::to_string(options.observableIndex):observableName)+"\"";
    } else throw std::invalid_argument("unknown decomposition stage");
    output.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
        metadata(sourcePath,method,label,details));
    const auto path=destination(root,sourcePath,method);
    dynsyslab::core::writeNpzAtomic(path,output); return path;
}

}  // namespace dynsyslab::workflow
