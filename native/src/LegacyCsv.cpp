#include <dynsyslab/core/LegacyCsv.hpp>

#include <dynsyslab/core/FileDigest.hpp>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using Row=std::vector<std::string>;

std::string trim(std::string value) {
    const auto visible=[](unsigned char item){return !std::isspace(item);};
    value.erase(value.begin(),std::find_if(value.begin(),value.end(),visible));
    value.erase(std::find_if(value.rbegin(),value.rend(),visible).base(),value.end());
    return value;
}

std::vector<Row> csvRows(const std::filesystem::path& path) {
    std::ifstream stream(path,std::ios::binary);if(!stream)throw std::runtime_error("legacy CSV does not exist: "+path.u8string());
    std::ostringstream buffer;buffer<<stream.rdbuf();std::string text=buffer.str();
    if(text.size()>=3&&static_cast<unsigned char>(text[0])==0xef&&static_cast<unsigned char>(text[1])==0xbb&&static_cast<unsigned char>(text[2])==0xbf)text.erase(0,3);
    std::vector<Row> rows;Row row;std::string cell;bool quoted=false;
    for(std::size_t index=0;index<text.size();++index){const char item=text[index];
        if(quoted){if(item=='"'&&index+1<text.size()&&text[index+1]=='"'){cell.push_back('"');++index;}else if(item=='"')quoted=false;else cell.push_back(item);continue;}
        if(item=='"'&&cell.empty()){quoted=true;continue;}if(item==','){row.push_back(cell);cell.clear();continue;}
        if(item=='\r'||item=='\n'){row.push_back(cell);cell.clear();rows.push_back(std::move(row));row.clear();if(item=='\r'&&index+1<text.size()&&text[index+1]=='\n')++index;continue;}cell.push_back(item);}
    if(quoted)throw std::invalid_argument("legacy CSV ends inside a quoted field");
    if(!cell.empty()||!row.empty()){row.push_back(cell);rows.push_back(std::move(row));}
    return rows;
}

std::vector<std::string> textPrefix(const Row& row) {
    std::vector<std::string> result;for(std::size_t index=1;index<row.size();++index){auto value=trim(row[index]);if(value.empty())break;result.push_back(std::move(value));}return result;
}

std::vector<double> numericPrefix(const Row& row,std::size_t line,const char* record) {
    std::vector<double> result;bool empty=false;for(std::size_t index=1;index<row.size();++index){const auto token=trim(row[index]);
        if(token.empty()){empty=true;continue;}if(empty)throw std::invalid_argument(std::string(record)+" record at line "+std::to_string(line)+" has data after an empty column");
        std::size_t used=0;const double value=std::stod(token,&used);if(used!=token.size()||!std::isfinite(value))throw std::invalid_argument(std::string(record)+" record at line "+std::to_string(line)+" contains non-numeric data");result.push_back(value);}return result;
}

std::size_t integer(const std::string& value,const char* field) {
    const double number=std::stod(value);if(!std::isfinite(number)||number<0||number!=std::floor(number))throw std::invalid_argument(std::string(field)+" must be a non-negative integer");return static_cast<std::size_t>(number);
}

void increasing(const std::vector<double>& time,std::size_t line) {
    if(time.size()<2)throw std::invalid_argument("time record at line "+std::to_string(line)+" needs at least two samples");
    for(std::size_t index=1;index<time.size();++index)if(!(time[index]>time[index-1]))throw std::invalid_argument("time record at line "+std::to_string(line)+" must be strictly increasing");
}

struct Builder {
    LegacyCsvSystem system;std::optional<std::size_t> cluster,trajectory;
    std::vector<double> trajectoryTime;
};

LegacyCsvSystem finish(Builder&& builder) {
    auto& system=builder.system;if(system.trajectories.empty())throw std::invalid_argument("legacy system "+std::to_string(system.index)+" contains no trajectories");
    if(!(system.declaredStep>0)||system.declaredTransitions<1)throw std::invalid_argument("legacy system "+std::to_string(system.index)+" has invalid declared time geometry");
    if(!system.time.empty()&&system.time.size()!=system.declaredTransitions+1)throw std::invalid_argument("legacy system "+std::to_string(system.index)+" declared transition count differs from stored time");
    if(!system.controlValues.empty()&&system.controlTime.size()!=system.controlValues.size())throw std::invalid_argument("legacy system "+std::to_string(system.index)+" has mismatched control time and values");
    return std::move(system);
}

json::Value strings(const std::vector<std::string>& values) {
    json::Value::Array output;for(const auto& value:values)output.emplace_back(value);return json::Value(std::move(output));
}

}  // namespace

std::size_t LegacyCsvBatch::trajectoryCount() const {
    std::size_t result=0;for(const auto& system:systems)result+=system.trajectories.size();return result;
}

json::Value LegacyCsvBatch::summary() const {
    std::set<std::string> names;std::set<std::size_t> dimensions;std::size_t controlled=0;
    for(const auto& system:systems){names.insert(system.dynamicsName);if(!system.controlValues.empty())++controlled;for(const auto& trajectory:system.trajectories)dimensions.insert(trajectory.stateDimension);}
    json::Value::Array dynamics;for(const auto& name:names)dynamics.emplace_back(name);json::Value::Array stateDimensions;for(auto dimension:dimensions)stateDimensions.emplace_back(static_cast<double>(dimension));
    json::Value::Object result;result["schema"]="dynsyslab.legacy-matlab-dynbatch-csv/v1";result["source_path"]=source.generic_u8string();result["source_sha256"]=sourceSha256;
    result["declared_system_count"]=declaredSystemCount?json::Value(static_cast<double>(*declaredSystemCount)):json::Value(nullptr);result["decoded_system_count"]=static_cast<double>(systems.size());
    result["trajectory_count"]=static_cast<double>(trajectoryCount());result["controlled_system_count"]=static_cast<double>(controlled);result["dynamics_names"]=json::Value(std::move(dynamics));result["state_dimensions"]=json::Value(std::move(stateDimensions));result["warnings"]=strings(warnings);return json::Value(std::move(result));
}

LegacyCsvBatch readLegacyDynBatchCsv(const std::filesystem::path& path,bool strict) {
    const auto source=std::filesystem::absolute(path);const auto rows=csvRows(source);if(rows.size()<2)throw std::invalid_argument("legacy CSV is empty or contains only a header");
    const std::set<std::string> tags{"dynbatch","dynsys","dynclst","inst","time","sol","ctrl","u"};
    LegacyCsvBatch result;result.source=source;result.sourceSha256=sha256File(source);result.header=rows.front();std::optional<Builder> current;
    for(std::size_t index=1;index<rows.size();++index){const Row& row=rows[index];const std::size_t line=index+1;const std::string tag=row.empty()?"":trim(row.front());
        if(tag.empty()){if(strict)throw std::invalid_argument("Unexpected continuation row at line "+std::to_string(line));result.warnings.push_back("Ignored unattached continuation row "+std::to_string(line));continue;}
        if(!tags.count(tag)){const std::string warning="Unknown legacy record tag '"+tag+"' at line "+std::to_string(line);if(strict)throw std::invalid_argument(warning);result.warnings.push_back(warning);continue;}
        if(tag=="dynbatch"){const auto fields=textPrefix(row);if(!fields.empty())result.declaredSystemCount=integer(fields[0],"dynbatch system count");continue;}
        if(tag=="dynsys"){if(current)result.systems.push_back(finish(std::move(*current)));const auto fields=textPrefix(row);if(fields.size()<4)throw std::invalid_argument("dynsys record at line "+std::to_string(line)+" needs index, dynamics, dt, and transitions");
            Builder builder;builder.system.index=integer(fields[0],"dynsys index");builder.system.dynamicsName=fields[1];builder.system.declaredStep=std::stod(fields[2]);builder.system.declaredTransitions=integer(fields[3],"dynsys transitions");builder.system.trailingFields.assign(fields.begin()+4,fields.end());builder.system.sourceLine=line;current=std::move(builder);continue;}
        if(!current)throw std::invalid_argument(tag+" record at line "+std::to_string(line)+" precedes the first dynsys record");auto& builder=*current;
        if(tag=="dynclst"){const auto fields=textPrefix(row);if(fields.empty())throw std::invalid_argument("dynclst record has no index");builder.cluster=integer(fields[0],"cluster index");builder.trajectory.reset();builder.trajectoryTime.clear();}
        else if(tag=="inst"){const auto fields=textPrefix(row);if(!builder.cluster||fields.empty())throw std::invalid_argument("inst record has no active cluster/index");builder.trajectory=integer(fields[0],"trajectory index");builder.trajectoryTime.clear();}
        else if(tag=="time"){auto values=numericPrefix(row,line,"time");increasing(values,line);if(builder.trajectory)builder.trajectoryTime=std::move(values);else builder.system.time=std::move(values);}
        else if(tag=="ctrl"){const auto fields=textPrefix(row);builder.system.forcing=fields.empty()?"off":fields[0];builder.system.subforcing=fields.size()>1?fields[1]:"off";}
        else if(tag=="u"){builder.system.controlValues=numericPrefix(row,line,"u");builder.system.controlTime=builder.system.time;}
        else if(tag=="sol"){if(!builder.cluster||!builder.trajectory)throw std::invalid_argument("sol record at line "+std::to_string(line)+" has no active cluster/trajectory");
            std::vector<std::vector<double>> components{numericPrefix(row,line,"sol")};std::size_t endLine=line,lookahead=index+1;
            while(lookahead<rows.size()){const auto& next=rows[lookahead];const std::string nextTag=next.empty()?"":trim(next.front());if(!nextTag.empty())break;components.push_back(numericPrefix(next,lookahead+1,"sol-continuation"));endLine=lookahead+1;++lookahead;}
            const std::size_t samples=components.front().size();if(samples<2)throw std::invalid_argument("sol block has fewer than two samples");for(const auto& component:components)if(component.size()!=samples)throw std::invalid_argument("sol block has inconsistent component lengths");
            const auto& time=builder.trajectoryTime.empty()?builder.system.time:builder.trajectoryTime;if(time.size()!=samples)throw std::invalid_argument("sol block has no matching time record");
            LegacyCsvTrajectory trajectory;trajectory.systemIndex=builder.system.index;trajectory.clusterIndex=*builder.cluster;trajectory.trajectoryIndex=*builder.trajectory;trajectory.time=time;trajectory.samples=samples;trajectory.stateDimension=components.size();trajectory.firstSourceLine=line;trajectory.lastSourceLine=endLine;trajectory.states.resize(samples*components.size());
            for(std::size_t sample=0;sample<samples;++sample)for(std::size_t axis=0;axis<components.size();++axis)trajectory.states[sample*components.size()+axis]=components[axis][sample];builder.system.trajectories.push_back(std::move(trajectory));index=lookahead-1;}
    }
    if(current)result.systems.push_back(finish(std::move(*current)));if(result.systems.empty())throw std::invalid_argument("no dynsys records were decoded from the legacy CSV");
    if(result.declaredSystemCount&&*result.declaredSystemCount!=result.systems.size()){const std::string warning="dynbatch declares "+std::to_string(*result.declaredSystemCount)+" systems but "+std::to_string(result.systems.size())+" were decoded";if(strict)throw std::invalid_argument(warning);result.warnings.push_back(warning);}
    return result;
}

}  // namespace dynsyslab::core
