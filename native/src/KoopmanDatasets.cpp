#include <dynsyslab/core/KoopmanDatasets.hpp>
#include <algorithm>
#include <stdexcept>

namespace dynsyslab::core {
void mergeKoopmanDatasetReference(std::vector<KoopmanDatasetReference>& refs,
                                  KoopmanDatasetReference ref) {
    if(ref.path.empty() || (ref.role!="training" && ref.role!="validation" && ref.role!="test")
        || (ref.dynamics!="forced" && ref.dynamics!="unforced" && ref.dynamics!="unknown"))
        throw std::invalid_argument("KROM dataset reference requires a path, training/validation/test role, and forced/unforced/unknown dynamics");
    const auto found=std::find_if(refs.begin(),refs.end(),[&](const auto& item){
        return item.path==ref.path && item.role==ref.role && item.inputPath==ref.inputPath;});
    if(found==refs.end()) refs.push_back(std::move(ref));
    else if(found->dynamics=="unknown") found->dynamics=ref.dynamics;
    else if(ref.dynamics!="unknown" && found->dynamics!=ref.dynamics)
        throw std::invalid_argument("Conflicting forcing classifications for a KROM dataset reference");
}

json::Value koopmanDatasetReferencesJson(const std::vector<KoopmanDatasetReference>& refs) {
    std::vector<KoopmanDatasetReference> checked;
    for(const auto& ref:refs) mergeKoopmanDatasetReference(checked,ref);
    json::Value::Array rows;
    for(const auto& ref:checked) {
        json::Value::Object row{{"path",ref.path},{"role",ref.role},{"dynamics",ref.dynamics}};
        if(!ref.inputPath.empty()) row["input_archive"]=ref.inputPath;
        rows.emplace_back(std::move(row));
    }
    return json::Value::Object{{"schema","dynsyslab.krom-datasets/v1"},{"references",std::move(rows)}};
}

std::vector<KoopmanDatasetReference> koopmanDatasetReferencesFromMetadata(const json::Value& meta) {
    std::vector<KoopmanDatasetReference> result;
    const auto& encoded=meta.at("dataset_references");
    if(!encoded.isNull()) {
        if(encoded.at("schema").string()!="dynsyslab.krom-datasets/v1")
            throw std::invalid_argument("Unsupported KROM dataset-reference schema");
        for(const auto& row:encoded.at("references").array())
            mergeKoopmanDatasetReference(result,{row.at("path").string(),row.at("role").string(),
                row.at("dynamics").string("unknown"),row.at("input_archive").string()});
        return result; // The explicit contract is authoritative, not stale legacy keys.
    }
    for(const auto* role:{"training","validation","test"}) {
        for(const auto& row:meta.at(std::string(role)+"_archives").array())
            mergeKoopmanDatasetReference(result,{row.string(),role,"unknown",{}});
        for(const auto& row:meta.at(std::string(role)+"_sources").array()) {
            const auto path=row.string(row.at("path").string(row.at("states").string()));
            if(!path.empty()) mergeKoopmanDatasetReference(result,{path,role,
                row.at("dynamics").string("unknown"),row.at("inputs").string()});
        }
    }
    const auto method=meta.at("method").string();
    const bool trainingMethod=method=="hdmd" || method=="hankel-dmd" || method=="hankel_dmd"
        || method=="koopman-with-sindy" || method=="sindy" || method=="edmdc" || method=="dmdc"
        || method=="echo-history-koopman" || method=="boundary-eigenfunctions";
    if(trainingMethod && !meta.at("source_archive").string().empty()
        && std::none_of(result.begin(),result.end(),[](const auto& ref){return ref.role=="training";}))
        mergeKoopmanDatasetReference(result,{meta.at("source_archive").string(),"training","unknown",{}});
    return result;
}
} // namespace dynsyslab::core
