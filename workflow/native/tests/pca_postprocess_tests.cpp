#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace core=dynsyslab::core;
namespace json=core::json;

void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::filesystem::path sourceArchive(const std::filesystem::path& root) {
    constexpr std::size_t samples=80,states=2;
    std::vector<double> time(samples),values(samples*states);
    for (std::size_t index=0;index<samples;++index) {
        time[index]=.01*static_cast<double>(index);
        values[index*states]=std::sin(2.0*time[index]);
        values[index*states+1]=2.0*values[index*states];
    }
    values[21]+=6.0;
    core::NpzRecord archive;
    archive.arrays["time"]=core::NpyArray::fromDoubles({samples},time);
    archive.arrays["values"]=core::NpyArray::fromDoubles({samples,states},values);
    archive.arrays["metadata_json"]=core::NpyArray::fromString(
        "{\"system\":\"duffing\",\"system_preset\":\"gibson_damped_stabilization_2022\"}");
    const auto path=root/"source.npz"; core::writeNpzAtomic(path,archive); return path;
}

json::Value job(const std::filesystem::path& source) {
    json::Value::Array methods;
    methods.emplace_back(json::Value::Object{{"name","bagged_pca"},{"enabled",true},
        {"settings",json::Value::Object{{"rank",1.0},{"bags",8.0},{"seed",42.0}}}});
    methods.emplace_back(json::Value::Object{{"name","robust_pca"},{"enabled",true},
        {"settings",json::Value::Object{{"tolerance",1.0e-7},{"max_iterations",1000.0}}}});
    methods.emplace_back(json::Value::Object{{"name","hankel_pod"},{"enabled",true},
        {"settings",json::Value::Object{{"matrix_size",8.0},{"delay_stride",2.0},
                                         {"rank",3.0},{"center",true}}}});
    return json::Value::Object{
        {"schema","dynsyslab.compute-job/v1"},{"kind","postprocess"},
        {"label","native-pca-family-test"},
        {"sources",json::Value::Array{source.generic_u8string()}},
        {"methods",std::move(methods)},
        {"archive",json::Value::Object{{"enabled",true}}}};
}

}  // namespace

int main() {
    std::filesystem::path root;
    try {
        root=std::filesystem::temp_directory_path()/(
            "dynsyslab-pca-postprocess-"+std::to_string(
                std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(root);
        const auto result=dynsyslab::workflow::executeJob(
            job(sourceArchive(root)),root/"archives",{});
        require(result.archives.size()==3,"PCA-family job did not publish three archives");
        const auto bagged=core::readNpz(result.archives[0]);
        require(bagged.contains("components") && bagged.contains("scores")
            && bagged.at("components").shape==std::vector<std::size_t>({1,2}),
            "bagged PCA archive contract is incomplete");
        const auto bagMetadata=json::parse(bagged.metadataJson());
        require(bagMetadata.at("seed").number()==42.0
            && !bagMetadata.at("execution_backend").at("python_required").boolean(true),
            "bagged PCA provenance is incomplete");
        const auto robust=core::readNpz(result.archives[1]);
        require(robust.at("low_rank").shape==std::vector<std::size_t>({80,2})
            && robust.at("sparse").shape==std::vector<std::size_t>({80,2}),
            "robust PCA archive contract is incomplete");
        const auto robustMetadata=json::parse(robust.metadataJson());
        require(!robustMetadata.at("converged").boolean(true)
            && robustMetadata.at("iterations").number()==1000.0
            && robustMetadata.at("relative_residual").number()>0.0,
                "robust PCA nonconvergence diagnostics were not retained");
        const auto hankel=core::readNpz(result.archives[2]);
        require(hankel.at("left_singular_vectors").shape
                    ==std::vector<std::size_t>({66,3})
            && hankel.at("right_singular_vectors").shape
                    ==std::vector<std::size_t>({3,16}),
            "Hankel POD archive geometry is wrong");
        std::filesystem::remove_all(root);
        std::cout<<"dynsyslab_pca_postprocess_tests: PASS\n"; return 0;
    } catch (const std::exception& error) {
        if (!root.empty()) { std::error_code ignored; std::filesystem::remove_all(root,ignored); }
        std::cerr<<"dynsyslab_pca_postprocess_tests: FAIL: "<<error.what()<<'\n'; return 1;
    }
}
