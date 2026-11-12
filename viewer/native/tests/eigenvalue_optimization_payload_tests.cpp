#include "dynsyslab_viewer/NativePayload.hpp"
#include "../src/ArchiveCatalog.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition,const char* message){
    if(!condition)throw std::runtime_error(message);
}
}

int main(){const auto path=std::filesystem::temp_directory_path()
        /"dynsyslab-eigenvalue-optimization-payload.npz";
    try{dynsyslab::core::NpzRecord archive;
        archive.arrays["initial_discrete_eigenvalues"]=
            dynsyslab::core::NpyArray::fromComplex({2},{{.7,.1},{.7,-.1}});
        archive.arrays["discrete_eigenvalues"]=
            dynsyslab::core::NpyArray::fromComplex({2},{{.8,.2},{.8,-.2}});
        archive.arrays["optimization_metrics"]=
            dynsyslab::core::NpyArray::fromDoubles({1,5},{.01,.1,12.,30.,1.});
        archive.arrays["metadata_json"]=dynsyslab::core::NpyArray::fromString(
            R"({"format":"dynsyslab.derived/v1","label":"Optimized spectrum","method":"trajectory-eigenvalue-optimization"})");
        dynsyslab::core::writeNpzAtomic(path,archive);
        const auto payload=dsl::buildNativePayload(path,"eigenvalue_optimization");
        require(payload.figures.size()==3,"optimization dashboard is incomplete");
        require(payload.figures[0].square&&payload.figures[0].traces.size()==3,
            "initial/optimized spectrum comparison is malformed");
        require(payload.figures[1].traces[0].y[0]==.1,
            "relative residual was not displayed");
        require(payload.warnings.size()==1,
            "spectrum-only warning was not displayed");
        dsl::CatalogItem item;item.source=path.u8string();
        dsl::enrichArchiveCatalogItem(item,path);
        require(item.kind=="koopman_spectrum"&&!item.views.empty()
            &&item.views.front()=="eigenvalue_optimization",
            "optimizer archive lacks its native catalog route");
        require(item.properties.at("eigenfunctions constructed")=="no",
            "optimizer spectrum could be mistaken for a KROM");
        std::error_code ignored;std::filesystem::remove(path,ignored);
        std::cout<<"dynsyslab_eigenvalue_optimization_payload_tests: PASS\n";
        return 0;
    }catch(const std::exception& error){std::error_code ignored;
        std::filesystem::remove(path,ignored);std::cerr
        <<"dynsyslab_eigenvalue_optimization_payload_tests: FAIL: "
        <<error.what()<<'\n';return 1;}
}
