#include "dynsyslab_compute/Workbench.hpp"

#include "dynsyslab_viewer/RuntimePaths.hpp"

#include <shellapi.h>

namespace dsl::compute {

Arguments parseArguments() {
    int count=0;
    wchar_t** values=CommandLineToArgvW(GetCommandLineW(),&count);
    Arguments result;
    for(int index=1;index<count;++index){
        const std::wstring key=values[index];
        if(key==L"--postprocessor"){result.openPostprocessor=true;continue;}
        if(key==L"--smoke-test"){result.smokeTest=true;continue;}
        if(key==L"--background-check"){result.backgroundCheck=true;continue;}
        if(index+1>=count)break;
        const std::filesystem::path value=values[++index];
        if(key==L"--project-root")result.backend.projectRoot=value;
        else if(key==L"--archive-root")result.backend.archiveRoot=value;
        else if(key==L"--registry")result.backend.registry=value;
        else if(key==L"--catalog")result.backend.catalog=value;
        else if(key==L"--cache")result.backend.cache=value;
        else if(key==L"--source")result.initialSource=value;
    }
    LocalFree(values);
    const auto paths=dsl::discoverRuntimePaths(L"Compute",result.backend.projectRoot,
        result.backend.archiveRoot,result.backend.cache);
    result.backend.projectRoot=paths.projectRoot;
    result.backend.archiveRoot=paths.archiveRoot;
    result.backend.cache=paths.cacheRoot;
    return result;
}

}  // namespace dsl::compute
