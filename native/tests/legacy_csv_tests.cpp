#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/LegacyCsv.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}

void write(const std::filesystem::path& path,const std::string& text){std::ofstream stream(path,std::ios::binary|std::ios::trunc);stream<<text;if(!stream)throw std::runtime_error("fixture write failed");}
}

int main(){try{
    const auto root=std::filesystem::temp_directory_path()/"dynsyslab-legacy-csv-test";std::error_code ignored;std::filesystem::remove_all(root,ignored);std::filesystem::create_directories(root);
    const auto digest=root/"digest.txt";write(digest,"abc");require(dynsyslab::core::sha256File(digest)=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad","SHA-256 file identity");
    const auto source=root/"legacy.csv";write(source,
        "Var1,Var2_1,Var2_2,Var2_3\r\n"
        "dynbatch,2,,\r\n"
        "dynsys,1,Duffing,0.1,2,2\r\n"
        "time,0,0.1,0.2\r\nctrl,random,uniform,\r\nu,0,1,1.5\r\n"
        "dynclst,1,,\r\ninst,1,,\r\ntime,0,0.1,0.2\r\n"
        "sol,0.2,0.3,0.4\r\n,0.9,0.8,0.7\r\n"
        "dynsys,2,Duffing,0.1,2,2\r\n"
        "time,0,0.1,0.2\r\nctrl,random,uniform,\r\nu,0,2,2.5\r\n"
        "dynclst,1,,\r\ninst,1,,\r\ntime,0,0.1,0.2\r\n"
        "sol,-0.4,-0.3,-0.2\r\n,-0.1,-0.2,-0.3\r\n");
    const auto decoded=dynsyslab::core::readLegacyDynBatchCsv(source);require(decoded.systems.size()==2&&decoded.trajectoryCount()==2,"legacy hierarchy decode");
    require(decoded.systems[0].trajectories[0].states==std::vector<double>({0.2,0.9,0.3,0.8,0.4,0.7}),"component rows become sample-major states");
    require(decoded.systems[0].controlValues==std::vector<double>({0,1,1.5}),"control history decode");
    dynsyslab::core::LegacyCsvImportOptions options;options.label="native legacy Duffing";options.inputScale=0.5;options.coordinateSystem="dimensionless";options.mappedSystem.identifier="duffing";options.mappedSystem.preset="double_well_reference";
    const auto imported=dynsyslab::core::importLegacyDynBatchCsv(source,root/"archives",options);require(imported.systemCount==1&&imported.archives.size()==2,"identical mapped dynamics grouped in one system run");
    const auto archive=dynsyslab::core::readNpz(imported.archives.front());require(archive.at("values").shape==std::vector<std::size_t>({3,2}),"canonical imported state shape");
    require(archive.at("control_values").doubles()==std::vector<double>({0,0.5,0.75}),"scaled sample control history");require(archive.at("control_interval_values").doubles()==std::vector<double>({0.5,0.75}),"legacy initializer removed only for regression inputs");
    const auto metadata=dynsyslab::core::json::parse(archive.metadataJson());require(metadata.at("system").string()=="duffing"&&metadata.at("execution_backend").at("python_required").boolean(true)==false,"explicit native mapping and provenance");
    const auto hierarchy=dynsyslab::core::json::load(imported.hierarchyManifest.wstring());require(hierarchy.at("systems").array().size()==1&&hierarchy.at("systems").array()[0].at("clusters").array().size()==2,"batch/system/cluster/trajectory import hierarchy");
    require(std::filesystem::is_regular_file(imported.importManifest),"import manifest");
    const auto malformed=root/"malformed.csv";write(malformed,"header,x\ndynbatch,1\n,1\n");bool rejected=false;try{(void)dynsyslab::core::readLegacyDynBatchCsv(malformed);}catch(const std::invalid_argument& error){rejected=std::string(error.what()).find("Unexpected continuation")!=std::string::npos;}require(rejected,"strict unattached continuation rejection");
    std::filesystem::remove_all(root,ignored);std::cout<<"dynsyslab_legacy_csv_tests: PASS\n";return 0;
}catch(const std::exception& error){std::cerr<<"dynsyslab_legacy_csv_tests: FAIL: "<<error.what()<<'\n';return 1;}}
