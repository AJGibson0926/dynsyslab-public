cmake_minimum_required(VERSION 3.24)
project(RpeForcedEdmdcValidation LANGUAGES CXX)
add_executable(rpe_edmdc_validation experiment.cpp)
target_compile_features(rpe_edmdc_validation PRIVATE cxx_std_17)
target_include_directories(rpe_edmdc_validation PRIVATE
  "${DYNSYSLAB_SOURCE_DIR}/native/include"
  "${DYNSYSLAB_SOURCE_DIR}/workflow/native/include")
foreach(module foundation eigenfunctions koopman spectral operations workflow)
  target_link_libraries(rpe_edmdc_validation PRIVATE "${DYNSYSLAB_SUITE_BUILD}/lib/dynsyslab_${module}.lib")
endforeach()
