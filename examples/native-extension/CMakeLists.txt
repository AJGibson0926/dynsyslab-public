cmake_minimum_required(VERSION 3.24)
project(DynSysLabOscillatorExtension LANGUAGES CXX)

if(NOT DYNSYSLAB_SDK_DIR)
  message(FATAL_ERROR "Set DYNSYSLAB_SDK_DIR to the installed DynSysLab root")
endif()

include("${DYNSYSLAB_SDK_DIR}/lib/cmake/DynSysLab/DynSysLabExtensions.cmake")
dynsyslab_add_extension(dynsyslab_example_oscillator
  OUTPUT_NAME "dynsyslab-example-oscillator"
  SOURCES oscillator_extension.cpp)
